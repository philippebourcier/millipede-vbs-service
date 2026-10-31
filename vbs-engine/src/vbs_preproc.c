#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_preproc.h"

/*
 * Cycle-slip detection, spec S4.2. See the header for why all three detectors are needed.
 */

/*
 * Geometry-free jump threshold, in metres: a noise floor plus an allowance that grows with
 * the time since the previous epoch.
 *
 *   threshold(dt) = GF_BASE_M + GF_RATE_MPS * dt
 *
 * The rate term is the point. Between epochs GF moves by the change in slant ionosphere, so
 * the allowance has to scale with how long the ionosphere had to move. RTKLIB's conventional
 * 0.05 m fixed threshold assumes roughly 1 Hz data; this dataset is sampled at 30 s, which
 * gives the ionosphere thirty times as long, and applying 0.05 m to it flags the ordinary
 * ionospheric tail as slips. Measured on an F9P base's quiet day: |dGF| between consecutive 30 s
 * epochs is a continuous distribution with 94% of samples under 0.1 m and the unambiguous
 * slip population above 2 m, so a 0.05 m cut declared 5.7% of all satellite-epochs slipped.
 *
 * At dt = 30 s this gives 0.15 m, which sits above that tail and below one cycle of L1
 * (0.19 m) or L2 (0.24 m). At dt = 1 s it gives 0.034 m, close to the conventional value --
 * which is the right behaviour, because the live engine will see 1 Hz base streams where GF
 * detection is genuinely more sensitive than it can be on this 30 s validation data.
 *
 * Consequence worth stating plainly: at 30 s sampling, GF cannot reliably resolve sub-cycle
 * slips, so Melbourne-Wubbena carries more of the load here than it would live.
 */
#define GF_BASE_M		0.03
#define GF_RATE_MPS		0.004

/*
 * Melbourne-Wubbena jump threshold, in wide-lane cycles, applied to the deviation from the
 * running mean rather than to a single epoch-to-epoch difference.
 *
 * MW carries code noise: tens of centimetres against a 0.86 m wide-lane wavelength, so a
 * single epoch is worth roughly half a cycle of noise. Averaging beats that down, and a slip
 * shows up as a step in the mean. The threshold is the larger of a fixed floor and a multiple
 * of the running standard deviation, so a quiet satellite gets a tight test and a noisy one
 * does not produce a false slip every other epoch.
 */
#define MW_THRESHOLD_CYC	1.5
#define MW_SIGMA_FACTOR		4.0
#define MW_MIN_SAMPLES		4	/* before the running mean is worth testing against */
/*
 * While the running mean is still warming up after a slip, test against the last value with
 * a fixed threshold instead of not testing at all. T30's Unicore stations produce phase
 * jumps of tens to hundreds of metres, equal on every band, 30 s apart, on low satellites;
 * the geometry-free test cannot see a band-common jump, and with the MW test silent for
 * four samples after the first jump the second went unflagged -- the filter then rejected
 * that satellite's phase for the rest of the arc while holding its fix (D49). 5 cycles is
 * 4.3 m of wide lane: far above MW noise on any receiver, far below the jumps in question.
 */
#define MW_WARMUP_CYC		5.0

#define DEFAULT_GAP_FACTOR	3.0

/* Per-satellite continuity state. */
struct sat_state {
	int valid;
	long arc_id;			/* advances on every slip; see vbs_preproc_arc_ids */
	gtime_t last;
	double gf;			/* last geometry-free value, metres */
	int gf_valid;
	int gf_b1, gf_b2;		/* the band pair GF was formed from */
	double mw_mean, mw_m2;		/* running mean and sum of squared deviations */
	long mw_n;
};

struct vbs_preproc {
	double sampling_s;
	double gap_factor;
	struct sat_state sat[MAXSAT];
	struct vbs_preproc_stats stats;
};

struct vbs_preproc *vbs_preproc_new(double sampling_s, double gap_factor) {
	struct vbs_preproc *this = (struct vbs_preproc *)calloc(1, sizeof *this);

	if (this == NULL)
		return NULL;
	this->sampling_s = sampling_s > 0.0 ? sampling_s : 30.0;
	this->gap_factor = gap_factor > 0.0 ? gap_factor : DEFAULT_GAP_FACTOR;
	return this;
}

void vbs_preproc_free(struct vbs_preproc *this) {
	free(this);
}

/* Slot carrying a given band, or -1. */
static int slot_of_band(const obsd_t *o, int band) {
	int j;

	for (j = 0; j < NFREQ + NEXOBS; j++)
		if (o->code[j] != CODE_NONE && vbs_signal_band(o->code[j]) == band)
			return j;
	return -1;
}

/*
 * The two lowest bands this satellite carries with both phase and code, which is what GF and
 * MW both need. Returns 1 if two were found.
 */
static int usable_band_pair(const obsd_t *o, const nav_t *nav, int *b1, int *b2,
	int *s1, int *s2) {
	int band, found = 0, j;

	for (band = 1; band < VBS_NBAND_MAX && found < 2; band++) {
		if ((j = slot_of_band(o, band)) < 0)
			continue;
		if (o->L[j] == 0.0 || o->P[j] == 0.0)
			continue;
		if (sat2freq(o->sat, o->code[j], nav) == 0.0)
			continue;
		if (found == 0) {
			*b1 = band;
			*s1 = j;
		} else {
			*b2 = band;
			*s2 = j;
		}
		found++;
	}
	return found == 2;
}

int vbs_preproc_epoch(struct vbs_preproc *this, const obsd_t *obs, int n, const nav_t *nav,
	struct vbs_slip *out, int maxout) {
	int i, nout = 0;

	if (this == NULL || obs == NULL)
		return -1;
	this->stats.epochs++;

	for (i = 0; i < n; i++) {
		struct sat_state *st;
		double f1, f2, lam_wl, gf, mw, dt, thresh, dev, sigma;
		int b1 = 0, b2 = 0, s1 = 0, s2 = 0, dual, j;
		unsigned reason = 0;

		if (obs[i].sat <= 0 || obs[i].sat > MAXSAT)
			continue;
		st = &this->sat[obs[i].sat - 1];
		this->stats.sat_epochs++;

		dual = usable_band_pair(&obs[i], nav, &b1, &b2, &s1, &s2);
		if (dual)
			this->stats.dual_band++;
		else
			this->stats.single_band++;

		/* Lock time: taken as a trigger, never as the only detector (see the header). */
		for (j = 0; j < NFREQ + NEXOBS; j++)
			if (obs[i].code[j] != CODE_NONE && (obs[i].LLI[j] & 1))
				reason |= VBS_SLIP_LLI;

		/* A long enough gap ends the arc regardless of what the combinations say. */
		dt = st->valid ? fabs(timediff(obs[i].time, st->last)) : 0.0;
		if (!st->valid)
			reason |= VBS_SLIP_NEW;
		else if (dt > this->gap_factor * this->sampling_s)
			reason |= VBS_SLIP_GAP;

		if (dual) {
			f1 = sat2freq(obs[i].sat, obs[i].code[s1], nav);
			f2 = sat2freq(obs[i].sat, obs[i].code[s2], nav);

			/* Geometry-free, in metres. */
			gf = obs[i].L[s1] * CLIGHT / f1 - obs[i].L[s2] * CLIGHT / f2;

			/*
			 * Melbourne-Wubbena, in wide-lane cycles: wide-lane phase minus
			 * narrow-lane code. Geometry, clocks, troposphere and ionosphere all
			 * cancel, leaving the wide-lane ambiguity plus code noise.
			 *
			 *   phi_wl = lambda_wl * (L1 - L2)          [L in cycles]
			 *   p_nl   = (f1 P1 + f2 P2) / (f1 + f2)
			 */
			lam_wl = CLIGHT / (f1 - f2);
			mw = (lam_wl * (obs[i].L[s1] - obs[i].L[s2])
			     - (f1 * obs[i].P[s1] + f2 * obs[i].P[s2]) / (f1 + f2)) / lam_wl;

			if (st->gf_valid && !(reason & (VBS_SLIP_GAP | VBS_SLIP_NEW)) &&
			    st->gf_b1 == b1 && st->gf_b2 == b2) {
				/* Noise floor plus an ionospheric allowance that grows with dt. */
				thresh = GF_BASE_M + GF_RATE_MPS *
				    (dt > 0.0 ? dt : this->sampling_s);
				if (fabs(gf - st->gf) > thresh)
					reason |= VBS_SLIP_GF;
			}

			if (st->mw_n >= 1 && !(reason & (VBS_SLIP_GAP | VBS_SLIP_NEW))) {
				if (st->mw_n >= MW_MIN_SAMPLES) {
					sigma = sqrt(st->mw_m2 / (st->mw_n - 1));
					thresh = MW_SIGMA_FACTOR * sigma;
					if (thresh < MW_THRESHOLD_CYC)
						thresh = MW_THRESHOLD_CYC;
				} else {
					thresh = MW_WARMUP_CYC;
				}
				dev = fabs(mw - st->mw_mean);
				if (dev > thresh)
					reason |= VBS_SLIP_MW;
			}

			/*
			 * Update continuity. A slip restarts both detectors: keeping the old
			 * running mean across a step is exactly how a slip detector goes blind
			 * to the next one.
			 */
			st->gf = gf;
			st->gf_valid = 1;
			st->gf_b1 = b1;
			st->gf_b2 = b2;
			if (reason & (VBS_SLIP_GF | VBS_SLIP_MW | VBS_SLIP_LLI |
			    VBS_SLIP_GAP | VBS_SLIP_NEW)) {
				st->mw_mean = mw;
				st->mw_m2 = 0.0;
				st->mw_n = 1;
			} else {
				double d = mw - st->mw_mean;

				st->mw_n++;
				st->mw_mean += d / st->mw_n;
				st->mw_m2 += d * (mw - st->mw_mean);
			}
		} else {
			st->gf_valid = 0;
			st->mw_n = 0;
			st->mw_mean = st->mw_m2 = 0.0;
		}

		st->valid = 1;
		st->last = obs[i].time;

		if (reason == 0)
			continue;

		st->arc_id++;
		this->stats.slips++;
		if (reason & VBS_SLIP_LLI) this->stats.slip_lli++;
		if (reason & VBS_SLIP_GF) this->stats.slip_gf++;
		if (reason & VBS_SLIP_MW) this->stats.slip_mw++;
		if (reason & VBS_SLIP_GAP) this->stats.slip_gap++;
		if (reason & VBS_SLIP_NEW) this->stats.slip_new++;
		/*
		 * The comparison spec S4.2 exists for: a slip the combinations found and the
		 * receiver's own flag did not. Arc starts and gaps are excluded -- nobody claims
		 * LLI should flag those.
		 */
		if ((reason & (VBS_SLIP_GF | VBS_SLIP_MW)) && !(reason & VBS_SLIP_LLI) &&
		    !(reason & (VBS_SLIP_NEW | VBS_SLIP_GAP)))
			this->stats.slip_missed_by_lli++;
		if ((reason & VBS_SLIP_LLI) && !(reason & (VBS_SLIP_GF | VBS_SLIP_MW)) &&
		    !(reason & (VBS_SLIP_NEW | VBS_SLIP_GAP)))
			this->stats.slip_lli_only++;

		if (nout < maxout) {
			out[nout].sat = obs[i].sat;
			out[nout].band = 0;
			out[nout].reason = reason;
			nout++;
		}
	}
	return nout;
}

void vbs_preproc_get_stats(const struct vbs_preproc *this, struct vbs_preproc_stats *out) {
	*out = this->stats;
}

void vbs_preproc_arc_ids(const struct vbs_preproc *this, long *out) {
	int i;

	for (i = 0; i < MAXSAT; i++)
		out[i] = this != NULL ? this->sat[i].arc_id : 0;
}

double vbs_preproc_mw_sigma(const struct vbs_preproc *this, int sat) {
	const struct sat_state *st;

	if (this == NULL || sat < 1 || sat > MAXSAT)
		return NAN;
	st = &this->sat[sat - 1];
	if (st->mw_n < MW_MIN_SAMPLES)
		return NAN;
	return sqrt(st->mw_m2 / (st->mw_n - 1));
}

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_ddres.h"

/*
 * Double-difference residuals. See the header for what is measured and why the arc mean is
 * removed.
 */

#define MAX_SAMPLES	400000		/* ~2880 epochs x 30 sats x 4 groups, generous */
#define MAX_ARCS	8192

/* Constellations that can carry a reference satellite of their own. */
#define VBS_NSYSREF	4

static int sysref_index(int sys) {
	switch (sys) {
	case SYS_GPS: return 0;
	case SYS_GAL: return 1;
	case SYS_CMP: return 2;
	case SYS_GLO: return 3;
	}
	return -1;
}

/* One double-difference sample, tagged with the arc it belongs to. */
struct dd_sample {
	int arc;
	double v;		/* metres */
	float elev;		/* degrees, of the non-reference satellite */
	/*
	 * The reference satellite's elevation, also in degrees. Needed because these are
	 * *double* differences: the elevation-dependent shape of a residual is
	 * f(el_sat) - f(el_ref), and the reference satellite's own elevation varies within an
	 * arc. Regressing on the non-reference elevation alone leaves that variation
	 * unmodelled, which is what made the first attempt at D35's antenna test
	 * under-specified and its r² uninterpretable (D36).
	 */
	float elev_ref;
};

/*
 * An arc is one (reference satellite, other satellite, group) combination, alive as long as
 * both satellites keep being observed without a slip. A break starts a new arc, because the
 * ambiguity that the de-meaning removes is only constant within an arc.
 */
struct dd_arc {
	int sat, refsat;
	int sys, band;
	enum vbs_dd_comb comb;
	gtime_t last;
	double last_v;		/* previous double difference, for the data-driven slip test */
	int n;
	double sum;
};

/*
 * Largest change in a double difference, between consecutive epochs of one arc, that is
 * accepted as real rather than a cycle slip.
 *
 * At a zero-length baseline the geometry cancels exactly, so from one 30 s epoch to the next
 * the double difference should move only by the change in the differential atmosphere --
 * millimetres to a couple of centimetres. One cycle is 0.19 m at L1 and 0.24 m at L2, so
 * 0.15 m separates "a slip happened" from "the ionosphere moved" with room on both sides.
 *
 * This exists because LLI is not trustworthy here. The dataset's loss-of-lock flags are
 * derived from MSM lock-time indicators (data/README.md S4), and measurement shows jumps of
 * 50-120 m between consecutive epochs with LLI clear on all four terms. Spec S4.2 already
 * requires slip detection from the data itself -- geometry-free and Melbourne-Wubbena
 * combinations *plus* lock time -- and this is the same lesson arriving early: an arc that
 * trusts the flag alone silently averages across a slip, and every sample in it is then
 * wrong by the slip size.
 */
#define DD_SLIP_THRESHOLD_M	0.15

struct vbs_ddres {
	double elev_mask;
	struct dd_sample *sample;
	int nsample;
	struct dd_arc *arc;
	int narc;
	int nepoch, nskipped;
};

struct vbs_ddres *vbs_ddres_new(double elev_mask_deg) {
	struct vbs_ddres *this = (struct vbs_ddres *)calloc(1, sizeof *this);

	if (this == NULL)
		return NULL;
	this->elev_mask = elev_mask_deg;
	this->sample = (struct dd_sample *)malloc(MAX_SAMPLES * sizeof(struct dd_sample));
	this->arc = (struct dd_arc *)malloc(MAX_ARCS * sizeof(struct dd_arc));
	if (this->sample == NULL || this->arc == NULL) {
		vbs_ddres_free(this);
		return NULL;
	}
	return this;
}

void vbs_ddres_free(struct vbs_ddres *this) {
	if (this == NULL)
		return;
	free(this->sample);
	free(this->arc);
	free(this);
}

/* Find the record for a satellite, or NULL. */
static const obsd_t *find_sat(const obsd_t *o, int n, int sat) {
	int i;

	for (i = 0; i < n; i++)
		if (o[i].sat == sat)
			return &o[i];
	return NULL;
}

/* Slot carrying a given band for this record, or -1. */
static int slot_of_band(const obsd_t *o, int band) {
	int j;

	for (j = 0; j < NFREQ + NEXOBS; j++)
		if (o->code[j] != CODE_NONE && vbs_signal_band(o->code[j]) == band)
			return j;
	return -1;
}

/*
 * Phase in metres for one band, or NAN if unusable. Converting to metres here is what makes
 * the geometry-free and ionosphere-free combinations below expressible at all, and it makes
 * every reported number a distance rather than a count of cycles.
 *
 * NAN, not 0.0, for "no observation": a single difference built from a missing term is not a
 * small number, it is the *other* station's raw range -- tens of millions of metres, which
 * then sails through any "is the difference zero" check and lands in the statistics. That bug
 * was real here (GPS L2C is not transmitted by the whole constellation, so the target
 * frequently has a band the master does not, or vice versa). NAN propagates through the
 * arithmetic and is caught once, at the end, where the sample is formed.
 */
static double phase_m(const obsd_t *o, int band, const nav_t *nav, int *slot_out) {
	double f;
	int j;

	if ((j = slot_of_band(o, band)) < 0)
		return NAN;
	if (o->L[j] == 0.0)
		return NAN;
	if ((f = sat2freq(o->sat, o->code[j], nav)) == 0.0)
		return NAN;
	if (slot_out != NULL)
		*slot_out = j;
	return o->L[j] * CLIGHT / f;
}

/*
 * Locate or open the arc for this (refsat, sat, group). A time gap, a reported slip, or an
 * unexplained jump in the value itself starts a new one.
 */
static int get_arc(struct vbs_ddres *this, int refsat, int sat, int sys, int band,
	enum vbs_dd_comb comb, gtime_t t, int slip, double v) {
	int i;

	for (i = this->narc - 1; i >= 0; i--) {
		struct dd_arc *a = &this->arc[i];

		if (a->sat != sat || a->refsat != refsat || a->comb != comb || a->band != band)
			continue;
		/*
		 * Same pair and group. Continue the arc only if the observation is contiguous
		 * with the last one, neither end reported a slip, and the value itself did not
		 * jump. 90 s tolerates a single missing 30 s epoch without pretending the
		 * ambiguity survived a long outage.
		 */
		if (!slip && fabs(timediff(t, a->last)) <= 90.0 &&
		    (a->n == 0 || fabs(v - a->last_v) <= DD_SLIP_THRESHOLD_M)) {
			a->last = t;
			a->last_v = v;
			return i;
		}
		break;		/* stale, flagged or jumped: fall through and open a new arc */
	}
	if (this->narc >= MAX_ARCS)
		return -1;
	struct dd_arc *a = &this->arc[this->narc];
	a->sat = sat;
	a->refsat = refsat;
	a->sys = sys;
	a->band = band;
	a->comb = comb;
	a->last = t;
	a->last_v = v;
	a->n = 0;
	a->sum = 0.0;
	return this->narc++;
}

static void add_sample(struct vbs_ddres *this, int arc, double v, double elev,
	double elev_ref) {
	if (arc < 0 || this->nsample >= MAX_SAMPLES)
		return;
	this->sample[this->nsample].arc = arc;
	this->sample[this->nsample].v = v;
	this->sample[this->nsample].elev = (float)elev;
	this->sample[this->nsample].elev_ref = (float)elev_ref;
	this->nsample++;
	this->arc[arc].n++;
	this->arc[arc].sum += v;
}

/* Did either end flag a loss of lock on this band? */
static int band_slip(const obsd_t *a, const obsd_t *b, int band) {
	int ja = slot_of_band(a, band), jb = slot_of_band(b, band);

	return (ja >= 0 && (a->LLI[ja] & 1)) || (jb >= 0 && (b->LLI[jb] & 1));
}

int vbs_ddres_add_epoch(struct vbs_ddres *this,
	const obsd_t *a, int na, const obsd_t *b, int nb,
	const double *rs, const double *rr, const nav_t *nav) {
	double pos[3], e[3], azel[2];
	double elev[MAXOBS], refelev[VBS_NSYSREF];
	int refidx[VBS_NSYSREF], i, k, used = 0, anyref = 0;
	gtime_t t;

	if (na <= 0 || nb <= 0)
		return 0;
	t = a[0].time;
	ecef2pos(rr, pos);

	for (k = 0; k < VBS_NSYSREF; k++) {
		refidx[k] = -1;
		refelev[k] = -1.0;
	}

	/*
	 * Elevation at the shared position, and a reference satellite *per constellation*: the
	 * highest one above the mask that both sets observe. A high reference keeps the double
	 * difference's noise down.
	 *
	 * Per constellation, not one overall: a double difference only cancels the receiver
	 * clock if both satellites are in the same system, so a single global reference would
	 * silently discard every satellite of every *other* system in that epoch -- which is
	 * exactly what an earlier version of this did, starving GPS in the epochs where a
	 * Galileo satellite happened to be highest.
	 */
	for (i = 0; i < na && i < MAXOBS; i++) {
		elev[i] = -1.0;
		if (rs == NULL || (rs[6 * i] == 0.0 && rs[6 * i + 1] == 0.0))
			continue;
		if (geodist(rs + 6 * i, rr, e) <= 0.0)
			continue;
		satazel(pos, e, azel);
		elev[i] = azel[1] * R2D;
		if (elev[i] < this->elev_mask)
			continue;
		if (find_sat(b, nb, a[i].sat) == NULL)
			continue;
		if ((k = sysref_index(satsys(a[i].sat, NULL))) < 0)
			continue;
		if (elev[i] > refelev[k]) {
			refelev[k] = elev[i];
			refidx[k] = i;
			anyref = 1;
		}
	}
	if (!anyref) {
		this->nskipped++;
		return 0;
	}

	for (i = 0; i < na && i < MAXOBS; i++) {
		const obsd_t *bi, *aref, *bref;
		int sys, bands[VBS_NBAND_MAX], nband = 0, ri;

		if (elev[i] < this->elev_mask)
			continue;
		sys = satsys(a[i].sat, NULL);
		if ((ri = sysref_index(sys)) < 0 || refidx[ri] < 0 || refidx[ri] == i)
			continue;
		if ((bi = find_sat(b, nb, a[i].sat)) == NULL)
			continue;
		aref = &a[refidx[ri]];
		bref = find_sat(b, nb, aref->sat);
		if (bref == NULL)
			continue;

		/* Bands this satellite has at both ends, and at the reference too. */
		for (int bd = 1; bd < VBS_NBAND_MAX; bd++) {
			if (slot_of_band(&a[i], bd) < 0 || slot_of_band(bi, bd) < 0)
				continue;
			if (slot_of_band(aref, bd) < 0 || slot_of_band(bref, bd) < 0)
				continue;
			if (nband < VBS_NBAND_MAX)
				bands[nband++] = bd;
		}
		if (nband == 0)
			continue;

		/* Per-band double differences. */
		for (int bk = 0; bk < nband; bk++) {
			double sd_i, sd_r, dd;
			int band = bands[bk], arc, slip;

			sd_i = phase_m(&a[i], band, nav, NULL) - phase_m(bi, band, nav, NULL);
			sd_r = phase_m(aref, band, nav, NULL) - phase_m(bref, band, nav, NULL);
			if (isnan(sd_i) || isnan(sd_r))
				continue;
			dd = sd_i - sd_r;
			slip = band_slip(&a[i], bi, band) || band_slip(aref, bref, band);
			arc = get_arc(this, aref->sat, a[i].sat, sys, band, VBS_DD_BAND, t, slip, dd);
			add_sample(this, arc, dd, elev[i], elev[refidx[ri]]);
			used = 1;
		}

		/*
		 * Two-band combinations, on the lowest two bands this satellite pair shares.
		 * Each is built from the four single-difference terms (this satellite and the
		 * reference, at each end), so the combination is formed per station and then
		 * double-differenced -- not the other way round, which would be equivalent here
		 * but is easier to get wrong.
		 */
		if (nband >= 2) {
			int b1 = bands[0], b2 = bands[1];
			double f1, f2, c1, c2, sd1_i, sd2_i, sd1_r, sd2_r, dd;
			int arc, slip, s1, s2;

			s1 = slot_of_band(&a[i], b1);
			s2 = slot_of_band(&a[i], b2);
			f1 = sat2freq(a[i].sat, a[i].code[s1], nav);
			f2 = sat2freq(a[i].sat, a[i].code[s2], nav);
			if (f1 == 0.0 || f2 == 0.0 || f1 == f2)
				continue;

			/* Single differences (VBS - real base) per band, this sat and the ref. */
			sd1_i = phase_m(&a[i], b1, nav, NULL) - phase_m(bi, b1, nav, NULL);
			sd2_i = phase_m(&a[i], b2, nav, NULL) - phase_m(bi, b2, nav, NULL);
			sd1_r = phase_m(aref, b1, nav, NULL) - phase_m(bref, b1, nav, NULL);
			sd2_r = phase_m(aref, b2, nav, NULL) - phase_m(bref, b2, nav, NULL);
			if (isnan(sd1_i) || isnan(sd2_i) || isnan(sd1_r) || isnan(sd2_r))
				continue;

			slip = band_slip(&a[i], bi, b1) || band_slip(&a[i], bi, b2) ||
			    band_slip(aref, bref, b1) || band_slip(aref, bref, b2);

			/* Geometry-free: L(b1) - L(b2). Geometry and clocks cancel. */
			dd = (sd1_i - sd2_i) - (sd1_r - sd2_r);
			arc = get_arc(this, aref->sat, a[i].sat, sys, 0, VBS_DD_GEOMFREE, t, slip, dd);
			add_sample(this, arc, dd, elev[i], elev[refidx[ri]]);

			/* Ionosphere-free: (f1^2 L1 - f2^2 L2) / (f1^2 - f2^2). */
			c1 = f1 * f1 / (f1 * f1 - f2 * f2);
			c2 = -f2 * f2 / (f1 * f1 - f2 * f2);
			dd = (c1 * sd1_i + c2 * sd2_i) - (c1 * sd1_r + c2 * sd2_r);
			arc = get_arc(this, aref->sat, a[i].sat, sys, 0, VBS_DD_IONOFREE, t, slip, dd);
			add_sample(this, arc, dd, elev[i], elev[refidx[ri]]);
		}
	}
	if (used)
		this->nepoch++;
	else
		this->nskipped++;
	return used;
}

void vbs_ddres_dump_outliers(struct vbs_ddres *this, double threshold, int max_lines,
	FILE *fp) {
	char sid[8], rid[8], ts[40];
	int i, printed = 0;

	for (i = 0; i < this->nsample && printed < max_lines; i++) {
		struct dd_arc *a = &this->arc[this->sample[i].arc];
		double resid;

		if (a->n < 2)
			continue;
		resid = this->sample[i].v - a->sum / a->n;
		if (fabs(resid) <= threshold)
			continue;
		satno2id(a->sat, sid);
		satno2id(a->refsat, rid);
		time2str(a->last, ts, 0);
		fprintf(fp, "# OUTLIER %10.3f m  sat=%s ref=%s band=%d comb=%d "
		    "arc_n=%d elev=%.1f arc_end=%s\n",
		    resid, sid, rid, a->band, (int)a->comb, a->n, this->sample[i].elev, ts);
		printed++;
	}
}

static int cmp_double(const void *a, const void *b) {
	double x = *(const double *)a, y = *(const double *)b;

	return (x > y) - (x < y);
}

int vbs_ddres_finish(struct vbs_ddres *this, struct vbs_dd_result *out) {
	double *abs_v;
	int i, g;

	memset(out, 0, sizeof *out);
	out->nepoch = this->nepoch;
	out->nskipped = this->nskipped;
	if (this->nsample <= 0)
		return 0;

	if ((abs_v = (double *)malloc(this->nsample * sizeof(double))) == NULL)
		return -1;

	/* Group samples by (sys, band, comb), de-meaning each arc as we go. */
	for (i = 0; i < this->nsample; i++) {
		struct dd_arc *a = &this->arc[this->sample[i].arc];
		double resid;
		int found = -1;

		/*
		 * An arc of one sample carries no information once its mean is removed -- the
		 * residual is identically zero and would only dilute the statistics.
		 */
		if (a->n < 2)
			continue;
		resid = this->sample[i].v - a->sum / a->n;

		for (g = 0; g < out->nstat; g++)
			if (out->stat[g].sys == a->sys && out->stat[g].band == a->band &&
			    out->stat[g].comb == a->comb) {
				found = g;
				break;
			}
		if (found < 0) {
			if (out->nstat >= VBS_DD_MAX_STATS)
				continue;
			found = out->nstat++;
			out->stat[found].sys = a->sys;
			out->stat[found].band = a->band;
			out->stat[found].comb = a->comb;
			out->stat[found].nsample = 0;
			out->stat[found].rms = 0.0;
			out->stat[found].max = 0.0;
		}
		out->stat[found].rms += resid * resid;		/* sum of squares for now */
		out->stat[found].nsample++;
		if (fabs(resid) > out->stat[found].max)
			out->stat[found].max = fabs(resid);
	}

	/* Arc counts per group. */
	for (i = 0; i < this->narc; i++) {
		if (this->arc[i].n < 2)
			continue;
		for (g = 0; g < out->nstat; g++)
			if (out->stat[g].sys == this->arc[i].sys &&
			    out->stat[g].band == this->arc[i].band &&
			    out->stat[g].comb == this->arc[i].comb)
				out->stat[g].narc++;
	}

	/* RMS, and the 95th percentile of |residual| per group. */
	for (g = 0; g < out->nstat; g++) {
		int n = 0;

		if (out->stat[g].nsample > 0)
			out->stat[g].rms = sqrt(out->stat[g].rms / out->stat[g].nsample);
		for (i = 0; i < this->nsample; i++) {
			struct dd_arc *a = &this->arc[this->sample[i].arc];

			if (a->n < 2 || a->sys != out->stat[g].sys ||
			    a->band != out->stat[g].band || a->comb != out->stat[g].comb)
				continue;
			abs_v[n++] = fabs(this->sample[i].v - a->sum / a->n);
		}
		if (n > 0) {
			qsort(abs_v, n, sizeof(double), cmp_double);
			out->stat[g].p95 = abs_v[(int)(0.95 * (n - 1))];
		}
	}
	free(abs_v);
	return out->nstat;
}

static int shape_terms(const struct dd_sample *s, double *sin_dd, double *map_dd);

int vbs_ddres_pair(const struct vbs_ddres *ref, const struct vbs_ddres *test,
	struct vbs_dd_pair *out, int maxout) {
	struct shpacc {
		double sy, syy, sr, srr, syr, sm, smm, sym;
		long n;
	} *shp;
	int i, g, nout = 0;
	double *sr, *sc, *src, *ashp_r = NULL, *ashp_m = NULL;
	long *ashp_n = NULL;
	int *cnt;

	if (ref == NULL || test == NULL || ref->nsample != test->nsample)
		return -1;
	if ((shp = (struct shpacc *)calloc(maxout, sizeof *shp)) == NULL)
		return -1;
	/* Arc means of the regressors, so they can be de-meaned the same way c is. */
	ashp_r = (double *)calloc(ref->narc, sizeof *ashp_r);
	ashp_m = (double *)calloc(ref->narc, sizeof *ashp_m);
	ashp_n = (long *)calloc(ref->narc, sizeof *ashp_n);
	if (ashp_r == NULL || ashp_m == NULL || ashp_n == NULL) {
		free(shp);
		free(ashp_r);
		free(ashp_m);
		free(ashp_n);
		return -1;
	}
	for (i = 0; i < ref->nsample; i++) {
		double se, sm;

		if (!shape_terms(&ref->sample[i], &se, &sm))
			continue;
		ashp_r[ref->sample[i].arc] += se;
		ashp_m[ref->sample[i].arc] += sm;
		ashp_n[ref->sample[i].arc]++;
	}
	if ((sr = (double *)calloc(maxout, sizeof *sr)) == NULL)
		return -1;
	if ((sc = (double *)calloc(maxout, sizeof *sc)) == NULL ||
	    (src = (double *)calloc(maxout, sizeof *src)) == NULL ||
	    (cnt = (int *)calloc(maxout, sizeof *cnt)) == NULL) {
		free(sr);
		return -1;
	}
	memset(out, 0, maxout * sizeof *out);

	for (i = 0; i < ref->nsample; i++) {
		const struct dd_arc *a = &ref->arc[ref->sample[i].arc];
		const struct dd_arc *b = &test->arc[test->sample[i].arc];
		double r, t, c;
		int found = -1;

		if (a->n < 2 || b->n < 2)
			continue;
		/*
		 * Both series are de-meaned per arc, as the statistics are. The correction's
		 * own constant part is removed along with the ambiguity, which is correct: a
		 * constant offset in the correction is absorbed by the rover's ambiguity exactly
		 * as the ambiguity itself is, and only the time-varying part can hurt it.
		 */
		r = ref->sample[i].v - a->sum / a->n;
		t = test->sample[i].v - b->sum / b->n;
		c = r - t;

		for (g = 0; g < nout; g++)
			if (out[g].sys == a->sys && out[g].band == a->band &&
			    out[g].comb == a->comb) {
				found = g;
				break;
			}
		if (found < 0) {
			if (nout >= maxout)
				continue;
			found = nout++;
			out[found].sys = a->sys;
			out[found].band = a->band;
			out[found].comb = a->comb;
		}
		out[found].rms_ref += r * r;
		out[found].rms_test += t * t;
		sr[found] += c * c;			/* sum c^2 -> var(c) */
		sc[found] += r * c;			/* sum r*c -> cov      */
		src[found] += 0.0;
		cnt[found]++;

		/*
		 * The correction's own shape. Both c and the regressors are de-meaned *per arc*,
		 * exactly as vbs_ddres_shape() does to the residual -- c is a difference of two
		 * arc-de-meaned series, so correlating it against a raw regressor would compare a
		 * within-arc signal against a variance dominated by between-arc spread and dilute
		 * the correlation to nothing. (It did: 0.001 before this was fixed, on a series
		 * that is proportional to the regressor by construction.)
		 */
		{
			double se, sm;

			if (shape_terms(&ref->sample[i], &se, &sm)) {
				long am = ashp_n[ref->sample[i].arc];

				if (am >= 2) {
					double sem = se - ashp_r[ref->sample[i].arc] / am;
					double smm2 = sm - ashp_m[ref->sample[i].arc] / am;

					shp[found].n++;
					shp[found].sy += c;
					shp[found].syy += c * c;
					shp[found].sr += sem;
					shp[found].srr += sem * sem;
					shp[found].syr += c * sem;
					shp[found].sm += smm2;
					shp[found].smm += smm2 * smm2;
					shp[found].sym += c * smm2;
				}
			}
		}
	}

	for (g = 0; g < nout; g++) {
		double n = cnt[g] ? (double)cnt[g] : 1.0, varc, cov;

		if (shp[g].n > 1) {
			double m = (double)shp[g].n;
			double vy = shp[g].syy / m - (shp[g].sy / m) * (shp[g].sy / m);
			double vr = shp[g].srr / m - (shp[g].sr / m) * (shp[g].sr / m);
			double vm = shp[g].smm / m - (shp[g].sm / m) * (shp[g].sm / m);
			double cr = shp[g].syr / m - (shp[g].sy / m) * (shp[g].sr / m);
			double cm = shp[g].sym / m - (shp[g].sy / m) * (shp[g].sm / m);

			out[g].r2_sin_corr = (vy > 0.0 && vr > 0.0) ? cr * cr / (vr * vy) : 0.0;
			out[g].r2_map_corr = (vy > 0.0 && vm > 0.0) ? cm * cm / (vm * vy) : 0.0;
		}

		out[g].nsample = cnt[g];
		out[g].rms_ref = sqrt(out[g].rms_ref / n);
		out[g].rms_test = sqrt(out[g].rms_test / n);
		varc = sr[g] / n;
		cov = sc[g] / n;
		out[g].rms_corr = sqrt(varc);
		out[g].cov = cov;
		out[g].alpha = varc > 0.0 ? cov / varc : 0.0;
		/*
		 * Residual variance at the optimal scale: var(ref) - cov^2/var(c). It cannot
		 * exceed var(ref) -- the worst the optimal scale can do is 0, i.e. not applying
		 * the correction at all -- so this is the floor the component could reach if it
		 * were damped rather than applied raw.
		 */
		{
			double vr = out[g].rms_ref * out[g].rms_ref;
			double best = varc > 0.0 ? vr - cov * cov / varc : vr;

			out[g].rms_best = sqrt(best > 0.0 ? best : 0.0);
		}
	}
	free(sr);
	free(sc);
	free(src);
	free(cnt);
	free(shp);
	free(ashp_r);
	free(ashp_m);
	free(ashp_n);
	return nout;
}

/*
 * The two elevation shapes a double-differenced residual can carry, in DD form.
 *
 * A vertical phase-centre or height error projects on the line of sight as sin(el); a
 * troposphere error as the mapping function m(el) ~ 1/sin(el). These are *double* differences,
 * so what appears in the data is the between-satellite difference of each, and the reference
 * satellite's elevation varies within an arc just as the other's does. Regressing on
 * f(el_sat) alone leaves that unmodelled -- the defect that made D36's first attempt
 * uninterpretable.
 *
 * Returns 0 if either elevation is unusable.
 */
static int shape_terms(const struct dd_sample *s, double *sin_dd, double *map_dd) {
	double es, er, ss, sr;

	if (s->elev <= 0.0 || s->elev_ref <= 0.0)
		return 0;
	es = s->elev * D2R;			/* stored in DEGREES */
	er = s->elev_ref * D2R;
	ss = sin(es);
	sr = sin(er);
	if (ss < 0.05)
		ss = 0.05;
	if (sr < 0.05)
		sr = 0.05;
	*sin_dd = ss - sr;
	*map_dd = 1.0 / ss - 1.0 / sr;
	return 1;
}

int vbs_ddres_shape(const struct vbs_ddres *this, struct vbs_dd_shape *out, int maxout) {
	/* Per group: sums for two single-regressor fits, both de-meaned within the arc. */
	struct acc {
		double sr, sm, srr, smm, srm, sy, syy, syr, sym;
		long n;
	};
	struct acc *a;
	int i, g, nout = 0, pass;

	if (this == NULL || this->nsample <= 0)
		return 0;
	if ((a = (struct acc *)calloc(maxout, sizeof *a)) == NULL)
		return -1;
	memset(out, 0, maxout * sizeof *out);

	/*
	 * Two passes: the first accumulates each arc's mean elevation terms so the regressors
	 * can be de-meaned the same way the residual already is, the second accumulates the
	 * cross products. Doing it in one pass would regress against raw sin(el), whose
	 * between-arc variation has nothing to do with the within-arc shape being tested.
	 */
	for (pass = 0; pass < 2; pass++) {
		double *arc_sin = NULL, *arc_map = NULL;
		long *arc_n = NULL;

		arc_sin = (double *)calloc(this->narc, sizeof *arc_sin);
		arc_map = (double *)calloc(this->narc, sizeof *arc_map);
		arc_n = (long *)calloc(this->narc, sizeof *arc_n);
		if (arc_sin == NULL || arc_map == NULL || arc_n == NULL) {
			free(arc_sin);
			free(arc_map);
			free(arc_n);
			free(a);
			return -1;
		}
		for (i = 0; i < this->nsample; i++) {
			double se, sr;

			if (!shape_terms(&this->sample[i], &se, &sr))
				continue;
			arc_sin[this->sample[i].arc] += se;
			arc_map[this->sample[i].arc] += sr;
			arc_n[this->sample[i].arc]++;
		}
		if (pass == 0) {
			/* Stash the arc means by reusing the accumulator pass structure. */
			for (i = 0; i < this->nsample; i++) {
				const struct dd_arc *arc = &this->arc[this->sample[i].arc];
				double se, sr, sv, mv, y;
				int found = -1;
				long an = arc_n[this->sample[i].arc];

				if (arc->n < 2 || an < 2)
					continue;
				if (!shape_terms(&this->sample[i], &se, &sr))
					continue;
				sv = se - arc_sin[this->sample[i].arc] / an;
				mv = sr - arc_map[this->sample[i].arc] / an;
				y = this->sample[i].v - arc->sum / arc->n;

				for (g = 0; g < nout; g++)
					if (out[g].sys == arc->sys && out[g].comb == arc->comb) {
						found = g;
						break;
					}
				if (found < 0) {
					if (nout >= maxout)
						continue;
					found = nout++;
					out[found].sys = arc->sys;
					out[found].comb = arc->comb;
				}
				a[found].sr += sv;
				a[found].sm += mv;
				a[found].srr += sv * sv;
				a[found].smm += mv * mv;
				a[found].srm += sv * mv;
				a[found].sy += y;
				a[found].syy += y * y;
				a[found].syr += y * sv;
				a[found].sym += y * mv;
				a[found].n++;
			}
		}
		free(arc_sin);
		free(arc_map);
		free(arc_n);
		if (pass == 0)
			break;			/* one accumulation pass is enough */
	}

	for (g = 0; g < nout; g++) {
		double n = a[g].n ? (double)a[g].n : 1.0;
		double vy = a[g].syy / n, vr = a[g].srr / n, vm = a[g].smm / n;
		double cyr = a[g].syr / n, cym = a[g].sym / n, crm = a[g].srm / n;

		out[g].nsample = (int)a[g].n;
		out[g].coef_sin = vr > 0.0 ? cyr / vr : 0.0;
		out[g].coef_map = vm > 0.0 ? cym / vm : 0.0;
		out[g].r2_sin = (vy > 0.0 && vr > 0.0) ? (cyr * cyr) / (vr * vy) : 0.0;
		out[g].r2_map = (vy > 0.0 && vm > 0.0) ? (cym * cym) / (vm * vy) : 0.0;
		out[g].regressor_corr = (vr > 0.0 && vm > 0.0) ?
		    crm / sqrt(vr * vm) : 0.0;
	}
	free(a);
	return nout;
}

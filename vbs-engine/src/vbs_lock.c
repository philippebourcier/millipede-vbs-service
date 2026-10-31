#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_lock.h"

/*
 * Lock-time tracking for the VBS output -- spec §5.3. See the header.
 */

#define NSLOT	(NFREQ + NEXOBS)

struct slot_state {
	int valid;
	int seen;			/* present in this epoch's output set */
	gtime_t last;
	double lock_s;			/* seconds of continuous phase at the VBS */
	double last_corr_dI, last_corr_dT;
	int have_corr;
};

struct vbs_lock {
	struct vbs_lock_opt opt;
	struct slot_state sat[MAXSAT][NSLOT];
	struct vbs_lock_stats stats;
};

void vbs_lock_default_opt(struct vbs_lock_opt *opt) {
	memset(opt, 0, sizeof *opt);
	/*
	 * 0.05 m. Above it sits the smallest step a re-fix can put into the correction -- one
	 * cycle of the widest lane in use is 0.86 m, and even a narrow-lane cycle is 0.19 m --
	 * and below it sits the real epoch-to-epoch atmospheric change, which is far smaller
	 * than the 1.4-12.5 cm `03-error-budget.md` measures across a whole *baseline*.
	 *
	 * Erring low is the safe direction: a spurious reset costs the rover a re-convergence,
	 * while a missed one leaves it holding a wrong fix with full confidence.
	 */
	opt->jump_threshold_m = 0.05;
	opt->max_gap_s = 10.0;
}

struct vbs_lock *vbs_lock_new(const struct vbs_lock_opt *opt) {
	struct vbs_lock *this = (struct vbs_lock *)calloc(1, sizeof *this);

	if (this == NULL)
		return NULL;
	if (opt != NULL)
		this->opt = *opt;
	else
		vbs_lock_default_opt(&this->opt);
	return this;
}

void vbs_lock_free(struct vbs_lock *this) {
	free(this);
}

static const struct vbs_corr *find_corr(const struct vbs_corr *corr, int ncorr, int sat) {
	int i;

	for (i = 0; i < ncorr; i++)
		if (corr[i].sat == sat)
			return &corr[i];
	return NULL;
}

static const obsd_t *find_obs(const obsd_t *obs, int n, int sat) {
	int i;

	for (i = 0; i < n; i++)
		if (obs[i].sat == sat)
			return &obs[i];
	return NULL;
}

static int in_list(const int *list, int n, int sat) {
	int i;

	for (i = 0; i < n; i++)
		if (list[i] == sat)
			return 1;
	return 0;
}

static void count_reason(struct vbs_lock_stats *st, unsigned reason) {
	int b;

	for (b = 0; b < 5; b++)
		if (reason & (1u << b))
			st->resets[b]++;
}

int vbs_lock_update(struct vbs_lock *this, gtime_t time,
	const obsd_t *master, int nmaster, obsd_t *synth, int nsynth,
	const struct vbs_corr *corr, int ncorr,
	const int *master_slip, int nmaster_slip,
	const int *refixed, int nrefixed) {
	int i, j, nreset = 0;

	/*
	 * An epoch with *no* output is not a no-op: it is precisely when satellites leave the
	 * output set, and the sweep at the end of this function is what records that. Returning
	 * early on synth == NULL skipped it, so a satellite withheld for a few epochs came back
	 * with its state untouched and was reported to the rover as continuously locked across
	 * the gap. The loop below is bounded by nsynth, so a NULL/0 output is safe to fall
	 * through.
	 */
	if (this == NULL)
		return 0;
	this->stats.epochs++;

	/*
	 * Clear the presence marks first. A satellite dropped from the output -- §4.6 forbids
	 * synthesizing one whose network correction is unusable -- is simply *absent* from
	 * `synth`, so the per-slot loop below never visits it and its state would otherwise
	 * survive untouched. It would then come back looking continuous, and because
	 * MRTKLIB's locktime() keeps its own lltime per satellite-signal and only restarts on
	 * LLI, the rover would be told the phase had been locked throughout a gap in which the
	 * network ambiguity may well have changed. The elapsed-time test alone does not cover
	 * it: a one- or two-epoch dropout at 1 Hz is shorter than any sane gap threshold.
	 */
	for (i = 0; i < MAXSAT; i++)
		for (j = 0; j < NSLOT; j++)
			this->sat[i][j].seen = 0;

	for (i = 0; i < nsynth; i++) {
		int sat = synth[i].sat;
		const obsd_t *m;
		const struct vbs_corr *c;
		unsigned sat_reason = 0;

		if (sat <= 0 || sat > MAXSAT)
			continue;
		m = find_obs(master, nmaster, sat);
		c = find_corr(corr, ncorr, sat);

		/*
		 * §5.3 case 1: the master slipped. Taken from what vbs_preproc detected, not
		 * from the master's LLI -- D21/D22 measured that flag reporting zero slips on the
		 * storm day while the combinations found 1 991.
		 */
		if (master_slip != NULL && in_list(master_slip, nmaster_slip, sat))
			sat_reason |= VBS_LOCK_MASTER_SLIP;

		/* §5.3 case 2: the network ambiguity was re-fixed or lost. Per satellite. */
		if (refixed != NULL && in_list(refixed, nrefixed, sat))
			sat_reason |= VBS_LOCK_AMB_REFIX;

		/*
		 * §5.3 case 3: the interpolated correction jumped. Per satellite, because the
		 * correction is per satellite -- an auxiliary dropping out of the cluster moves
		 * the plane fit for every satellite at once, but by different amounts.
		 */
		if (c != NULL) {
			struct slot_state *any = &this->sat[sat - 1][0];

			if (any->have_corr &&
			    (fabs(c->dI - any->last_corr_dI) > this->opt.jump_threshold_m ||
			     fabs(c->dT - any->last_corr_dT) > this->opt.jump_threshold_m))
				sat_reason |= VBS_LOCK_CORR_JUMP;
		}

		for (j = 0; j < NSLOT; j++) {
			struct slot_state *s = &this->sat[sat - 1][j];
			unsigned reason = sat_reason;
			double dt;

			if (synth[i].code[j] == CODE_NONE || synth[i].L[j] == 0.0) {
				s->valid = 0;
				continue;
			}
			this->stats.tracked++;

			/* The master's own flag, as an additional trigger only (see the header). */
			if (m != NULL && (m->LLI[j] & 1))
				reason |= VBS_LOCK_MASTER_SLIP;

			/*
			 * §5.3 case 4 is folded in here: a satellite that was excluded by FDE
			 * simply stops being synthesized, so it comes back as an untracked slot
			 * or after a gap. Either way continuity cannot be assumed.
			 */
			dt = s->valid ? fabs(timediff(time, s->last)) : 0.0;
			if (!s->valid)
				reason |= VBS_LOCK_NEW;
			else if (dt > this->opt.max_gap_s)
				reason |= VBS_LOCK_FDE_RETURN;

			if (reason != 0) {
				s->lock_s = 0.0;
				nreset++;
				count_reason(&this->stats, reason);
			} else {
				s->lock_s += dt;
			}
			s->valid = 1;
			s->seen = 1;
			s->last = time;

			/*
			 * Stamp the synthesized observation. Setting LLI bit 0 is what makes
			 * MRTKLIB's encoder restart the MSM Lock Time Indicator, which is the only
			 * way to tell the rover that continuity is broken (§5.3).
			 */
			if (reason != 0)
				synth[i].LLI[j] |= 1;
			else
				synth[i].LLI[j] &= (uint8_t)~1;
		}

		if (c != NULL) {
			struct slot_state *any = &this->sat[sat - 1][0];

			any->last_corr_dI = c->dI;
			any->last_corr_dT = c->dT;
			any->have_corr = 1;
		}
	}

	/*
	 * Anything tracked but absent from this epoch's output has left the output set. Mark it
	 * unusable so that its return raises VBS_LOCK_NEW on the first epoch back, whatever the
	 * gap length.
	 */
	for (i = 0; i < MAXSAT; i++)
		for (j = 0; j < NSLOT; j++)
			if (this->sat[i][j].valid && !this->sat[i][j].seen)
				this->sat[i][j].valid = 0;

	return nreset;
}

void vbs_lock_get_stats(const struct vbs_lock *this, struct vbs_lock_stats *out) {
	*out = this->stats;
}

double vbs_lock_seconds(const struct vbs_lock *this, int sat, int slot) {
	if (sat <= 0 || sat > MAXSAT || slot < 0 || slot >= NSLOT)
		return -1.0;
	if (!this->sat[sat - 1][slot].valid)
		return -1.0;
	return this->sat[sat - 1][slot].lock_s;
}

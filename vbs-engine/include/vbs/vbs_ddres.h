#ifndef __VBS_DDRES_H__
#define __VBS_DDRES_H__

#include <stdio.h>

#include "rtklib.h"

#include "vbs/vbs_signal.h"

/*
 * Double-difference residuals between two observation sets taken at the same position.
 *
 * This is the measurement spec S8 step 1 and data/README.md S5 step 1 ask for: exclude a
 * base, synthesize a VBS *at that base's own coordinates*, and double-difference the
 * synthesized observations against the base's real ones. Because both sets belong to the same
 * point, the geometric range cancels exactly, so what is left is the part the network
 * corrections are supposed to supply -- plus an integer ambiguity per satellite pair, plus
 * noise.
 *
 * With corrections forced to zero this measures the *uncorrected* error budget: how far the
 * master's observations are from the target's once geometry is removed, as a function of
 * baseline length. That is the quantity the network EKF has to explain, and it is worth
 * knowing before designing the filter that has to remove it. Run again later with real
 * corrections, the same code reports how much of it the engine actually removed.
 *
 * Residual convention: each (satellite pair, signal) arc carries an unknown constant integer
 * ambiguity, so the arc mean is removed before any statistic is taken. What is reported is
 * therefore the *time-varying* part of the double difference -- which is what matters, since
 * a rover's own ambiguity resolution absorbs the constant part. Code observables have no
 * ambiguity, but they are reported the same way for comparability; their scatter is dominated
 * by multipath and noise, not by the atmosphere.
 */

/* Combination reported alongside the per-band residuals (spec S8 step 1). */
enum vbs_dd_comb {
	VBS_DD_BAND,		/* one carrier band, as observed */
	VBS_DD_GEOMFREE,	/* L_a - L_b: geometry and clocks cancel, ionosphere remains */
	VBS_DD_IONOFREE		/* first-order ionosphere cancels, non-dispersive remains */
};

/* Statistics for one (system, band-or-combination) group. */
struct vbs_dd_stat {
	int sys;
	int band;		/* VBS_BAND_* for VBS_DD_BAND, else 0 */
	enum vbs_dd_comb comb;
	int nsample;
	int narc;
	double rms;		/* metres, after per-arc mean removal */
	double p95;		/* metres, absolute value */
	double max;		/* metres, absolute value */
};

#define VBS_DD_MAX_STATS	32

struct vbs_dd_result {
	int nstat;
	struct vbs_dd_stat stat[VBS_DD_MAX_STATS];
	int nepoch;		/* epochs with a usable reference satellite */
	int nskipped;		/* epochs dropped for want of common satellites */
};

/*
 * Accumulator. Opaque: it holds every DD sample until the arcs can be de-meaned at the end,
 * which is not something the caller should have to manage.
 */
struct vbs_ddres;

/*
 * elev_mask_deg applies at the *target* position; the reference satellite is the highest
 * one above the mask that both sets observe.
 */
struct vbs_ddres *vbs_ddres_new(double elev_mask_deg);
void vbs_ddres_free(struct vbs_ddres *this);

/*
 * Add one epoch. a/b are the two sets of observations for the same epoch and the same
 * position (a = synthesized VBS, b = the real base). rs holds satellite positions for set a,
 * 6 doubles per record, as satposs() produces them. rr is the shared ECEF position.
 *
 * Returns 1 if the epoch contributed samples, 0 if it was skipped.
 */
int vbs_ddres_add_epoch(struct vbs_ddres *this,
	const obsd_t *a, int na, const obsd_t *b, int nb,
	const double *rs, const double *rr, const nav_t *nav);

/*
 * Finish: de-mean every arc and reduce to per-group statistics.
 */
int vbs_ddres_finish(struct vbs_ddres *this, struct vbs_dd_result *out);

/*
 * Diagnostic: print every residual whose magnitude exceeds `threshold` metres, with the
 * satellite pair, band and arc it came from. A double difference at a fixed point has no
 * business being metres large, so when one is, the question is always "which pair, and was
 * the arc really continuous" -- this answers it directly instead of by inference.
 * Call after vbs_ddres_finish(). Prints at most `max_lines`.
 */
void vbs_ddres_dump_outliers(struct vbs_ddres *this, double threshold, int max_lines,
	FILE *fp);

/*
 * Paired comparison of two accumulators that saw identical epochs and satellites -- the
 * corrected run against a reference run with the corrections zeroed downstream of every
 * decision that selects samples (see tools/vbs_pipeline --no-corrections).
 *
 * Comparing two RMS numbers only says which is smaller. Because the two runs differ by exactly
 * the applied correction, `c = reference - test` sample by sample, and the RMS decomposes:
 *
 *     RMS_test² = RMS_ref² - 2·cov(ref, c) + var(c)
 *
 * which separates the two things a correction does: `2·cov` is the error it removes, `var(c)`
 * is the noise it adds. A correction can be perfectly correct on average and still make the
 * residual worse, if it is noisier than the signal it is removing -- and the two RMS figures
 * alone cannot tell that case from a correction that is simply wrong.
 *
 * `alpha = cov(ref, c) / var(c)` is the scale the correction should be multiplied by to
 * minimise the residual. alpha ≈ 1 means the correction is the right size; alpha ≈ 0.3 means
 * it carries three times more noise than signal and should be damped or withheld. That number
 * is what D33 item 1 asks for when it says to apply a component "only where it measurably
 * reduces the error", and its variation with baseline length is the crossover distance.
 */
struct vbs_dd_pair {
	int sys;
	int band;
	enum vbs_dd_comb comb;
	int nsample;
	double rms_ref;			/* residual RMS with the correction zeroed */
	double rms_test;		/* residual RMS with it applied */
	double rms_corr;		/* RMS of the correction itself */
	double cov;			/* cov(reference residual, correction) */
	double alpha;			/* cov/var(c): the optimal scale for the correction */
	double rms_best;		/* residual RMS the correction would give at that scale */
	/*
	 * The correction's *own* elevation shape, as a check on the shape regression rather
	 * than on the engine. The troposphere model term is ΔZHD·(m(el_sat) − m(el_ref)) by
	 * construction, so `r2_map_corr` must come out near 1 whenever that term dominates the
	 * correction. If it does not, the shape regression is wrong and its verdict on the
	 * residual means nothing -- which is the standing rule (D37) applied to the tool that
	 * most needs it.
	 */
	double r2_sin_corr, r2_map_corr;
};

/*
 * Both accumulators must already have been through vbs_ddres_finish(). Samples are paired by
 * index, which is valid exactly when the two runs processed the same epochs and satellites --
 * the function refuses (returns -1) if the sample counts differ, rather than silently
 * comparing misaligned series.
 */
int vbs_ddres_pair(const struct vbs_ddres *ref, const struct vbs_ddres *test,
	struct vbs_dd_pair *out, int maxout);

/*
 * Which elevation shape does a residual actually have? -- D35 item 1.
 *
 * A vertical phase-centre or height error projects onto the line of sight as **sin(el)**. A
 * troposphere error maps as **m(el) ≈ 1/sin(el)**. Within a single arc those two are strongly
 * anti-correlated -- as a satellite rises, sin(el) grows while m(el) shrinks -- which is the
 * classic height/ZTD correlation that makes the two hard to separate anywhere in GNSS, and
 * which D35 proposes as the reason a correctly signed troposphere correction can show ρ ≈ -0.75
 * against a residual that is really antenna error.
 *
 * Because the two regressors are collinear, fitting them jointly would produce unstable
 * coefficients that say more about the conditioning than the physics. Each is therefore fitted
 * *alone* and reported by the variance it explains, with their mutual correlation given so the
 * reader can see how much the comparison is worth: at |corr| near 1 the two are barely
 * distinguishable and neither r² should be over-read.
 */
struct vbs_dd_shape {
	int sys;
	enum vbs_dd_comb comb;
	int nsample;
	double r2_sin;			/* variance explained by sin(el) alone: antenna-like */
	double r2_map;			/* by m(el) alone: troposphere-like */
	double coef_sin;		/* metres per unit sin(el) */
	double coef_map;		/* metres per unit m(el) */
	double regressor_corr;		/* corr(sin(el), m(el)) within arcs */
};

/* Call after vbs_ddres_finish(). Returns the number of groups written. */
int vbs_ddres_shape(const struct vbs_ddres *this, struct vbs_dd_shape *out, int maxout);

#endif /* __VBS_DDRES_H__ */

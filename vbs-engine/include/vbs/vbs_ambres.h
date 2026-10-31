#ifndef __VBS_AMBRES_H__
#define __VBS_AMBRES_H__

#include "vbs/vbs_netekf.h"

/*
 * Network ambiguity resolution -- spec §4.3.
 *
 * D24 measured why this is load-bearing rather than an optimisation: with float ambiguities
 * the per-satellite slant ionosphere converges only to 30-370 mm, because code alone pins the
 * *level* of ΔI while phase pins only how it changes. That is an order of magnitude above the
 * 1.4-12.5 cm of ionospheric interpolation error `03-error-budget.md` measured, so float
 * corrections cannot meet the budget they exist to remove. Fixing the ambiguities is what
 * turns ΔI into a phase-determined, millimetre-level quantity.
 *
 * **The cascade.** §4.3 is explicit that "direct L1/L2 fixing is not realistic" at 50-150 km,
 * and the reason is wavelength: an L1 ambiguity is 19 cm, so fixing it needs the float
 * solution good to a few centimetres, which at 150 km it is not. The way through is to fix a
 * long-wavelength combination first, feed that back, and let it tighten the rest:
 *
 *   extra-wide lane   E5a-E5b, B2a-B2b     λ ≈ 9.8 m    trivially fixable, needs 3 bands
 *   wide lane         L1-L2, E1-E5b        λ ≈ 0.75-0.86 m
 *   narrow lane       the individual bands  λ ≈ 0.19-0.25 m
 *
 * Each step is fixed, then applied to the filter as a constraint (`vbs_netekf_constrain_dd`),
 * so the next step starts from a float solution the previous one has already sharpened. On a
 * dual-frequency baseline the extra-wide step simply finds nothing to do and the cascade
 * starts at the wide lane -- which is §4.3's "otherwise" branch, reached by the same code
 * rather than by a separate path.
 *
 * Measured on this network (`tools/vbs_rtcm_inv` on a live Unicore UM980): Galileo carries
 * four bands on 100% of satellite-epochs and BeiDou four on 90%, so the triple-frequency
 * branch is the normal case between Unicore bases, not an exception. It is the F9P baselines
 * that are dual-frequency.
 *
 * **Partial resolution.** §4.3 asks for it, and it matters most exactly where the cascade is
 * most needed: one satellite with a bad float ambiguity should cost that satellite, not the
 * epoch. Validation is the ratio test on the best and second-best integer candidates; when it
 * fails, the worst-conditioned ambiguity is dropped and the subset is retried.
 */

#define VBS_AMB_MAXDD		(VBS_NET_MAXSAT - 1)

/* Which cascade step fixed an ambiguity -- reported, because "how" is the diagnostic. */
#define VBS_AMB_STEP_EWL	0
#define VBS_AMB_STEP_WL		1
#define VBS_AMB_STEP_NL		2
#define VBS_AMB_NSTEP		3

struct vbs_ambres_opt {
	double min_ratio;		/* ratio-test threshold */
	double min_elev_deg;		/* exclude low satellites from resolution entirely */
	double max_float_sigma;		/* skip a step whose float solution is hopeless, cycles */
	/*
	 * Multiplier on the float covariance handed to the search. The filter's sigmas describe
	 * independent observations; at 1 Hz consecutive residuals are not, so the formal
	 * ambiguity covariance is optimistic by roughly (correlation time / interval), and the
	 * ratio test admits fixes at 1 Hz it refuses at 30 s (D48). Inflating Q at the search
	 * alone -- not in the filter -- keeps the float converging at full rate while the
	 * search sees the information the data really holds. 1.0 = off (30 s data).
	 */
	double q_inflate;
	double fix_var;			/* variance of the constraint applied to a fix, cycles² */
	int min_fixed;			/* fewer than this in a subset is not worth applying */
	int partial;			/* allow dropping ambiguities to rescue a subset */
	/*
	 * Highest cascade step the resolver is allowed to attempt: VBS_AMB_STEP_EWL, _WL or _NL
	 * (the default). A diagnostic first -- on the Unicore networks the wide lanes close and
	 * the narrow lane misdoses by whole cycles in bulk, and stopping the cascade short is how
	 * that is attributed rather than argued -- and a policy knob if stopping short turns out
	 * to serve a rover better than a wrong narrow lane does.
	 */
	int max_step;
	/*
	 * D67, the ionosphere-free narrow lane. With the wide-lane integer known, the remaining
	 * integer is resolved on the ionosphere-free combination instead of on the raw carrier,
	 * which removes the ionosphere from the model entirely and with it the baseline-length
	 * limit. Its float is averaged over the arc because one epoch of an IF observable is
	 * about a third of a narrow-lane cycle.
	 */
	int if_narrow_lane;		/* use it in place of the raw-carrier NL steps */
	int if_min_epochs;		/* epochs of averaging before a fix is attempted */
	double if_max_sigma;		/* accept only below this, cycles of the narrow lane */
	double if_sys_m;		/* systematic floor that averaging cannot beat down, metres */
	double if_max_frac;		/* refuse a mean further than this from its integer, cycles */
	/*
	 * Order the cascade by each combination's actual wavelength on each constellation,
	 * rather than by its row in the table (D74). The default, 1, is the correct behaviour;
	 * 0 restores the pre-D74 row order so the change can be measured from one binary.
	 */
	int cascade_by_lambda;
	/*
	 * Use the band-6 combinations -- Galileo E6, BeiDou B3I (D76). The default, 1, uses
	 * whatever the stations actually stream; 0 drops them, to measure what they are worth.
	 */
	int use_band6;
	/*
	 * Constellations for which the L1/L2 wide lane (CASCADE row "WL 1-2") is not attempted.
	 *
	 * That row is the only one built from the L1/L2 Melbourne-Wubbena combination, so it is
	 * the only one a degraded L2 *code* reaches. Measured on day 210 (D129): at the moment the
	 * cascade rounds, a UM980 cluster's WL 1-2 floats sit 0.055-0.111 cycles from integer with
	 * 10-27 % of decisions beyond 0.20 of the 0.25-cycle tolerance, while every combination on
	 * the same baselines that avoids band 2 -- WL 1-5, EWL 5-6, WL 1-7 -- sits at 0.000-0.004.
	 * A Mosaic cluster's WL 1-2 sits at 0.018-0.034. BeiDou's band 2 is B1I, not GPS L2, and
	 * measures clean, which is why this is a system mask and not a row toggle.
	 *
	 * 0 disables nothing (the default).
	 */
	int skip_wl12_sys;
	/*
	 * Constellations for which the L1/L2 wide lane is DERIVED from two already-fixed
	 * combinations instead of being resolved from its own float.
	 *
	 * N(1-2) = N1 - N2 = (N1 - N5) - (N2 - N5) = N(1-5) - N(2-5), exactly, as integers.
	 * D129 measured WL 1-2's float sitting 0.055-0.111 cycles from integer on UM980 clusters
	 * at the moment the cascade rounds -- 10-27 % of decisions within 0.05 cycles of the 0.25
	 * tolerance -- because it is the only cascade row built from the L1/L2 Melbourne-Wubbena
	 * combination, and a UM980 reports L2C under the L2W label so that code carries two
	 * different satellite biases. WL 1-5 and EWL 2-5, on the same baselines and the same
	 * epochs, measure 0.000. The integer is therefore already known once those two are fixed:
	 * take the difference rather than round a biased float.
	 *
	 * Needs band 5, so a satellite that does not transmit L5 simply gets no WL 1-2 from this
	 * path. Implies skipping the ordinary WL 1-2 row for those constellations.
	 */
	int wl12_derive_sys;
	/*
	 * Constellations for which the cascade holds ONE reference satellite across all of its
	 * steps, instead of letting each combination pick its own highest-elevation candidate.
	 *
	 * `vbs_ambres_ion_from_fixes` can only use satellites whose wide lane and narrow lane are
	 * fixed against a single common reference -- mixing references gives each satellite a
	 * different unknown ionosphere level, which D66 measured taking M30 from 98.4 % fixed to
	 * 13.2 % with a 12.6 cm vertical bias. With per-combination references, 64 % of corrections
	 * had no usable pair and the feature replaced nothing at all (D131). Holding one reference
	 * makes the pair exist by construction.
	 *
	 * A combination whose bands the held reference lacks falls back to its own choice, so this
	 * can only add pairs, never remove them. 0 disables it (the default).
	 */
	int one_ref_sys;
};

void vbs_ambres_default_opt(struct vbs_ambres_opt *opt);

struct vbs_ambres_stats {
	long attempts[VBS_AMB_NSTEP];
	long fixed[VBS_AMB_NSTEP];		/* ambiguities constrained */
	long ratio_fail[VBS_AMB_NSTEP];
	long partial_used[VBS_AMB_NSTEP];
	double last_ratio[VBS_AMB_NSTEP];
	int last_ndd[VBS_AMB_NSTEP];
	int last_nfix[VBS_AMB_NSTEP];
};

/*
 * One ambiguity this pass fixed, for the network loop-closure check below.
 */
struct vbs_amb_fix {
	int sat, refsat, sys;
	int comb_id;			/* index into the cascade table: identifies the combination */
	int step;			/* VBS_AMB_STEP_* */
	double integer;			/* cycles of that combination */
	double value;			/* the float it was fixed from, cycles */
	double sigma;			/* its float standard deviation, cycles */
};

/*
 * Run the cascade once, for the filter's current epoch, over every constellation in
 * `sys_mask`. Fixes are applied to the filter as constraints before returning.
 *
 * `fixed`/`maxfixed` optionally receive what was fixed (pass NULL to ignore); `*nfixed` gets
 * the count. Returns the number of ambiguities constrained across all steps.
 */
/*
 * Recompute the dispersive corrections from the resolved integers instead of from the filter's
 * ionosphere states (D66). For every satellite whose double-differenced ambiguities are fixed
 * against its reference, the geometry-free combination gives the DD ionosphere exactly; the
 * reference's own level is taken from the filter, which D55 showed is constellation-common and
 * cancels at the rover. Satellites without the needed integers are left untouched.
 *
 * Returns how many corrections were replaced. `sd_I` may be NULL; when given it receives the
 * propagated phase-noise sigma for those entries.
 */
int vbs_ambres_ion_from_fixes(const struct vbs_netekf *ekf, const struct vbs_amb_fix *fix,
	int nfix, struct vbs_corr *corr, int ncorr, double *sd_I);

/*
 * Why ion_from_fixes replaced nothing. It is inert in the shipped configuration -- 0 of
 * 1 750 839 corrections -- and "implemented but not reached" is the defect class this engine
 * keeps producing (D90, six instances). Counts are cumulative over the run: corrections
 * skipped for want of a usable reference satellite, skipped for want of a wide-lane /
 * narrow-lane pair fixed against that reference, and actually replaced.
 */
void vbs_ambres_ion_diag(long *no_ref, long *no_pair, long *replaced);

int vbs_ambres_resolve(struct vbs_netekf *ekf, const struct vbs_ambres_opt *opt,
	int sys_mask, struct vbs_ambres_stats *stats,
	struct vbs_amb_fix *fixed, int maxfixed, int *nfixed);

/*
 * --- Network loop closure (D33 item 2) ---------------------------------------------------
 *
 * A ratio threshold is a guess that has to be calibrated per case, and D32 caught it letting
 * wrong fixes through at 290 km while a tighter value rescued one constellation on one case --
 * exactly the evidence that invites over-fitting a number. Loop closure is an *identity*
 * instead: it either holds or it does not.
 *
 * For a triangle of master M and auxiliaries A and B, the double-difference ambiguities of the
 * three baselines are not independent. With the same satellite s, the same reference satellite
 * r, and the same combination:
 *
 *     N_AB = N_MB − N_MA
 *
 * because each is a difference of the same four per-station ambiguities. A non-zero misclosure
 * proves at least one of the three fixes is wrong. It does *not* say which -- a single triangle
 * cannot isolate it -- so all three are unfixed, which is the conservative reading and the only
 * one the evidence supports.
 *
 * The three baselines must share a reference satellite for the identity to apply. Each filter
 * picks its own (highest elevation), and they normally agree because all three see nearly the
 * same sky; pairs where they disagree are skipped rather than compared, since a mismatched
 * reference makes the identity simply not hold.
 */
/*
 * Misclosure bins. The distinction that matters is **integer versus not**: a genuine wrong fix
 * moves an ambiguity by a whole number of cycles, so the misclosure lands on ±1, ±2, ... A
 * misclosure that is not near an integer cannot be a wrong fix at all -- it is a comparison
 * bug, and reporting a "wrong-fix rate" that includes those would be measuring the code rather
 * than the network.
 */
#define VBS_LOOP_BIN_CLOSED	0	/* |mis| < tolerance */
#define VBS_LOOP_BIN_1CYC	1
#define VBS_LOOP_BIN_2CYC	2
#define VBS_LOOP_BIN_LARGE	3	/* > 2.5 cycles but still near an integer */
#define VBS_LOOP_BIN_NONINT	4	/* not near any integer: not a wrong fix */
#define VBS_LOOP_NBIN		5

struct vbs_amb_loop_stats {
	long checked;			/* triangles × ambiguities where all three were fixed */
	long misclosed;			/* of those, failing the identity */
	long skipped_ref;		/* skipped: the three baselines used different references */
	double worst_misclosure;	/* cycles */
	long hist[VBS_AMB_NSTEP][VBS_LOOP_NBIN];
	/*
	 * Cascade propagation. A wide-lane fix wrong by one cycle (λ ≈ 0.86 m) constrains the
	 * narrow-lane search to a wrong N1 − N2, and the narrow lane then lands 4–5 L1 cycles
	 * off to satisfy it. That predicts exactly the shape seen in D44: a hole at 2 cycles
	 * and a tail at 3–8. These count narrow-lane misclosures by whether the same satellite's
	 * wide lane *also* misclosed in the same triangle -- if the tail lives almost entirely
	 * in `with_wl`, it is one wrong fix counted twice more, not a separate population.
	 */
	long nl_tail_with_wl, nl_tail_without_wl;	/* narrow-lane |mis| > 2.5 */
	long nl_one_with_wl, nl_one_without_wl;		/* narrow-lane |mis| ~ 1 */
	/*
	 * The same question one level up. The cascade fixes the extra-wide lane *first* and feeds
	 * it back before the wide lane runs, so a wrong extra-wide fix -- one cycle is 6-10 metres
	 * -- propagates downward into both the wide and the narrow lane. On the triple-frequency
	 * case the extra-wide lane misclosed 21 % (D45), which no ratio test lets through for a
	 * 10 m error; if the wide- and narrow-lane tails there live almost entirely under an
	 * extra-wide misclosure of the same satellite, the whole population is one defect at the
	 * top of the cascade, not three.
	 */
	long wl_tail_with_ewl, wl_tail_without_ewl;	/* wide-lane |mis| > 2.5 */
	long nl_tail_with_ewl, nl_tail_without_ewl;	/* narrow-lane |mis| > 2.5 */
};

/*
 * Check one triangle. Returns the number of misclosures found, and appends the (sat, refsat,
 * comb) triples that failed to `bad` so the caller can unfix them on all three baselines.
 */
int vbs_ambres_loop_check(const struct vbs_amb_fix *ma, int nma,
	const struct vbs_amb_fix *mb, int nmb,
	const struct vbs_amb_fix *ab, int nab,
	double tol_cycles, struct vbs_amb_loop_stats *stats,
	struct vbs_amb_fix *bad, int maxbad);

/*
 * One step of the cascade as it will actually be run, for tests and diagnostics.
 *
 * `lam` is the combination's wavelength on that constellation, in metres, and is what orders
 * the cascade -- the table of combinations is keyed on RINEX band digits, which mean a
 * different physical combination on each constellation (D74).
 */
struct vbs_amb_step_info {
	int sys;			/* SYS_* */
	int step;			/* VBS_AMB_STEP_*, derived from lam, not from table order */
	double lam;			/* wavelength magnitude, metres */
	const char *name;		/* the combination's name, e.g. "WL 1-2" */
};

/*
 * Report the cascade order for a set of constellations, longest wavelength first. Writes at
 * most `maxout` entries; returns how many were written.
 */
int vbs_ambres_cascade_order(int sys_mask, int by_lambda, int use_band6,
	struct vbs_amb_step_info *out, int maxout);

#endif /* __VBS_AMBRES_H__ */

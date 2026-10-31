#ifndef __VBS_NETEKF_H__
#define __VBS_NETEKF_H__

#include "rtklib.h"

#include "vbs/vbs_signal.h"
#include "vbs/vbs_synth.h"

/*
 * Network EKF for one master->auxiliary baseline -- spec §4.3, and §4.4's correction
 * extraction, which falls out of the state vector rather than needing a second step.
 *
 * The spec's state vector is "residual ZTD difference, per-satellite slant ionosphere
 * differences (STEC), and double-difference (DD) ambiguities per frequency". This filter
 * estimates the first two exactly as written and keeps the ambiguities in *single*-difference
 * form, forming DD combinations only when ambiguity resolution needs them.
 *
 * Why single difference. §4.4 requires that "corrections stay continuous when the DD
 * reference satellite changes". In a DD state vector that is a datum problem: every state is
 * expressed against one satellite, so when that satellite sets, the whole vector has to be
 * transformed, and any error in the transform is a step in the corrections -- which §5.3 would
 * then have to signal as a slip on every satellite at once. Estimating between stations only
 * (aux - master, per satellite) has no reference satellite in it at all, so the datum cannot
 * move and the corrections are continuous by construction. The price is one extra state per
 * constellation for the receiver clock difference, and single-difference ambiguities that are
 * not integers because they absorb the two receivers' phase biases. Neither costs anything
 * here: the corrections §4.4 asks for (ΔI and ΔT per satellite, relative to the master) are
 * themselves single-difference quantities, and integer structure is only ever needed in the
 * double differences formed at resolution time, where the receiver biases cancel.
 *
 * Observation model, per satellite s and band i, after removing the computed geometry and the
 * a priori troposphere (everything in metres, phase converted from cycles):
 *
 *   z_P = (P_aux - P_master) - Δρ - ΔT_model = Δclk + m_w·ΔZTD + μ_i·ΔI + ε
 *   z_L = (Φ_aux - Φ_master) - Δρ - ΔT_model = Δclk + m_w·ΔZTD - μ_i·ΔI + λ_i·ΔN_i + ε
 *
 * with μ_i = (f_ref/f_i)² and ΔI on the constellation's reference frequency (the same
 * convention `vbs_synth_ref_freq()` applies on the way back out, §4.6). The opposite
 * ionosphere sign on code and phase is what makes ΔI separable from ΔN; the μ_i spread across
 * bands is what makes it separable from the non-dispersive terms.
 *
 * Δρ is computed by `vbs_synth_drho()` -- the same transmission-time-corrected, Sagnac-aware
 * range difference the synthesis uses, so the filter and the synthesis cannot disagree about
 * geometry. The satellite clock cancels in the difference and is not modelled.
 *
 * What this filter does NOT estimate, deliberately: a per-satellite non-dispersive term. The
 * spec's state vector does not have one, so orbit error and any azimuthally asymmetric
 * troposphere are absorbed into the single-difference ambiguities instead of appearing in ΔT.
 * §4.7 argues that is acceptable (~2.5 mm per metre of orbit error at 50 km). Whether it
 * really is, is a measurement, not an assumption -- `tools/vbs_loo` is what settles it.
 */

#define VBS_NET_MAXSAT		48	/* satellites tracked simultaneously on one baseline */
/*
 * 5, not 3. Unicore bases carry four bands on every Galileo and BeiDou satellite (D26), and
 * with three slots band_index() silently dropped whichever band was seen fourth -- for Galileo
 * that was E5b, the one every F9P streams and the partner of E5a in the only true extra-wide
 * lane. Which band got dropped also depended on first-sight order, so two filters could hold
 * different triples for the same satellite. NX grows from 198 to 294 states; the covariance is
 * still under 700 kB per baseline.
 */
#define VBS_NET_MAXBAND		5	/* bands per satellite carried in the state */

/* Why an observation was rejected -- reported, never silently dropped (spec §4.3 FDE). */
#define VBS_NET_REJ_ELEV	0	/* below the elevation mask */
#define VBS_NET_REJ_NOSLOT	1	/* state vector full */
#define VBS_NET_REJ_CODE	2	/* code innovation outside the FDE gate */
#define VBS_NET_REJ_PHASE	3	/* phase innovation outside the FDE gate */
#define VBS_NET_REJ_GEOM	4	/* no usable satellite position */
#define VBS_NET_REJ_SNR		5	/* below the C/N0 mask (D51) */
#define VBS_NET_NREJ		6

struct vbs_netekf_opt {
	int sys_mask;			/* OR of SYS_*; D5 runs GPS+Galileo, §4.6 adds BeiDou */
	double elev_mask_deg;
	double humi;			/* relative humidity for the a priori Saastamoinen */

	double sig_code;		/* zenith code sigma, m */
	double sig_phase;		/* zenith phase sigma, m */

	double q_ztd;			/* ZTD difference random walk, m²/s */
	double q_ion;			/* slant ionosphere random walk, m²/s */
	/*
	 * D168: the same two figures for GPS satellites alone, 0 = as every constellation. A UM980
	 * cluster's GPS ionosphere tracks the truth far worse than its Galileo one (D167), and the
	 * code it is levelled on is the obvious per-constellation difference.
	 */
	double sig_code_gps, q_ion_gps;

	double var_ztd0;		/* initial ZTD difference variance, m² */
	/*
	 * D68: horizontal troposphere gradients of the station difference, north and east, with
	 * Chen-Herring mapping m_g(e) = 1/(sin e·tan e + 0.0032). A single zenith scalar cannot
	 * represent a tilted atmosphere; what it leaves behind is azimuth- and elevation-
	 * dependent, is worth several centimetres on low satellites, and biases every baseline
	 * of a cluster alike -- so loop closure never sees it while the corrections carry it.
	 */
	double var_grad0;		/* initial gradient variance, m² */
	double q_grad;			/* gradient random walk, m²/s */
	int est_gradients;		/* estimate them at all */
	double max_mapg;		/* cap on the gradient mapping function */
	double var_ion0;		/* initial slant ionosphere variance, m² */
	/*
	 * Variance of the per-epoch *seed* for the receiver clock difference, not of the clock
	 * itself. The state is re-seeded every epoch from the median of that epoch's code
	 * single differences, so this only has to cover the median's own error.
	 */
	double var_clk;
	double var_amb0;		/* initial ambiguity variance, cycles² */

	/*
	 * FDE gate: reject an observation whose innovation exceeds this many standard
	 * deviations of its own predicted innovation variance. §4.3 requires that one badly
	 * calibrated base cannot corrupt the state; this is the per-observation half of that,
	 * and it works because the filter updates observations one at a time, so a bad one is
	 * tested against a state it has not yet touched.
	 */
	double max_innov_sigma;
	/* Do not gate until the state has had this many epochs to converge from its priors. */
	int gate_after_epochs;
	/*
	 * After this many consecutive epochs with every phase of a satellite gated out, its
	 * ambiguities are reset as if a slip had been seen (D49): a fix that the phase no longer
	 * touches is not a fix, and the slip detector will not catch everything.
	 */
	int reject_streak_reset;
	/*
	 * C/N0 mask, dB-Hz; 0 = off. D51: T30's Unicore stations produce measurements on low
	 * satellites where phase AND code move together by ~100 m and come back -- invisible to
	 * geometry-free (equal on bands) and to Melbourne-Wübbena (equal on L and P). Their
	 * median C/N0 is 33 dB-Hz against 46 for the rest; a mask is the generic safety net, and
	 * its coverage cost is reported rather than assumed. Applied to both stations' signal.
	 */
	double min_snr_dbhz;
	/*
	 * D53: when a GPS record carries both semi-codeless L2 (2W/2P/2Y) and L2C (2S/2L/2X),
	 * take L2C. The filter otherwise takes the first usable slot of the band, which on the
	 * Unicore stations is L2W (listed first); the F9P network, which works, has no L2W at
	 * all. Diagnostic knob until measured.
	 */
	int prefer_l2c;
	int drop_band5;			/* D53 diagnostic: leave band 5 (L5/E5a/B2a) out of the signal set */
	int drop_band2;			/* diagnostic: leave GPS band 2 out, so the GPS ionosphere state is observed
					 * through L1/L5 only -- is the UM980's L2 the corrupted input (D105)? */
};

void vbs_netekf_default_opt(struct vbs_netekf_opt *opt);

struct vbs_netekf_stats {
	long epochs;
	long obs_used;
	long obs_rejected[VBS_NET_NREJ];
	long slips_master, slips_aux;	/* ambiguity states reset because of a cycle slip */
	long reject_resets;		/* ambiguity resets after a phase-rejection streak (D49) */
	long slots_evicted;		/* satellites dropped because the state vector was full */
	/*
	 * Innovation scale of the accepted observations. `norm` is the innovation divided by
	 * its own predicted standard deviation, so a correctly specified filter gives ~1: it
	 * says directly whether sig_code/sig_phase describe the errors the data actually has,
	 * which is not something the residuals alone can tell you.
	 */
	double innov_rms_code, innov_rms_phase;		/* m */
	double innov_norm_code, innov_norm_phase;	/* dimensionless */
	int nsat_active;		/* at the last epoch */
	int nband_max;			/* most bands carried for any one satellite */
	double ztd;			/* last estimate, m */
	double ztd_sigma;
};

struct vbs_netekf;

/*
 * master_xyz / aux_xyz are the two bases' known ECEF coordinates. Unlike rover RTK there is no
 * position state: in network processing both ends are surveyed points, which is what makes a
 * filter with this few states well conditioned over a 150 km baseline.
 *
 * Pass opt = NULL for the defaults.
 */
struct vbs_netekf *vbs_netekf_new(const double *master_xyz, const double *aux_xyz,
	const struct vbs_netekf_opt *opt);
void vbs_netekf_free(struct vbs_netekf *this);

/*
 * Process one epoch. `master`/`aux` are the two ends' observations for the same epoch; `rs`
 * holds satellite positions and velocities for the master's records, 6 doubles each, as
 * satposs() produces them, and `dts` the satellite clocks (2 per record) -- both are passed in
 * rather than computed here so that one satposs() call can serve every baseline of a cluster.
 *
 * `arc_master`/`arc_aux` are the two stations' arc ids (MAXSAT each, from
 * `vbs_preproc_arc_ids()`). The filter does **not** detect slips itself: it resets a
 * satellite's ambiguities when either station's arc id for it has changed since the previous
 * epoch. Slip detection is a property of a station's data stream, so it is done once per
 * station and shared -- D45 records why a per-baseline detector was wrong: the same slip landed
 * in different filters at different epochs, the corrections mixed inconsistent arcs, closure
 * differenced across arc boundaries, and no single detector could drive the output LLI.
 *
 * Returns the number of observations accepted, or -1 on error.
 */
/*
 * `rs` is indexed by SATELLITE NUMBER -- position and velocity of satellite `sat` at
 * rs + (sat - 1) * 6, zero when unknown -- not by observation record. Every filter of a
 * network is then handed the same geometry for the same satellite, whatever order its own
 * station's records come in (D42's 293 000-cycle misclosures were record-order indexing
 * meeting a differently ordered array). Build it with vbs_satpos_by_sat(). `dts` is unused.
 */
int vbs_netekf_epoch(struct vbs_netekf *this, gtime_t time,
	const obsd_t *master, int nm, const double *rs, const double *dts,
	const obsd_t *aux, int na, const nav_t *nav,
	const long *arc_master, const long *arc_aux);

/*
 * The §4.4 corrections for this auxiliary relative to the master, as of the last epoch:
 * ΔI on the constellation's reference frequency and ΔT non-dispersive, both in metres, with
 * their standard deviations. Only satellites whose states have converged past
 * `max_sigma` are marked valid -- §4.6 forbids synthesizing a satellite without a usable
 * correction, and an unconverged state is not a usable correction.
 *
 * `max_sigma` MUST be below sqrt(var_ion0), or the call returns 0 and marks nothing valid: a
 * threshold at or above the prior cannot tell an informed state from an untouched one. Pick
 * it from the error the correction has to beat rather than from the filter -- 0.02 m is the
 * sensible operating point on this network, against the 1.4-12.5 cm of ionospheric
 * interpolation error `docs/vbs/03-error-budget.md` measured across quiet and storm days.
 *
 * Writes at most maxout records; returns how many were written.
 */
/*
 * The de-geometried single-difference phase of one satellite on one band, as the last epoch
 * saw it, in metres, together with that band's ionosphere coefficient and wavelength (D66).
 *
 * This is what makes the ionosphere an *output* of ambiguity resolution rather than an input
 * to it. Between two bands of one satellite the clock and the troposphere cancel, leaving the
 * geometry-free combination; differenced again against a reference satellite the receiver
 * phase biases cancel too, and what remains is (mu2 - 1) times the double-differenced
 * ionosphere plus the two DD integer ambiguities. Once those integers are fixed the ionosphere
 * follows exactly, to phase noise -- no state, no prior, no convergence time, and no
 * dependence on baseline length.
 *
 * Returns 1 and fills whichever of zL/mu/lam are non-NULL, or 0 if that satellite-band was not
 * observed this epoch.
 */
int vbs_netekf_phase_sd(const struct vbs_netekf *this, int sat, int band, double *zL,
	double *mu, double *lam);

/*
 * Accumulate one epoch's ionosphere-free narrow-lane float for a satellite and return the
 * running mean and sample count (D67). The float is constant over an arc, so averaging is what
 * turns an observable too noisy to round in one epoch into one that rounds reliably in a few
 * minutes -- the standard way network ambiguity resolution reaches long baselines. The
 * accumulator restarts whenever the reference satellite, the band pair, or either station's
 * arc for that satellite changes.
 */
int vbs_netekf_if_average(struct vbs_netekf *this, int sat, int refsat, int band, double value,
	double *mean, int *count);

/* The wet mapping function and elevation this filter last used for a satellite. */
int vbs_netekf_slot_mapfw(const struct vbs_netekf *this, int sat, double *mapfw, double *el);

/* The estimated residual ZTD difference, metres. */
double vbs_netekf_ztd(const struct vbs_netekf *this);

/* The filter's zenith phase sigma, for callers propagating measurement noise. */
double vbs_netekf_sig_phase(const struct vbs_netekf *this);

int vbs_netekf_corrections(const struct vbs_netekf *this, struct vbs_corr *out, int maxout,
	double max_sigma);

/*
 * Standard deviations of the last correction set, parallel to vbs_netekf_corrections()'s
 * output: sd_I[k] and sd_T[k] belong to out[k]. §4.4 requires the variances to come out with
 * the corrections, and §4.5 needs them to weight the interpolation.
 */
int vbs_netekf_correction_sigmas(const struct vbs_netekf *this, const struct vbs_corr *corr,
	int ncorr, double *sd_I, double *sd_T);

void vbs_netekf_get_stats(const struct vbs_netekf *this, struct vbs_netekf_stats *out);

/*
 * --- Ambiguity resolution support (§4.3) ------------------------------------------------
 *
 * The filter's own ambiguities are single-difference and therefore not integers: they absorb
 * the two receivers' phase biases. Integer structure lives in the double differences, so
 * resolution works on DD combinations built from the SD state, and a fix is fed back as an
 * ordinary constrained update rather than by restructuring the state vector. That is what
 * makes "held once fixed" (§4.3) cheap: the constraint is simply reapplied.
 */

/*
 * A linear combination of bands, in cycles: sum(coef[i] * N_band[i]).
 *
 * The cascade is expressed entirely through these. {+1,-1} on L1/L2 is the wide lane
 * (λ ≈ 0.86 m); {+1,-1} on E5a/E5b is the extra-wide lane (λ ≈ 9.8 m); {+1,0} is a plain
 * band. Generalising over combinations rather than special-casing each cascade step is what
 * lets one code path serve dual- and triple-frequency baselines (§4.3 requires both, chosen
 * per baseline from the data).
 */
struct vbs_amb_comb {
	const char *name;
	int nband;
	int band[VBS_NET_MAXBAND];	/* RINEX band digits */
	double coef[VBS_NET_MAXBAND];
};

/* One double-difference ambiguity of some combination. */
struct vbs_dd_amb {
	int sat, refsat, sys;
	double value;			/* float DD ambiguity, cycles of the combination */
	double lam;			/* the combination's wavelength, m */
	double el;			/* the lower of the two satellites' elevations, rad */
};

/*
 * Build the DD set of `comb` for one constellation, against the highest-elevation satellite
 * that carries every band the combination needs. Writes the float values to `out` and their
 * covariance to `Q` (n×n, column-major, caller-allocated at maxout²).
 *
 * Returns the number of DD ambiguities, 0 if this constellation cannot form the combination
 * on this baseline (which is the normal answer for a dual-frequency baseline asked for an
 * extra-wide lane, and must not be treated as an error).
 */
int vbs_netekf_dd_set(const struct vbs_netekf *this, int sys, const struct vbs_amb_comb *comb,
	struct vbs_dd_amb *out, int maxout, double *Q);

/*
 * As `vbs_netekf_dd_set`, with `prefer_sat` used as the reference satellite when it is a valid
 * candidate for the combination. 0 means "choose as before", and reproduces it exactly.
 */
int vbs_netekf_dd_set_ref(const struct vbs_netekf *this, int sys, const struct vbs_amb_comb *comb,
	int prefer_sat, struct vbs_dd_amb *out, int maxout, double *Q);

/*
 * The satellite of `sys` carrying the most of `bands` (highest elevation breaking ties), or 0.
 * Lets the cascade hold one reference across all of a constellation's steps, which is what
 * `vbs_ambres_ion_from_fixes` needs and did not have (D131).
 */
int vbs_netekf_best_ref(const struct vbs_netekf *this, int sys, const int *bands, int nband);

/*
 * Constrain one DD combination to an integer: a pseudo-observation with variance `var`.
 *
 * Feeding a fix back this way rather than collapsing the state is deliberate. The constraint
 * propagates through the covariance into ΔI, ΔZTD and every other ambiguity exactly as the
 * information content says it should, and "hold" costs nothing but reapplying it. It also
 * degrades safely: a wrong fix is a wrong pseudo-observation, which the FDE gate and the
 * residuals can still see, whereas a wrong fix baked into the state structure is invisible.
 *
 * Returns 1 if applied, 0 if the ambiguity could not be located (satellite set, band gone).
 */
int vbs_netekf_constrain_dd(struct vbs_netekf *this, const struct vbs_dd_amb *a,
	const struct vbs_amb_comb *comb, double integer, double var);

/*
 * Standard deviation of the *double-differenced* slant ionosphere, ΔI^s - ΔI^ref, in metres,
 * for the highest-elevation satellite of each constellation as reference.
 *
 * This, not the single-difference sigma, is what decides whether the corrections are good
 * enough. A correction error common to every satellite enters the synthesized phase as
 * -μ_i·δ on all of them alike, and the rover differences satellites against each other, so
 * the common part cancels before it can reach a position. The single-difference sigma is
 * dominated by exactly that harmless common term -- it carries the reference satellite's own
 * code-determined level -- and reporting it overstates the error by an order of magnitude.
 *
 * Writes up to maxout values; returns how many. Satellites with no usable pair are skipped.
 */
int vbs_netekf_dd_ion_sigma(const struct vbs_netekf *this, int sys, double *out, int maxout);

/*
 * Per-satellite state, for diagnosing *why* a correction is or is not usable. The DD
 * ionosphere sigma alone says a satellite is bad; this says whether that is because it just
 * rose, sits too low for resolution to be attempted, carries too few bands, or keeps failing
 * the ratio test -- four problems with four different answers.
 */
struct vbs_netekf_satinfo {
	int sat, sys;
	int nband;
	long epochs_tracked;		/* consecutive epochs in the state */
	double el_deg;
	double sd_ion;			/* single-difference ΔI sigma, m */
	double amb_sigma_max;		/* worst ambiguity sigma of this satellite, cycles */
};

/*
 * Scatter satposs() output (one 6-vector per observation record) into a satellite-indexed
 * array of MAXSAT * 6 doubles, zero for satellites without a position.
 */
void vbs_satpos_by_sat(const obsd_t *obs, int n, const double *rs_rec, double *rs_sat);

/*
 * Report every satellite in the state at the last epoch. Returns how many were written.
 */
int vbs_netekf_satinfo(const struct vbs_netekf *this, struct vbs_netekf_satinfo *out,
	int maxout);

/*
 * Undo every fix on one satellite: its ambiguity states go back to their prior.
 *
 * This is what loop closure does with a misclosure (D33 item 2). Because a fix is held by
 * re-applying a tight constraint each epoch, "unfixing" has to widen the state again -- simply
 * not re-applying the constraint would leave the ambiguity pinned by the covariance from the
 * epochs when it *was* applied, and the next cascade pass would refix it to the same wrong
 * integer with a high ratio, because the filter is now confident in it.
 */
void vbs_netekf_unfix_sat(struct vbs_netekf *this, int sat);

/*
 * The station arc ids this filter's ambiguities for `sat` currently rest on, as last seen.
 *
 * With arcs shared per station (D45), every filter containing a station sees the same
 * boundary in the same epoch, so three baselines of a triangle *must* agree on these. Loop
 * closure asserts that rather than assuming it: a disagreement here would mean the sharing is
 * broken, not that a fix is wrong. Returns 0 if the satellite is not in the state.
 */
int vbs_netekf_station_arcs(const struct vbs_netekf *this, int sat, long *arc_master,
	long *arc_aux);


/*
 * The common signal set the filter settled on for this baseline, and how many bands the
 * weakest constellation has -- §4.3 requires the per-baseline signal set to be determined
 * from the data and logged, because it is what selects dual- vs triple-frequency resolution.
 */
const struct vbs_sigset *vbs_netekf_sigset(const struct vbs_netekf *this);

/*
 * Which RINEX code the filter would actually take for one band of one satellite.
 *
 * The band-to-slot rule lives in one static function inside the filter, and until now nothing
 * outside could ask what it chose. That is how a whole investigation went astray: a UM980
 * publishes L2C twice, once honestly as `2X` and once as `2W`, and whether a baseline ends up
 * differencing 2W against 2W or 2X against 2L was invisible from the outside. Reports that
 * guess at it are worse than none, so this exposes the decision itself rather than a copy of
 * the rule.
 *
 * Returns CODE_NONE when the satellite carries no usable observation on that band.
 */
uint8_t vbs_netekf_band_code(const obsd_t *o, int band, int prefer_l2c);

#endif /* __VBS_NETEKF_H__ */

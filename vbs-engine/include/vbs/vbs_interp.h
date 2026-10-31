#ifndef __VBS_INTERP_H__
#define __VBS_INTERP_H__

#include "rtklib.h"

#include "vbs/vbs_synth.h"

/*
 * Spatial interpolation of the network corrections at the VBS point -- spec §4.5.
 *
 * Each auxiliary k contributes, per satellite, a correction relative to the master: ΔI_k^s and
 * ΔT_k^s from `vbs_netekf`. Those are values of a field sampled at the auxiliaries' positions,
 * and the VBS point needs that field evaluated where it sits. The baseline method §4.5 asks
 * for is a weighted plane fit -- Wanninger's linear interpolation method -- which is what this
 * is: for each satellite, fit
 *
 *     Δ(E, N) = a + b·E + c·N
 *
 * in a local frame centred on the *master* (so a is the master-relative correction at the
 * master, which is zero by construction and is therefore not forced), and evaluate at the VBS.
 *
 * Three things this gets right that a naive plane fit does not:
 *
 *  - **Weighting by the network's own variances.** §4.4 delivers a variance with every
 *    correction and §4.5 needs them; an auxiliary whose ambiguities never fixed should not
 *    pull the plane as hard as one whose did. Weights are 1/σ², so an unconverged correction
 *    fades out rather than having to be excluded by a threshold.
 *  - **The troposphere height term is modelled, not interpolated.** §4.5 is explicit: the
 *    a priori model difference between master and VBS *heights* is computed from the
 *    troposphere model, and only the residual is interpolated. Interpolating the whole
 *    non-dispersive term across stations at different altitudes would smear a height-driven
 *    signal across a horizontal plane, which is exactly the error the model already knows how
 *    to avoid.
 *  - **An honest extrapolation metric.** The fitted value's own variance inflation at the VBS
 *    (the leverage, A⁻¹[0,0]) is reported alongside the correction. Inside a well-spread
 *    cluster it is ~1/n; outside the hull it grows without bound. `04-geophysical-corrections.md`
 *    measured it at 0.25-0.35 for enclosed clusters against 6-12 for the extrapolated cases,
 *    which is the difference §4.5's `vbs.max_extrapolation_km` has to be calibrated against.
 */

#define VBS_INTERP_MAXAUX	16

/*
 * One auxiliary's contribution: its position and its per-satellite corrections.
 *
 * Zero it before filling it in (memset, or a static array). Fields are added to this struct
 * over time -- `backup` was, on 2026-10-02 -- and a caller that sets only the fields it knows
 * leaves the new ones holding whatever was on the stack: a stray non-zero `backup` silently
 * turns a primary into a reserve. The unit tests filled it field by field and broke exactly so.
 */
struct vbs_interp_aux {
	double xyz[3];
	const struct vbs_corr *corr;
	const double *sd_I;		/* per-correction σ, metres; NULL for unweighted */
	const double *sd_T;
	int ncorr;
	int backup;			/* held in reserve: fitted only when the primaries fall short */
};

struct vbs_interp_opt {
	double max_sigma;		/* ignore a contribution worse than this, metres */
	/*
	 * Apply the *estimated* non-dispersive residual (the interpolated m_w·ΔZTD), as
	 * distinct from the troposphere model term, which is deterministic geometry and is
	 * always applied.
	 *
	 * These two halves deserve separate switches because they have completely different
	 * standing. The model term is computed from station heights and a standard atmosphere:
	 * it is as reliable as the coordinates. The residual is a single scalar per baseline
	 * that the filter estimates, and it absorbs whatever else the filter cannot explain --
	 * so it is only worth applying where it measurably reduces the error, which D33 item 1
	 * requires to be demonstrated rather than assumed.
	 */
	int apply_nondispersive;
	int min_aux;			/* fewer usable auxiliaries than this -> no output */
	double humi;			/* for the a priori troposphere height term */
	int datum;			/* VBS_INTERP_DATUM_*, see below */
	/*
	 * Diagnostic scale on the interpolated dispersive correction. The independent truth
	 * check reports the applied correction as consistently ~0.75 of the real differential
	 * ionosphere; if that under-correction is what leaves a horizontal bias growing with the
	 * master-to-VBS displacement, scaling it should move the bias and nothing else. 1 = off.
	 */
	double iono_scale;
	/*
	 * Ionosphere consistency gate, metres; 0 = off (the default, and what is served).
	 *
	 * The network fixes wide lanes only (D132), so each auxiliary's per-satellite ionosphere
	 * takes its level from code, and one station's code error can sit in one satellite's
	 * estimate as a constant of tens of centimetres (D164: 15-56 cm on GB85's cluster). The
	 * fit then hands it to the VBS, where a rover's float ambiguity absorbs it and the integer
	 * never comes. With the gate on, each auxiliary is predicted from the others plus the
	 * master anchor (a plane, leave-one-out); while the worst misfit exceeds this value and
	 * more than three auxiliaries remain, that one is dropped for this satellite and the rest
	 * refitted. If three remain and they still disagree, the satellite gets no dispersive
	 * correction at all -- what a UM980 VBS always gets (D106) -- rather than a dragged one.
	 */
	double iono_gate_m;
	/*
	 * Include the master itself in the plane fit, as a point of value zero (D69). The
	 * corrections are differences against the master, so the field is zero there by
	 * definition; a free plane through three auxiliaries discards that and lets a spurious
	 * offset and tilt through, which reaches the rover as a horizontal bias growing with the
	 * master-to-VBS displacement.
	 */
	int anchor_master;
	int method;			/* VBS_INTERP_METHOD_* */
	double corr_len_km;		/* collocation: e-folding distance of the field */
	double sig_field_m;		/* collocation: field standard deviation, metres */
};

/*
 * Remove each auxiliary's own ionosphere datum before fitting the plane (D55).
 *
 * A single-difference slant ionosphere state is not observable on its own: it absorbs the two
 * receivers' differential code biases, so `dI` for baseline master->aux_j is the true
 * differential ionosphere plus a constant of that *pair*. Measured on the corrections
 * themselves, that constant is 0.5-2.5 m, differs per auxiliary and per constellation, and
 * drifts 16-95 cm over a day while stepping less than a centimetre between epochs -- too slow
 * for any slip test to see, far too large to interpolate.
 *
 * Fitting a plane through three such values is fitting three receiver biases: with three
 * auxiliaries and three plane parameters the fit is exact, so the VBS receives a weighted
 * average of them (weights summing to one, and outside the hull not even bounded by the
 * inputs). The comment at the fit -- "the field is zero at the master by construction" -- is
 * true of the ionosphere and false of the state that carries it.
 *
 * With this set, each auxiliary's constellation datum is estimated from the satellites that
 * every participating auxiliary carries this epoch and subtracted before the fit. What
 * reaches the VBS is then a double difference in all but name, and the constellation-common
 * term it no longer carries is the one term a rover can absorb.
 */
/*
 * How the correction field is carried from the auxiliaries to the VBS point (D70).
 *
 * PLANE is §4.5's baseline: Wanninger's weighted linear fit. It is exact for a field that
 * really is planar, and the measured residual bias says the field is not -- the error grows at
 * about 0.19 cm per km of master-to-VBS displacement, survives anchoring, and disappears when
 * more samples are added (at a cost in noise that is worse than the bias).
 *
 * COLLOCATION is least-squares collocation: predict the field at the VBS as
 *
 *     s0 = c0' (C + N)^-1 s
 *
 * with C(d) an isotropic covariance of the field between samples, c0 its value between each
 * sample and the VBS, and N the samples' own noise. It carries curvature without needing more
 * parameters than samples, and it weighs a distant or noisy auxiliary by how much the field at
 * the VBS is actually expected to resemble it -- the two failures of the plane fit in one
 * estimator. Distance correlation is exponential, which is the usual choice for tropospheric
 * and ionospheric fields at these scales.
 */
#define VBS_INTERP_METHOD_PLANE		0
#define VBS_INTERP_METHOD_COLLOCATION	1

#define VBS_INTERP_DATUM_OFF	0
#define VBS_INTERP_DATUM_COMMON	1	/* weighted mean over the commonly-held satellites */

void vbs_interp_default_opt(struct vbs_interp_opt *opt);

/*
 * What the consistency gate did since the last reset (D90: a path with no visible count is a
 * path nobody can tell is running): satellite fits it examined, auxiliaries it dropped from a
 * satellite's fit, satellites whose dispersive correction it withheld.
 */
void vbs_interp_gate_stats(long *checked, long *dropped, long *withheld);
void vbs_interp_gate_reset(void);

/* Per-satellite diagnostics, parallel to the corrections the call produces. */
struct vbs_interp_quality {
	int sat;
	int naux;			/* auxiliaries that actually contributed */
	int backup_used;		/* the reserve was needed for this satellite */
	double leverage;		/* geometry-only variance inflation at the VBS: ~1/n
					   inside the cluster, unbounded outside it */
	double sd_I, sd_T;		/* propagated standard deviations, metres */
	double dT_height;		/* the modelled height term, metres, for the record */
};

/*
 * Interpolate to `vbs_xyz`. `master_xyz` anchors the local frame and supplies the
 * troposphere height reference.
 *
 * `azel` is an array of 2*MAXSAT doubles indexed by (satellite number - 1), giving each
 * satellite's azimuth and elevation in radians; pass NULL to omit the troposphere model term
 * entirely. The *same* angles are used for both master and VBS deliberately: the difference
 * this term must carry is the one driven by the height and position of the two points, and
 * evaluating it at a common elevation isolates exactly that. Master and VBS are tens of km
 * apart against a satellite 20 000 km away, so their true elevations differ by about a tenth
 * of a degree, and letting that difference in would add an angular effect on top of the
 * height effect the term is for. Writes corrections to `out` (and diagnostics to `qual`, if
 * given), at most `maxout`; returns how many.
 *
 * A satellite is emitted only where enough auxiliaries carry a valid correction for it --
 * §4.6 forbids synthesizing a satellite without one.
 */
int vbs_interp_epoch(gtime_t time, const double *master_xyz, const double *vbs_xyz,
	const struct vbs_interp_aux *aux, int naux, const struct vbs_interp_opt *opt,
	const double *azel, struct vbs_corr *out, struct vbs_interp_quality *qual,
	int maxout);

#endif /* __VBS_INTERP_H__ */

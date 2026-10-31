#ifndef __VBS_C2WCAL_H__
#define __VBS_C2WCAL_H__

#include <stdint.h>

#include "rtklib.h"

/*
 * Per-base GPS L2 P(Y) (C2W) code-bias correction for UM980 bases.
 *
 * A UM980 publishes two GPS L2 code columns, `2W` (P(Y)) and `2X` (L2C). The engine reads `2W`
 * (prefer_l2c off, D140), and on a UM980 `2W` sits a per-unit offset away from `2X`: 3.05 m of
 * spread across the fleet (D159/D161). At a given time the offset is one number per unit, the same
 * on every satellite including the L2C-less ones (a UM980 against a Mosaic on a shared
 * antenna, 2026-10-17). But it is not a calibration constant: D161 measured units drifting by
 * metres over days without losing lock (one by 3.9 m in five days) -- the P(Y) code moves, L2C does
 * not -- and a 0.3-0.4 m diurnal term locked to the satellite geometry (differential multipath
 * between the two codes). D160 removed a static table of these offsets (one 27 Sep snapshot applied
 * to days 64 days apart) and lost. So this module never takes the offset as a configured constant:
 * each base measures its own, continuously, from its own stream.
 *
 * Two corrections, chosen by the caller:
 *   - vbs_c2wcal_apply(): subtract the base's current 2W-2X estimate from every GPS P(Y) code,
 *     bringing `2W` to that unit's L2C level on all satellites;
 *   - vbs_c2wcal_apply_l2c(): on a satellite carrying both codes, put the L2C code itself into the
 *     P(Y) slot (no estimate needed, so no drift and no diurnal term); only on an L2C-less satellite
 *     subtract the estimate. The carrier phase is untouched either way: the two columns' phases are
 *     the same L2 carrier to under a centimetre (D159).
 *
 * The gate. Vendor bit 26 says the P(Y) slot can carry L2C (D159), in which case `2W` is not a
 * second signal at all and there is nothing to correct. That case is recognised by what l2pair.py
 * defined it as -- one signal twice, bit-identical columns -- not by the size of the offset: about
 * 40 % of genuine units sit within +-0.2 m of zero (D161), and a Mosaic's genuine semi-codeless `2W`
 * sits 0.025 m from its L2C, so a gate on |mean| refuses real P(Y). Whether a receiver relabels is a
 * property of its firmware build (behaviour is coherent within a build), so
 * the verdict is pooled over every base of a build and applies to all of them; a base whose build is
 * unknown is judged on its own samples.
 */

/* GPS L2 P(Y) code labels (CODE_L2W/L2P/L2Y) and L2C labels (CODE_L2S/L2L/L2X). CODE_L2P is also
 * GLONASS G2P, so every function below checks the satellite is GPS. */
int vbs_c2wcal_is_l2_py_code(uint8_t code);
int vbs_c2wcal_is_l2c_code(uint8_t code);

/* --- the gate --------------------------------------------------------------------------- */

#define VBS_C2WCAL_MIN_SAMPLES	300	/* P(Y)+L2C satellite-epochs before a verdict */
#define VBS_C2WCAL_IDENT_M	0.0005	/* |2W-2X| below this is the same number (RINEX has 1 mm) */
#define VBS_C2WCAL_RELABEL_FRAC	0.5	/* identical in at least this share: one signal twice */
#define VBS_C2WCAL_BLUNDER_M	20.0	/* a 2W-2X beyond this is a slip or a wrong fix, not a bias */

struct vbs_c2wcal_probe {
	double sum_m;		/* sum of 2W-2X over the `both` samples, metres */
	long n_both;		/* GPS satellite-epochs with both a P(Y) and an L2C code */
	long n_ident;		/* ... of which the two codes are the same number */
	long n_py_only;		/* GPS satellite-epochs with P(Y) and no L2C */
};

void vbs_c2wcal_probe_reset(struct vbs_c2wcal_probe *p);
void vbs_c2wcal_probe_add_epoch(struct vbs_c2wcal_probe *p, const obsd_t *obs, int nobs);

/* 1: `2W` is a distinct signal, correct it. 0: `2W` is L2C again, leave it. -1: not enough samples.
 * *mean_m, if non-NULL, receives the mean 2W-2X (NAN when n_both == 0). */
int vbs_c2wcal_probe_verdict(const struct vbs_c2wcal_probe *p, double *mean_m);

/* --- pooling the gate per firmware build ---------------------------------------------------- */

#define VBS_C2WCAL_MAXGROUP	64

struct vbs_c2wcal_groups {
	int n;
	long build[VBS_C2WCAL_MAXGROUP];
	struct vbs_c2wcal_probe probe[VBS_C2WCAL_MAXGROUP];
};

void vbs_c2wcal_groups_reset(struct vbs_c2wcal_groups *g);

/* The pooled probe of a build (vbs_policy_fw_build()), created on first use. NULL when the build
 * is unknown (< 0) or the table is full: the caller then judges that base on its own probe. */
struct vbs_c2wcal_probe *vbs_c2wcal_group(struct vbs_c2wcal_groups *g, long build);

/* --- the estimate: one base's own 2W-2X, causal ------------------------------------------- */

#define VBS_C2WCAL_BUCKET_S	60	/* samples are binned per minute */
#define VBS_C2WCAL_MAXBUCKET	1440	/* so the longest trailing window is a day */
#define VBS_C2WCAL_EST_MIN	60	/* samples in the window before the estimate is used */

struct vbs_c2wcal_est {
	int nbucket;			/* trailing window in buckets; 0 = everything seen so far */
	long id[VBS_C2WCAL_MAXBUCKET];	/* minute number held by each ring slot, -1 empty */
	double sum[VBS_C2WCAL_MAXBUCKET];
	long n[VBS_C2WCAL_MAXBUCKET];
	double cum_sum;
	long cum_n;
};

/* window_s <= 0: cumulative (the estimate is the mean of everything added). Otherwise rounded up to
 * whole minutes and capped at a day. */
void vbs_c2wcal_est_init(struct vbs_c2wcal_est *e, double window_s);
void vbs_c2wcal_est_add_epoch(struct vbs_c2wcal_est *e, gtime_t t, const obsd_t *obs, int nobs);

/* The mean 2W-2X over the window ending at t (the minute holding t included). Returns 1 and sets
 * *mean_m when the window holds at least VBS_C2WCAL_EST_MIN samples, else 0. *n, if non-NULL,
 * receives the sample count either way. */
int vbs_c2wcal_est_value(const struct vbs_c2wcal_est *e, gtime_t t, double *mean_m, long *n);

/* --- the corrections (the caller applies the gate first) ---------------------------------- */

/* Subtract bias_m from every GPS P(Y) code in the epoch. No-op at 0. Returns codes changed. */
int vbs_c2wcal_apply(obsd_t *obs, int nobs, double bias_m);

/* Per GPS satellite: with both codes, the P(Y) code becomes the L2C code; with P(Y) only, it has
 * bias_m subtracted. *n_l2c / *n_py (either may be NULL) count the two cases. Returns the total. */
int vbs_c2wcal_apply_l2c(obsd_t *obs, int nobs, double bias_m, int *n_l2c, int *n_py);

/* --- one live base: gate + estimate + correction, as the server runs it ------------------ */

#define VBS_C2WCAL_CONST	1	/* subtract the trailing estimate from every GPS P(Y) code */
#define VBS_C2WCAL_L2C		2	/* L2C code where both are tracked, the estimate elsewhere */

struct vbs_c2wcal_base {
	int mode;			/* VBS_C2WCAL_CONST or VBS_C2WCAL_L2C */
	struct vbs_c2wcal_probe probe;	/* this base's own samples, the gate when its build is unknown */
	struct vbs_c2wcal_est est;
};

void vbs_c2wcal_base_init(struct vbs_c2wcal_base *b, int mode, double window_s);

/*
 * One epoch of one base, in arrival order. Feeds the base's own probe, its build's pooled probe
 * (g may be NULL; build from vbs_policy_fw_build(), < 0 unknown) and its estimate, then -- when
 * the gate says "distinct" and the estimate has settled -- writes the corrected epoch to out
 * (at most nmax records) and returns its record count; otherwise returns 0 and out is untouched,
 * and the caller uses the epoch as it came. *verdict and *offset_m (either may be NULL) receive
 * the gate verdict used (vbs_c2wcal_probe_verdict) and the estimate applied (NAN when none).
 */
int vbs_c2wcal_base_epoch(struct vbs_c2wcal_base *b, struct vbs_c2wcal_groups *g, long build,
	const obsd_t *obs, int n, obsd_t *out, int nmax, int *verdict, double *offset_m);

#endif

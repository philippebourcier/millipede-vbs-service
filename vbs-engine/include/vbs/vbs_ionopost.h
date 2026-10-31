#ifndef __VBS_IONOPOST_H__
#define __VBS_IONOPOST_H__

#include "rtklib.h"
#include "vbs/vbs_synth.h"

/*
 * What a constellation's dispersive correction becomes after interpolation, before synthesis.
 *
 * D167 serves the network ionosphere on Galileo and BeiDou only: a UM980 cluster's GPS
 * ionosphere is estimated through its L2 P(Y) observable and is erratic (correlation with the
 * truth +0.10 to +0.77 where Galileo's is +0.71 to +0.80), so GPS gets none. Whether GPS could be given a good one
 * instead was asked on 2026-10-19. This module holds the candidates,
 * every one off by default:
 *
 *   need_band5_sys  keep a constellation's own correction only on satellites the master tracks
 *                   on band 5 (GPS L5: the IIF/III satellites). With the network filter's GPS L2
 *                   dropped (`--no-band2`) those satellites' ionosphere rests on L1/L5 alone;
 *   scale_gps       GPS's own correction scaled (it tracks the truth partly);
 *   model_sys       the correction of these constellations predicted from the others': a thin
 *                   shell (350 km) fitted per epoch to the kept satellites' corrections -- one
 *                   datum per constellation plus a mapped vertical term and its north/east
 *                   gradient over the pierce points (model_mode 0), plus an inverse-distance
 *                   interpolation of the fit's residuals between pierce points (1), or a
 *                   quadratic surface instead of the plane (2). Fills satellites whose own
 *                   correction was withheld; with model_check_m > 0 it also replaces an own
 *                   correction that departs from the prediction by more than that, after the
 *                   constellation's median offset to the prediction is removed;
 *   drop_sys        these constellations are not output at all;
 *   shrink_m        each kept correction's departure d from its constellation's median (on L1)
 *                   scaled by d^2 / (d^2 + shrink_m^2): a large one passes nearly whole, one near
 *                   the estimate's noise fades toward the median. Day 251's leave-one-out (2026-10-19) showed E+C losing where the served
 *                   arm already fixes > 90 % and winning where it struggles: a small correction is
 *                   mostly noise. Continuous in dI, so it makes no jump at the rover.
 *
 * The datum per constellation is free because a rover double-differences within a
 * constellation; a model prediction therefore carries none.
 */

struct vbs_ionopost_opt {
	int sys_mask;			/* constellations keeping their own correction; 0 = every one */
	int need_band5_sys;
	double scale_gps;		/* 0 or 1: unscaled */
	int model_sys;
	int model_mode;			/* VBS_IONOPOST_PLANE, _IDW, _QUAD */
	double model_check_m;
	int drop_sys;
	double shrink_m;
};

#define VBS_IONOPOST_PLANE	0
#define VBS_IONOPOST_IDW	1
#define VBS_IONOPOST_QUAD	2

struct vbs_ionopost_stats {
	long fits, fit_refused;		/* epochs fitted, epochs with too few satellites to fit */
	long filled, replaced, kept;	/* model-system corrections filled, replaced, kept own */
	long band5_withheld;
	long shrunk;			/* corrections scaled by shrink_m */
	double shrink_w_sum;		/* sum of their weights */
	double fit_rms_sum;		/* sum over fits of the fit's residual RMS, m */
};

/* 1 when any option differs from the D167 behaviour of `sys_mask` alone. */
int vbs_ionopost_active(const struct vbs_ionopost_opt *opt);

/*
 * Apply the options to one epoch's corrections in place (dI and valid only). azel is indexed
 * (sat - 1) * 2 at the VBS, vpos its geodetic position; mep/nm the master's epoch, for band 5.
 * Returns the number of corrections changed.
 */
int vbs_ionopost_apply(struct vbs_corr *corr, int n, const double *azel, const double *vpos,
	const obsd_t *mep, int nm, const nav_t *nav, const struct vbs_ionopost_opt *opt,
	struct vbs_ionopost_stats *stats);

/* Parse a constellation letter string (GRECJ) into a SYS_* mask; -1 on a bad letter. */
int vbs_ionopost_parse_sys(const char *s);

#endif

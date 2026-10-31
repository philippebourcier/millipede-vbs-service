#ifndef __VBS_SYNTH_H__
#define __VBS_SYNTH_H__

#include "rtklib.h"

/*
 * Observation synthesis -- spec S4.6.
 *
 * Translates the master's observations to the VBS point:
 *
 *   drho^s      = |X^s - X_VBS| - |X^s - X_Master|          (with Sagnac correction)
 *   Phi_VBS,i^s = Phi_Master,i^s + drho^s + dT^s - mu_i dI^s
 *   P_VBS,i^s   = P_Master,i^s   + drho^s + dT^s + mu_i dI^s
 *
 * with mu_i = (f1/fi)^2 and dI expressed on the constellation's reference frequency f1.
 * Note the opposite ionosphere sign on phase and code.
 *
 * The critical requirement (spec S4.6) is that the translation preserves *exactly* the
 * integer nature of the master's ambiguity: every term added to the phase is a known
 * geometric or atmospheric quantity in metres, converted to cycles by the signal's own
 * wavelength. Nothing here introduces a non-integer bias, and nothing here introduces a
 * discontinuity -- discontinuities are signalled separately through the lock time (S5.3).
 *
 * Phases are carried in cycles (as obsd_t does), so the metre-domain correction is divided
 * by that signal's wavelength. Code observations stay in metres.
 */

/*
 * Interpolated network correction for one satellite at the VBS point, relative to the
 * master (spec S4.4/S4.5). Both terms are in metres.
 */
struct vbs_corr {
	int sat;		/* satellite number (RTKLIB numbering) */
	double dI;		/* slant ionosphere on the reference frequency f1 */
	double dT;		/* non-dispersive: troposphere residual + orbit */
	int valid;		/* 0: no usable correction, satellite must not be output */
};

/*
 * Reference frequency f1 on which dI is expressed, for one constellation.
 *
 * This is a convention internal to the engine: as long as extraction (S4.4) and application
 * (S4.6) agree on it, any consistent choice works. We use each constellation's L1-equivalent
 * carrier. Returns 0.0 for an unsupported system.
 */
double vbs_synth_ref_freq(int sys);

/*
 * Geometric range difference drho for one satellite, VBS relative to master.
 *
 * rs/vs are the satellite position and velocity in ECEF at the transmission time that
 * corresponds to the *master's* reception of this signal (i.e. straight out of satposs()).
 * The VBS sees the same signal at the same reception epoch but over a different range, so it
 * was transmitted at a slightly different time: over a 30 km displacement that is ~100 us,
 * during which the satellite moves ~0.4 m and the range to it can change by several
 * centimetres. Ignoring that is not acceptable at the mm level this engine targets, so the
 * satellite position is propagated to the VBS's own transmission time and the range is
 * re-converged. Both ranges include the Sagnac (earth-rotation) correction via geodist().
 */
double vbs_synth_drho(const double *rs, const double *vs,
	const double *master_xyz, const double *vbs_xyz);

/*
 * Synthesize one epoch of VBS observations from the master's epoch.
 *
 * master/n   master's observations for one epoch
 * rs         satellite positions+velocities, 6 doubles per record, from satposs()
 * corr/ncorr per-satellite interpolated corrections; a satellite with no valid correction is
 *            not output at all (spec S4.6: output only satellites with valid corrections).
 *            Pass ncorr = 0 to synthesize with zero corrections, which is what the
 *            zero-baseline test needs.
 * out        caller-provided array of at least n records
 *
 * Returns the number of records written, or -1 on error.
 */
int vbs_synth_epoch(const obsd_t *master, int n, const double *rs,
	const double *master_xyz, const double *vbs_xyz,
	const struct vbs_corr *corr, int ncorr,
	const nav_t *nav, obsd_t *out);

#endif /* __VBS_SYNTH_H__ */

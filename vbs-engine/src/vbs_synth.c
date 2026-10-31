#include <math.h>
#include <string.h>

#include "vbs/vbs_signal.h"
#include "vbs/vbs_synth.h"

/*
 * Observation synthesis, spec S4.6. See the header for the model.
 */

double vbs_synth_ref_freq(int sys) {
	switch (sys) {
	case SYS_GPS:
	case SYS_GAL:
	case SYS_QZS:
		return FREQ1;		/* L1 / E1 */
	case SYS_CMP:
		return FREQ1_CMP;	/* B1I */
	case SYS_GLO:
		/*
		 * GLONASS is FDMA, so there is no single carrier frequency for the constellation.
		 * What this function supplies is not a carrier but a *datum*: the ionosphere state
		 * is a slant delay expressed at this frequency, and each observation's coefficient
		 * is mu = (fref/f)^2 with f taken per satellite from sat2freq(). Any consistent
		 * choice of datum works, provided the filter and the synthesis make the same one --
		 * and they do, because both call this function.
		 *
		 * The k = 0 channel is the natural choice: it is the constellation's nominal G1 and
		 * it keeps GLONASS mu values in the same range as the CDMA systems', so nothing
		 * downstream has to special-case their magnitude.
		 *
		 * This used to return 0 to refuse GLONASS outright, deliberately, so that a caller
		 * feeding it would notice rather than silently receive a wrong mu. That guard has
		 * been replaced by a real convention now that vbs_signal_sat_freq() exists to give
		 * the per-satellite frequency this datum is scaled against (D90).
		 */
		return FREQ1_GLO;
	}
	return 0.0;
}

double vbs_synth_drho(const double *rs, const double *vs,
	const double *master_xyz, const double *vbs_xyz) {
	double e[3], rs_vbs[3];
	double r_master, r_vbs, dt;
	int i, it;

	/*
	 * geodist() returns -1 when the satellite position is unavailable (satposs() leaves it
	 * zeroed if there is no usable ephemeris at this epoch). Returning 0.0 for that would
	 * be a trap: the caller would emit the master's phase *untranslated*, which looks like a
	 * perfectly ordinary observation but carries the entire master-to-VBS geometry as
	 * error. Returning NAN forces the caller to deal with it.
	 */
	r_master = geodist(rs, master_xyz, e);
	if (r_master <= 0.0)
		return NAN;

	/*
	 * First pass with the master's satellite position, then propagate the satellite to the
	 * VBS's own transmission time and re-converge. t_tx_vbs - t_tx_master = -drho/c, so a
	 * VBS further from the satellite than the master sees an *earlier* transmission.
	 * Two iterations are far more than enough: the correction to drho is second order and
	 * already sub-micrometre after the first.
	 */
	memcpy(rs_vbs, rs, sizeof rs_vbs);
	r_vbs = geodist(rs_vbs, vbs_xyz, e);
	if (r_vbs <= 0.0)
		return NAN;

	for (it = 0; it < 2; it++) {
		dt = (r_vbs - r_master) / CLIGHT;
		for (i = 0; i < 3; i++)
			rs_vbs[i] = rs[i] - dt * vs[i];
		r_vbs = geodist(rs_vbs, vbs_xyz, e);
		if (r_vbs <= 0.0)
			return NAN;
	}
	return r_vbs - r_master;
}

/* Look up the correction for one satellite; NULL if it has none. */
static const struct vbs_corr *find_corr(const struct vbs_corr *corr, int ncorr, int sat) {
	int i;

	for (i = 0; i < ncorr; i++)
		if (corr[i].sat == sat)
			return &corr[i];
	return NULL;
}

int vbs_synth_epoch(const obsd_t *master, int n, const double *rs,
	const double *master_xyz, const double *vbs_xyz,
	const struct vbs_corr *corr, int ncorr,
	const nav_t *nav, obsd_t *out) {
	const struct vbs_corr *c;
	double drho, f1, fi, lam, mu, dI, dT, dphi, dp;
	int i, j, sys, nout = 0;

	if (master == NULL || out == NULL || master_xyz == NULL || vbs_xyz == NULL)
		return -1;

	for (i = 0; i < n; i++) {
		sys = satsys(master[i].sat, NULL);
		if ((f1 = vbs_synth_ref_freq(sys)) == 0.0)
			continue;

		/*
		 * Spec S4.6: output only satellites that are tracked by the master AND have valid
		 * interpolated corrections. With ncorr == 0 the caller is explicitly asking for
		 * zero corrections (the zero-baseline test), so every satellite passes.
		 */
		dI = dT = 0.0;
		if (ncorr > 0) {
			if ((c = find_corr(corr, ncorr, master[i].sat)) == NULL || !c->valid)
				continue;
			dI = c->dI;
			dT = c->dT;
		}

		/*
		 * No usable satellite position means no translation is possible, so the satellite
		 * is not output at all. Spec S4.6 already requires exactly this for a satellite
		 * without valid corrections; the same reasoning applies to one without geometry.
		 * Emitting it untranslated would hand the rover an observation that silently
		 * belongs to the master's position instead of the VBS's.
		 */
		drho = (rs != NULL) ? vbs_synth_drho(rs + 6 * i, rs + 6 * i + 3,
		    master_xyz, vbs_xyz) : 0.0;
		if (isnan(drho))
			continue;

		out[nout] = master[i];		/* carries time, sat, SNR, LLI, code */

		for (j = 0; j < NFREQ + NEXOBS; j++) {
			if (master[i].code[j] == CODE_NONE)
				continue;
			if ((fi = sat2freq(master[i].sat, master[i].code[j], nav)) == 0.0) {
				/* Unknown frequency: we cannot translate it, so do not emit it. */
				out[nout].code[j] = CODE_NONE;
				out[nout].L[j] = out[nout].P[j] = 0.0;
				out[nout].D[j] = 0.0f;
				out[nout].SNR[j] = 0;
				out[nout].LLI[j] = 0;
				continue;
			}
			mu = (f1 / fi) * (f1 / fi);
			lam = CLIGHT / fi;

			dphi = drho + dT - mu * dI;	/* metres */
			dp   = drho + dT + mu * dI;	/* metres */

			/*
			 * Phase is in cycles: convert the metre-domain correction with this
			 * signal's own wavelength. The master's integer ambiguity rides through
			 * untouched, which is the S4.6 requirement.
			 */
			if (master[i].L[j] != 0.0)
				out[nout].L[j] = master[i].L[j] + dphi / lam;
			if (master[i].P[j] != 0.0)
				out[nout].P[j] = master[i].P[j] + dp;
		}
		nout++;
	}
	return nout;
}

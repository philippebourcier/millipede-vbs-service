/*
 * vbs_netrun -- run the network EKF and the §4.3 ambiguity cascade over one baseline.
 *
 * The per-baseline counterpart of tools/vbs_loo: where vbs_loo measures what the corrections
 * *must* remove, this measures what the filter actually achieves on that baseline -- fix rates
 * per cascade step, and the double-differenced ionosphere sigma, which D25 establishes as the
 * metric that decides whether the corrections are usable (the single-difference figure is
 * dominated by a term that cancels at the rover).
 *
 * Takes coordinates and paths as arguments rather than reading a manifest, for the same reason
 * vbs_loo does: manifests are validation-dataset metadata, and the engine's real basemeta
 * module reads the network catalogue. A Python driver iterates cases.
 *
 * Usage:
 *   vbs_netrun <master.rnx> <mx> <my> <mz> <aux.rnx> <ax> <ay> <az> <nav.rnx> [label]
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_ambres.h"
#include "vbs/vbs_netekf.h"
#include "vbs/vbs_preproc.h"

static int epoch_slice(const obs_t *obs, int *i, const obsd_t **out) {
	int start = *i, n = 0;

	if (start >= obs->n)
		return 0;
	*out = &obs->data[start];
	while (start + n < obs->n &&
	    fabs(timediff(obs->data[start + n].time, obs->data[start].time)) <= 1e-3)
		n++;
	*i = start + n;
	return n;
}

static void span(const double *v, int n, double *lo, double *hi) {
	int i;

	*lo = 1e9;
	*hi = 0.0;
	for (i = 0; i < n; i++) {
		if (v[i] < *lo)
			*lo = v[i];
		if (v[i] > *hi)
			*hi = v[i];
	}
	if (n == 0)
		*lo = *hi = 0.0;
}

int main(int argc, char **argv) {
	static obs_t mobs, aobs;
	static nav_t nav;
	static double rs[MAXOBS * 6], dts[MAXOBS * 2], var[MAXOBS], rs_sat[MAXSAT * 6];
	static int svh[MAXOBS];
	static struct vbs_ambres_stats ar;
	static const char *SN[VBS_AMB_NSTEP] = {"extra-wide", "wide", "narrow"};
	static const int SL[3] = {SYS_GPS, SYS_GAL, SYS_CMP};
	static const char *SLN[3] = {"GPS", "GAL", "BDS"};
	struct vbs_netekf *ekf;
	struct vbs_preproc *ppm, *ppa;
	static long arc_m[MAXSAT], arc_a[MAXSAT];
	struct vbs_netekf_opt opt;
	struct vbs_netekf_stats st;
	const struct vbs_sigset *ss;
	gtime_t t0 = {0, 0.0};
	double mxyz[3], axyz[3], sd[MAXOBS], lo, hi;
	const char *label;
	const obsd_t *ep, *ap;
	int i = 0, j = 0, n, na, k;

	if (argc < 10) {
		fprintf(stderr, "usage: %s <master.rnx> <mx> <my> <mz> <aux.rnx> <ax> <ay> <az> "
		    "<nav.rnx> [label]\n", argv[0]);
		return 1;
	}
	mxyz[0] = atof(argv[2]); mxyz[1] = atof(argv[3]); mxyz[2] = atof(argv[4]);
	axyz[0] = atof(argv[6]); axyz[1] = atof(argv[7]); axyz[2] = atof(argv[8]);
	label = argc > 10 ? argv[10] : "";

	memset(&mobs, 0, sizeof mobs);
	memset(&aobs, 0, sizeof aobs);
	memset(&nav, 0, sizeof nav);
	if (readrnxt(argv[1], 1, t0, t0, 0.0, "", &mobs, &nav, NULL) < 0 || mobs.n <= 0) {
		fprintf(stderr, "cannot read master %s\n", argv[1]);
		return 1;
	}
	if (readrnxt(argv[5], 2, t0, t0, 0.0, "", &aobs, &nav, NULL) < 0 || aobs.n <= 0) {
		fprintf(stderr, "cannot read aux %s\n", argv[5]);
		return 1;
	}
	if (readrnxt(argv[9], 1, t0, t0, 0.0, "", NULL, &nav, NULL) < 0) {
		fprintf(stderr, "cannot read nav %s\n", argv[9]);
		return 1;
	}

	vbs_netekf_default_opt(&opt);
	if ((ekf = vbs_netekf_new(mxyz, axyz, &opt)) == NULL) {
		fprintf(stderr, "cannot create filter\n");
		return 1;
	}
	ppm = vbs_preproc_new(30.0, 3.0);
	ppa = vbs_preproc_new(30.0, 3.0);

	while ((n = epoch_slice(&mobs, &i, &ep)) > 0) {
		int jj = j;

		while ((na = epoch_slice(&aobs, &jj, &ap)) > 0 &&
		    timediff(ap[0].time, ep[0].time) < -1e-3)
			j = jj;
		if (na <= 0)
			break;
		if (fabs(timediff(ap[0].time, ep[0].time)) > 1e-3)
			continue;
		j = jj;
		if (n > MAXOBS)
			n = MAXOBS;
		satposs(ep[0].time, ep, n, &nav, EPHOPT_BRDC, rs, dts, var, svh);
		vbs_satpos_by_sat(ep, n, rs, rs_sat);
		vbs_preproc_epoch(ppm, ep, n, &nav, NULL, 0);
		vbs_preproc_epoch(ppa, ap, na, &nav, NULL, 0);
		vbs_preproc_arc_ids(ppm, arc_m);
		vbs_preproc_arc_ids(ppa, arc_a);
		vbs_netekf_epoch(ekf, ep[0].time, ep, n, rs_sat, dts, ap, na, &nav, arc_m, arc_a);
		vbs_ambres_resolve(ekf, NULL, opt.sys_mask, &ar, NULL, 0, NULL);
	}
	vbs_netekf_get_stats(ekf, &st);
	ss = vbs_netekf_sigset(ekf);

	printf("=== %s: %ld epochs, %ld observations, %d satellites, %d band(s) max\n",
	    label, st.epochs, st.obs_used, st.nsat_active, st.nband_max);
	printf("    signals:");
	for (k = 0; k < ss->n; k++)
		printf(" %s%s", vbs_signal_str(ss->sig[k].sys, ss->sig[k].band),
		    ss->sig[k].attr_mismatch ? "*" : "");
	printf("\n    ZTD difference %.4f +- %.4f m, phase rejected %ld of %ld\n",
	    st.ztd, st.ztd_sigma, st.obs_rejected[VBS_NET_REJ_PHASE], st.obs_used);
	for (k = 0; k < VBS_AMB_NSTEP; k++)
		printf("    %-11s %7ld attempt(s) %8ld fixed %7ld ratio-fail %6ld partial\n",
		    SN[k], ar.attempts[k], ar.fixed[k], ar.ratio_fail[k], ar.partial_used[k]);
	for (k = 0; k < 3; k++) {
		int m = vbs_netekf_dd_ion_sigma(ekf, SL[k], sd, MAXOBS);

		if (m <= 0)
			continue;
		span(sd, m, &lo, &hi);
		printf("    DD ionosphere sigma %s: %6.1f to %6.1f mm over %d pair(s)\n",
		    SLN[k], lo * 1000.0, hi * 1000.0, m);
	}

	{
		struct vbs_netekf_satinfo si[VBS_NET_MAXSAT];
		int m = vbs_netekf_satinfo(ekf, si, VBS_NET_MAXSAT);

		printf("    per-satellite at the last epoch (why a correction is or is not usable):\n");
		printf("      %-4s %-4s %5s %8s %7s %10s %12s\n",
		    "sat", "sys", "bands", "epochs", "elev", "sd_ion_mm", "amb_sigma_cyc");
		for (k = 0; k < m; k++) {
			char id[8];

			satno2id(si[k].sat, id);
			printf("      %-4s %-4s %5d %8ld %7.1f %10.1f %12.3f\n",
			    id, SLN[si[k].sys == SYS_GPS ? 0 : si[k].sys == SYS_GAL ? 1 : 2],
			    si[k].nband, si[k].epochs_tracked, si[k].el_deg,
			    si[k].sd_ion * 1000.0, si[k].amb_sigma_max);
		}
	}

	vbs_netekf_free(ekf);
	vbs_preproc_free(ppm);
	vbs_preproc_free(ppa);
	free(mobs.data);
	free(aobs.data);
	freenav(&nav, 0xFF);
	return 0;
}

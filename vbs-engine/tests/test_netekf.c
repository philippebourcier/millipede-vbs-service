/*
 * vbs_netekf tests -- the network EKF, spec §4.3/§4.4.
 *
 * Two layers:
 *
 *   1. closed loop. Take a real master station-day, synthesize a second station at a known
 *      offset carrying *known* corrections (via vbs_synth_epoch, which test_zero_baseline
 *      already validates), and require the filter to recover those corrections. This is the
 *      test that pins the observation model down: a sign error on the ionosphere, a wrong
 *      mu_i, or a mapping function applied to the wrong station all survive a "does it
 *      converge" check and all fail here.
 *
 *   2. a real baseline. Run STA_M->STA_T and report what the filter does with data nobody
 *      constructed. There is no ground truth, so the checks are the ones that can be made
 *      without it: the states must converge, the clock state must absorb the clock, and the
 *      FDE gate must not be throwing away most of the data.
 */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_ambres.h"
#include "vbs/vbs_netekf.h"
#include "vbs/vbs_preproc.h"
#include "vbs/vbs_synth.h"

#define DATA_DIR	"testdata"
#define MASTER_FILE	DATA_DIR "/sta_m_2026210_30s.rnx"
#define TARGET_FILE	DATA_DIR "/sta_t_2026210_30s.rnx"
#define NAV_FILE	DATA_DIR "/nav__2026__210__BRDC00WRD_S_20262100000_01D_MN.rnx"

/* STA_M and STA_T, the §3.2 acceptance baseline (29.1 km). */
static const double MASTER_XYZ[3] = {4582081.2, 55183.6, 4421736.6};
static const double TARGET_XYZ[3] = {4603034.1, 55436.0, 4400265.6};

static int failures;

static void check(int ok, const char *fmt, ...) {
	va_list ap;
	char buf[200];

	va_start(ap, fmt);
	vsnprintf(buf, sizeof buf, fmt, ap);
	va_end(ap);
	printf("  %-72s %s\n", buf, ok ? "ok" : "FAIL");
	if (!ok)
		failures++;
}

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

static int load(const char *file, int rcv, obs_t *obs, nav_t *nav) {
	gtime_t t0 = {0, 0.0};

	if (readrnxt(file, rcv, t0, t0, 0.0, "", obs, nav, NULL) < 0 || obs->n <= 0)
		return 0;
	return 1;
}

/* --- layer 1: closed loop against known corrections ------------------------------------ */

/*
 * The corrections injected into the synthetic station. dI is a per-satellite constant so the
 * filter has something satellite-specific to separate; dT is built from the same a priori
 * model and mapping function the filter uses, plus a known residual ZTD difference, so that
 * the truth is exactly representable by the filter's state vector. Injecting a dT the state
 * vector *cannot* represent would test the wrong thing -- that is what layer 2 is for.
 */
#define TRUE_ZTD	0.042		/* m of residual ZTD difference */

static double true_dI(int sat) {
	return 0.03 + 0.002 * (sat % 17);	/* 3.0 to 6.2 cm, satellite-dependent */
}

static void test_closed_loop(void) {
	static obs_t obs;
	static nav_t nav;
	static double rs[MAXOBS * 6], dts[MAXOBS * 2], var[MAXOBS], rs_sat[MAXSAT * 6];
	static int svh[MAXOBS];
	struct vbs_netekf *ekf;
	struct vbs_netekf_opt opt;
	struct vbs_netekf_stats st;
	struct vbs_corr inj[MAXOBS], got[MAXOBS];
	static double inj_model[MAXOBS];
	struct vbs_preproc *ppm, *ppa;
	static long arc_m[MAXSAT], arc_a[MAXSAT];
	static struct vbs_ambres_stats ar;
	obsd_t synth[MAXOBS];
	double mpos[3], apos[3];
	double worst_dI = 0.0, worst_dT = 0.0;
	int i = 0, k, n, nepoch = 0, ncheck = 0;
	const obsd_t *ep;

	printf("closed loop: synthesize a station with known corrections, recover them\n");

	memset(&obs, 0, sizeof obs);
	memset(&nav, 0, sizeof nav);
	if (!load(MASTER_FILE, 1, &obs, &nav)) {
		printf("  SKIPPED: cannot read %s\n", MASTER_FILE);
		return;
	}
	if (readrnxt(NAV_FILE, 1, (gtime_t){0, 0.0}, (gtime_t){0, 0.0}, 0.0, "",
	    NULL, &nav, NULL) < 0) {
		printf("  SKIPPED: cannot read %s\n", NAV_FILE);
		return;
	}
	ecef2pos(MASTER_XYZ, mpos);
	ecef2pos(TARGET_XYZ, apos);

	vbs_netekf_default_opt(&opt);
	if ((ekf = vbs_netekf_new(MASTER_XYZ, TARGET_XYZ, &opt)) == NULL) {
		check(0, "filter created");
		return;
	}
	ppm = vbs_preproc_new(30.0, 3.0);
	ppa = vbs_preproc_new(30.0, 3.0);

	while ((n = epoch_slice(&obs, &i, &ep)) > 0) {
		int ninj = 0, nsynth, ngot;

		if (n > MAXOBS)
			n = MAXOBS;
		satposs(ep[0].time, ep, n, &nav, EPHOPT_BRDC, rs, dts, var, svh);
		vbs_satpos_by_sat(ep, n, rs, rs_sat);

		/* Build the truth for this epoch, in exactly the filter's parameterisation. */
		for (k = 0; k < n; k++) {
			double e[3], azel_m[2], azel_a[2], mapfw_m, mapfw_a, mapfw;

			if (geodist(rs + k * 6, MASTER_XYZ, e) <= 0.0)
				continue;
			satazel(mpos, e, azel_m);
			if (azel_m[1] < opt.elev_mask_deg * D2R)
				continue;
			if (geodist(rs + k * 6, TARGET_XYZ, e) <= 0.0)
				continue;
			satazel(apos, e, azel_a);

			tropmapf(ep[0].time, mpos, azel_m, &mapfw_m);
			tropmapf(ep[0].time, apos, azel_a, &mapfw_a);
			mapfw = 0.5 * (mapfw_m + mapfw_a);

			inj[ninj].sat = ep[k].sat;
			inj[ninj].dI = true_dI(ep[k].sat);
			inj_model[ninj] = tropmodel(ep[0].time, apos, azel_a, opt.humi)
					- tropmodel(ep[0].time, mpos, azel_m, opt.humi);
			inj[ninj].dT = inj_model[ninj] + mapfw * TRUE_ZTD;
			inj[ninj].valid = 1;
			ninj++;
		}
		if (ninj == 0)
			continue;

		nsynth = vbs_synth_epoch(ep, n, rs, MASTER_XYZ, TARGET_XYZ,
		    inj, ninj, &nav, synth);
		if (nsynth <= 0)
			continue;

		/*
		 * Reverse the auxiliary's record order.
		 *
		 * `vbs_netekf_epoch` indexes `rs` by the *master's* record order, and looks the
		 * auxiliary's satellites up by number. D42 found a caller that violated that --
		 * it handed an auxiliary's records alongside the master's satellite positions,
		 * pairing each satellite with another's position, and produced misclosures of
		 * ~293 000 cycles. The master↔auxiliary path is correct by construction, but
		 * "by construction" is what the other caller looked like too, so it is asserted
		 * here instead: with the auxiliary's records in a different order from the
		 * master's, the injected corrections must still come back exactly.
		 */
		for (k = 0; k < nsynth / 2; k++) {
			obsd_t tmp = synth[k];

			synth[k] = synth[nsynth - 1 - k];
			synth[nsynth - 1 - k] = tmp;
		}

		vbs_preproc_epoch(ppm, ep, n, &nav, NULL, 0);
		vbs_preproc_epoch(ppa, synth, nsynth, &nav, NULL, 0);
		vbs_preproc_arc_ids(ppm, arc_m);
		vbs_preproc_arc_ids(ppa, arc_a);
		vbs_netekf_epoch(ekf, ep[0].time, ep, n, rs_sat, dts, synth, nsynth, &nav,
		    arc_m, arc_a);
		vbs_ambres_resolve(ekf, NULL, opt.sys_mask, &ar, NULL, 0, NULL);
		nepoch++;

		/* Let the filter converge before judging it. */
		if (nepoch < 120)
			continue;
		ngot = vbs_netekf_corrections(ekf, got, MAXOBS, 0.05);
		for (k = 0; k < ngot; k++) {
			int j;

			if (!got[k].valid)
				continue;
			for (j = 0; j < ninj; j++) {
				double eI, eT;

				if (inj[j].sat != got[k].sat)
					continue;
				eI = fabs(got[k].dI - inj[j].dI);
				/*
				 * vbs_netekf reports the *residual* non-dispersive term only; the
				 * a priori model difference is §4.5's job, computed at the VBS
				 * rather than interpolated (see the comment on dT in
				 * vbs_netekf_corrections). The injected correction is the full
				 * thing, as a real correction is, so the model part is removed
				 * here before comparing.
				 */
				eT = fabs(got[k].dT - (inj[j].dT - inj_model[j]));
				if (eI > worst_dI)
					worst_dI = eI;
				if (eT > worst_dT)
					worst_dT = eT;
				ncheck++;
			}
		}
	}
	vbs_netekf_get_stats(ekf, &st);

	printf("  %ld epochs, %ld observations used, %d corrections compared\n",
	    st.epochs, st.obs_used, ncheck);
	printf("  recovered ZTD difference %.4f m (injected %.4f, sigma %.4f)\n",
	    st.ztd, TRUE_ZTD, st.ztd_sigma);
	printf("  worst dI error %.2f mm, worst dT error %.2f mm\n",
	    worst_dI * 1000.0, worst_dT * 1000.0);
	{
		struct vbs_corr c[MAXOBS];
		double sI[MAXOBS], sT[MAXOBS], lo = 1e9, hi = 0.0;
		int nc = vbs_netekf_corrections(ekf, c, MAXOBS, 0.49), q;

		vbs_netekf_correction_sigmas(ekf, c, nc, sI, sT);
		for (q = 0; q < nc; q++) {
			if (sI[q] < lo) lo = sI[q];
			if (sI[q] > hi) hi = sI[q];
		}
		printf("  float dI sigma across %d satellites: %.1f to %.1f mm\n",
		    nc, lo * 1000.0, hi * 1000.0);
		/*
		 * This is the number that says why §4.3's ambiguity resolution is not an
		 * optimisation. Without fixed ambiguities the *level* of the slant ionosphere
		 * is pinned only by code -- metres of noise -- while phase pins only how it
		 * changes. Decimetre sigmas are the honest float answer, and they are an order
		 * of magnitude above the 1.4-12.5 cm the corrections have to beat
		 * (docs/vbs/03-error-budget.md). Fixing is what turns this into a mm-level
		 * quantity.
		 */
		check(lo < 0.05, "at least one satellite's float dI beat 50 mm (best %.1f mm)",
		    lo * 1000.0);
	}

	check(st.epochs > 100, "epochs processed (%ld)", st.epochs);
	check(ncheck > 1000, "corrections compared against truth (%d)", ncheck);
	check(fabs(st.ztd - TRUE_ZTD) < 0.005,
	    "ZTD difference recovered within 5 mm (%.2f mm off)",
	    fabs(st.ztd - TRUE_ZTD) * 1000.0);
	check(worst_dI < 0.01, "every slant ionosphere within 10 mm (worst %.2f mm)",
	    worst_dI * 1000.0);
	check(worst_dT < 0.01, "every non-dispersive term within 10 mm (worst %.2f mm)",
	    worst_dT * 1000.0);

	/*
	 * A validity threshold at or above the prior is not a gate: it would pass a satellite
	 * whose ionosphere state is still untouched. The filter refuses that outright, and this
	 * is the check that it does -- because the failure mode it prevents is silent, and was
	 * measured here at 61.8 mm of error reported as a valid correction.
	 */
	check(vbs_netekf_corrections(ekf, got, MAXOBS, sqrt(opt.var_ion0)) == 0,
	    "a validity threshold at the prior is refused, not honoured");

	vbs_netekf_free(ekf);
	vbs_preproc_free(ppm);
	vbs_preproc_free(ppa);
	free(obs.data);
	freenav(&nav, 0xFF);
}

/* --- layer 2: a real baseline ---------------------------------------------------------- */

static void test_real_baseline(void) {
	static obs_t mobs, aobs;
	static nav_t nav;
	static double rs[MAXOBS * 6], dts[MAXOBS * 2], var[MAXOBS], rs_sat[MAXSAT * 6];
	static int svh[MAXOBS];
	struct vbs_netekf *ekf;
	struct vbs_netekf_opt opt;
	struct vbs_netekf_stats st;
	static struct vbs_ambres_stats ar;
	const struct vbs_sigset *ss;
	struct vbs_preproc *ppm, *ppa;
	static long arc_m[MAXSAT], arc_a[MAXSAT];
	int i = 0, j = 0, n, na, k;
	const obsd_t *ep, *ap;

	printf("real baseline: STA_M -> STA_T\n");

	memset(&mobs, 0, sizeof mobs);
	memset(&aobs, 0, sizeof aobs);
	memset(&nav, 0, sizeof nav);
	if (!load(MASTER_FILE, 1, &mobs, &nav) || !load(TARGET_FILE, 2, &aobs, &nav)) {
		printf("  SKIPPED: cannot read the baseline's observation files\n");
		return;
	}
	if (readrnxt(NAV_FILE, 1, (gtime_t){0, 0.0}, (gtime_t){0, 0.0}, 0.0, "",
	    NULL, &nav, NULL) < 0) {
		printf("  SKIPPED: cannot read %s\n", NAV_FILE);
		return;
	}

	vbs_netekf_default_opt(&opt);
	if ((ekf = vbs_netekf_new(MASTER_XYZ, TARGET_XYZ, &opt)) == NULL) {
		check(0, "filter created");
		return;
	}
	ppm = vbs_preproc_new(30.0, 3.0);
	ppa = vbs_preproc_new(30.0, 3.0);

	while ((n = epoch_slice(&mobs, &i, &ep)) > 0) {
		int jj = j;

		/* Walk the auxiliary to the master's epoch. */
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
		if (getenv("VBS_NO_AR") == NULL)
			vbs_ambres_resolve(ekf, NULL, opt.sys_mask, &ar, NULL, 0, NULL);
	}
	vbs_netekf_get_stats(ekf, &st);
	ss = vbs_netekf_sigset(ekf);

	printf("  %ld epochs, %ld observations used, %d satellites active at the end\n",
	    st.epochs, st.obs_used, st.nsat_active);
	printf("  common signals (%d):", ss->n);
	for (k = 0; k < ss->n; k++)
		printf(" %s%s", vbs_signal_str(ss->sig[k].sys, ss->sig[k].band),
		    ss->sig[k].attr_mismatch ? "*" : "");
	printf("   (* = tracking attribute differs between ends, D4)\n");
	printf("  rejected: %ld elevation, %ld code gate, %ld phase gate, %ld geometry, %ld no slot\n",
	    st.obs_rejected[VBS_NET_REJ_ELEV], st.obs_rejected[VBS_NET_REJ_CODE],
	    st.obs_rejected[VBS_NET_REJ_PHASE], st.obs_rejected[VBS_NET_REJ_GEOM],
	    st.obs_rejected[VBS_NET_REJ_NOSLOT]);
	printf("  ambiguity resets from slips: %ld master, %ld auxiliary\n",
	    st.slips_master, st.slips_aux);
	printf("  residual ZTD difference %.4f +- %.4f m\n", st.ztd, st.ztd_sigma);
	printf("  innovation rms: code %.3f m (normalised %.2f), phase %.4f m (normalised %.2f)\n",
	    st.innov_rms_code, st.innov_norm_code, st.innov_rms_phase, st.innov_norm_phase);

	{
		static const char *SN[VBS_AMB_NSTEP] = {"extra-wide", "wide", "narrow"};
		struct vbs_corr c[MAXOBS];
		double sI[MAXOBS], sT[MAXOBS], lo = 1e9, hi = 0.0;
		int nc = vbs_netekf_corrections(ekf, c, MAXOBS, 0.49), q;

		printf("  ambiguity cascade:\n");
		for (q = 0; q < VBS_AMB_NSTEP; q++)
			printf("    %-11s %7ld attempt(s) %8ld fixed %7ld ratio-fail "
			    "%6ld partial\n", SN[q], ar.attempts[q], ar.fixed[q],
			    ar.ratio_fail[q], ar.partial_used[q]);
		vbs_netekf_correction_sigmas(ekf, c, nc, sI, sT);
		for (q = 0; q < nc; q++) {
			if (sI[q] < lo) lo = sI[q];
			if (sI[q] > hi) hi = sI[q];
		}
		printf("  single-difference dI sigma, %d satellites: %.1f to %.1f mm"
		    " (dominated by the common term)\n", nc, lo * 1000.0, hi * 1000.0);
		{
			static const int SL[3] = {SYS_GPS, SYS_GAL, SYS_CMP};
			static const char *SLN[3] = {"GPS", "GAL", "BDS"};
			double sd[MAXOBS];
			int y, m2, w;

			for (y = 0; y < 3; y++) {
				if ((m2 = vbs_netekf_dd_ion_sigma(ekf, SL[y], sd, MAXOBS)) <= 0)
					continue;
				lo = 1e9; hi = 0.0;
				for (w = 0; w < m2; w++) {
					if (sd[w] < lo) lo = sd[w];
					if (sd[w] > hi) hi = sd[w];
				}
				printf("  DOUBLE-difference dI sigma %s, %d pair(s): "
				    "%.1f to %.1f mm\n", SLN[y], m2, lo * 1000.0, hi * 1000.0);
			}
		}
	}

	check(st.epochs > 1000, "epochs processed (%ld)", st.epochs);
	check(st.obs_used > 10000, "observations used (%ld)", st.obs_used);
	check(st.obs_rejected[VBS_NET_REJ_NOSLOT] == 0,
	    "state vector never ran out of satellite slots");
	check(st.ztd_sigma < 0.05, "ZTD difference converged (sigma %.1f mm)",
	    st.ztd_sigma * 1000.0);
	check(fabs(st.ztd) < 0.5, "ZTD difference is physically plausible (%.3f m)", st.ztd);
	/*
	 * The FDE gate has to reject the genuinely bad and keep almost everything else. A gate
	 * throwing away a large fraction is not robustness, it is a mis-tuned filter quietly
	 * discarding the network.
	 */
	check(st.obs_rejected[VBS_NET_REJ_PHASE] < st.obs_used / 20,
	    "FDE gate rejected under 5%% of phase observations (%ld of %ld)",
	    st.obs_rejected[VBS_NET_REJ_PHASE], st.obs_used);

	vbs_netekf_free(ekf);
	vbs_preproc_free(ppm);
	vbs_preproc_free(ppa);
	free(mobs.data);
	free(aobs.data);
	freenav(&nav, 0xFF);
}

/*
 * Network loop closure (D33 item 2). The identity N_AB = N_MB - N_MA is exact, so the test can
 * be exact too: build three consistent fix sets and require zero misclosures, then move one
 * integer by one cycle and require it to be caught. A wrong fix is always at least a whole
 * cycle, so there is no borderline case to tune a tolerance against.
 */
static void test_loop_closure(void) {
	struct vbs_amb_fix ma[4], mb[4], ab[4], bad[8];
	struct vbs_amb_loop_stats st;
	int i, n;

	printf("network loop closure (D33 item 2)\n");

	for (i = 0; i < 4; i++) {
		ma[i].sat = 10 + i; ma[i].refsat = 3; ma[i].sys = SYS_GPS;
		ma[i].comb_id = 2;  ma[i].step = VBS_AMB_STEP_WL; ma[i].integer = 100.0 + i;
		mb[i] = ma[i];      mb[i].integer = 250.0 - i;
		ab[i] = ma[i];      ab[i].integer = mb[i].integer - ma[i].integer;
	}

	memset(&st, 0, sizeof st);
	n = vbs_ambres_loop_check(ma, 4, mb, 4, ab, 4, 0.25, &st, bad, 8);
	printf("  consistent triangle: %ld checked, %ld misclosed, worst %.3f cycles\n",
	    st.checked, st.misclosed, st.worst_misclosure);
	check(st.checked == 4, "all four ambiguities checked (%ld)", st.checked);
	check(n == 0 && st.misclosed == 0, "a consistent triangle closes exactly");
	check(st.worst_misclosure < 1e-12,
	    "misclosure is zero, not merely small (%.3g)", st.worst_misclosure);

	/* One wrong fix, by exactly one cycle -- the smallest error a wrong fix can be. */
	ab[2].integer += 1.0;
	memset(&st, 0, sizeof st);
	n = vbs_ambres_loop_check(ma, 4, mb, 4, ab, 4, 0.25, &st, bad, 8);
	printf("  one ambiguity off by a cycle: %ld misclosed, worst %.3f\n",
	    st.misclosed, st.worst_misclosure);
	check(n == 1 && st.misclosed == 1, "a single wrong fix is caught (%d)", n);
	check(n == 1 && bad[0].sat == 12,
	    "and the offending satellite pair is reported (%d)", n ? bad[0].sat : -1);

	/*
	 * A reference-satellite mismatch is re-referenced, not skipped: N(s,r) = N(s,k) - N(r,k)
	 * is exact integer arithmetic, and each filter picks its own reference, so skipping
	 * would discard perfectly well-defined checks whenever they disagreed.
	 *
	 * Here AB is expressed against satellite 10 instead of 3. Re-referencing needs AB's fix
	 * for satellite 3 as well, so that is added; the identity must then still close.
	 */
	ab[2].integer -= 1.0;			/* undo the injected error */
	{
		struct vbs_amb_fix ab2[5];
		double n10 = ab[0].integer;	/* AB's value for satellite 10 against ref 3 */

		for (i = 0; i < 4; i++) {
			ab2[i] = ab[i];
			ab2[i].refsat = 10;
			ab2[i].integer = ab[i].integer - n10;	/* N(s,10) = N(s,3) - N(10,3) */
		}
		/* AB's own entry for the old reference, against the new one. */
		ab2[4] = ab[0];
		ab2[4].sat = 3;
		ab2[4].refsat = 10;
		ab2[4].integer = -n10;

		memset(&st, 0, sizeof st);
		n = vbs_ambres_loop_check(ma, 4, mb, 4, ab2, 5, 0.25, &st, bad, 8);
		printf("  AB re-referenced to another satellite: %ld checked, %ld misclosed, "
		    "%ld skipped\n", st.checked, st.misclosed, st.skipped_ref);
		check(st.checked >= 3 && n == 0 && st.misclosed == 0,
		    "re-referencing keeps the triangle checkable and it still closes (%ld checked)",
		    st.checked);
	}

	/* An ambiguity fixed on only two of the three baselines says nothing either way. */
	memset(&st, 0, sizeof st);
	n = vbs_ambres_loop_check(ma, 4, mb, 4, ab, 0, 0.25, &st, bad, 8);
	check(n == 0 && st.checked == 0,
	    "an ambiguity not fixed on all three baselines is not judged");
}

int main(void) {
	test_closed_loop();
	test_loop_closure();
	test_real_baseline();
	printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
	return failures ? 1 : 0;
}

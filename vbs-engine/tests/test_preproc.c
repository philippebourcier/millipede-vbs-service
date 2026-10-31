/*
 * vbs_preproc tests -- cycle-slip detection, spec S4.2.
 *
 * Two layers:
 *   1. synthetic: inject a known slip into a clean arc and check each detector fires where it
 *      should, including the case each one is individually blind to;
 *   2. real: run a full station-day from the validation dataset and report what each detector
 *      found. The number that matters there is `slip_missed_by_lli` -- slips the combinations
 *      caught and the receiver's own flag did not. docs/vbs/03-error-budget.md measured the
 *      consequence of missing those (GPS L2 residual RMS 25.5 m instead of 7.3 cm); this
 *      counts them directly.
 */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_preproc.h"

#define DATA_DIR	"testdata"
#define QUIET_FILE	DATA_DIR "/sta_m_2026210_30s.rnx"
#define STORM_FILE	DATA_DIR "/sta_m_2024284_30s.rnx"
#define NAV_FILE	DATA_DIR "/nav__2026__210__BRDC00WRD_S_20262100000_01D_MN.rnx"

static int failures;

static void check(int ok, const char *fmt, ...) {
	va_list ap;
	char buf[160];

	va_start(ap, fmt);
	vsnprintf(buf, sizeof buf, fmt, ap);
	va_end(ap);
	printf("  %-66s %s\n", buf, ok ? "ok" : "FAIL");
	if (!ok)
		failures++;
}

/* --- layer 1: synthetic slips --------------------------------------------------------- */

/*
 * Build a clean GPS L1/L2 arc for one satellite: ranges that drift smoothly, no slips.
 * Values need not be physically exact -- the detectors only look at continuity.
 */
static void make_epoch(obsd_t *o, gtime_t t, int sat, double range, double dn1, double dn2) {
	double f1 = FREQ1, f2 = FREQ2;

	memset(o, 0, sizeof *o);
	o->time = t;
	o->sat = sat;
	o->code[0] = CODE_L1C;
	o->code[1] = CODE_L2W;
	o->P[0] = range;
	o->P[1] = range;
	o->L[0] = range / (CLIGHT / f1) + dn1;
	o->L[1] = range / (CLIGHT / f2) + dn2;
	o->SNR[0] = o->SNR[1] = 45 * 1000 / 4;
}

static unsigned run_synthetic(double dn1_at_10, double dn2_at_10) {
	struct vbs_preproc *pp = vbs_preproc_new(30.0, 3.0);
	static nav_t nav;
	struct vbs_slip slip[8];
	gtime_t t = epoch2time((double[]){2026, 7, 29, 0, 0, 0});
	unsigned reason_at_10 = 0;
	int k, ns;

	memset(&nav, 0, sizeof nav);
	for (k = 0; k < 20; k++) {
		obsd_t o;
		double dn1 = (k >= 10) ? dn1_at_10 : 0.0;
		double dn2 = (k >= 10) ? dn2_at_10 : 0.0;

		/* 20 000 km, drifting 300 m per epoch: smooth, slip-free geometry. */
		make_epoch(&o, timeadd(t, k * 30.0), satno(SYS_GPS, 1),
		    20000000.0 + k * 300.0, dn1, dn2);
		ns = vbs_preproc_epoch(pp, &o, 1, &nav, slip, 8);
		if (k == 10 && ns > 0)
			reason_at_10 = slip[0].reason;
	}
	vbs_preproc_free(pp);
	return reason_at_10;
}

/*
 * Two band-common jumps 30 s apart (D49): a phase step equal in metres on L1 and L2 leaves
 * the geometry-free combination unchanged, so only the MW test can see it, and the second
 * step lands inside the MW warm-up after the first. Both must be flagged.
 */
static void test_two_step_jump(void) {
	struct vbs_preproc *pp = vbs_preproc_new(30.0, 3.0);
	static nav_t nav;
	struct vbs_slip slip[8];
	gtime_t t = epoch2time((double[]){2026, 7, 29, 0, 0, 0});
	const double lam1 = CLIGHT / FREQ1, lam2 = CLIGHT / FREQ2;
	unsigned r10 = 0, r11 = 0, r12 = 0;
	int k, ns;

	memset(&nav, 0, sizeof nav);
	for (k = 0; k < 20; k++) {
		obsd_t o;
		double jump_m = (k >= 11) ? 55.6 + 155.7 : (k >= 10) ? 55.6 : 0.0;

		make_epoch(&o, timeadd(t, k * 30.0), satno(SYS_GPS, 1),
		    20000000.0 + k * 300.0, jump_m / lam1, jump_m / lam2);
		ns = vbs_preproc_epoch(pp, &o, 1, &nav, slip, 8);
		if (k == 10 && ns > 0) r10 = slip[0].reason;
		if (k == 11 && ns > 0) r11 = slip[0].reason;
		if (k == 12 && ns > 0) r12 = slip[0].reason;
	}
	vbs_preproc_free(pp);
	printf("two band-common jumps 30 s apart (T30's pattern, D49)\n");
	check((r10 & VBS_SLIP_MW) != 0, "first 55.6 m band-common jump flagged by MW");
	check((r11 & VBS_SLIP_MW) != 0, "second 155.7 m jump, 30 s later, flagged by MW during warm-up");
	check(r12 == 0, "no slip on the epoch after, once the phase is steady again");
}

static void test_synthetic(void) {
	unsigned r;

	printf("synthetic slips (GPS L1/L2, clean arc with one injected step)\n");

	r = run_synthetic(0.0, 0.0);
	check(r == 0, "no slip declared on a clean arc");

	/* One cycle on L1 only: geometry-free shifts by a full L1 wavelength. */
	r = run_synthetic(1.0, 0.0);
	check((r & VBS_SLIP_GF) != 0, "1 cycle on L1 alone is caught by geometry-free");

	/* One cycle on L2 only. */
	r = run_synthetic(0.0, 1.0);
	check((r & VBS_SLIP_GF) != 0, "1 cycle on L2 alone is caught by geometry-free");

	/*
	 * The classic geometry-free blind spot: a 9:7 slip moves GF by well under a
	 * millimetre, so GF cannot see it -- but it changes the wide-lane count by 2, which is
	 * exactly what Melbourne-Wubbena is for. This is the case that justifies running both.
	 */
	r = run_synthetic(9.0, 7.0);
	check((r & VBS_SLIP_MW) != 0, "9:7 slip (geometry-free blind spot) is caught by MW");
	printf("      9:7 slip reason bits: GF=%d MW=%d\n",
	    (r & VBS_SLIP_GF) != 0, (r & VBS_SLIP_MW) != 0);

	/* A large simultaneous slip should be caught by at least one detector. */
	r = run_synthetic(50.0, 39.0);
	check((r & (VBS_SLIP_GF | VBS_SLIP_MW)) != 0, "50:39 slip is caught");
}

/* --- layer 2: a real station-day ------------------------------------------------------- */

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

static void test_real_day(const char *file, const char *label) {
	static obs_t obs;
	static nav_t nav;
	struct vbs_preproc *pp;
	struct vbs_preproc_stats st;
	struct vbs_slip slip[MAXOBS];
	gtime_t t0 = {0, 0.0};
	const obsd_t *ep;
	int i = 0, n;

	printf("real station-day: %s\n", label);

	memset(&obs, 0, sizeof obs);
	memset(&nav, 0, sizeof nav);
	if (readrnxt(file, 1, t0, t0, 0.0, "", &obs, &nav, NULL) < 0 || obs.n <= 0) {
		printf("  SKIPPED: cannot read %s\n", file);
		return;
	}
	if (readrnxt(NAV_FILE, 1, t0, t0, 0.0, "", NULL, &nav, NULL) < 0)
		printf("  (no nav file; GLONASS frequencies unavailable, others unaffected)\n");

	if ((pp = vbs_preproc_new(30.0, 3.0)) == NULL) {
		check(0, "preprocessor created");
		return;
	}
	while ((n = epoch_slice(&obs, &i, &ep)) > 0)
		vbs_preproc_epoch(pp, ep, n, &nav, slip, MAXOBS);

	vbs_preproc_get_stats(pp, &st);
	printf("  %ld epochs, %ld satellite-epochs (%ld dual-band, %ld single-band)\n",
	    st.epochs, st.sat_epochs, st.dual_band, st.single_band);
	printf("  slips: %ld total -- LLI %ld, geom-free %ld, MW %ld, gap %ld, arc-start %ld\n",
	    st.slips, st.slip_lli, st.slip_gf, st.slip_mw, st.slip_gap, st.slip_new);
	printf("  caught by GF/MW but NOT flagged by LLI: %ld\n", st.slip_missed_by_lli);
	printf("  flagged by LLI but neither GF nor MW agreed: %ld\n", st.slip_lli_only);

	check(st.epochs > 0 && st.sat_epochs > 0, "observations processed");
	check(st.dual_band + st.single_band == st.sat_epochs,
	    "every satellite-epoch is accounted for as dual- or single-band");
	check(st.single_band > 0,
	    "single-band satellites exist and are reported (BeiDou here, S3.1)");
	check(st.slip_missed_by_lli > 0,
	    "LLI alone would have missed %ld slip(s) -- spec S4.2's point, measured",
	    st.slip_missed_by_lli);

	vbs_preproc_free(pp);
	free(obs.data);
	freenav(&nav, 0xFF);
}

int main(void) {
	test_synthetic();
	test_two_step_jump();
	test_real_day(QUIET_FILE, "STA_M 2026-07-29 (quiet)");
	test_real_day(STORM_FILE, "STA_M 2024-10-10 (storm, Kp 8.667)");
	printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
	return failures ? 1 : 0;
}

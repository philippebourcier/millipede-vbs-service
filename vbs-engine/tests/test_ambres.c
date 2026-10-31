/*
 * vbs_ambres cascade-ordering tests (D74).
 *
 * The cascade's whole mechanism is "longest wavelength first": each step's integers are applied
 * to the filter before the next runs, so a step fixed out of order commits integers on a short
 * wavelength while the long, nearly-free ones are still float. The combination table is keyed
 * on RINEX band digits so one set of rows covers every constellation -- which is what made the
 * order wrong, because the same digits are a different physical combination on each. {2,5} is
 * GPS L2-L5 at 5.86 m and BeiDou B1I-B2a at 0.78 m.
 *
 * These tests go through `vbs_ambres_cascade_order`, which the resolver itself calls to build
 * the sequence it runs. That matters: a test that recomputed the expected order from its own
 * copy of the table would keep passing if the resolver stopped using it.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "rtklib.h"
#include "vbs/vbs_ambres.h"
#include "vbs/vbs_signal.h"

#define MAXSTEP		64

static int failures;

static void check(int ok, const char *what) {
	printf("  %-62s %s\n", what, ok ? "ok" : "FAIL");
	if (!ok)
		failures++;
}

static const char *sysname(int sys) {
	return sys == SYS_GPS ? "GPS" : sys == SYS_GAL ? "GAL" : sys == SYS_CMP ? "BDS" : "?";
}

/* First step of one constellation, or NULL. */
static const struct vbs_amb_step_info *first_of(const struct vbs_amb_step_info *o, int n,
	int sys) {
	int i;

	for (i = 0; i < n; i++)
		if (o[i].sys == sys)
			return &o[i];
	return NULL;
}

static void test_descending(void) {
	struct vbs_amb_step_info o[MAXSTEP];
	int n = vbs_ambres_cascade_order(SYS_GPS | SYS_GAL | SYS_CMP, 1, 1, o, MAXSTEP), i;
	int mono = 1, per_sys = 1, s;

	printf("wavelength ordering, GPS+Galileo+BeiDou\n");
	for (i = 0; i < n; i++)
		printf("    %2d. %s %-8s lambda %8.4f m  step %d\n", i + 1, sysname(o[i].sys),
		    o[i].name, o[i].lam, o[i].step);

	check(n > 0, "the cascade has steps to run");
	for (i = 1; i < n; i++)
		if (o[i].lam > o[i - 1].lam)
			mono = 0;
	check(mono, "every step's wavelength <= the previous step's");

	/* And, separately, within each constellation -- the property the cascade depends on. */
	for (s = 0; s < 3; s++) {
		int sys = s == 0 ? SYS_GPS : s == 1 ? SYS_GAL : SYS_CMP;
		double prev = 1e30;

		for (i = 0; i < n; i++) {
			if (o[i].sys != sys)
				continue;
			if (o[i].lam > prev)
				per_sys = 0;
			prev = o[i].lam;
		}
	}
	check(per_sys, "and within each constellation taken on its own");
}

static void test_beidou_starts_long(void) {
	struct vbs_amb_step_info o[MAXSTEP];
	int n = vbs_ambres_cascade_order(SYS_CMP, 1, 1, o, MAXSTEP), i;
	const struct vbs_amb_step_info *f = first_of(o, n, SYS_CMP);
	int found_short_as_ewl = 0;

	printf("BeiDou: B1C-B1I (20.93 m) is the first step, not the third\n");
	check(f != NULL && fabs(f->lam - 20.9323) < 0.01,
	    "BeiDou's first step is the 20.93 m B1C-B1I wide lane");
	/*
	 * And the 0.78 m B1I-B2a combination must not be bucketed as an extra-wide lane. It sat
	 * in the EWL row of the table, so `max_step` admitted it as an essentially free integer
	 * and the cascade attempted it before everything else.
	 */
	for (i = 0; i < n; i++)
		if (fabs(o[i].lam - 0.7794) < 0.01 && o[i].step == VBS_AMB_STEP_EWL)
			found_short_as_ewl = 1;
	check(!found_short_as_ewl,
	    "BeiDou's 0.78 m combination is classified WL, not extra-wide");
	for (i = 0; i < n; i++)
		if (o[i].step == VBS_AMB_STEP_EWL)
			check(o[i].lam >= 2.0, "every step in the EWL bucket is >= 2 m");
}

static void test_gps_unchanged(void) {
	struct vbs_amb_step_info o[MAXSTEP];
	int n = vbs_ambres_cascade_order(SYS_GPS, 1, 1, o, MAXSTEP);
	const struct vbs_amb_step_info *f = first_of(o, n, SYS_GPS);

	/*
	 * GPS is the constellation the table's row order was written for, so reordering must not
	 * move it: L2-L5 at 5.86 m first, then L1-L2, then L1-L5.
	 */
	printf("GPS is unaffected -- the table's row order was already correct for it\n");
	check(f != NULL && fabs(f->lam - 5.8610) < 0.01, "GPS starts at the 5.86 m L2-L5 EWL");
	check(n >= 3 && fabs(o[1].lam - 0.8619) < 0.01, "then the 0.86 m L1-L2 wide lane");
	check(n >= 3 && fabs(o[2].lam - 0.7514) < 0.01, "then the 0.75 m L1-L5 wide lane");
}

static void test_table_arm_differs(void) {
	struct vbs_amb_step_info lam[MAXSTEP], tab[MAXSTEP];
	int nl = vbs_ambres_cascade_order(SYS_CMP, 1, 1, lam, MAXSTEP);
	int nt = vbs_ambres_cascade_order(SYS_CMP, 0, 1, tab, MAXSTEP);

	/*
	 * `--cascade-order table` exists so the change can be measured against the behaviour it
	 * replaced from one binary. If it silently produced the new order the A/B would compare
	 * a configuration against itself and report no difference -- which is indistinguishable
	 * from a change that does nothing.
	 */
	printf("the pre-D74 arm really is the old behaviour\n");
	check(nl == nt, "both arms offer the same steps");
	check(nt > 0 && fabs(tab[0].lam - 0.7794) < 0.01,
	    "table order starts BeiDou on the 0.78 m combination, as it used to");
	check(nt > 0 && tab[0].step == VBS_AMB_STEP_EWL,
	    "and still mis-buckets it as extra-wide, as it used to");
	check(nl > 0 && fabs(lam[0].lam - tab[0].lam) > 1.0, "the two arms differ");
}

static void test_no_phantom_bands(void) {
	struct vbs_amb_step_info o[MAXSTEP];
	int n = vbs_ambres_cascade_order(SYS_GPS | SYS_GAL | SYS_CMP, 1, 1, o, MAXSTEP), i;
	int gps_band6 = 0;

	/*
	 * A constellation must contribute steps only for bands it actually transmits. GPS shared
	 * its frequency-table case with QZSS, so it answered FREQ6 for band 6; that was inert
	 * until a band-6 combination existed, and then produced four GPS steps on a signal no
	 * GPS satellite sends.
	 */
	printf("no steps on bands a constellation does not transmit\n");
	check(vbs_signal_band_freq(SYS_GPS, 6) == 0.0, "GPS has no band 6 frequency");
	check(vbs_signal_band_freq(SYS_QZS, 6) != 0.0, "QZSS still does");
	check(vbs_signal_band_freq(SYS_GAL, 6) != 0.0, "Galileo E6 is present");
	check(vbs_signal_band_freq(SYS_CMP, 6) != 0.0, "BeiDou B3I is present");
	for (i = 0; i < n; i++)
		if (o[i].sys == SYS_GPS && strstr(o[i].name, "6") != NULL)
			gps_band6 = 1;
	check(!gps_band6, "and the cascade offers GPS no band-6 step");
}

static void test_band6_present(void) {
	struct vbs_amb_step_info on[MAXSTEP], off[MAXSTEP];
	int n_on = vbs_ambres_cascade_order(SYS_GAL | SYS_CMP, 1, 1, on, MAXSTEP);
	int n_off = vbs_ambres_cascade_order(SYS_GAL | SYS_CMP, 1, 0, off, MAXSTEP), i;
	int ewl6 = 0;

	/*
	 * Band 6 was in the filter's state all along and in no combination, so those ambiguities
	 * were estimated and never fixed. The combinations it brings are worth having on their
	 * wavelengths: 4.19 m (E6-E5b) and 4.88 m (B3I-B2b) are extra-wide lanes, and the
	 * band-6 wide lanes at 0.98-1.02 m are longer than any wide lane the table had.
	 */
	printf("band 6 (Galileo E6, BeiDou B3I) contributes real combinations\n");
	check(n_on > n_off, "enabling band 6 adds steps");
	for (i = 0; i < n_on; i++)
		if (strstr(on[i].name, "6") != NULL && on[i].step == VBS_AMB_STEP_EWL &&
		    on[i].lam > 2.0)
			ewl6++;
	check(ewl6 >= 4, "including extra-wide lanes above 2 m on both constellations");
	for (i = 0; i < n_off; i++)
		check(strstr(off[i].name, "6") == NULL, "--no-band6 really removes them");
}

static void test_sys_mask(void) {
	struct vbs_amb_step_info o[MAXSTEP];
	int n = vbs_ambres_cascade_order(SYS_GAL, 1, 1, o, MAXSTEP), i, only_gal = 1;

	printf("a constellation not in the mask contributes no steps\n");
	for (i = 0; i < n; i++)
		if (o[i].sys != SYS_GAL)
			only_gal = 0;
	check(n > 0 && only_gal, "Galileo-only mask yields Galileo steps only");
	/* Galileo has no band 2, so every row needing it must be absent rather than mis-sized. */
	for (i = 0; i < n; i++)
		if (o[i].lam == 0.0)
			check(0, "no step has a zero wavelength");
}

int main(void) {
	test_descending();
	test_beidou_starts_long();
	test_gps_unchanged();
	test_table_arm_differs();
	test_no_phantom_bands();
	test_band6_present();
	test_sys_mask();
	printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
	return failures ? 1 : 0;
}

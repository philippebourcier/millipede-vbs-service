/*
 * vbs_signal tests.
 *
 * Two layers:
 *   1. pure unit checks of the band mapping;
 *   2. the acceptance case named in docs/vbs/02-architecture.md S4 -- run the band matcher
 *      over the real STA_M/STA_T baseline from the validation dataset, the pair that
 *      docs/vbs/01-design-study.md S3.2 showed has NO common GPS L2 and NO common Galileo
 *      under exact-code matching (STA_M reports L2X/L1X/L7X, STA_T reports L2L,L2S/L1C/L7Q).
 *      Band matching must find GPS L1+L2 and Galileo E1+E5b on it anyway, and must flag the
 *      attribute mismatch rather than hide it.
 *
 * The second layer also answers a question the engine depends on and that no amount of
 * reading settles: whether MRTKLIB's RINEX reader actually delivers Galileo E5b observations
 * at all, given its positioning-side frequency-slot convention prefers E5a (D19).
 */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "vbs/vbs_signal.h"

#define DATA_DIR	"testdata"
#define NAV_FILE	DATA_DIR "/nav__2026__210__BRDC00WRD_S_20262100000_01D_MN.rnx"
#define DEFLAU_FILE	DATA_DIR "/sta_m_2026210_30s.rnx"
#define LAPO_FILE	DATA_DIR "/sta_t_2026210_30s.rnx"

static int failures;

static void check(int ok, const char *what) {
	printf("  %-62s %s\n", what, ok ? "ok" : "FAIL");
	if (!ok)
		failures++;
}

static void test_band_mapping(void) {
	printf("band mapping (code -> RINEX band digit)\n");
	check(vbs_signal_band(obs2code("1C")) == VBS_BAND_L1, "1C -> band 1");
	check(vbs_signal_band(obs2code("2X")) == VBS_BAND_L2, "2X -> band 2 (GPS L2C, M+L)");
	check(vbs_signal_band(obs2code("2L")) == VBS_BAND_L2, "2L -> band 2 (GPS L2C, L)");
	check(vbs_signal_band(obs2code("2S")) == VBS_BAND_L2, "2S -> band 2 (GPS L2C, M)");
	check(vbs_signal_band(obs2code("2W")) == VBS_BAND_L2, "2W -> band 2 (GPS L2 P(Y))");
	check(vbs_signal_band(obs2code("7X")) == VBS_BAND_L7, "7X -> band 7 (GAL E5b, I+Q)");
	check(vbs_signal_band(obs2code("7Q")) == VBS_BAND_L7, "7Q -> band 7 (GAL E5b, Q)");
	check(vbs_signal_band(obs2code("5Q")) == VBS_BAND_L5, "5Q -> band 5 (GAL E5a)");
	check(vbs_signal_band(CODE_NONE) == 0, "CODE_NONE -> 0");

	/* The whole point of D4: different attributes, same band. */
	check(vbs_signal_band(obs2code("2X")) == vbs_signal_band(obs2code("2L")),
	    "2X and 2L are the same band");
	check(vbs_signal_band(obs2code("7X")) == vbs_signal_band(obs2code("7Q")),
	    "7X and 7Q are the same band");
	/* ...and different bands stay different. */
	check(vbs_signal_band(obs2code("5Q")) != vbs_signal_band(obs2code("7Q")),
	    "E5a and E5b are NOT the same band");
}

/*
 * Read a RINEX obs file and return how many records belong to its first epoch.
 *
 * nav_t and obs_t are far too large for the stack once MRTKLIB is built with NFREQ=5/
 * NEXOBS=5 (nav_t alone carries MAXSAT-sized arrays), so callers pass heap or static
 * storage.
 */
static int read_first_epoch(const char *file, obs_t *obs, nav_t *nav) {
	gtime_t t0 = {0, 0.0};
	int n;

	memset(obs, 0, sizeof *obs);
	memset(nav, 0, sizeof *nav);

	/* tint=30 s matches the dataset's sampling; we only need the leading epoch. */
	if (readrnxt(file, 1, t0, t0, 30.0, "", obs, nav, NULL) < 0) {
		printf("  cannot read %s\n", file);
		return -1;
	}
	if (obs->n <= 0) {
		printf("  no observations in %s\n", file);
		return -1;
	}
	/*
	 * Records arrive grouped by epoch. Compare with a tolerance rather than for exact
	 * equality: epoch times come from parsed decimal seconds, so identical epochs are not
	 * guaranteed to be bit-identical doubles.
	 */
	for (n = 1; n < obs->n; n++)
		if (fabs(timediff(obs->data[n].time, obs->data[0].time)) > 1e-3)
			break;
	return n;
}

static void dump_bands(const char *name, const obsd_t *o, int n) {
	int i, j, seen[8][10];
	int si, band, sys;

	memset(seen, 0, sizeof seen);
	for (i = 0; i < n; i++) {
		sys = satsys(o[i].sat, NULL);
		si = (sys == SYS_GPS) ? 0 : (sys == SYS_GAL) ? 1 : (sys == SYS_CMP) ? 2 :
		     (sys == SYS_GLO) ? 3 : 4;
		for (j = 0; j < NFREQ + NEXOBS; j++) {
			if (o[i].code[j] == CODE_NONE)
				continue;
			band = vbs_signal_band(o[i].code[j]);
			if (band > 0 && band < 10)
				seen[si][band]++;
		}
	}
	printf("  %s bands present (sat-slots): GPS", name);
	for (band = 1; band < 10; band++)
		if (seen[0][band])
			printf(" L%d:%d", band, seen[0][band]);
	printf(" | GAL");
	for (band = 1; band < 10; band++)
		if (seen[1][band])
			printf(" L%d:%d", band, seen[1][band]);
	printf(" | BDS");
	for (band = 1; band < 10; band++)
		if (seen[2][band])
			printf(" L%d:%d", band, seen[2][band]);
	printf("\n");
}

static void test_real_baseline(void) {
	/* static: obs_t/nav_t are too large for the stack at NFREQ=5/NEXOBS=5. */
	static obs_t oa, ob;
	static nav_t na, nb;
	struct vbs_sigset set;
	int n_a, n_b, i, rc;
	int have_g1 = 0, have_g2 = 0, have_e1 = 0, have_e7 = 0, mismatches = 0;

	printf("real baseline STA_M <-> STA_T (the S3.2 acceptance case)\n");

	if ((n_a = read_first_epoch(DEFLAU_FILE, &oa, &na)) < 0 ||
	    (n_b = read_first_epoch(LAPO_FILE, &ob, &nb)) < 0) {
		printf("  SKIPPED: dataset files not available\n");
		return;
	}
	printf("  first epoch: STA_M %d sats, STA_T %d sats\n", n_a, n_b);
	dump_bands("STA_M", oa.data, n_a);
	dump_bands("STA_T  ", ob.data, n_b);

	rc = vbs_sigset_common(oa.data, n_a, ob.data, n_b, SYS_GPS | SYS_GAL | SYS_CMP, &set);
	check(rc >= 0, "common signal set built");

	printf("  common signals:");
	for (i = 0; i < set.n; i++) {
		printf(" %s(%s/%s)%s", vbs_signal_str(set.sig[i].sys, set.sig[i].band),
		    code2obs(set.sig[i].code_a), code2obs(set.sig[i].code_b),
		    set.sig[i].attr_mismatch ? "*" : "");
		if (set.sig[i].sys == SYS_GPS && set.sig[i].band == VBS_BAND_L1) have_g1 = 1;
		if (set.sig[i].sys == SYS_GPS && set.sig[i].band == VBS_BAND_L2) have_g2 = 1;
		if (set.sig[i].sys == SYS_GAL && set.sig[i].band == VBS_BAND_L1) have_e1 = 1;
		if (set.sig[i].sys == SYS_GAL && set.sig[i].band == VBS_BAND_L7) have_e7 = 1;
		if (set.sig[i].attr_mismatch) mismatches++;
	}
	printf("   (* = same band, different tracking attribute)\n");

	check(have_g1 && have_g2, "GPS L1+L2 common despite L2X vs L2L/L2S labelling");
	check(vbs_sigset_nbands(&set, SYS_GPS) >= 2, "GPS is dual-frequency on this baseline");
	check(have_e1 && have_e7, "Galileo E1+E5b common despite L1X/L7X vs L1C/L7Q labelling");
	check(vbs_sigset_nbands(&set, SYS_GAL) >= 2, "Galileo is dual-frequency on this baseline");
	check(mismatches > 0, "attribute mismatches are recorded, not hidden (D4 open risk)");

	free(oa.data);
	free(ob.data);
	freenav(&na, 0xFF);
	freenav(&nb, 0xFF);
}

/*
 * GLONASS per-satellite frequencies (D90).
 *
 * The point of vbs_signal_sat_freq is that GLONASS satellites do *not* share a frequency, so
 * the test that matters is that different satellites come back different, that each matches the
 * channel-number formula, and that the CDMA systems are untouched. A function that returned a
 * constant would pass a weaker test and break every GLONASS ambiguity.
 */
static void test_glonass_freq(void) {
	nav_t nav = {0};
	int i, nglo = 0, distinct = 0, mismatch = 0;
	double seen[64];

	printf("GLONASS per-satellite frequency (FDMA)\n");
	if (access(NAV_FILE, R_OK) != 0) {
		printf("  SKIPPED: navigation file not available\n");
		return;
	}
	if (readrnxt(NAV_FILE, 1, (gtime_t){0}, (gtime_t){0}, 0.0, "", NULL, &nav, NULL) < 0) {
		check(0, "nav file read for GLONASS ephemerides");
		return;
	}
	check(nav.ng > 0, "nav file carries GLONASS ephemerides");

	for (i = 0; i < nav.ng; i++) {
		int sat = nav.geph[i].sat, k = nav.geph[i].frq, j;
		double f1 = vbs_signal_sat_freq(sat, VBS_BAND_L1, &nav);
		double f2 = vbs_signal_sat_freq(sat, VBS_BAND_L2, &nav);
		double want1 = FREQ1_GLO + k * DFRQ1_GLO;
		double want2 = FREQ2_GLO + k * DFRQ2_GLO;

		if (satsys(sat, NULL) != SYS_GLO)
			continue;
		nglo++;
		if (fabs(f1 - want1) > 1.0 || fabs(f2 - want2) > 1.0)
			mismatch++;
		for (j = 0; j < distinct; j++)
			if (fabs(seen[j] - f1) < 1.0)
				break;
		if (j == distinct && distinct < 64)
			seen[distinct++] = f1;
	}
	printf("  %d GLONASS ephemerides, %d distinct G1 frequencies\n", nglo, distinct);
	check(nglo > 0, "GLONASS satellites present");
	check(mismatch == 0, "every G1/G2 matches 1602+k*0.5625 / 1246+k*0.4375 MHz");
	check(distinct > 4, "satellites really are on different frequencies (not a constant)");

	/* A CDMA system must be unaffected, and must not need nav at all. */
	check(fabs(vbs_signal_sat_freq(satno(SYS_GPS, 1), VBS_BAND_L1, &nav) - FREQ1) < 1.0,
	    "GPS L1 unchanged by the new path");
	check(fabs(vbs_signal_sat_freq(satno(SYS_GPS, 1), VBS_BAND_L1, NULL) - FREQ1) < 1.0,
	    "and needs no nav data");
	/* GLONASS without nav must report "no frequency", not a plausible wrong number. */
	check(vbs_signal_sat_freq(satno(SYS_GLO, 1), VBS_BAND_L1, NULL) == 0.0,
	    "GLONASS without ephemerides returns 0, not a nominal");
	freenav(&nav, 0xFF);
}

int main(void) {
	test_band_mapping();
	test_glonass_freq();
	test_real_baseline();
	printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
	return failures ? 1 : 0;
}

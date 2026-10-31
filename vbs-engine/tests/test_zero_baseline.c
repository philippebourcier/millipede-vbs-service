/*
 * Zero-baseline and RTCM round-trip test -- the Phase 1 acceptance target from
 * docs/vbs/DECISIONS.md D17 item 1:
 *
 *   "Zero-baseline test. Not bit-exact: Rough Range is recomputed from the VBS position, and
 *    LSQ plane-fit interpolation is not exactly zero at the master. With corrections forced
 *    to zero, decoded observables must match the master's within 1 MSM LSB. Test
 *    interpolation separately."
 *
 * Three layers, narrow to wide:
 *   1. synthesis identity: with the VBS *at* the master and corrections zero, the synthesized
 *      observations must equal the master's exactly (this part genuinely is bit-exact -- no
 *      encoding is involved yet, and every added term is identically zero);
 *   2. geometry: drho must be exactly zero at the master, and behave sanely a few km away;
 *   3. the full path: synthesize -> encode MSM7 -> decode with an *independent* decoder ->
 *      compare. Spec S9 Phase 1 asks for "encoding round-trip with an independent decoder";
 *      we encode with MRTKLIB and decode with RTKLIB-EX (tests/helpers/rtklib_rtcm_dump), a
 *      separate codebase in a separate process, so a shared bug cannot cancel itself out.
 *      RTKLIB-EX's convbin would have been simpler but can only emit RINEX, whose ~1 mm
 *      printed resolution is coarser than the 0.14 mm phase LSB this test has to resolve.
 */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "vbs/vbs_rtcm_out.h"
#include "vbs/vbs_signal.h"
#include "vbs/vbs_synth.h"

#define DATA_DIR	"testdata"
#define MASTER_FILE	DATA_DIR "/sta_m_2026210_30s.rnx"
#define NAV_FILE	DATA_DIR "/nav__2026__210__BRDC00WRD_S_20262100000_01D_MN.rnx"

/* STA_M, rounded catalogue ECEF. */
static const double MASTER_XYZ[3] = { 4582081.2, 55183.6, 4421736.6 };

/* Independent decoder helper, built from tests/helpers/ against RTKLIB-EX. */
#ifndef RTKLIB_DUMP_PATH
#define RTKLIB_DUMP_PATH	"tests/helpers/rtklib_rtcm_dump"
#endif

/*
 * MSM7 quantisation, in metres. DF405 (fine pseudorange, extended) has a resolution of
 * 2^-29 ms of range and DF406 (fine phaserange, extended) 2^-31 ms; one millisecond of range
 * is CLIGHT/1000 metres.
 */
#define MS_RANGE	(CLIGHT * 1E-3)
#define MSM7_PR_LSB	(MS_RANGE / 536870912.0)	/* 2^29  -> ~0.56 mm */
#define MSM7_CP_LSB	(MS_RANGE / 2147483648.0)	/* 2^31  -> ~0.14 mm */

#define VBS_STAID	4000			/* from the reserved range, S5.1 */

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

/* Read a RINEX obs (and optionally nav) file; returns records in the first epoch. */
static int read_first_epoch(const char *obsfile, const char *navfile,
	obs_t *obs, nav_t *nav) {
	gtime_t t0 = {0, 0.0};
	int n;

	memset(obs, 0, sizeof *obs);
	memset(nav, 0, sizeof *nav);

	if (readrnxt(obsfile, 1, t0, t0, 30.0, "", obs, nav, NULL) < 0)
		return -1;
	if (navfile != NULL && readrnxt(navfile, 1, t0, t0, 0.0, "", NULL, nav, NULL) < 0)
		return -1;
	if (obs->n <= 0)
		return -1;
	for (n = 1; n < obs->n; n++)
		if (fabs(timediff(obs->data[n].time, obs->data[0].time)) > 1e-3)
			break;
	return n;
}

/* --- layer 1: synthesis identity ------------------------------------------------------ */

/*
 * Systems the synthesis can translate: GPS, Galileo, BeiDou and now GLONASS (D90).
 *
 * GLONASS joined once vbs_signal_sat_freq() existed to give a satellite's own frequency -- the
 * synthesis translates phase at that frequency, so FDMA was the only thing keeping it out.
 * Whether GLONASS is *used* is a cluster decision (matching receivers, for inter-channel bias),
 * but the synthesis itself is not the place that decides, and this test is about the synthesis.
 */
static int is_supported(int sat) {
	int sys = satsys(sat, NULL);

	return sys == SYS_GPS || sys == SYS_GLO || sys == SYS_GAL || sys == SYS_CMP;
}

static void test_identity(const obsd_t *master, int n, const nav_t *nav, obsd_t *vbs) {
	int i, j, k, nout, diffs = 0, compared = 0, nsupported = 0;

	printf("synthesis identity (VBS at the master, corrections zero)\n");

	for (i = 0; i < n; i++)
		if (is_supported(master[i].sat))
			nsupported++;

	nout = vbs_synth_epoch(master, n, NULL, MASTER_XYZ, MASTER_XYZ, NULL, 0, nav, vbs);

	/*
	 * The synthesized set is the translatable subset of the master's: a satellite whose
	 * frequency cannot be determined is not emitted at all, which for GLONASS means one
	 * without ephemerides. Asserting "nout == n" would assert that the engine emits
	 * satellites it cannot translate.
	 */
	check(nout == nsupported, "every supported satellite synthesized (%d of %d; %d "
	    "untranslatable)", nout, nsupported, n - nsupported);

	/* Match by satellite number: the output is compacted, so indices do not line up. */
	for (i = 0; i < nout; i++) {
		for (k = 0; k < n; k++)
			if (master[k].sat == vbs[i].sat)
				break;
		if (k >= n) {
			diffs++;
			continue;
		}
		for (j = 0; j < NFREQ + NEXOBS; j++) {
			if (master[k].code[j] == CODE_NONE)
				continue;
			compared++;
			if (vbs[i].L[j] != master[k].L[j] || vbs[i].P[j] != master[k].P[j] ||
			    vbs[i].code[j] != master[k].code[j])
				diffs++;
		}
	}
	check(compared > 0, "observables compared: %d", compared);
	check(diffs == 0, "synthesized == master exactly (%d difference(s))", diffs);
}

/* --- layer 2: geometry ---------------------------------------------------------------- */

static void test_geometry(const obsd_t *master, int n, const nav_t *nav) {
	static double rs[6 * MAXOBS], dts[2 * MAXOBS], var[MAXOBS];
	static int svh[MAXOBS];
	double vbs_xyz[3], drho, worst_zero = 0.0, maxabs = 0.0;
	int i, nvalid = 0;

	printf("geometric range difference (spec S4.6)\n");

	satposs(master[0].time, master, n, nav, EPHOPT_BRDC, rs, dts, var, svh);

	/* At the master itself, drho must be identically zero for every satellite. */
	for (i = 0; i < n; i++) {
		if (rs[6 * i] == 0.0 && rs[6 * i + 1] == 0.0)
			continue;
		nvalid++;
		drho = vbs_synth_drho(rs + 6 * i, rs + 6 * i + 3, MASTER_XYZ, MASTER_XYZ);
		if (fabs(drho) > worst_zero)
			worst_zero = fabs(drho);
	}
	check(nvalid > 0, "satellite positions computed for %d satellites", nvalid);
	check(worst_zero == 0.0, "drho is exactly 0 at the master (worst |drho| = %.3e m)",
	    worst_zero);

	/* 5 km east of the master: |drho| must be bounded by the displacement. */
	memcpy(vbs_xyz, MASTER_XYZ, sizeof vbs_xyz);
	vbs_xyz[1] += 5000.0;
	for (i = 0; i < n; i++) {
		if (rs[6 * i] == 0.0 && rs[6 * i + 1] == 0.0)
			continue;
		drho = vbs_synth_drho(rs + 6 * i, rs + 6 * i + 3, MASTER_XYZ, vbs_xyz);
		if (fabs(drho) > maxabs)
			maxabs = fabs(drho);
	}
	check(maxabs > 1.0, "drho is non-zero 5 km away (max |drho| = %.1f m)", maxabs);
	check(maxabs <= 5000.0 + 1e-6,
	    "drho never exceeds the 5 km displacement (max |drho| = %.1f m)", maxabs);
}

/* --- layer 3: RTCM round-trip through an independent decoder --------------------------- */

/* One observable as the independent decoder reported it. */
struct dec_obs {
	char id[8];		/* satellite id, e.g. "G25" -- compared as a string because the
				   two libraries need not number satellites identically */
	char code[4];		/* RINEX observation code, e.g. "2X" */
	double P, L;
	int lli;
};

#define MAX_DEC	1024

/* Read the helper's output: "<satid> <obscode> <P> <L>" per line. */
static int read_dec(const char *path, struct dec_obs *out, int max) {
	char line[128];
	FILE *fp;
	int n = 0;

	if ((fp = fopen(path, "r")) == NULL)
		return -1;
	while (n < max && fgets(line, sizeof line, fp) != NULL) {
		out[n].lli = 0;
		if (sscanf(line, "%7s %3s %lf %lf %d", out[n].id, out[n].code,
		    &out[n].P, &out[n].L, &out[n].lli) >= 4)
			n++;
	}
	fclose(fp);
	return n;
}

static const struct dec_obs *find_dec(const struct dec_obs *d, int n,
	const char *id, const char *code) {
	int i;

	for (i = 0; i < n; i++)
		if (!strcmp(d[i].id, id) && !strcmp(d[i].code, code))
			return &d[i];
	return NULL;
}

static void test_roundtrip(const obsd_t *vbs, int nvbs, const obsd_t *master, int nmaster,
	const nav_t *nav) {
	static uint8_t buf[16384];
	static struct dec_obs dec[MAX_DEC];
	struct vbs_rtcm_out *enc;
	char dir[] = "/tmp/vbs_rt_XXXXXX";
	char rtcm_path[256], dump_path[256], cmd[1024], id[8];
	FILE *fp;
	int len, len1005, i, j, ndec, rc;
	static obsd_t flagged[MAXOBS];
	int half_ok = 0, bad_half = 0;
	int compared = 0, missing = 0, bad_p = 0, bad_l = 0;
	double worst_p = 0.0, worst_l = 0.0, lam, fi, dp, dl;

	printf("RTCM round-trip: MRTKLIB encode -> RTKLIB-EX decode (independent)\n");

	if (access(RTKLIB_DUMP_PATH, X_OK) != 0) {
		printf("  SKIPPED: no independent decoder at %s\n", RTKLIB_DUMP_PATH);
		return;
	}
	if (mkdtemp(dir) == NULL) {
		check(0, "temp dir created");
		return;
	}
	snprintf(rtcm_path, sizeof rtcm_path, "%s/vbs.rtcm3", dir);
	snprintf(dump_path, sizeof dump_path, "%s/vbs.txt", dir);

	if ((enc = vbs_rtcm_out_new(VBS_STAID, MASTER_XYZ)) == NULL) {
		check(0, "encoder created");
		return;
	}
	/*
	 * Without the GLONASS channel numbers the encoder drops every GLONASS satellite from the
	 * MSM message without complaint, so the round-trip must install them before encoding.
	 */
	printf("  GLONASS channel numbers installed: %d\n",
	    vbs_rtcm_out_set_glo_fcn(enc, nav));

	/*
	 * The guard, verified rather than assumed. An encoder that never had channel numbers
	 * installed must REFUSE an epoch carrying GLONASS, not emit it minus the constellation.
	 * This is the fifth instance in this engine of a capability that was implemented,
	 * documented and silently not connected (D90), so the protection against the sixth is
	 * itself tested: a safety net nobody has seen catch anything is not a safety net.
	 */
	{
		struct vbs_rtcm_out *bare = vbs_rtcm_out_new(VBS_STAID, MASTER_XYZ);
		uint8_t gbuf[8192];
		int nglo = 0, gi, grc;

		for (gi = 0; gi < nvbs; gi++)
			if (satsys(vbs[gi].sat, NULL) == SYS_GLO)
				nglo++;
		if (bare == NULL) {
			check(0, "bare encoder created");
		} else if (nglo == 0) {
			printf("  (no GLONASS in this epoch; guard not exercised)\n");
			vbs_rtcm_out_free(bare);
		} else {
			grc = vbs_rtcm_out_msm(bare, vbs[0].time, vbs, nvbs, gbuf, sizeof gbuf);
			check(grc == -1, "encoder without channel numbers refuses %d GLONASS "
			    "satellite(s) instead of dropping them silently (rc %d)", nglo, grc);
			vbs_rtcm_out_free(bare);
		}
	}
	check(vbs_rtcm_out_glo_dropped(enc) == 0,
	    "no GLONASS satellite skipped for a missing channel number (%ld)",
	    vbs_rtcm_out_glo_dropped(enc));
	len1005 = vbs_rtcm_out_1005(enc, buf, sizeof buf);
	check(len1005 > 0, "1005 encoded (%d bytes)", len1005);

	/*
	 * Raise the half-cycle flag on every second synthesized phase before encoding, so the
	 * round-trip has both states to carry and a decoder that simply cleared the field could
	 * not pass by accident.
	 */
	{
		int hi, hj;

		memcpy(flagged, vbs, (size_t)nvbs * sizeof *flagged);
		for (hi = 0; hi < nvbs; hi++)
			for (hj = 0; hj < NFREQ + NEXOBS; hj++)
				if (flagged[hi].code[hj] != CODE_NONE && flagged[hi].L[hj] != 0.0 && (hi % 2) == 0)
					flagged[hi].LLI[hj] |= 2;
		vbs = flagged;			/* everything below compares against what was encoded */
	}
	len = vbs_rtcm_out_msm(enc, vbs[0].time, vbs, nvbs, buf + len1005,
	    (int)sizeof buf - len1005);
	check(len > 0, "MSM7 encoded (%d bytes for %d satellites)", len, nvbs);
	if (len <= 0 || len1005 <= 0) {
		vbs_rtcm_out_free(enc);
		return;
	}

	if ((fp = fopen(rtcm_path, "wb")) == NULL) {
		check(0, "rtcm file written");
		vbs_rtcm_out_free(enc);
		return;
	}
	fwrite(buf, 1, len1005 + len, fp);
	fclose(fp);
	vbs_rtcm_out_free(enc);

	snprintf(cmd, sizeof cmd, "%s %s 2026/07/29 00:00:00 > %s 2>/dev/null",
	    RTKLIB_DUMP_PATH, rtcm_path, dump_path);
	rc = system(cmd);
	check(rc == 0, "independent decoder ran (exit %d)", rc);

	if ((ndec = read_dec(dump_path, dec, MAX_DEC)) <= 0) {
		check(0, "decoder produced observables");
		return;
	}
	check(ndec > 0, "decoded %d observables", ndec);

	/*
	 * Compare against the *synthesized* set, which is what was actually encoded. Comparing
	 * against the raw master would count the GLONASS satellites the engine deliberately
	 * never emitted as losses.
	 */
	/*
	 * Half-cycle ambiguity (LLI bit 1) must survive to the rover.
	 *
	 * A receiver raises it when it has locked the carrier but cannot yet tell a whole cycle
	 * from a half one; a rover that fixes integers on such a phase can land half a cycle out,
	 * and a half-cycle error is not detectable downstream. It is live on this network -- a
	 * BeiDou satellite raised it on every epoch of a Mosaic master during a Phase 2 run -- and
	 * the flag has to travel master -> synthesis -> MSM (DF420) -> rover intact. Nothing else
	 * in the suite would notice if it were dropped: the phase values would be identical.
	 */
	for (i = 0; i < nvbs; i++) {
		satno2id(vbs[i].sat, id);
		for (j = 0; j < NFREQ + NEXOBS; j++) {
			const struct dec_obs *d;

			if (vbs[i].code[j] == CODE_NONE)
				continue;
			d = find_dec(dec, ndec, id, code2obs(vbs[i].code[j]));
			if (d == NULL) {
				missing++;
				printf("  MISSING after decode: %s %s\n", id,
				    code2obs(vbs[i].code[j]));
				continue;
			}
			if ((fi = sat2freq(vbs[i].sat, vbs[i].code[j], nav)) == 0.0)
				continue;
			lam = CLIGHT / fi;
			compared++;

			if (vbs[i].P[j] != 0.0 && d->P != 0.0) {
				dp = fabs(d->P - vbs[i].P[j]);
				if (dp > worst_p)
					worst_p = dp;
				if (dp > MSM7_PR_LSB)
					bad_p++;
			}
			if (vbs[i].L[j] != 0.0 && d->L != 0.0) {
				dl = fabs(d->L - vbs[i].L[j]) * lam;
				if (dl > worst_l)
					worst_l = dl;
				if (dl > MSM7_CP_LSB)
					bad_l++;
				if (((vbs[i].LLI[j] & 2) != 0) != ((d->lli & 2) != 0)) {
					if (bad_half < 5)
						printf("  HALF-CYCLE flag lost/invented: %s %s sent %d got %d\n",
						    id, code2obs(vbs[i].code[j]), vbs[i].LLI[j], d->lli);
					bad_half++;
				} else if (vbs[i].LLI[j] & 2)
					half_ok++;
			}
		}
	}
	printf("  compared %d observables, %d missing after decode\n", compared, missing);
	printf("  half-cycle flags carried through encode/decode: %d (%d wrong)\n", half_ok, bad_half);
	printf("  worst pseudorange   error %.6f mm (1 MSM7 LSB = %.6f mm)\n",
	    worst_p * 1000.0, MSM7_PR_LSB * 1000.0);
	printf("  worst carrier-phase error %.6f mm (1 MSM7 LSB = %.6f mm)\n",
	    worst_l * 1000.0, MSM7_CP_LSB * 1000.0);

	check(compared > 0, "observables survived the round-trip");
	check(missing == 0, "no observable lost in encode/decode (%d missing)", missing);
	check(bad_p == 0, "every pseudorange within 1 MSM7 LSB (%d over)", bad_p);
	check(bad_l == 0, "every carrier phase within 1 MSM7 LSB (%d over)", bad_l);
	check(half_ok > 0 && bad_half == 0,
	    "half-cycle ambiguity flag survives to the rover (%d carried, %d wrong)", half_ok, bad_half);

	unlink(rtcm_path);
	unlink(dump_path);
	rmdir(dir);
}

int main(void) {
	static obs_t obs;
	static nav_t nav;
	static obsd_t vbs[MAXOBS];
	int n, nout;

	if ((n = read_first_epoch(MASTER_FILE, NAV_FILE, &obs, &nav)) < 0) {
		printf("SKIPPED: dataset not available (%s)\n", MASTER_FILE);
		return 0;
	}
	if (n > MAXOBS)
		n = MAXOBS;
	printf("master STA_M, first epoch: %d satellites\n", n);

	test_identity(obs.data, n, &nav, vbs);
	test_geometry(obs.data, n, &nav);

	nout = vbs_synth_epoch(obs.data, n, NULL, MASTER_XYZ, MASTER_XYZ, NULL, 0, &nav, vbs);
	test_roundtrip(vbs, nout, obs.data, n, &nav);

	printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
	return failures ? 1 : 0;
}

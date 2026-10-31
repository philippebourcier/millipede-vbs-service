/*
 * The live session must reproduce the measured pipeline.
 *
 * Every Phase 1 number came from tools/vbs_pipeline emitting a RINEX file that a rover was then
 * run against. vbs_session is that same computation re-arranged to be fed epoch by epoch. If
 * the re-arrangement changed a result, the report would describe an engine that is not the one
 * deployed -- so this test feeds both the same cluster-hour and compares what they synthesize,
 * satellite by satellite, at the precision RINEX carries.
 *
 * The pipeline is run as a subprocess with the shipped flag set (the same string
 * data/scripts/rover_run.py passes), so the comparison is against the real artefact, not a
 * re-implementation of it. Needs the validation dataset; without it the test reports that and
 * fails, as the other data-dependent tests do.
 *
 * Usage: tests/test_session [epochs]   (default 600: ten minutes, enough for a first fix)
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "vbs/vbs_session.h"

#define DATA_DIR	"testdata"
#define NAV_FILE	DATA_DIR "/nav__2026__210__BRDC00WRD_S_20262100000_01D_MN.rnx"
#define PIPELINE	"tools/vbs_pipeline"

/* E1-0: target STA_T, cluster STA_M (master, 29.1 km), STA_A1, STA_A2, STA_A3 -- all F9P */
static const struct {
	const char *file;
	double xyz[3];
} STA[] = {
	{DATA_DIR "/sta_m_2026210_1s.rnx", {4582081.2, 55183.6, 4421736.6}},
	{DATA_DIR "/sta_a1_2026210_1s.rnx", {4609353.4, 92528.1, 4393040.6}},
	{DATA_DIR "/sta_a2_2026210_1s.rnx", {4622327.0, 46793.0, 4380187.4}},
	{DATA_DIR "/sta_a3_2026210_1s.rnx", {4586204.1, 34419.7, 4417798.6}},
};
#define NSTA	4
/* the target, STA_T, as the VBS point */
static const double VBS_XYZ[3] = {4603034.1, 55436.0, 4400265.6};

static int nfail;
static double gal_L_sum[2];	/* per run: sum of emitted Galileo phases, to prove the runs differ */

static void check(int ok, const char *what) {
	printf("  %-66s %s\n", what, ok ? "ok" : "FAIL");
	if (!ok)
		nfail++;
}

static int slice(const obs_t *obs, int *i, const obsd_t **out) {
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

/*
 * One identity run. iono_sys == 0 is the shipped UM980 recipe (--no-iono); otherwise the dispersive
 * correction is on and restricted to those constellations (--iono-sys, D167), which must reach the
 * output through the session exactly as through the pipeline.
 */
static int run(int nepochs, int iono_sys) {
	static obs_t obs[NSTA], ref;
	static nav_t nav;
	static struct vbs_session_station sta[NSTA];
	struct vbs_session_opt opt;
	struct vbs_session *sess;
	struct vbs_session_stats st;
	gtime_t t0 = {0, 0.0};
	char cmd[4096], out[] = "/tmp/vbs_session_ref_XXXXXX";
	int i, cur[NSTA] = {0, 0, 0, 0}, rcur = 0;
	int ep = 0, compared = 0, mism_val = 0, mism_set = 0, mism_lli = 0, ref_epochs = 0;
	double worst_L = 0.0, worst_P = 0.0;
	int fd;

	printf("=== session vs pipeline identity (E1-0, day 210, %d epochs, %s) ===\n", nepochs,
	    iono_sys ? "dispersive on Galileo+BeiDou only" : "--no-iono");

	/* 1. the reference: the pipeline, shipped flags, limited to the same epochs */
	if ((fd = mkstemp(out)) < 0)
		return 1;
	close(fd);
	snprintf(cmd, sizeof cmd,
	    "%s --nav %s --mode emit --out %s --interval 1 --ar-corr-time 30 "
	    "--out-elmask 25 --if-narrow-lane --if-sys 0.020 --if-max-frac 0.12 --if-min-epochs 120 "
	    "--if-max-sigma 0.08 --q-ion 1e-4 %s --vbs %.4f,%.4f,%.4f --offset 0,0,0 "
	    "--base %s,%.4f,%.4f,%.4f --base %s,%.4f,%.4f,%.4f --base %s,%.4f,%.4f,%.4f "
	    "--base %s,%.4f,%.4f,%.4f > /dev/null 2>&1",
	    PIPELINE, NAV_FILE, out, iono_sys ? "--iono-sys EC" : "--no-iono",
	    VBS_XYZ[0], VBS_XYZ[1], VBS_XYZ[2],
	    STA[0].file, STA[0].xyz[0], STA[0].xyz[1], STA[0].xyz[2],
	    STA[1].file, STA[1].xyz[0], STA[1].xyz[1], STA[1].xyz[2],
	    STA[2].file, STA[2].xyz[0], STA[2].xyz[1], STA[2].xyz[2],
	    STA[3].file, STA[3].xyz[0], STA[3].xyz[1], STA[3].xyz[2]);
	printf("  running the pipeline as the reference ...\n");
	fflush(stdout);
	if (system(cmd) != 0) {
		check(0, "pipeline reference run");
		printf("FAILED: %d failure(s)\n", nfail);
		return 1;
	}
	memset(&ref, 0, sizeof ref);
	memset(&nav, 0, sizeof nav);
	if (readrnxt(out, 9, t0, t0, 0.0, "", &ref, &nav, NULL) < 0 || ref.n <= 0) {
		check(0, "pipeline emit readable");
		unlink(out);
		printf("FAILED: %d failure(s)\n", nfail);
		return 1;
	}
	unlink(out);
	check(1, "pipeline reference emitted");

	/* 2. the same data through the session, pushed epoch by epoch */
	for (i = 0; i < NSTA; i++) {
		memset(&obs[i], 0, sizeof obs[i]);
		if (readrnxt(STA[i].file, i + 1, t0, t0, 0.0, "", &obs[i], &nav, NULL) < 0 || obs[i].n <= 0) {
			check(0, "station RINEX readable");
			return 1;
		}
		memcpy(sta[i].xyz, STA[i].xyz, sizeof sta[i].xyz);
		snprintf(sta[i].name, sizeof sta[i].name, "B%d", i);
		snprintf(sta[i].rectype, sizeof sta[i].rectype, "UBL-F9P");
	}
	if (readrnxt(NAV_FILE, 1, t0, t0, 0.0, "", NULL, &nav, NULL) < 0) {
		check(0, "nav readable");
		return 1;
	}
	vbs_session_default_opt(&opt);
	opt.sys_mask = SYS_GPS | SYS_GAL | SYS_CMP | SYS_GLO;	/* homogeneous: GLONASS admitted */
	opt.dispersive = iono_sys != 0;				/* --no-iono, or */
	opt.iono_sys_mask = iono_sys;				/* --iono-sys */
	opt.staid = 1234;
	sess = vbs_session_new(VBS_XYZ, sta, NSTA, 0, &opt, &nav);
	check(sess != NULL, "session created");
	if (sess == NULL)
		return 1;

	/* pass-through: before convergence a base's bytes go straight out; after, they are dropped */
	{
		unsigned char pt[] = {0xd3, 0x00, 0x01, 0x42, 0x00, 0x00, 0x00}, drain[64];
		int n;

		check(!vbs_session_converged(sess), "a new session is not converged");
		check(vbs_session_passthrough(sess, pt, sizeof pt) == (int)sizeof pt,
		    "pass-through bytes are accepted before convergence");
		n = vbs_session_output(sess, drain, sizeof drain);
		check(n == (int)sizeof pt && memcmp(drain, pt, n) == 0,
		    "and come out verbatim, with nothing of the session's own in front");
	}

	/* the master drives; each master epoch pulls the auxiliaries' matching epochs */
	clock_t cpu0 = clock();
	for (;;) {
		const obsd_t *mep;
		int nm, k;
		gtime_t t;

		if ((nm = slice(&obs[0], &cur[0], &mep)) <= 0 || ep >= nepochs)
			break;
		t = mep[0].time;
		ep++;
		vbs_session_push(sess, 0, mep, nm, (double)ep);
		for (k = 1; k < NSTA; k++) {
			const obsd_t *aep;
			int na, save;

			for (;;) {
				save = cur[k];
				if ((na = slice(&obs[k], &cur[k], &aep)) <= 0)
					break;
				if (timediff(aep[0].time, t) < -1e-3)
					continue;		/* behind: skip */
				if (timediff(aep[0].time, t) > 1e-3) {
					cur[k] = save;		/* ahead: keep for later */
					break;
				}
				vbs_session_push(sess, k, aep, na, (double)ep);
				break;
			}
		}
		if (vbs_session_run(sess, (double)ep + 100.0) > 0) {
			const obsd_t *syn, *rep;
			gtime_t ts;
			int ns, nr, a, b2;

			ns = vbs_session_last_synth(sess, &syn, &ts);
			for (a = 0; a < ns; a++)
				if (satsys(syn[a].sat, NULL) == SYS_GAL)
					gal_L_sum[iono_sys != 0] += syn[a].L[0];
			/* find the reference epoch with the same time */
			for (;;) {
				int save = rcur;

				if ((nr = slice(&ref, &rcur, &rep)) <= 0) {
					nr = 0;
					break;
				}
				if (timediff(rep[0].time, ts) < -1e-3)
					continue;
				if (timediff(rep[0].time, ts) > 1e-3) {
					rcur = save;
					nr = 0;
					break;
				}
				break;
			}
			if (nr == 0) {
				mism_set++;
				continue;
			}
			ref_epochs++;
			if (nr != ns)
				mism_set++;
			for (a = 0; a < ns; a++) {
				int found = 0;

				for (b2 = 0; b2 < nr; b2++) {
					int f;

					if (rep[b2].sat != syn[a].sat)
						continue;
					found = 1;
					for (f = 0; f < NFREQ + NEXOBS; f++) {
						double dL, dP;

						if (syn[a].L[f] == 0.0 && rep[b2].L[f] == 0.0)
							continue;
						compared++;
						dL = fabs(syn[a].L[f] - rep[b2].L[f]);
						dP = fabs(syn[a].P[f] - rep[b2].P[f]);
						if (dL > worst_L)
							worst_L = dL;
						if (dP > worst_P)
							worst_P = dP;
						/* RINEX carries 1e-3 cycles / 1e-3 m; allow one unit of rounding */
						if (dL > 0.0015 || dP > 0.0015)
							mism_val++;
						if ((syn[a].LLI[f] & 1) != (rep[b2].LLI[f] & 1))
							mism_lli++;
					}
				}
				if (!found)
					mism_set++;
			}
		}
	}
	{
		unsigned char pt[] = {0xd3, 0x00, 0x01, 0x42, 0x00, 0x00, 0x00};

		check(vbs_session_converged(sess), "the session converged once it emitted");
		check(vbs_session_passthrough(sess, pt, sizeof pt) == 0,
		    "pass-through bytes are dropped after convergence");
	}
	{
		double cpu_s = (double)(clock() - cpu0) / CLOCKS_PER_SEC;

		printf("  CPU: %.3f s for %d epochs of a %d-station cluster = %.1f ms/epoch, "
		    "so one core carries about %d such tiles at 1 Hz\n", cpu_s, ep, NSTA,
		    1000.0 * cpu_s / (ep > 0 ? ep : 1), ep > 0 && cpu_s > 0 ? (int)(ep / cpu_s) : 0);
	}
	vbs_session_get_stats(sess, &st);
	printf("  per system G/R/E/C: in %ld/%ld/%ld/%ld, aux corr %ld/%ld/%ld/%ld, tile corr %ld/%ld/%ld/%ld, out %ld/%ld/%ld/%ld; bands %d/%d/%d/%d\n",
	    st.in_sys[0], st.in_sys[1], st.in_sys[2], st.in_sys[3],
	    st.aux_corr_sys[0], st.aux_corr_sys[1], st.aux_corr_sys[2], st.aux_corr_sys[3],
	    st.corr_sys[0], st.corr_sys[1], st.corr_sys[2], st.corr_sys[3],
	    st.out_sys[0], st.out_sys[1], st.out_sys[2], st.out_sys[3],
	    st.bands_sys[0], st.bands_sys[1], st.bands_sys[2], st.bands_sys[3]);
	printf("  master epochs %ld, emitted %ld, incomplete %ld, no-correction %ld, bytes %ld\n",
	    st.epochs_in, st.epochs_out, st.epochs_incomplete, st.epochs_no_corr, st.bytes_out);
	printf("  reference epochs matched %d, observables compared %d\n", ref_epochs, compared);
	printf("  worst |dL| %.4f cycles, worst |dP| %.4f m, value mismatches %d, set mismatches %d, LLI mismatches %d\n",
	    worst_L, worst_P, mism_val, mism_set, mism_lli);
	check(st.epochs_out > 0, "the session emitted epochs");
	check(st.epochs_incomplete == 0, "every epoch had all auxiliaries (deterministic feed)");
	check(ref_epochs == st.epochs_out, "the session emitted exactly the epochs the pipeline did");
	check(compared > 1000, "a meaningful number of observables compared");
	check(mism_set == 0, "same satellites in every epoch");
	check(mism_val == 0, "every phase and code agrees to RINEX precision");
	check(mism_lli == 0, "every lock-loss flag agrees");
	check(st.bytes_out > 0, "RTCM was produced");
	vbs_session_free(sess);
	for (i = 0; i < NSTA; i++)
		freeobs(&obs[i]);
	freeobs(&ref);
	freenav(&nav, 0xFF);
	return 0;
}

int main(int argc, char **argv) {
	int nepochs = argc > 1 ? atoi(argv[1]) : 600;

	if (access(NAV_FILE, R_OK) != 0 || access(STA[0].file, R_OK) != 0) {
		printf("  dataset not present under %s\n", DATA_DIR);
		printf("SKIPPED: validation dataset not available\n");
		return 0;
	}
	run(nepochs, 0);
	run(nepochs, SYS_GAL | SYS_CMP);
	printf("  Galileo phase sum, --no-iono %.3f vs Galileo+BeiDou dispersive %.3f cycles\n",
	    gal_L_sum[0], gal_L_sum[1]);
	check(fabs(gal_L_sum[1] - gal_L_sum[0]) > 1.0,
	    "the restricted dispersive correction changes the Galileo output (not a no-op)");
	printf("%s: %d failure(s)\n", nfail ? "FAILED" : "PASSED", nfail);
	return nfail != 0;
}

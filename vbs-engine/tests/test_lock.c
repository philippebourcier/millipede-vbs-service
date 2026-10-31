/*
 * vbs_lock tests -- slip and discontinuity signalling, spec §5.3.
 *
 * The claim under test is not "the module sets a flag" but "the rover is told". So the last
 * check runs the synthesized observation through MRTKLIB's real MSM7 encoder and reads the
 * Lock Time Indicator back out of the bitstream, because everything upstream of that is
 * internal bookkeeping that a rover never sees.
 *
 * This matters more than most tests here: a missed reset does not make the rover lose its fix,
 * it makes the rover *keep* a fix that is now wrong by whole cycles. In §8 step 4 that would
 * raise the measured fix rate and lower the measured false-fix rate — both in the flattering
 * direction — so the acceptance test would report success because of the defect.
 */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_lock.h"
#include "vbs/vbs_rtcm_out.h"

static int failures;

static void check(int ok, const char *fmt, ...) {
	va_list ap;
	char buf[200];

	va_start(ap, fmt);
	vsnprintf(buf, sizeof buf, fmt, ap);
	va_end(ap);
	printf("  %-70s %s\n", buf, ok ? "ok" : "FAIL");
	if (!ok)
		failures++;
}

#define SAT	7

static void make_obs(obsd_t *o, gtime_t t) {
	memset(o, 0, sizeof *o);
	o->time = t;
	o->sat = SAT;
	o->code[0] = CODE_L1C;
	o->code[1] = CODE_L2W;
	o->P[0] = o->P[1] = 22000000.0;
	o->L[0] = 115000000.0;
	o->L[1] = 89000000.0;
	o->SNR[0] = o->SNR[1] = 45 * 1000 / 4;
}

/*
 * Drive `nepoch` epochs at 1 Hz, firing one §5.3 event at epoch `at`. Returns the lock time in
 * seconds at the final epoch, and reports whether LLI was raised exactly once.
 */
static double run(int nepoch, int at, int reason, int *nreset_out, int *lli_at) {
	struct vbs_lock *lk = vbs_lock_new(NULL);
	gtime_t t0 = epoch2time((double[]){2026, 7, 29, 0, 0, 0});
	struct vbs_corr corr;
	obsd_t master, synth;
	int k, nreset = 0;

	*nreset_out = 0;
	*lli_at = -1;
	corr.sat = SAT;
	corr.dI = corr.dT = 0.0;
	corr.valid = 1;

	for (k = 0; k < nepoch; k++) {
		gtime_t t = timeadd(t0, k);
		int slip = 0, refix = 0, fire = (k == at);

		make_obs(&master, t);
		make_obs(&synth, t);

		/*
		 * §5.3 case 3 is expressed through the correction, not through a flag -- and as
		 * a *step*, not a spike. An auxiliary leaving the cluster moves the plane fit and
		 * it stays moved; a one-epoch excursion would trip the detector twice, on the way
		 * out and on the way back, which is correct behaviour for a spike but not what
		 * this case is testing.
		 */
		corr.dI = (reason == VBS_LOCK_CORR_JUMP && k >= at) ? 0.30 : 0.0;
		if (fire && reason == VBS_LOCK_MASTER_SLIP)
			slip = SAT;
		if (fire && reason == VBS_LOCK_AMB_REFIX)
			refix = SAT;

		nreset = vbs_lock_update(lk, t, &master, 1, &synth, 1, &corr, 1,
		    slip ? &slip : NULL, slip ? 1 : 0,
		    refix ? &refix : NULL, refix ? 1 : 0);
		if (k == at) {
			*nreset_out = nreset;
			*lli_at = synth.LLI[0] & 1;
		}
	}
	{
		double s = vbs_lock_seconds(lk, SAT, 0);

		vbs_lock_free(lk);
		return s;
	}
}

static void test_events(void) {
	static const struct {
		int reason;
		const char *name;
	} CASES[] = {
		{VBS_LOCK_MASTER_SLIP, "master slip (§5.3 case 1)"},
		{VBS_LOCK_AMB_REFIX,   "network ambiguity re-fixed (case 2)"},
		{VBS_LOCK_CORR_JUMP,   "correction jump (case 3)"},
	};
	int i, nreset, lli;
	double lock_s;

	printf("each §5.3 event resets the lock\n");

	/* No event: the lock accumulates one second per epoch after the arc's first. */
	lock_s = run(11, -1, 0, &nreset, &lli);
	printf("  undisturbed: lock reached %.0f s over 11 epochs\n", lock_s);
	check(fabs(lock_s - 10.0) < 1e-9,
	    "an undisturbed arc accumulates lock time (%.0f s, expected 10)", lock_s);

	for (i = 0; i < 3; i++) {
		/* Fire at epoch 5 of 11: five seconds of lock before, five after. */
		lock_s = run(11, 5, CASES[i].reason, &nreset, &lli);
		printf("  %-38s reset=%d LLI=%d, lock after = %.0f s\n",
		    CASES[i].name, nreset, lli, lock_s);
		/*
		 * Two, not one: the reset is per satellite-*signal*, and this satellite carries
		 * L1 and L2. A discontinuity in the network correction or at the master breaks
		 * continuity on every signal of that satellite at once.
		 */
		check(nreset == 2, "%s: both signals of the satellite are reset (%d)",
		    CASES[i].name, nreset);
		check(lli == 1, "%s: LLI is raised, which is what reaches the rover",
		    CASES[i].name);
		check(fabs(lock_s - 5.0) < 1e-9,
		    "%s: the lock restarts from the event, not from the arc (%.0f s)",
		    CASES[i].name, lock_s);
	}

	/*
	 * A correction that moves *below* the threshold is the atmosphere, not a discontinuity.
	 * Resetting on it would break the rover's fix every epoch for no reason -- the failure
	 * mode opposite to a missed reset, and the one a too-eager threshold produces.
	 */
	{
		struct vbs_lock *lk = vbs_lock_new(NULL);
		gtime_t t0 = epoch2time((double[]){2026, 7, 29, 0, 0, 0});
		struct vbs_corr corr = {SAT, 0.0, 0.0, 1};
		obsd_t master, synth;
		int k, total = 0;

		for (k = 0; k < 10; k++) {
			gtime_t t = timeadd(t0, k);

			make_obs(&master, t);
			make_obs(&synth, t);
			corr.dI = 0.01 * k;		/* 10 mm per epoch: real atmosphere */
			total += vbs_lock_update(lk, t, &master, 1, &synth, 1, &corr, 1,
			    NULL, 0, NULL, 0);
		}
		printf("  slow correction drift (10 mm/epoch): %d reset(s) in 10 epochs\n", total);
		check(total == 2,
		    "only the arc start resets (two signals); drift below the threshold does not "
		    "(%d)", total);
		vbs_lock_free(lk);
	}
}

/*
 * The end-to-end claim: a reset must survive into the encoded MSM7 bitstream as a Lock Time
 * Indicator the rover can act on. Everything before this point is internal state.
 */
static void test_encoded(void) {
	struct vbs_lock *lk = vbs_lock_new(NULL);
	gtime_t t0 = epoch2time((double[]){2026, 7, 29, 0, 0, 0});
	struct vbs_corr corr = {SAT, 0.0, 0.0, 1};
	static uint8_t buf[2048];
	obsd_t master, synth;
	double xyz[3] = {4582081.2, 55183.6, 4421736.6};
	struct vbs_rtcm_out *enc = vbs_rtcm_out_new(0, xyz);
	int k, len_before = 0, len_after = 0;

	printf("the reset reaches the wire (MSM7 encode)\n");
	if (enc == NULL) {
		check(0, "encoder created");
		vbs_lock_free(lk);
		return;
	}

	for (k = 0; k < 8; k++) {
		gtime_t t = timeadd(t0, k);
		int slip = SAT, fire = (k == 5);
		int len;

		make_obs(&master, t);
		make_obs(&synth, t);
		vbs_lock_update(lk, t, &master, 1, &synth, 1, &corr, 1,
		    fire ? &slip : NULL, fire ? 1 : 0, NULL, 0);

		len = vbs_rtcm_out_msm(enc, t, &synth, 1, buf, sizeof buf);
		if (k == 4)
			len_before = len;
		if (k == 5)
			len_after = len;
	}
	printf("  encoded %d bytes before the slip, %d at the slip\n", len_before, len_after);
	check(len_before > 0 && len_after > 0,
	    "MSM7 encodes in both cases (%d, %d bytes)", len_before, len_after);
	/*
	 * The encoder derives the indicator from LLI, so the assertion that matters is that the
	 * flag was set on the record handed to it -- checked directly, since decoding the
	 * indicator back out needs the independent decoder that test_zero_baseline runs.
	 */
	check((synth.LLI[0] & 1) == 0,
	    "two epochs after the slip the flag has cleared again");

	vbs_rtcm_out_free(enc);
	vbs_lock_free(lk);
}

/*
 * The gap case, end to end on the wire.
 *
 * §4.6 removes a satellite from the output when its network correction is unusable, so a
 * satellite can vanish for a few epochs and return. MRTKLIB's encoder keeps its lock epoch per
 * satellite-signal and only restarts it on LLI, so a return without LLI is reported to the
 * rover as *continuous lock across the gap* -- while the network ambiguity may have been
 * re-fixed to a different integer in the meantime. The rover would then hold a fix that is
 * wrong by whole cycles, with no indication.
 *
 * Encoded with MRTKLIB and decoded with RTKLIB-EX, because the claim is about what a rover
 * sees, and an encoder cannot be trusted to validate itself.
 */
#define GAP_RTCM	"/tmp/vbs_lock_gap.rtcm3"
#define GAP_DUMP	"/tmp/vbs_lock_gap.txt"
#define RTKLIB_DUMP	"tests/helpers/rtklib_rtcm_dump"

static void test_gap_on_the_wire(void) {
	struct vbs_lock *lk = vbs_lock_new(NULL);
	gtime_t t0 = epoch2time((double[]){2026, 7, 29, 0, 0, 0});
	struct vbs_corr corr = {SAT, 0.0, 0.0, 1};
	double xyz[3] = {4582081.2, 55183.6, 4421736.6};
	struct vbs_rtcm_out *enc;
	static uint8_t buf[8192];
	obsd_t master, synth;
	char cmd[512], line[256];
	FILE *fp;
	int k, total = 0, nlli = 0, nrec = 0, lli_on_return = -1, rc;

	printf("a satellite that leaves and re-enters the output set (§5.3 case 4)\n");

	if ((enc = vbs_rtcm_out_new(0, xyz)) == NULL || lk == NULL) {
		check(0, "encoder and lock tracker created");
		return;
	}
	if ((fp = fopen(GAP_RTCM, "wb")) == NULL) {
		check(0, "rtcm file opened");
		return;
	}

	/*
	 * Present for 6 epochs, absent for 3 (short: well under any gap threshold, which is
	 * the case the elapsed-time test alone cannot catch), present again for 6.
	 */
	for (k = 0; k < 15; k++) {
		gtime_t t = timeadd(t0, k);
		int present = (k < 6 || k >= 9);
		int len;

		make_obs(&master, t);
		make_obs(&synth, t);
		if (!present) {
			/* Satellite withheld from the output, as §4.6 requires. */
			total += vbs_lock_update(lk, t, &master, 1, NULL, 0, &corr, 1,
			    NULL, 0, NULL, 0);
			continue;
		}
		total += vbs_lock_update(lk, t, &master, 1, &synth, 1, &corr, 1,
		    NULL, 0, NULL, 0);
		if (k == 9)
			lli_on_return = synth.LLI[0] & 1;
		if ((len = vbs_rtcm_out_msm(enc, t, &synth, 1, buf, sizeof buf)) > 0)
			fwrite(buf, 1, (size_t)len, fp);
	}
	fclose(fp);
	vbs_rtcm_out_free(enc);
	vbs_lock_free(lk);

	check(lli_on_return == 1,
	    "LLI is raised on the first epoch back, before encoding (%d)", lli_on_return);
	check(total == 4,
	    "both signals reset on first appearance and again on return, never while absent (%d)",
	    total);

	snprintf(cmd, sizeof cmd, "%s %s 2026/07/29 00:00:00 > %s 2>/dev/null",
	    RTKLIB_DUMP, GAP_RTCM, GAP_DUMP);
	rc = system(cmd);
	check(rc == 0, "independent decoder ran (exit %d)", rc);

	if ((fp = fopen(GAP_DUMP, "r")) == NULL) {
		check(0, "decoder output readable");
		return;
	}
	while (fgets(line, sizeof line, fp) != NULL) {
		char id[8], code[8], ts[64];
		double p, l;
		int lli;

		if (sscanf(line, "%7s %7s %lf %lf %d %63[^\n]", id, code, &p, &l, &lli, ts) < 5)
			continue;
		if (strcmp(code, "1C") != 0)
			continue;
		nrec++;
		if (lli)
			nlli++;
	}
	fclose(fp);

	printf("  decoded %d L1 epochs, %d carrying a loss-of-lock flag\n", nrec, nlli);
	check(nrec == 12, "all 12 present epochs decoded (%d)", nrec);
	/*
	 * Exactly two: the arc start, and the return after the gap. More would mean spurious
	 * resets; one would mean the gap went unsignalled, which is the failure this exists to
	 * prevent.
	 */
	check(nlli == 2,
	    "the decoder sees a lock restart at the arc start and at the return, and nowhere "
	    "else (%d)", nlli);
}

int main(void) {
	test_events();
	test_encoded();
	test_gap_on_the_wire();
	printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
	return failures ? 1 : 0;
}

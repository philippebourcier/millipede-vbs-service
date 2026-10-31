/*
 * vbs_c2wcal tests -- per-base GPS L2 P(Y) (C2W) code-bias correction.
 *
 * Pure synthetic observations, no dataset files:
 *   1. which code labels are GPS P(Y) and which are L2C;
 *   2. the gate judges "one signal twice" by bit-identity, not by the size of the offset: a
 *      genuine unit whose offset happens to be near zero is still a distinct signal (D161: ~40 %
 *      of units within +-0.2 m), and a relabelling one is refused whatever its numbers;
 *   3. the gate pools per firmware build, and a base of unknown build gets no pool;
 *   4. the estimate is causal and windowed: it follows a drift (D161) and forgets what left the
 *      window, and withholds a value until it has enough samples;
 *   5. both corrections touch the GPS P(Y) code only, never L2C, other bands, GLONASS or phase.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "vbs/vbs_c2wcal.h"
#include "vbs/vbs_policy.h"

static int failures;

static void check(int ok, const char *what) {
	printf("  %-66s %s\n", what, ok ? "ok" : "FAIL");
	if (!ok)
		failures++;
}

/* A GPS record: P(Y) 2W in slot 0, L2C 2X in slot 1, L1 C/A in slot 2 (CODE_NONE drops a slot). */
static obsd_t gps_obs(int prn, uint8_t c2w, double p2w, uint8_t c2x, double p2x,
	uint8_t c1, double p1) {
	obsd_t o;

	memset(&o, 0, sizeof o);
	o.sat = (uint8_t)satno(SYS_GPS, prn);
	o.code[0] = c2w; o.P[0] = p2w; o.L[0] = 1.0e8 + prn;
	o.code[1] = c2x; o.P[1] = p2x; o.L[1] = 1.0e8 + prn;
	o.code[2] = c1;  o.P[2] = p1;
	return o;
}

static gtime_t at(double s) {
	gtime_t t;

	t.time = (time_t)(1780000000 + (long)floor(s));
	t.sec = s - floor(s);
	return t;
}

static void test_codes(void) {
	printf("code identification\n");
	check(vbs_c2wcal_is_l2_py_code(CODE_L2W) && vbs_c2wcal_is_l2_py_code(CODE_L2P) &&
	    vbs_c2wcal_is_l2_py_code(CODE_L2Y), "L2W, L2P, L2Y are P(Y)");
	check(vbs_c2wcal_is_l2c_code(CODE_L2X) && vbs_c2wcal_is_l2c_code(CODE_L2L) &&
	    vbs_c2wcal_is_l2c_code(CODE_L2S), "L2X, L2L, L2S are L2C");
	check(!vbs_c2wcal_is_l2_py_code(CODE_L2X) && !vbs_c2wcal_is_l2c_code(CODE_L2W),
	    "the two sets are disjoint");
	check(!vbs_c2wcal_is_l2_py_code(CODE_L1C) && !vbs_c2wcal_is_l2c_code(CODE_NONE),
	    "L1C and NONE are neither");
}

/* epochs x nsat records with 2W-2X = offset (+ an alternating +-noise), into p. */
static void feed(struct vbs_c2wcal_probe *p, int epochs, int nsat, double offset, double noise) {
	obsd_t obs[16];
	int e, k;

	for (e = 0; e < epochs; e++) {
		for (k = 0; k < nsat && k < 16; k++) {
			double d = offset + (((e + k) & 1) ? noise : -noise);

			obs[k] = gps_obs(k + 1, CODE_L2W, 24000000.0 + d, CODE_L2X, 24000000.0,
			    CODE_NONE, 0.0);
		}
		vbs_c2wcal_probe_add_epoch(p, obs, nsat);
	}
}

static void test_gate(void) {
	struct vbs_c2wcal_probe p;
	double mean;

	printf("gate: one signal twice is bit-identity, not a small offset\n");
	vbs_c2wcal_probe_reset(&p);
	feed(&p, 50, 8, -0.71, 0.3);
	check(vbs_c2wcal_probe_verdict(&p, &mean) == 1, "UM980-like -0.71 m, noisy -> distinct");
	check(fabs(mean + 0.71) < 1e-6, "mean reported -0.71 m");

	vbs_c2wcal_probe_reset(&p);
	feed(&p, 50, 8, 0.02, 0.3);
	check(vbs_c2wcal_probe_verdict(&p, NULL) == 1,
	    "a genuine unit at +0.02 m is still distinct (the old |mean| gate refused it)");

	vbs_c2wcal_probe_reset(&p);
	feed(&p, 50, 8, 0.0, 0.0);
	check(vbs_c2wcal_probe_verdict(&p, &mean) == 0, "bit-identical columns -> relabelled, refuse");
	check(p.n_ident == p.n_both, "every identical pair counted");

	vbs_c2wcal_probe_reset(&p);
	feed(&p, 2, 8, -0.71, 0.3);
	check(vbs_c2wcal_probe_verdict(&p, NULL) == -1, "16 samples -> no verdict yet");

	printf("gate book-keeping\n");
	{
		obsd_t o[4];

		vbs_c2wcal_probe_reset(&p);
		o[0] = gps_obs(1, CODE_L2W, 24000000.0 - 0.71, CODE_L2X, 24000000.0, CODE_NONE, 0.0);
		o[1] = gps_obs(2, CODE_L2W, 24000050.0, CODE_L2X, 24000000.0, CODE_NONE, 0.0);
		o[2] = gps_obs(3, CODE_L2W, 24000000.0, CODE_NONE, 0.0, CODE_NONE, 0.0);
		memset(&o[3], 0, sizeof o[3]);
		o[3].sat = (uint8_t)satno(SYS_GLO, 1);
		o[3].code[0] = CODE_L2P; o[3].P[0] = 21000000.0;
		o[3].code[1] = CODE_L2X; o[3].P[1] = 21000000.0;
		vbs_c2wcal_probe_add_epoch(&p, o, 4);
		check(p.n_both == 1 && fabs(p.sum_m + 0.71) < 1e-6,
		    "a 50 m blunder is kept out; GLONASS G2P is not GPS P(Y)");
		check(p.n_py_only == 1, "the L2C-less satellite is counted");
	}
}

static void test_groups(void) {
	struct vbs_c2wcal_groups g;
	struct vbs_c2wcal_probe *a, *b, *a2;

	printf("pooling per firmware build\n");
	vbs_c2wcal_groups_reset(&g);
	a = vbs_c2wcal_group(&g, vbs_policy_fw_build("2.7.0 R4.10Build17548"));
	b = vbs_c2wcal_group(&g, vbs_policy_fw_build("2.7.0 R4.10Build11833"));
	a2 = vbs_c2wcal_group(&g, vbs_policy_fw_build("2.7.0 R4.10Build17548"));
	check(a != NULL && b != NULL && a != b, "two builds, two pools");
	check(a == a2 && g.n == 2, "the same build lands in the same pool");
	check(vbs_c2wcal_group(&g, vbs_policy_fw_build("")) == NULL &&
	    vbs_c2wcal_group(&g, vbs_policy_fw_build("4.14.2")) == NULL,
	    "unknown build -> no pool (judge the base alone)");
	check(vbs_policy_fw_build("2.7.0 R4.10Build1754") == 1754 &&
	    vbs_policy_fw_build("2.7.0 R4.10Build1754") != vbs_policy_fw_build("2.7.0 R4.10Build1183"),
	    "RINEX-truncated builds still separate");

	/* pooled verdict: two units of one build, opposite signs, both genuine */
	feed(a, 40, 8, -1.5, 0.3);
	feed(a, 40, 8, +1.2, 0.3);
	check(vbs_c2wcal_probe_verdict(a, NULL) == 1,
	    "units at -1.5 and +1.2 m pool to a distinct build (no cancellation)");
}

static void test_estimate(void) {
	struct vbs_c2wcal_est e;
	obsd_t obs[8];
	double m;
	long n;
	int s, k;

	printf("estimate: causal, windowed\n");
	vbs_c2wcal_est_init(&e, 600.0);
	for (s = 0; s < 3600; s++) {			/* 2W-2X drifts +1 m over the hour */
		double off = (double)s / 3600.0;

		for (k = 0; k < 8; k++)
			obs[k] = gps_obs(k + 1, CODE_L2W, 22000000.0 + off, CODE_L2X, 22000000.0,
			    CODE_NONE, 0.0);
		vbs_c2wcal_est_add_epoch(&e, at(s), obs, 8);
		if (s == 3)
			check(!vbs_c2wcal_est_value(&e, at(s), &m, &n) && n == 32,
			    "32 samples < minimum -> no estimate yet");
	}
	check(vbs_c2wcal_est_value(&e, at(3599), &m, &n), "estimate available");
	check(fabs(m - (3000.0 + 3599.0) / 2.0 / 3600.0) < 0.01,
	    "10-min window tracks the drift (~0.92 m, not the hour mean 0.5)");
	check(!vbs_c2wcal_est_value(&e, at(3599 + 1200), &m, &n) && n == 0,
	    "20 min later the window has emptied");

	vbs_c2wcal_est_init(&e, 0.0);
	for (s = 0; s < 100; s++) {
		for (k = 0; k < 4; k++)
			obs[k] = gps_obs(k + 1, CODE_L2W, 22000000.0 + (s < 50 ? -1.0 : 1.0),
			    CODE_L2X, 22000000.0, CODE_NONE, 0.0);
		vbs_c2wcal_est_add_epoch(&e, at(s * 3600.0), obs, 4);
	}
	check(vbs_c2wcal_est_value(&e, at(1e7), &m, NULL) && fabs(m) < 1e-9,
	    "window 0 = cumulative mean over everything");
}

static void test_apply(void) {
	obsd_t obs[4];
	int a, b;

	printf("apply: GPS P(Y) code only\n");
	obs[0] = gps_obs(1, CODE_L2W, 24000001.0, CODE_L2X, 24000002.0, CODE_L1C, 24000000.0);
	obs[1] = gps_obs(5, CODE_L2W, 22000001.5, CODE_NONE, 0.0, CODE_L1C, 22000000.0);
	obs[2] = gps_obs(10, CODE_L2W, 0.0, CODE_NONE, 0.0, CODE_L1C, 23000000.0);
	memset(&obs[3], 0, sizeof obs[3]);
	obs[3].sat = (uint8_t)satno(SYS_GLO, 3);
	obs[3].code[0] = CODE_L2P; obs[3].P[0] = 21000000.0;

	check(vbs_c2wcal_apply(obs, 4, 0.5) == 2, "constant: the two populated GPS P(Y) codes");
	check(fabs(obs[0].P[0] - 24000000.5) < 1e-9 && fabs(obs[1].P[0] - 22000001.0) < 1e-9,
	    "constant: each reduced by 0.5 m");
	check(fabs(obs[0].P[1] - 24000002.0) < 1e-9 && fabs(obs[0].P[2] - 24000000.0) < 1e-9,
	    "constant: L2C and L1 untouched");
	check(obs[2].P[0] == 0.0 && fabs(obs[3].P[0] - 21000000.0) < 1e-9,
	    "constant: empty slot stays empty; GLONASS G2P untouched");
	check(vbs_c2wcal_apply(obs, 4, 0.0) == 0, "constant: zero is a no-op");

	obs[0] = gps_obs(1, CODE_L2W, 24000001.0, CODE_L2X, 24000002.0, CODE_L1C, 24000000.0);
	obs[1] = gps_obs(5, CODE_L2W, 22000001.5, CODE_NONE, 0.0, CODE_L1C, 22000000.0);
	check(vbs_c2wcal_apply_l2c(obs, 4, 0.5, &a, &b) == 2 && a == 1 && b == 1,
	    "l2c: one satellite takes L2C, one the estimate");
	check(fabs(obs[0].P[0] - 24000002.0) < 1e-9 && fabs(obs[0].P[1] - 24000002.0) < 1e-9,
	    "l2c: P(Y) slot now holds the L2C code, L2C slot unchanged");
	check(fabs(obs[1].P[0] - 22000001.0) < 1e-9, "l2c: L2C-less satellite gets the estimate");
	check(obs[0].L[0] == 1.0e8 + 1, "l2c: carrier phase untouched");
	check(fabs(obs[3].P[0] - 21000000.0) < 1e-9, "l2c: GLONASS untouched");
}

/* The live path vbsd runs: gate pooled per build, causal estimate, corrected copy or nothing. */
static void test_live_base(void) {
	struct vbs_c2wcal_groups g;
	struct vbs_c2wcal_base a, r;
	obsd_t in[8], out[8];
	double off;
	int s, k, v, n, first = -1;
	long build = vbs_policy_fw_build("2.7.0 R4.10Build17548");

	printf("live base (vbsd): gate, estimate, correction\n");
	vbs_c2wcal_groups_reset(&g);
	vbs_c2wcal_base_init(&a, VBS_C2WCAL_L2C, 1800.0);
	for (s = 0; s < 120; s++) {
		for (k = 0; k < 6; k++)		/* 5 L2C satellites, 2W-2X = -0.71 m +- noise */
			in[k] = gps_obs(k + 1, CODE_L2W, 22000000.0 - 0.71 + ((s + k) & 1 ? 0.2 : -0.2),
			    CODE_L2X, 22000000.0, CODE_NONE, 0.0);
		in[5] = gps_obs(20, CODE_L2W, 23000000.0, CODE_NONE, 0.0, CODE_NONE, 0.0);
		n = vbs_c2wcal_base_epoch(&a, &g, build, in, 6, out, 8, &v, &off);
		if (n > 0 && first < 0)
			first = s;
	}
	check(first == 59, "no correction until the gate has 300 samples (5 per epoch: epoch 59)");
	check(n == 6 && v == 1 && fabs(off + 0.71) < 0.01, "then a corrected copy, offset -0.71 m");
	check(fabs(out[0].P[0] - out[0].P[1]) < 1e-9 && fabs(in[0].P[0] - in[0].P[1]) > 0.4,
	    "l2c: the copy carries the L2C code in the P(Y) slot; the input is untouched");
	check(fabs(out[5].P[0] - (23000000.0 - off)) < 1e-9,
	    "l2c: the L2C-less satellite takes the estimate");

	vbs_c2wcal_base_init(&r, VBS_C2WCAL_CONST, 1800.0);
	for (s = 0; s < 120; s++) {
		for (k = 0; k < 6; k++)		/* a relabelling unit: 2W == 2X exactly */
			in[k] = gps_obs(k + 1, CODE_L2W, 22000000.0, CODE_L2X, 22000000.0, CODE_NONE, 0.0);
		n = vbs_c2wcal_base_epoch(&r, &g, vbs_policy_fw_build("x Build99999"), in, 6, out, 8,
		    &v, &off);
	}
	check(n == 0 && v == 0 && isnan(off), "a relabelling build is never corrected");
	check(g.n == 2, "one pool per build");

	vbs_c2wcal_base_init(&r, VBS_C2WCAL_CONST, 1800.0);
	for (s = 0; s < 120; s++) {	/* same build as `a`: its pooled verdict is already distinct */
		for (k = 0; k < 6; k++)
			in[k] = gps_obs(k + 1, CODE_L2W, 22000000.0 + 1.2, CODE_L2X, 22000000.0, CODE_NONE,
			    0.0);
		n = vbs_c2wcal_base_epoch(&r, &g, build, in, 6, out, 8, &v, &off);
		if (s == 10)
			check(n == 6 && fabs(off - 1.2) < 1e-6,
			    "a second unit of a known-good build corrects on its own estimate at once");
	}
	check(fabs(out[0].P[0] - 22000000.0) < 1e-6, "constant: P(Y) brought to the L2C level");
}

int main(void) {
	printf("vbs_c2wcal tests\n");
	test_codes();
	test_gate();
	test_groups();
	test_estimate();
	test_apply();
	test_live_base();
	printf("%s (%d failure%s)\n", failures ? "FAIL" : "PASS", failures,
	    failures == 1 ? "" : "s");
	return failures != 0;
}

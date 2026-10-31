/*
 * vbs_interp tests -- spatial interpolation at the VBS point, spec §4.5.
 *
 * The properties worth pinning down are algebraic, so they are tested against fields whose
 * right answer is known exactly rather than against recorded data:
 *
 *   - a field that really is planar must be reproduced *exactly* at the VBS, inside the
 *     cluster and outside it alike. That is the whole claim of a plane fit, and it is the one
 *     a sign error in the local-frame offset breaks;
 *   - at the master the correction must be zero, because the corrections are differences
 *     against the master and the field is zero there by construction;
 *   - the leverage must separate interpolation from extrapolation, since
 *     `vbs.max_extrapolation_km` is calibrated against exactly that number (D17 item 4);
 *   - a noisy auxiliary must be down-weighted, not silently treated as equal.
 */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_interp.h"

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

/* Master near Toulouse; auxiliaries placed by ENU offset around it. */
static const double MASTER[3] = {4582081.2, 55183.6, 4421736.6};

#define NAUX	4
#define NSAT	3
static const int SATS[NSAT] = {1, 5, 9};

/*
 * The planar truth: dI = gI·(E, N), dT = gT·(E, N), in metres with E/N in km.
 *
 * Zero at the master, and per satellite only through its *gradient*. These corrections are
 * differences against the master, so the field is zero there by definition -- a synthetic field
 * with a per-satellite constant added would be one the engine cannot produce, and testing
 * against it would forbid the fit from using the master as the sample it is (D69).
 */
static const double GI[2] = {0.0011, -0.0007};
static const double GT[2] = {-0.0004, 0.0009};
#define SGI(s)	(0.0002 * (s))		/* per-satellite tilt, so satellites stay distinguishable */
#define SGT(s)	(-0.0003 * (s))

static void enu_to_ecef(const double *mpos, double e, double n, double *xyz) {
	double denu[3] = {e, n, 0.0}, dr[3];
	int i;

	enu2ecef(mpos, denu, dr);
	for (i = 0; i < 3; i++)
		xyz[i] = MASTER[i] + dr[i];
}

/* Build a cluster and interpolate at (ve, vn) km from the master. Returns n corrections. */
static int run(double ve, double vn, const double *aux_e, const double *aux_n,
	const double *sigma_scale, struct vbs_corr *out, struct vbs_interp_quality *qual) {
	struct vbs_interp_aux aux[NAUX];
	static struct vbs_corr corr[NAUX][NSAT];
	static double sdI[NAUX][NSAT], sdT[NAUX][NSAT];
	struct vbs_interp_opt popt;
	double mpos[3], vbs[3];
	gtime_t t = epoch2time((double[]){2026, 7, 29, 12, 0, 0});
	int k, s;

	memset(aux, 0, sizeof aux);	/* a field added later must start at 0 */
	ecef2pos(MASTER, mpos);
	for (k = 0; k < NAUX; k++) {
		enu_to_ecef(mpos, aux_e[k] * 1000.0, aux_n[k] * 1000.0, aux[k].xyz);
		for (s = 0; s < NSAT; s++) {
			corr[k][s].sat = SATS[s];
			corr[k][s].dI = (GI[0] + SGI(s)) * aux_e[k] + GI[1] * aux_n[k];
			corr[k][s].dT = (GT[0] + SGT(s)) * aux_e[k] + GT[1] * aux_n[k];
			corr[k][s].valid = 1;
			sdI[k][s] = sdT[k][s] = 0.004 * (sigma_scale ? sigma_scale[k] : 1.0);
		}
		aux[k].corr = corr[k];
		aux[k].sd_I = sdI[k];
		aux[k].sd_T = sdT[k];
		aux[k].ncorr = NSAT;
	}
	enu_to_ecef(mpos, ve * 1000.0, vn * 1000.0, vbs);
	/*
	 * Explicitly the plane fit. Exact reproduction of a planar field is the plane's defining
	 * property and the thing a sign error in the local frame breaks; collocation deliberately
	 * does not have it, reverting toward the field mean where there is no data, so testing it
	 * through whichever method happens to be the default would test the wrong claim.
	 */
	vbs_interp_default_opt(&popt);
	popt.method = VBS_INTERP_METHOD_PLANE;
	return vbs_interp_epoch(t, MASTER, vbs, aux, NAUX, &popt, NULL, out, qual, NSAT * 2);
}

/* Truth at the VBS for satellite index s. */
static double truth_dI(double ve, double vn, int s) {
	return (GI[0] + SGI(s)) * ve + GI[1] * vn;
}

static void test_planar(const char *label, double ve, double vn, double tol_mm) {
	struct vbs_corr out[NSAT * 2];
	struct vbs_interp_quality q[NSAT * 2];
	/* A square cluster, 30 km on a side, centred on the master. */
	static const double AE[NAUX] = {-15.0, 15.0, 15.0, -15.0};
	static const double AN[NAUX] = {-15.0, -15.0, 15.0, 15.0};
	double worst = 0.0, lev = 0.0;
	int n, i, s;

	n = run(ve, vn, AE, AN, NULL, out, q);
	for (i = 0; i < n; i++) {
		for (s = 0; s < NSAT; s++) {
			double err;

			if (out[i].sat != SATS[s])
				continue;
			err = fabs(out[i].dI - truth_dI(ve, vn, s));
			if (err > worst)
				worst = err;
		}
		if (q[i].leverage > lev)
			lev = q[i].leverage;
	}
	printf("  %-28s %d satellite(s), worst dI error %.4f mm, leverage %.2f\n",
	    label, n, worst * 1000.0, lev);
	check(n == NSAT, "%s: every satellite interpolated (%d of %d)", label, n, NSAT);
	check(worst * 1000.0 < tol_mm, "%s: planar field reproduced within %.3g mm (%.4f)",
	    label, tol_mm, worst * 1000.0);
}

static void test_leverage(void) {
	struct vbs_corr out[NSAT * 2];
	struct vbs_interp_quality qin[NSAT * 2], qout[NSAT * 2];
	static const double AE[NAUX] = {-15.0, 15.0, 15.0, -15.0};
	static const double AN[NAUX] = {-15.0, -15.0, 15.0, 15.0};
	double lin = 0.0, lout = 0.0;
	int n, i;

	n = run(0.0, 0.0, AE, AN, NULL, out, qin);
	for (i = 0; i < n; i++)
		if (qin[i].leverage > lin)
			lin = qin[i].leverage;
	n = run(120.0, 0.0, AE, AN, NULL, out, qout);
	for (i = 0; i < n; i++)
		if (qout[i].leverage > lout)
			lout = qout[i].leverage;

	printf("  leverage at the cluster centre %.3f, 120 km outside the hull %.1f\n",
	    lin, lout);
	check(lin < 0.5, "inside a 4-station cluster the leverage is ~1/n (%.3f)", lin);
	check(lout > 5.0 * lin,
	    "extrapolation inflates the leverage, which is what D17's "
	    "max_extrapolation_km is calibrated on (%.1f vs %.3f)", lout, lin);
}

static void test_at_master(void) {
	struct vbs_corr out[NSAT * 2];
	struct vbs_interp_quality q[NSAT * 2];
	static const double AE[NAUX] = {-15.0, 15.0, 15.0, -15.0};
	static const double AN[NAUX] = {-15.0, -15.0, 15.0, 15.0};
	double worst = 0.0;
	int n, i, s;

	n = run(0.0, 0.0, AE, AN, NULL, out, q);
	for (i = 0; i < n; i++)
		for (s = 0; s < NSAT; s++) {
			if (out[i].sat != SATS[s])
				continue;
			/*
			 * Zero: the corrections are differences against the master, so the
			 * field is zero there by definition -- which is exactly what
			 * `anchor_master` now spends the fit's redundancy on (D69).
			 */
			if (fabs(out[i].dI) > worst)
				worst = fabs(out[i].dI);
		}
	printf("  at the master, worst departure from the master-relative value %.4f mm\n",
	    worst * 1000.0);
	check(worst < 1e-9, "a VBS at the master reproduces the master-relative field exactly");
}

/*
 * Weighting. One auxiliary is given a correction far off the plane and a large sigma; the fit
 * must lean on the other three. Without weighting the outlier drags the plane by roughly a
 * quarter of its error, which at 10 cm is 2.5 cm of correction error at the VBS.
 */
static void test_weighting(void) {
	struct vbs_interp_aux aux[NAUX];
	static struct vbs_corr corr[NAUX][1];
	static double sdI[NAUX][1], sdT[NAUX][1];
	struct vbs_corr out[4];
	struct vbs_interp_quality q[4];
	static const double AE[NAUX] = {-15.0, 15.0, 15.0, -15.0};
	static const double AN[NAUX] = {-15.0, -15.0, 15.0, 15.0};
	struct vbs_interp_opt opt;
	double mpos[3], vbs[3];
	gtime_t t = epoch2time((double[]){2026, 7, 29, 12, 0, 0});
	int k, n;

	memset(aux, 0, sizeof aux);	/* a field added later must start at 0 */
	vbs_interp_default_opt(&opt);
	opt.max_sigma = 0.0;			/* keep the outlier in, judge it by weight alone */
	ecef2pos(MASTER, mpos);
	for (k = 0; k < NAUX; k++) {
		enu_to_ecef(mpos, AE[k] * 1000.0, AN[k] * 1000.0, aux[k].xyz);
		corr[k][0].sat = SATS[0];
		corr[k][0].dI = GI[0] * AE[k] + GI[1] * AN[k];
		corr[k][0].dT = GT[0] * AE[k] + GT[1] * AN[k];
		corr[k][0].valid = 1;
		sdI[k][0] = sdT[k][0] = 0.004;
		aux[k].corr = corr[k];
		aux[k].sd_I = sdI[k];
		aux[k].sd_T = sdT[k];
		aux[k].ncorr = 1;
	}
	corr[2][0].dI += 0.10;			/* 10 cm off the plane */
	sdI[2][0] = 1.0;			/* and it says so: 1 m sigma */

	enu_to_ecef(mpos, 5000.0, 5000.0, vbs);
	n = vbs_interp_epoch(t, MASTER, vbs, aux, NAUX, &opt, NULL, out, q, 4);

	printf("  one auxiliary 10 cm off the plane at 1 m sigma -> dI error %.2f mm\n",
	    n > 0 ? fabs(out[0].dI - (GI[0] * 5.0 + GI[1] * 5.0)) * 1000.0 : -1.0);
	check(n == 1, "satellite still interpolated (%d)", n);
	check(n == 1 && fabs(out[0].dI - (GI[0] * 5.0 + GI[1] * 5.0)) < 0.005,
	    "the outlier is down-weighted, not averaged in");
}

/*
 * The troposphere model-term path, rewritten in D34 and — until this test — carrying no
 * verification of its sign at all. D34 drew a conclusion about the non-dispersive correction
 * while resting on it, which is precisely when an unverified sign is most dangerous: a flipped
 * model term would produce exactly the anti-correlation D34 reported and would have been read
 * as physics.
 *
 * The check is a closed loop. A VBS placed *above* the master sits in less atmosphere, so its
 * slant delay is smaller and the correction that converts a master observation into a VBS one
 * must be negative; below the master, positive. Magnitude is checked against the model
 * evaluated directly, so this pins the sign and the scale, not just the sign.
 */
static void test_model_term_sign(void) {
	struct vbs_interp_aux aux[NAUX];
	static struct vbs_corr corr[NAUX][1];
	static double sdI[NAUX][1], sdT[NAUX][1];
	struct vbs_corr out[4];
	static double azel[MAXSAT * 2];
	static const double AE[NAUX] = {-15.0, 15.0, 15.0, -15.0};
	static const double AN[NAUX] = {-15.0, -15.0, 15.0, 15.0};
	struct vbs_interp_opt opt;
	double mpos[3], vbs[3], vpos[3], expect, el = 30.0 * D2R;
	gtime_t t = epoch2time((double[]){2026, 7, 29, 12, 0, 0});
	int k, n, s;

	memset(aux, 0, sizeof aux);	/* a field added later must start at 0 */
	printf("troposphere model term: sign and scale (D34 rewrite, D35 item 3)\n");

	vbs_interp_default_opt(&opt);
	opt.apply_nondispersive = 0;		/* isolate the model term */
	ecef2pos(MASTER, mpos);

	/* Every correction zero, so anything that comes out is the model term alone. */
	for (k = 0; k < NAUX; k++) {
		enu_to_ecef(mpos, AE[k] * 1000.0, AN[k] * 1000.0, aux[k].xyz);
		corr[k][0].sat = SATS[0];
		corr[k][0].dI = corr[k][0].dT = 0.0;
		corr[k][0].valid = 1;
		sdI[k][0] = sdT[k][0] = 0.004;
		aux[k].corr = corr[k];
		aux[k].sd_I = sdI[k];
		aux[k].sd_T = sdT[k];
		aux[k].ncorr = 1;
	}
	azel[(SATS[0] - 1) * 2] = 0.0;
	azel[(SATS[0] - 1) * 2 + 1] = el;

	for (s = 0; s < 2; s++) {
		/*
		 * +-200 m keeps both points inside tropmodel()'s valid range (the master sits at
		 * 110 m, so -200 m lands at -90 m, just above the -100 m floor). The
		 * out-of-range case is tested separately below, because it has a different
		 * correct answer: omit the term, not compute it against a zero.
		 */
		double dh = s == 0 ? +200.0 : -200.0;
		double denu[3] = {0.0, 0.0, 0.0}, dr[3];
		int i;

		denu[2] = dh;
		enu2ecef(mpos, denu, dr);
		for (i = 0; i < 3; i++)
			vbs[i] = MASTER[i] + dr[i];
		ecef2pos(vbs, vpos);

		n = vbs_interp_epoch(t, MASTER, vbs, aux, NAUX, &opt, azel, out, NULL, 4);
		expect = tropmodel(t, vpos, (double[]){0.0, el}, opt.humi)
		       - tropmodel(t, mpos, (double[]){0.0, el}, opt.humi);

		printf("  VBS %+.0f m: dT = %+.4f m (model says %+.4f)\n",
		    dh, n > 0 ? out[0].dT : 0.0, expect);
		check(n == 1, "satellite emitted with the model term (%d)", n);
		check(n == 1 && fabs(out[0].dT - expect) < 1e-9,
		    "model term matches the model, sign and scale (%+.4f vs %+.4f)",
		    n > 0 ? out[0].dT : 0.0, expect);
		check(n == 1 && ((dh > 0.0) == (out[0].dT < 0.0)),
		    "a VBS %s the master gets a %s correction: less atmosphere above",
		    dh > 0 ? "above" : "below", dh > 0 ? "negative" : "positive");
	}

	/*
	 * A VBS outside tropmodel()'s valid height range. The model cannot answer, so the term
	 * must be omitted -- not differenced against the 0.0 the model returns to say so.
	 */
	{
		double denu[3] = {0.0, 0.0, -5000.0}, dr[3];
		int i;

		enu2ecef(mpos, denu, dr);
		for (i = 0; i < 3; i++)
			vbs[i] = MASTER[i] + dr[i];
		n = vbs_interp_epoch(t, MASTER, vbs, aux, NAUX, &opt, azel, out, NULL, 4);
		printf("  VBS 5 km below the master: dT = %+.4f m\n", n > 0 ? out[0].dT : 0.0);
		check(n == 1 && out[0].dT == 0.0,
		    "an unmodellable height omits the term rather than differencing against 0 "
		    "(%.4f)", n > 0 ? out[0].dT : -1.0);
	}

	/* With no azel supplied the term must be omitted entirely, not guessed at. */
	n = vbs_interp_epoch(t, MASTER, vbs, aux, NAUX, &opt, NULL, out, NULL, 4);
	check(n == 1 && out[0].dT == 0.0,
	    "without elevations the model term is omitted, not approximated (%.6f)",
	    n > 0 ? out[0].dT : -1.0);
}

/*
 * D55: each auxiliary's single-difference ionosphere state carries a constant of that
 * master-auxiliary pair -- the two receivers' differential code bias -- which is not
 * ionosphere and must not be interpolated.
 *
 * The property that matters is independence, not a particular value: whatever the biases are,
 * the correction handed to the VBS must be the same. So the same planar field is interpolated
 * twice with different per-auxiliary biases and the two results are required to agree. With
 * the datum off, the same comparison must *fail* -- a test that passes either way would not be
 * testing anything.
 *
 * What the datum version returns is the field minus its own constellation mean at the VBS
 * (a double difference in all but name), which is checked against that closed form rather than
 * against the raw field.
 */
static void run_biased(double ve, double vn, const double *bias, int datum,
	struct vbs_corr *out, int *nout) {
	struct vbs_interp_aux aux[NAUX];
	static struct vbs_corr corr[NAUX][NSAT];
	static double sdI[NAUX][NSAT], sdT[NAUX][NSAT];
	struct vbs_interp_opt opt;
	static const double AE[NAUX] = {-15.0, 15.0, 15.0, -15.0};
	static const double AN[NAUX] = {-15.0, -15.0, 15.0, 15.0};
	double mpos[3], vbs[3];
	gtime_t t = epoch2time((double[]){2026, 7, 29, 12, 0, 0});
	int k, s;

	memset(aux, 0, sizeof aux);	/* a field added later must start at 0 */
	vbs_interp_default_opt(&opt);
	opt.method = VBS_INTERP_METHOD_PLANE;	/* exact arithmetic: see run() */
	opt.datum = datum;
	ecef2pos(MASTER, mpos);
	for (k = 0; k < NAUX; k++) {
		enu_to_ecef(mpos, AE[k] * 1000.0, AN[k] * 1000.0, aux[k].xyz);
		for (s = 0; s < NSAT; s++) {
			corr[k][s].sat = SATS[s];
			corr[k][s].dI = (GI[0] + SGI(s)) * AE[k] + GI[1] * AN[k] + bias[k];
			corr[k][s].dT = (GT[0] + SGT(s)) * AE[k] + GT[1] * AN[k];
			corr[k][s].valid = 1;
			sdI[k][s] = sdT[k][s] = 0.004;
		}
		aux[k].corr = corr[k];
		aux[k].sd_I = sdI[k];
		aux[k].sd_T = sdT[k];
		aux[k].ncorr = NSAT;
	}
	enu_to_ecef(mpos, ve * 1000.0, vn * 1000.0, vbs);
	*nout = vbs_interp_epoch(t, MASTER, vbs, aux, NAUX, &opt, NULL, out, NULL, NSAT * 2);
}

static double worst_diff(const struct vbs_corr *a, int na, const struct vbs_corr *b, int nb) {
	double worst = 0.0;
	int i, j;

	if (na != nb)
		return 1e9;
	for (i = 0; i < na; i++)
		for (j = 0; j < nb; j++)
			if (a[i].sat == b[j].sat && fabs(a[i].dI - b[j].dI) > worst)
				worst = fabs(a[i].dI - b[j].dI);
	return worst;
}

static void test_datum(void) {
	/* Biases of the size measured on the real networks: 0.5-2.5 m, differing per station. */
	static const double B1[NAUX] = {-0.67, 0.21, 0.86, 1.55};
	static const double B2[NAUX] = {1.40, -0.30, 0.05, -1.10};
	static const double B0[NAUX] = {0.0, 0.0, 0.0, 0.0};
	struct vbs_corr o1[NSAT * 2], o2[NSAT * 2], o0[NSAT * 2];
	int n1, n2, n0, s;
	double mean = 0.0;

	run_biased(7.0, -4.0, B1, VBS_INTERP_DATUM_OFF, o1, &n1);
	run_biased(7.0, -4.0, B2, VBS_INTERP_DATUM_OFF, o2, &n2);
	printf("  datum off:  two bias sets differ by %.1f mm at the VBS\n",
	    1000.0 * worst_diff(o1, n1, o2, n2));
	check(worst_diff(o1, n1, o2, n2) > 0.1,
	    "without the datum the receiver bias reaches the VBS (it must, or nothing is tested)");

	run_biased(7.0, -4.0, B1, VBS_INTERP_DATUM_COMMON, o1, &n1);
	run_biased(7.0, -4.0, B2, VBS_INTERP_DATUM_COMMON, o2, &n2);
	printf("  datum on:   two bias sets differ by %.6f mm at the VBS\n",
	    1000.0 * worst_diff(o1, n1, o2, n2));
	check(n1 == NSAT && worst_diff(o1, n1, o2, n2) < 1e-6,
	    "with the datum removed the correction is independent of the receiver biases");

	/* And it equals the field minus its constellation mean, evaluated at the VBS. */
	for (s = 0; s < NSAT; s++)
		mean += truth_dI(7.0, -4.0, s) / NSAT;
	run_biased(7.0, -4.0, B0, VBS_INTERP_DATUM_COMMON, o0, &n0);
	{
		double worst = 0.0;
		int i;

		for (i = 0; i < n0; i++)
			for (s = 0; s < NSAT; s++)
				if (o0[i].sat == SATS[s]) {
					double want = truth_dI(7.0, -4.0, s) - mean;

					if (fabs(o0[i].dI - want) > worst)
						worst = fabs(o0[i].dI - want);
				}
		printf("  datum on:   worst departure from (field - constellation mean) %.6f mm\n",
		    1000.0 * worst);
		check(worst < 1e-9,
		    "the datum version returns the planar field referred to its own constellation mean");
	}
}

/*
 * Anchoring (D69). A correction field is zero at the master by definition, so three auxiliaries
 * that all sit a common amount above zero describe something the physics does not allow -- a
 * per-baseline bias, not a field. A free plane through three points reproduces that bias at the
 * VBS in full; anchored, the master pulls it back. What is tested is the direction and that a
 * genuinely planar field is still exact, not a particular fraction.
 */
static void test_anchor(void) {
	struct vbs_interp_aux aux[NAUX];
	static struct vbs_corr corr[NAUX][NSAT];
	static double sdI[NAUX][NSAT], sdT[NAUX][NSAT];
	static const double AE[NAUX] = {-15.0, 15.0, 15.0, -15.0};
	static const double AN[NAUX] = {-15.0, -15.0, 15.0, 15.0};
	struct vbs_interp_opt opt;
	struct vbs_corr out[NSAT * 2];
	double mpos[3], vbs[3], free_v = 0.0, anch_v = 0.0;
	gtime_t t = epoch2time((double[]){2026, 7, 29, 12, 0, 0});
	int k, s, n, i;

	memset(aux, 0, sizeof aux);	/* a field added later must start at 0 */
	printf("anchoring the fit at the master (D69)\n");
	ecef2pos(MASTER, mpos);
	for (k = 0; k < NAUX; k++) {
		enu_to_ecef(mpos, AE[k] * 1000.0, AN[k] * 1000.0, aux[k].xyz);
		for (s = 0; s < NSAT; s++) {
			corr[k][s].sat = SATS[s];
			corr[k][s].dI = (GI[0] + SGI(s)) * AE[k] + GI[1] * AN[k] + 0.05;
			corr[k][s].dT = 0.0;
			corr[k][s].valid = 1;
			sdI[k][s] = sdT[k][s] = 0.004;
		}
		aux[k].corr = corr[k];
		aux[k].sd_I = sdI[k];
		aux[k].sd_T = sdT[k];
		aux[k].ncorr = NSAT;
	}
	enu_to_ecef(mpos, 7000.0, -4000.0, vbs);
	vbs_interp_default_opt(&opt);
	opt.method = VBS_INTERP_METHOD_PLANE;
	opt.anchor_master = 0;
	n = vbs_interp_epoch(t, MASTER, vbs, aux, NAUX, &opt, NULL, out, NULL, NSAT * 2);
	for (i = 0; i < n; i++)
		if (out[i].sat == SATS[0])
			free_v = out[i].dI;
	opt.anchor_master = 1;
	n = vbs_interp_epoch(t, MASTER, vbs, aux, NAUX, &opt, NULL, out, NULL, NSAT * 2);
	for (i = 0; i < n; i++)
		if (out[i].sat == SATS[0])
			anch_v = out[i].dI;
	printf("  a 5 cm offset common to every auxiliary reaches the VBS as %.1f mm free, "
	    "%.1f mm anchored\n", 1000.0 * (free_v - truth_dI(7.0, -4.0, 0)),
	    1000.0 * (anch_v - truth_dI(7.0, -4.0, 0)));
	check(fabs(anch_v - truth_dI(7.0, -4.0, 0)) < fabs(free_v - truth_dI(7.0, -4.0, 0)),
	    "anchoring pulls a physically impossible common offset back toward zero");
}

/*
 * Collocation (D70). Two properties worth pinning, and they are not the plane's.
 *
 * Inside the cluster it must track a planar field closely -- not exactly, because it shrinks
 * toward the field's zero mean by design, but to well within the centimetre the corrections
 * are worth. Far outside it must *stop* extrapolating: a plane fitted to four stations and
 * evaluated 60 km beyond them runs away linearly, while collocation reverts to the mean, which
 * for a VBS is the safer failure and is why `max_extrapolation_km` means something different
 * under it.
 */
static void test_collocation(void) {
	struct vbs_interp_aux aux[NAUX];
	static struct vbs_corr corr[NAUX][NSAT];
	static double sdI[NAUX][NSAT], sdT[NAUX][NSAT];
	static const double AE[NAUX] = {-15.0, 15.0, 15.0, -15.0};
	static const double AN[NAUX] = {-15.0, -15.0, 15.0, 15.0};
	struct vbs_interp_opt opt;
	struct vbs_corr out[NSAT * 2];
	double mpos[3], vbs[3], inside = 0.0, far_plane = 0.0, far_coll = 0.0;
	gtime_t t = epoch2time((double[]){2026, 7, 29, 12, 0, 0});
	int k, s, n, i;

	memset(aux, 0, sizeof aux);	/* a field added later must start at 0 */
	printf("least-squares collocation (D70)\n");
	ecef2pos(MASTER, mpos);
	for (k = 0; k < NAUX; k++) {
		enu_to_ecef(mpos, AE[k] * 1000.0, AN[k] * 1000.0, aux[k].xyz);
		for (s = 0; s < NSAT; s++) {
			corr[k][s].sat = SATS[s];
			corr[k][s].dI = (GI[0] + SGI(s)) * AE[k] + GI[1] * AN[k];
			corr[k][s].dT = 0.0;
			corr[k][s].valid = 1;
			sdI[k][s] = sdT[k][s] = 0.002;
		}
		aux[k].corr = corr[k];
		aux[k].sd_I = sdI[k];
		aux[k].sd_T = sdT[k];
		aux[k].ncorr = NSAT;
	}
	vbs_interp_default_opt(&opt);
	opt.method = VBS_INTERP_METHOD_COLLOCATION;

	enu_to_ecef(mpos, 7000.0, -4000.0, vbs);
	n = vbs_interp_epoch(t, MASTER, vbs, aux, NAUX, &opt, NULL, out, NULL, NSAT * 2);
	for (i = 0; i < n; i++)
		if (out[i].sat == SATS[1])
			inside = fabs(out[i].dI - truth_dI(7.0, -4.0, 1));
	printf("  inside the cluster: %.2f mm from the planar truth\n", 1000.0 * inside);
	check(inside < 0.005, "collocation tracks a planar field inside the cluster (%.2f mm)",
	    1000.0 * inside);

	enu_to_ecef(mpos, 75000.0, 0.0, vbs);
	n = vbs_interp_epoch(t, MASTER, vbs, aux, NAUX, &opt, NULL, out, NULL, NSAT * 2);
	for (i = 0; i < n; i++)
		if (out[i].sat == SATS[1])
			far_coll = fabs(out[i].dI);
	opt.method = VBS_INTERP_METHOD_PLANE;
	n = vbs_interp_epoch(t, MASTER, vbs, aux, NAUX, &opt, NULL, out, NULL, NSAT * 2);
	for (i = 0; i < n; i++)
		if (out[i].sat == SATS[1])
			far_plane = fabs(out[i].dI);
	printf("  60 km outside the hull: plane %.1f mm, collocation %.1f mm\n",
	    1000.0 * far_plane, 1000.0 * far_coll);
	check(far_coll < far_plane,
	    "collocation stops extrapolating where the plane runs away");
}

/* --- "2 + 1 backup" and the fit floor ---------------------------------------------------- */

/*
 * Two primaries and a reserve around the master, on the planar field -- except the reserve,
 * whose corrections are deliberately 5 cm off the plane, so that any use of it shows in the
 * output. `drop` removes one primary's correction for satellite index `drop_s` (-1: none);
 * `with_backup` includes the reserve at all; `datum` selects the auxiliary datum mode.
 */
static int run_backup(int with_backup, int drop_s, int datum, int method, struct vbs_corr *out,
	struct vbs_interp_quality *qual) {
	static const double AE[3] = {-15.0, 15.0, 15.0}, AN[3] = {-15.0, -15.0, 15.0};
	struct vbs_interp_aux aux[3];
	static struct vbs_corr corr[3][NSAT];
	static double sdI[3][NSAT], sdT[3][NSAT];
	struct vbs_interp_opt popt;
	double mpos[3], vbs[3];
	gtime_t t = epoch2time((double[]){2026, 7, 29, 12, 0, 0});
	int k, s, n = with_backup ? 3 : 2;

	memset(aux, 0, sizeof aux);
	ecef2pos(MASTER, mpos);
	for (k = 0; k < 3; k++) {
		enu_to_ecef(mpos, AE[k] * 1000.0, AN[k] * 1000.0, aux[k].xyz);
		for (s = 0; s < NSAT; s++) {
			corr[k][s].sat = SATS[s];
			corr[k][s].dI = (GI[0] + SGI(s)) * AE[k] + GI[1] * AN[k] + (k == 2 ? 0.05 : 0.0);
			corr[k][s].dT = (GT[0] + SGT(s)) * AE[k] + GT[1] * AN[k];
			corr[k][s].valid = !(k == 1 && s == drop_s);
			sdI[k][s] = sdT[k][s] = 0.004;
		}
		aux[k].corr = corr[k];
		aux[k].sd_I = sdI[k];
		aux[k].sd_T = sdT[k];
		aux[k].ncorr = NSAT;
		aux[k].backup = k == 2;
	}
	enu_to_ecef(mpos, 3000.0, -2000.0, vbs);
	vbs_interp_default_opt(&popt);
	popt.method = method;
	popt.datum = datum;
	return vbs_interp_epoch(t, MASTER, vbs, aux, n, &popt, NULL, out, qual, NSAT * 2);
}

static const struct vbs_corr *find_sat(const struct vbs_corr *c, int n, int sat) {
	int i;

	for (i = 0; i < n; i++)
		if (c[i].sat == sat)
			return &c[i];
	return NULL;
}

static void test_backup(void) {
	struct vbs_corr with[NSAT * 2], without[NSAT * 2];
	struct vbs_interp_quality qw[NSAT * 2], qo[NSAT * 2];
	struct vbs_interp_opt def;
	int nw, no, i, same, used, d, m;
	static const int METHODS[2] = {VBS_INTERP_METHOD_PLANE, VBS_INTERP_METHOD_COLLOCATION};

	vbs_interp_default_opt(&def);
	check(def.min_aux == 2, "the fit floor is 2, separate from the admission floor of 3 (%d)",
	    def.min_aux);

	/* A reserve that is not needed must change nothing, in either datum mode or method. */
	for (m = 0; m < 2; m++)
		for (d = 0; d < 2; d++) {
			int datum = d ? VBS_INTERP_DATUM_COMMON : VBS_INTERP_DATUM_OFF;

			nw = run_backup(1, -1, datum, METHODS[m], with, qw);
			no = run_backup(0, -1, datum, METHODS[m], without, qo);
			same = nw == no;
			used = 0;
			for (i = 0; same && i < nw; i++) {
				const struct vbs_corr *o = find_sat(without, no, with[i].sat);

				same = o != NULL && o->dI == with[i].dI && o->dT == with[i].dT &&
				    o->valid == with[i].valid;
				used |= qw[i].backup_used;
			}
			check(same && !used, "%s, datum %s: an unneeded reserve changes no bit",
			    m ? "collocation" : "plane", d ? "COMMON" : "OFF");
		}

	/*
	 * One primary loses satellite index 1: it now has one primary, below the floor of 2, so
	 * the reserve is fitted for that satellite only -- and its 5 cm shows there, nowhere else.
	 */
	nw = run_backup(1, 1, VBS_INTERP_DATUM_OFF, VBS_INTERP_METHOD_PLANE, with, qw);
	no = run_backup(0, -1, VBS_INTERP_DATUM_OFF, VBS_INTERP_METHOD_PLANE, without, qo);
	{
		int s1 = -1, others_clean = 1;

		for (i = 0; i < nw; i++) {
			const struct vbs_corr *o = find_sat(without, no, with[i].sat);

			if (with[i].sat == SATS[1])
				s1 = i;
			else if (o == NULL || o->dI != with[i].dI || qw[i].backup_used)
				others_clean = 0;
		}
		check(s1 >= 0 && qw[s1].backup_used,
		    "the satellite one primary lost is still emitted, through the reserve");
		check(s1 >= 0 && fabs(with[s1].dI - find_sat(without, no, SATS[1])->dI) > 0.005,
		    "... and the reserve's offset reaches that satellite (%.1f mm)",
		    s1 >= 0 ? 1000.0 * (with[s1].dI - find_sat(without, no, SATS[1])->dI) : 0.0);
		check(others_clean, "the other satellites are untouched and did not use it");
	}

	/* Without a reserve, a satellite down to one primary is withheld (floor 2). */
	no = run_backup(0, 1, VBS_INTERP_DATUM_OFF, VBS_INTERP_METHOD_PLANE, without, qo);
	check(find_sat(without, no, SATS[1]) == NULL,
	    "with no reserve, a satellite carried by one primary is not synthesized");
}


/*
 * The ionosphere consistency gate (vbs_interp_opt.iono_gate_m, D165). One auxiliary's estimate
 * of one satellite sits 30 cm off an otherwise exact planar field -- the shape a code-levelled
 * ionosphere leaves (D132, D164). Without the gate the offset reaches the VBS; with it that
 * auxiliary leaves that satellite's fit and the VBS gets the field; with only three auxiliaries
 * there is nothing to drop and the satellite's dispersive correction is withheld instead; at 0
 * the gate changes no bit.
 */
static int run_gate(int naux, double gate, double outlier, struct vbs_corr *out) {
	static const double AE[5] = {-15.0, 15.0, 15.0, -15.0, 25.0};
	static const double AN[5] = {-15.0, -15.0, 15.0, 15.0, 2.0};
	struct vbs_interp_aux aux[5];
	static struct vbs_corr corr[5][NSAT];
	static double sdI[5][NSAT], sdT[5][NSAT];
	struct vbs_interp_opt opt;
	double mpos[3], vbs[3];
	gtime_t t = epoch2time((double[]){2026, 7, 29, 12, 0, 0});
	int k, s;

	memset(aux, 0, sizeof aux);
	ecef2pos(MASTER, mpos);
	for (k = 0; k < naux; k++) {
		enu_to_ecef(mpos, AE[k] * 1000.0, AN[k] * 1000.0, aux[k].xyz);
		for (s = 0; s < NSAT; s++) {
			corr[k][s].sat = SATS[s];
			corr[k][s].dI = (GI[0] + SGI(s)) * AE[k] + GI[1] * AN[k] +
			    (k == 2 && s == 0 ? outlier : 0.0);
			corr[k][s].dT = 0.0;
			corr[k][s].valid = 1;
			sdI[k][s] = sdT[k][s] = 0.004;
		}
		aux[k].corr = corr[k];
		aux[k].sd_I = sdI[k];
		aux[k].sd_T = sdT[k];
		aux[k].ncorr = NSAT;
	}
	enu_to_ecef(mpos, 7000.0, -4000.0, vbs);
	vbs_interp_default_opt(&opt);
	opt.method = VBS_INTERP_METHOD_PLANE;
	opt.iono_gate_m = gate;
	return vbs_interp_epoch(t, MASTER, vbs, aux, naux, &opt, NULL, out, NULL, NSAT * 2);
}

static void test_gate(void) {
	struct vbs_corr a[NSAT * 2], b[NSAT * 2];
	double truth = truth_dI(7.0, -4.0, 0);
	long chk, drop, wh;
	int na, nb, i, same;
	const struct vbs_corr *s0;

	printf("ionosphere consistency gate (D165)\n");
	{
		struct vbs_interp_opt def;

		vbs_interp_default_opt(&def);
		check(def.iono_gate_m == 0.0, "off by default");
	}
	na = run_gate(5, 0.0, 0.30, a);
	s0 = find_sat(a, na, SATS[0]);
	check(s0 != NULL && fabs(s0->dI - truth) > 0.02,
	    "without the gate a 30 cm outlier reaches the VBS (%.1f mm)",
	    s0 ? 1000.0 * (s0->dI - truth) : 0.0);

	vbs_interp_gate_reset();
	nb = run_gate(5, 0.05, 0.30, b);
	s0 = find_sat(b, nb, SATS[0]);
	vbs_interp_gate_stats(&chk, &drop, &wh);
	check(s0 != NULL && fabs(s0->dI - truth) < 1e-4,
	    "with it, five auxiliaries: the outlier is dropped and the VBS gets the field (%.3f mm)",
	    s0 ? 1000.0 * (s0->dI - truth) : 0.0);
	check(drop == 1 && wh == 0 && chk == NSAT,
	    "... one auxiliary dropped, nothing withheld, every satellite checked (%ld/%ld/%ld)",
	    drop, wh, chk);
	for (i = 1; i < NSAT; i++) {
		const struct vbs_corr *x = find_sat(b, nb, SATS[i]);

		check(x != NULL && fabs(x->dI - truth_dI(7.0, -4.0, i)) < 1e-6,
		    "... satellite %d, consistent, is untouched", SATS[i]);
	}

	vbs_interp_gate_reset();
	nb = run_gate(3, 0.05, 0.30, b);
	s0 = find_sat(b, nb, SATS[0]);
	vbs_interp_gate_stats(&chk, &drop, &wh);
	check(s0 != NULL && s0->dI == 0.0 && wh == 1 && drop == 0,
	    "three auxiliaries that disagree: emitted with no dispersive correction (%ld withheld)", wh);

	na = run_gate(5, 0.0, 0.0, a);
	nb = run_gate(5, 0.05, 0.0, b);
	same = na == nb;
	for (i = 0; same && i < na; i++)
		same = a[i].sat == b[i].sat && a[i].dI == b[i].dI && a[i].dT == b[i].dT;
	check(same, "a consistent field passes the gate bit-identical");
}

int main(void) {
	printf("planar field reproduction (spec §4.5 baseline method)\n");
	test_planar("VBS at cluster centre", 0.0, 0.0, 1e-6);
	test_planar("VBS inside the hull", 7.0, -4.0, 1e-6);
	test_planar("VBS 60 km outside", 75.0, 0.0, 1e-6);
	printf("leverage as the extrapolation metric\n");
	test_leverage();
	printf("the master is the datum\n");
	test_at_master();
	printf("variance weighting\n");
	test_weighting();
	test_model_term_sign();
	test_anchor();
	test_collocation();
	printf("the auxiliary ionosphere datum (D55)\n");
	test_datum();
	printf("\"2 + 1 backup\": the reserve and the fit floor\n");
	test_backup();
	test_gate();
	printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
	return failures ? 1 : 0;
}

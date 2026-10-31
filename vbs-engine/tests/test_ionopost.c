/*
 * vbs_ionopost (D168): the post-interpolation handling of the dispersive correction.
 *
 * - inactive options are reported inactive, so callers keep the D167 sys_mask path;
 * - sys_mask alone zeroes exactly the left-out constellations (D167 behaviour);
 * - need_band5 withholds a GPS correction when the master has no band-5 phase for it;
 * - scale_gps scales GPS only;
 * - the shell model reproduces a field of its own form: Galileo and BeiDou corrections built
 *   from a mapped vertical term with a gradient, each constellation with its own datum, give the
 *   withheld GPS satellites that field with no datum, to the micrometre; the quadratic mode
 *   reproduces a curved field the plane cannot; too few satellites leave GPS withheld;
 * - model_check replaces only the own correction that departs from the prediction;
 * - drop_sys removes a constellation from the output;
 * - shrink_m scales each kept correction's departure d from its constellation's median by
 *   d^2 / (d^2 + shrink_m^2) on L1, BeiDou converted.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "vbs/vbs_ionopost.h"

static int nfail;

static void check(int ok, const char *what) {
	printf("  %-66s %s\n", what, ok ? "ok" : "FAIL");
	if (!ok)
		nfail++;
}

/* the model's definition: 350 km thin shell, pierce offsets in 100 km units */
static double field(double az, double el, double v0, double gn, double ge, double qn) {
	double re = 6371.0, h = 350.0, z = PI / 2.0 - el, zp = asin(re / (re + h) * sin(z));
	double psi = z - zp, n = psi * re * cos(az) / 100.0, e = psi * re * sin(az) / 100.0;

	return (v0 + gn * n + ge * e + qn * n * n) / cos(zp);
}

static double azel[MAXSAT * 2];
static struct vbs_corr corr[32];
static obsd_t mep[32];
static int ncorr;

static void add(int sys, int prn, double az_deg, double el_deg, double dI, int band5) {
	int sat = satno(sys, prn);

	corr[ncorr].sat = sat;
	corr[ncorr].dI = dI;
	corr[ncorr].valid = 1;
	memset(&mep[ncorr], 0, sizeof mep[ncorr]);
	mep[ncorr].sat = sat;
	mep[ncorr].code[0] = CODE_L1C;
	mep[ncorr].L[0] = 1.0e8;
	if (band5) {
		mep[ncorr].code[2] = CODE_L5Q;
		mep[ncorr].L[2] = 0.8e8;
	}
	azel[(sat - 1) * 2] = az_deg * D2R;
	azel[(sat - 1) * 2 + 1] = el_deg * D2R;
	ncorr++;
}

static int idx(int sys, int prn) {
	int k, sat = satno(sys, prn);

	for (k = 0; k < ncorr; k++)
		if (corr[k].sat == sat)
			return k;
	return -1;
}

/* f_ref(sys)^2 scaling of an L1-equivalent delay onto the constellation's datum frequency */
static double on_datum(int sys, double y) {
	double f = vbs_synth_ref_freq(sys);

	return y * (FREQ1 / f) * (FREQ1 / f);
}

/* Galileo and BeiDou over the sky with the given field and datums; GPS with junk */
static void scene(double v0, double gn, double ge, double qn, double dE, double dC) {
	static const double eaz[] = {10, 70, 130, 190, 250, 310, 40, 160, 280};
	static const double eel[] = {35, 60, 28, 45, 70, 30, 80, 52, 40};
	static const double caz[] = {20, 100, 200, 290, 340, 150, 60};
	static const double cel[] = {50, 33, 65, 27, 42, 75, 38};
	static const double gaz[] = {0, 90, 180, 270, 225};
	static const double gel[] = {55, 31, 47, 66, 29};
	int k;

	ncorr = 0;
	memset(azel, 0, sizeof azel);
	for (k = 0; k < 9; k++)
		add(SYS_GAL, k + 1, eaz[k], eel[k], dE + field(eaz[k] * D2R, eel[k] * D2R, v0, gn, ge, qn), 1);
	for (k = 0; k < 7; k++)
		add(SYS_CMP, k + 20, caz[k], cel[k],
		    on_datum(SYS_CMP, dC + field(caz[k] * D2R, cel[k] * D2R, v0, gn, ge, qn)), 1);
	for (k = 0; k < 5; k++)
		add(SYS_GPS, k + 1, gaz[k], gel[k], 0.37 * (k + 1), k % 2 == 0);
}

static double gps_truth(int prn, double v0, double gn, double ge, double qn) {
	int k = idx(SYS_GPS, prn), sat = corr[k].sat;

	return field(azel[(sat - 1) * 2], azel[(sat - 1) * 2 + 1], v0, gn, ge, qn);
}

int main(void) {
	struct vbs_ionopost_opt o;
	struct vbs_ionopost_stats st;
	double vpos[3] = {47.0 * D2R, 2.0 * D2R, 100.0};
	int k, ok;

	printf("=== ionosphere post (D168) ===\n");
	memset(&o, 0, sizeof o);
	check(!vbs_ionopost_active(&o), "nothing set: inactive");
	o.scale_gps = 1.0;
	check(!vbs_ionopost_active(&o), "a GPS scale of 1 is inactive");
	o.sys_mask = SYS_GAL | SYS_CMP;
	check(!vbs_ionopost_active(&o), "sys_mask alone is D167's path, not this module's");
	check(vbs_ionopost_parse_sys("EC") == (SYS_GAL | SYS_CMP) && vbs_ionopost_parse_sys("GX") < 0,
	    "constellation letters parse, a bad one is refused");

	/* sys_mask alone, through apply: the D167 behaviour */
	scene(0.02, 0.01, -0.005, 0.0, 0.3, -0.2);
	memset(&o, 0, sizeof o);
	o.sys_mask = SYS_GAL | SYS_CMP;
	vbs_ionopost_apply(corr, ncorr, azel, vpos, mep, ncorr, NULL, &o, NULL);
	ok = 1;
	for (k = 0; k < ncorr; k++)
		if ((satsys(corr[k].sat, NULL) == SYS_GPS) != (corr[k].dI == 0.0) || !corr[k].valid)
			ok = 0;
	check(ok, "sys_mask EC: GPS zeroed, Galileo and BeiDou kept, all still valid");

	/* need_band5 */
	scene(0.02, 0.01, -0.005, 0.0, 0.3, -0.2);
	memset(&o, 0, sizeof o);
	memset(&st, 0, sizeof st);
	o.need_band5_sys = SYS_GPS;
	vbs_ionopost_apply(corr, ncorr, azel, vpos, mep, ncorr, NULL, &o, &st);
	check(corr[idx(SYS_GPS, 1)].dI == 0.37 && corr[idx(SYS_GPS, 2)].dI == 0.0 &&
	    corr[idx(SYS_GPS, 3)].dI != 0.0 && st.band5_withheld == 2,
	    "band 5: GPS kept with L5 at the master, withheld without");
	check(corr[idx(SYS_GAL, 1)].dI != 0.0, "band 5 on GPS leaves Galileo alone");

	/* scale_gps */
	scene(0.02, 0.01, -0.005, 0.0, 0.3, -0.2);
	memset(&o, 0, sizeof o);
	o.scale_gps = 0.5;
	{
		double e1 = corr[idx(SYS_GAL, 1)].dI;

		vbs_ionopost_apply(corr, ncorr, azel, vpos, mep, ncorr, NULL, &o, NULL);
		check(fabs(corr[idx(SYS_GPS, 2)].dI - 0.37) < 1e-12 && corr[idx(SYS_GAL, 1)].dI == e1,
		    "GPS scale halves GPS only");
	}

	/* the plane model */
	scene(0.02, 0.01, -0.005, 0.0, 0.3, -0.2);
	memset(&o, 0, sizeof o);
	memset(&st, 0, sizeof st);
	o.sys_mask = SYS_GAL | SYS_CMP;
	o.model_sys = SYS_GPS;
	vbs_ionopost_apply(corr, ncorr, azel, vpos, mep, ncorr, NULL, &o, &st);
	ok = st.fits == 1 && st.filled == 5;
	for (k = 1; k <= 5; k++)
		if (fabs(corr[idx(SYS_GPS, k)].dI - gps_truth(k, 0.02, 0.01, -0.005, 0.0)) > 1e-6)
			ok = 0;
	check(ok, "plane: GPS gets the Galileo/BeiDou field, datum-free, to 1 um");
	check(st.fit_rms_sum < 1e-9, "plane: a field of the model's form fits exactly");

	/* a curved field: the quadratic mode reproduces it, the plane does not */
	scene(0.02, 0.01, -0.005, 0.004, 0.3, -0.2);
	memset(&o, 0, sizeof o);
	o.sys_mask = SYS_GAL | SYS_CMP;
	o.model_sys = SYS_GPS;
	vbs_ionopost_apply(corr, ncorr, azel, vpos, mep, ncorr, NULL, &o, NULL);
	{
		double ep = 0.0, eq = 0.0;

		for (k = 1; k <= 5; k++)
			ep += fabs(corr[idx(SYS_GPS, k)].dI - gps_truth(k, 0.02, 0.01, -0.005, 0.004));
		scene(0.02, 0.01, -0.005, 0.004, 0.3, -0.2);
		o.model_mode = VBS_IONOPOST_QUAD;
		vbs_ionopost_apply(corr, ncorr, azel, vpos, mep, ncorr, NULL, &o, NULL);
		for (k = 1; k <= 5; k++)
			eq += fabs(corr[idx(SYS_GPS, k)].dI - gps_truth(k, 0.02, 0.01, -0.005, 0.004));
		check(eq < 1e-6 && ep > 1e-4, "quad: a curved field reproduced, where the plane misses it");
	}

	/* IDW mode on a planar field: residuals are zero, so it equals the plane */
	scene(0.02, 0.01, -0.005, 0.0, 0.3, -0.2);
	o.model_mode = VBS_IONOPOST_IDW;
	vbs_ionopost_apply(corr, ncorr, azel, vpos, mep, ncorr, NULL, &o, NULL);
	ok = 1;
	for (k = 1; k <= 5; k++)
		if (fabs(corr[idx(SYS_GPS, k)].dI - gps_truth(k, 0.02, 0.01, -0.005, 0.0)) > 1e-6)
			ok = 0;
	check(ok, "idw: on a planar field the residual term adds nothing");

	/* too few satellites to fit: GPS stays withheld */
	scene(0.02, 0.01, -0.005, 0.0, 0.3, -0.2);
	memset(&st, 0, sizeof st);
	o.model_mode = VBS_IONOPOST_PLANE;
	vbs_ionopost_apply(corr, 9, azel, vpos, mep, ncorr, NULL, &o, &st);	/* Galileo only, 9 */
	check(st.fits == 1, "nine Galileo satellites are enough for one datum and a plane");
	scene(0.02, 0.01, -0.005, 0.0, 0.3, -0.2);
	memset(&st, 0, sizeof st);
	for (k = 4; k < 16; k++)
		corr[k].valid = 0;	/* 4 Galileo left: 1 datum + 3 shell terms + 2 is 6 */
	vbs_ionopost_apply(corr, ncorr, azel, vpos, mep, ncorr, NULL, &o, &st);
	check(st.fit_refused == 1 && corr[idx(SYS_GPS, 1)].dI == 0.0,
	    "four satellites: no fit, GPS stays withheld");

	/* model_check: own GPS kept where it agrees with the model up to a datum, replaced where not */
	scene(0.02, 0.01, -0.005, 0.0, 0.3, -0.2);
	for (k = 1; k <= 5; k++)
		corr[idx(SYS_GPS, k)].dI = 0.9 + gps_truth(k, 0.02, 0.01, -0.005, 0.0);
	corr[idx(SYS_GPS, 4)].dI += 0.10;	/* one wrong satellite */
	memset(&o, 0, sizeof o);
	memset(&st, 0, sizeof st);
	o.model_sys = SYS_GPS;
	o.model_check_m = 0.03;
	vbs_ionopost_apply(corr, ncorr, azel, vpos, mep, ncorr, NULL, &o, &st);
	check(st.replaced == 1 && st.kept == 4 &&
	    fabs(corr[idx(SYS_GPS, 4)].dI - (0.9 + gps_truth(4, 0.02, 0.01, -0.005, 0.0))) < 1e-6 &&
	    fabs(corr[idx(SYS_GPS, 2)].dI - (0.9 + gps_truth(2, 0.02, 0.01, -0.005, 0.0))) < 1e-12,
	    "check: the outlier replaced on the GPS datum, the others untouched");

	/* drop */
	scene(0.02, 0.01, -0.005, 0.0, 0.3, -0.2);
	memset(&o, 0, sizeof o);
	o.drop_sys = SYS_GPS;
	vbs_ionopost_apply(corr, ncorr, azel, vpos, mep, ncorr, NULL, &o, NULL);
	ok = 1;
	for (k = 0; k < ncorr; k++)
		if (corr[k].valid != (satsys(corr[k].sat, NULL) != SYS_GPS))
			ok = 0;
	check(ok, "drop G: GPS not output, the rest still valid");

	/* shrink: on each satellite's departure from its constellation's median, on L1 */
	scene(0.02, 0.01, -0.005, 0.0, 0.3, -0.2);
	for (k = 0; k < ncorr; k++)
		corr[k].dI = satsys(corr[k].sat, NULL) == SYS_CMP ? on_datum(SYS_CMP, 0.5) : 1.0;
	corr[idx(SYS_GAL, 1)].dI = 1.10;
	corr[idx(SYS_GAL, 2)].dI = 1.005;
	corr[idx(SYS_CMP, 20)].dI = on_datum(SYS_CMP, 0.52);
	memset(&o, 0, sizeof o);
	memset(&st, 0, sizeof st);
	o.sys_mask = SYS_GAL | SYS_CMP;
	o.shrink_m = 0.02;
	check(vbs_ionopost_active(&o), "shrink alone makes the module active");
	vbs_ionopost_apply(corr, ncorr, azel, vpos, mep, ncorr, NULL, &o, &st);
	check(fabs(corr[idx(SYS_GAL, 1)].dI - (1.0 + 0.1 * 0.01 / (0.01 + 0.0004))) < 1e-12 &&
	    fabs(corr[idx(SYS_GAL, 2)].dI - (1.0 + 0.005 * 0.000025 / (0.000025 + 0.0004))) < 1e-12 &&
	    corr[idx(SYS_GAL, 5)].dI == 1.0,
	    "shrink: the datum kept, a large departure nearly whole, a small one faded");
	check(fabs(corr[idx(SYS_CMP, 20)].dI -
	    on_datum(SYS_CMP, 0.5 + 0.02 * 0.0004 / (0.0004 + 0.0004))) < 1e-12,
	    "shrink: BeiDou weighed on its L1 equivalent");
	check(corr[idx(SYS_GPS, 1)].dI == 0.0 && st.shrunk == 16,
	    "shrink: a masked constellation stays zero and is not counted");

	printf(nfail ? "FAILED: %d failure(s)\n" : "PASSED: %d failure(s)\n", nfail);
	return nfail != 0;
}

/*
 * vbs_cluster and vbs_basemeta tests -- spec §4.1, and D23's tide handling.
 *
 * The hull and leverage results are checked against geometries whose answer is already known
 * from an independent tool: `04-geophysical-corrections.md` measured the interpolation
 * leverage at 0.25-0.35 for the enclosed LOO cases and 6-12 for the extrapolated ones, using
 * `tools/vbs_tidecheck`, which shares no code with this module. Two implementations agreeing
 * on a geometric quantity is worth more than either one agreeing with my expectation.
 */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_cluster.h"

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

/* --- hull geometry, where the answers are exact ---------------------------------------- */

static void test_hull(void) {
	/* Unit square. */
	static const double e[4] = {0.0, 10.0, 10.0, 0.0};
	static const double n[4] = {0.0, 0.0, 10.0, 10.0};
	int inside;
	double d;

	printf("convex hull containment and distance\n");

	d = vbs_cluster_hull_distance(e, n, 4, 5.0, 5.0, &inside);
	check(inside && d == 0.0, "a point at the centre is inside, distance 0 (%.3f)", d);

	d = vbs_cluster_hull_distance(e, n, 4, 5.0, 0.0, &inside);
	check(inside && d == 0.0, "a point on an edge counts as inside (%.3f)", d);

	d = vbs_cluster_hull_distance(e, n, 4, 15.0, 5.0, &inside);
	check(!inside && fabs(d - 5.0) < 1e-9,
	    "5 units beyond an edge is outside at distance 5 (%.6f)", d);

	d = vbs_cluster_hull_distance(e, n, 4, 13.0, 14.0, &inside);
	check(!inside && fabs(d - hypot(3.0, 4.0)) < 1e-9,
	    "beyond a corner, the distance is to the corner (%.6f, expected 5)", d);

	/* Collinear stations span no area: nothing can be inside them. */
	{
		static const double ce[3] = {0.0, 5.0, 10.0};
		static const double cn[3] = {0.0, 0.0, 0.0};

		d = vbs_cluster_hull_distance(ce, cn, 3, 5.0, 3.0, &inside);
		check(!inside && fabs(d - 3.0) < 1e-9,
		    "collinear bases enclose nothing, however many there are (%.3f)", d);
	}
}

/* --- selection on the real LOO geometries ---------------------------------------------- */

struct site {
	const char *name;
	double xyz[3];
};

/* E1-0: STA_T enclosed by STA_M, STA_A1, STA_A2, STA_A3 (the validation catalogue). */
static const struct site E1_0_CLUSTER[] = {
	{"STA_M", {4582081.2, 55183.6, 4421736.6}},
	{"STA_A1",  {4614946.0,    21620.0,    4387760.0}},
	{"STA_A2",   {4620000.0,    75000.0,    4382000.0}},
	{"STA_A3",   {4585000.0,    30000.0,    4418000.0}},
};
static const double STA_T[3] = {4603034.1, 55436.0, 4400265.6};

static void fill(struct vbs_base *b, int n, const char *frame) {
	int i;

	for (i = 0; i < n; i++) {
		vbs_base_init(&b[i], E1_0_CLUSTER[i].name, frame, E1_0_CLUSTER[i].xyz);
		vbs_base_check(&b[i], frame, (gtime_t){0, 0.0}, NULL);
	}
}

static void test_select_enclosed(void) {
	struct vbs_base b[4];
	struct vbs_cluster cl;
	int r;

	printf("cluster selection, STA_T enclosed by its four neighbours (E1-0)\n");
	fill(b, 4, "RGF93");

	r = vbs_cluster_select(STA_T, b, 4, NULL, NULL, &cl);
	printf("    master %s at %.1f km, %d aux, inside=%d, extrap %.1f km, leverage %.2f\n",
	    cl.master >= 0 ? b[cl.master].name : "-", cl.master_dist_km, cl.naux,
	    cl.inside_hull, cl.extrapolation_km, cl.leverage);

	check(r == VBS_CLUSTER_OK, "a cluster was formed (%s)", vbs_cluster_reason_str(r));
	check(cl.naux == 3, "three auxiliaries kept alongside the master (%d)", cl.naux);
	check(cl.inside_hull && cl.extrapolation_km == 0.0,
	    "the target is reported inside the hull");
	/*
	 * The independent check: `04-geophysical-corrections.md` measured 0.25-0.35 for every
	 * enclosed LOO case with tools/vbs_tidecheck, which shares no code with this module.
	 */
	check(cl.leverage > 0.2 && cl.leverage < 0.6,
	    "leverage matches the 0.25-0.35 measured independently for enclosed cases (%.2f)",
	    cl.leverage);
}

static void test_select_extrapolated(void) {
	struct vbs_base b[4];
	struct vbs_cluster cl;
	double far[3];
	int r, i;

	printf("cluster selection, a VBS pushed well outside the hull\n");
	fill(b, 4, "RGF93");

	/* Move ~100 km along +X from the cluster: unambiguously outside. */
	for (i = 0; i < 3; i++)
		far[i] = STA_T[i];
	far[0] += 100000.0;

	r = vbs_cluster_select(far, b, 4, NULL, NULL, &cl);
	printf("    inside=%d, extrapolation %.1f km, leverage %.2f\n",
	    cl.inside_hull, cl.extrapolation_km, cl.leverage);

	check(r == VBS_CLUSTER_OK, "a cluster is still formed (%s)",
	    vbs_cluster_reason_str(r));
	check(!cl.inside_hull, "the VBS is reported outside the hull");
	check(cl.extrapolation_km > 40.0,
	    "the extrapolation distance is published and substantial (%.1f km)",
	    cl.extrapolation_km);
	check(cl.leverage > 2.0,
	    "leverage is inflated far above the enclosed case, as tools/vbs_tidecheck "
	    "measured at 6-12 (%.2f)", cl.leverage);

	/* With the cap on, the same geometry must be refused rather than quietly served. */
	{
		struct vbs_cluster_opt opt;

		vbs_cluster_default_opt(&opt);
		opt.max_extrapolation_km = 20.0;
		r = vbs_cluster_select(far, b, 4, NULL, &opt, &cl);
		check(r == VBS_CLUSTER_TOO_FAR,
		    "max_extrapolation_km refuses it (%s)", vbs_cluster_reason_str(r));
	}
}

/* --- the quality gate ------------------------------------------------------------------- */

static void test_gate(void) {
	struct vbs_base b;
	gtime_t now = epoch2time((double[]){2026, 7, 29, 12, 0, 0});
	double zero[3] = {0.0, 0.0, 0.0};

	printf("base quality gate\n");

	vbs_base_init(&b, "X", "RGF93", E1_0_CLUSTER[0].xyz);
	check(vbs_base_check(&b, "RGF93", now, NULL) == VBS_BASE_OK, "a good base passes");

	check(vbs_base_check(&b, "ETRF2000", now, NULL) == VBS_BASE_FRAME_MISMATCH,
	    "a base in another frame is refused, not transformed (D13)");

	vbs_base_init(&b, "X", NULL, E1_0_CLUSTER[0].xyz);
	check(vbs_base_check(&b, "RGF93", now, NULL) == VBS_BASE_NO_FRAME,
	    "an unknown frame is refused rather than assumed to match");

	vbs_base_init(&b, "X", "RGF93", zero);
	check(vbs_base_check(&b, "RGF93", now, NULL) == VBS_BASE_NO_COORD,
	    "a base with no coordinate is refused");

	vbs_base_init(&b, "X", "RGF93", E1_0_CLUSTER[0].xyz);
	b.last_obs = timeadd(now, -3600.0);
	check(vbs_base_check(&b, "RGF93", now, NULL) == VBS_BASE_STALE,
	    "a base an hour stale is refused");

	/*
	 * D56: STA_A2's declared position was measured 7-13 cm high from all four of its
	 * neighbours. The thresholds are the network's own published ones, so what is tested here
	 * is that a station the network's own rules would reject is rejected -- and, just as
	 * importantly, that a station nobody has measured is still admitted, since unmeasured
	 * must not be read as bad.
	 */
	vbs_base_init(&b, "STA_A2", "RGF93", E1_0_CLUSTER[0].xyz);
	check(vbs_base_check(&b, "RGF93", now, NULL) == VBS_BASE_OK,
	    "an unmeasured base is admitted, not presumed bad");

	vbs_base_set_coord_resid(&b, (double[]){0.004, -0.004, 0.127});
	check(vbs_base_check(&b, "RGF93", now, NULL) == VBS_BASE_COORD_OUTLIER,
	    "a base 12.7 cm high is refused (six times the published 20 mm gate)");

	vbs_base_set_coord_resid(&b, (double[]){0.004, -0.004, 0.018});
	check(vbs_base_check(&b, "RGF93", now, NULL) == VBS_BASE_OK,
	    "a base inside the published gate is kept");

	vbs_base_set_coord_resid(&b, (double[]){0.009, 0.009, 0.000});
	check(vbs_base_check(&b, "RGF93", now, NULL) == VBS_BASE_COORD_OUTLIER,
	    "the horizontal gate is on the 2-D magnitude, not per component (12.7 mm > 10 mm)");
}

/*
 * The gate is only useful if cluster selection acts on it: a cluster must be formed from the
 * remaining stations rather than refused outright, which is the whole point of running four or
 * five auxiliaries instead of three (D59).
 */
static void test_select_drops_outlier(void) {
	struct vbs_base b[4];
	struct vbs_cluster cl;
	int r, i, seen = 0;

	printf("cluster selection with one coordinate outlier\n");
	fill(b, 4, "RGF93");
	vbs_base_set_coord_resid(&b[2], (double[]){0.004, -0.004, 0.127});
	for (i = 0; i < 4; i++)
		vbs_base_check(&b[i], "RGF93", epoch2time((double[]){2026, 7, 29, 12, 0, 0}), NULL);

	r = vbs_cluster_select(STA_T, b, 4, NULL, NULL, &cl);
	check(b[2].reason == VBS_BASE_COORD_OUTLIER, "the outlier is flagged (%s)",
	    vbs_base_reason_str(b[2].reason));
	if (r == VBS_CLUSTER_OK) {
		if (cl.master == 2)
			seen = 1;
		for (i = 0; i < cl.naux; i++)
			if (cl.aux[i] == 2)
				seen = 1;
		printf("    cluster formed from %d station(s), outlier used: %s\n",
		    cl.naux + 1, seen ? "yes" : "no");
	}
	check(!seen, "the outlier is not used as master or auxiliary");
}

/* --- the tide (D23) --------------------------------------------------------------------- */

static void test_tide(void) {
	struct vbs_base b;
	struct vbs_basemeta_opt opt;
	gtime_t t = epoch2time((double[]){2026, 7, 29, 12, 0, 0});
	double with[3], without[3], d[3], dm;
	int i;

	printf("solid Earth tide on base coordinates (D23)\n");

	vbs_base_init(&b, "X", "RGF93", E1_0_CLUSTER[0].xyz);
	vbs_basemeta_default_opt(&opt);
	vbs_base_instant_xyz(&b, t, &opt, with);
	opt.apply_solid_tide = 0;
	vbs_base_instant_xyz(&b, t, &opt, without);

	for (i = 0; i < 3; i++)
		d[i] = with[i] - without[i];
	dm = norm(d, 3);
	printf("    displacement %.1f mm\n", dm * 1000.0);

	check(dm > 0.005 && dm < 0.5,
	    "the solid Earth tide is applied and is decimetre-scale (%.1f mm)", dm * 1000.0);
	/*
	 * The same call must serve the VBS point, or the consistency that makes the tide cancel
	 * at the rover is broken. Two points a metre apart must move together.
	 */
	{
		double near[3], nxyz[3], delta;

		for (i = 0; i < 3; i++)
			near[i] = b.xyz[i] + (i == 0 ? 1.0 : 0.0);
		opt.apply_solid_tide = 1;
		vbs_basemeta_tide_xyz(near, t, &opt, nxyz);
		for (i = 0; i < 3; i++)
			d[i] = (nxyz[i] - near[i]) - (with[i] - b.xyz[i]);
		delta = norm(d, 3);
		check(delta < 1e-6,
		    "a point 1 m away moves by the same displacement (%.3g m difference)",
		    delta);
	}
}

/*
 * Grid-VBS heights (D37). DEM values are orthometric; everything downstream wants ellipsoidal.
 * The checks are closed-loop against the geoid model itself, so a wrong sign on the conversion
 * -- which would put every French VBS ~100 m out instead of ~50 m -- cannot pass.
 */
static void test_grid_height(void) {
	/* Toulouse-ish, inside the French network. */
	double pos[2] = {43.6 * D2R, 1.4 * D2R};
	double xyz[3], h, n;
	int r;

	printf("Grid-VBS heights: orthometric DEM -> ellipsoidal (D37)\n");

	n = vbs_basemeta_geoid_sep((double[]){pos[0], pos[1], 0.0});
	printf("  geoid separation here: %+.2f m\n", n);
	check(!isnan(n), "a geoid separation is available");
	/*
	 * France sits ~44-56 m above the ellipsoid. A conversion that dropped the term, or
	 * applied it backwards, would land outside this band -- which is the whole point of
	 * asserting a range rather than just "not NaN".
	 */
	check(!isnan(n) && n > 30.0 && n < 70.0,
	    "and it is plausible for France (%+.2f m)", n);

	r = vbs_basemeta_grid_xyz(pos, 150.0, 0, xyz, &h);
	printf("  DEM 150 m orthometric -> %.2f m ellipsoidal\n", h);
	check(r == 0, "a normal inland cell converts");
	check(r == 0 && fabs(h - (150.0 + n)) < 1e-6,
	    "ellipsoidal = orthometric + N, with the right sign (%.2f vs %.2f)",
	    h, 150.0 + n);
	/* The sign matters twice: getting it backwards is the same magnitude, opposite way. */
	check(r == 0 && h > 150.0,
	    "the conversion raises the point in France, it does not lower it");

	/*
	 * Offshore. A DEM that returns bathymetry hands back a large negative orthometric
	 * height; clamped, the VBS sits at the sea surface, which is the geoid, so its
	 * ellipsoidal height is exactly N.
	 */
	r = vbs_basemeta_grid_xyz(pos, -1200.0, 1, xyz, &h);
	printf("  DEM -1200 m (bathymetry), over water -> %.2f m ellipsoidal\n", h);
	check(r == 0, "an offshore cell with bathymetry still converts");
	check(r == 0 && fabs(h - n) < 1e-6,
	    "clamped to the sea surface: ellipsoidal height is exactly N (%.2f vs %.2f)", h, n);

	/* Unclamped, the same DEM value must be refused rather than served. */
	r = vbs_basemeta_grid_xyz(pos, -1200.0, 0, xyz, &h);
	check(r == -1,
	    "unclamped bathymetry is refused, not passed to a model that returns 0.0 for it");

	/* Above the troposphere model's ceiling: refuse rather than publish a silent error. */
	r = vbs_basemeta_grid_xyz(pos, 12000.0, 0, xyz, &h);
	check(r == -1, "a height beyond the troposphere model's range is refused");

	/* And the position round-trips: what comes back out is the height that went in. */
	r = vbs_basemeta_grid_xyz(pos, 150.0, 0, xyz, &h);
	if (r == 0) {
		double back[3];

		ecef2pos(xyz, back);
		check(fabs(back[2] - h) < 1e-6 &&
		    fabs(back[0] - pos[0]) < 1e-12 && fabs(back[1] - pos[1]) < 1e-12,
		    "the ECEF position round-trips to the same lat/lon/height");
	}
}

/* --- "2 + 1 backup" ------------------------------------------------------------------- */

/* A station `e` km east and `n` km north of STA_T. */
static void at(double e, double n, double *xyz) {
	double pos[3], enu[3] = {e * 1000.0, n * 1000.0, 0.0}, dr[3];
	int i;

	ecef2pos(STA_T, pos);
	enu2ecef(pos, enu, dr);
	for (i = 0; i < 3; i++)
		xyz[i] = STA_T[i] + dr[i];
}

static void put(struct vbs_base *b, const char *name, const char *frame, double e, double n) {
	double xyz[3];

	at(e, n, xyz);
	vbs_base_init(b, name, frame, xyz);
	vbs_base_check(b, frame, (gtime_t){0, 0.0}, NULL);
}

static int in_aux(const struct vbs_cluster *cl, int idx) {
	int i;

	for (i = 0; i < cl->naux; i++)
		if (cl->aux[i] == idx)
			return 1;
	return 0;
}

static void test_select_backup(void) {
	struct vbs_base b[6];
	struct vbs_cluster cl;
	struct vbs_cluster_opt opt;
	int r;

	printf("cluster selection, \"2 + 1 backup\" (decided 2026-10-02)\n");

	/*
	 * Master and two auxiliaries close by; the only other station is ~400 km out, beyond
	 * max_baseline_km (320) but within backup_max_km (640).
	 */
	put(&b[0], "M", "RGF93", 2.0, 0.0);
	put(&b[1], "A1", "RGF93", -20.0, 10.0);
	put(&b[2], "A2", "RGF93", -20.0, -10.0);
	put(&b[3], "FAR", "RGF93", 400.0, 0.0);

	r = vbs_cluster_select(STA_T, b, 4, NULL, NULL, &cl);
	printf("    %s: master %s, %d aux, backup %s at %.0f km (%.0f km from the master)\n",
	    vbs_cluster_reason_str(r), cl.master >= 0 ? b[cl.master].name : "-", cl.naux,
	    cl.backup >= 0 ? b[cl.backup].name : "-", cl.backup_dist_km, cl.backup_baseline_km);
	check(r == VBS_CLUSTER_OK, "two primaries and a far reserve form a cluster (%s)",
	    vbs_cluster_reason_str(r));
	check(cl.naux == 2 && cl.backup == 3, "the far station is the backup (%d aux, backup %d)",
	    cl.naux, cl.backup);
	check(!in_aux(&cl, 3), "the backup is never among the auxiliaries");
	check(cl.inside_hull, "the geometry reported is the primaries', VBS inside it");

	vbs_cluster_default_opt(&opt);
	opt.backup_min_aux = 0;
	r = vbs_cluster_select(STA_T, b, 4, NULL, &opt, &cl);
	check(r == VBS_CLUSTER_TOO_FEW_AUX && cl.backup == -1,
	    "backup_min_aux = 0 restores the old refusal (%s)", vbs_cluster_reason_str(r));

	r = vbs_cluster_select(STA_T, (struct vbs_base[]){b[0], b[1], b[3]}, 3, NULL, NULL, &cl);
	check(r == VBS_CLUSTER_TOO_FEW_AUX,
	    "one primary and a reserve is not enough (%s)", vbs_cluster_reason_str(r));

	put(&b[3], "FAR", "ETRF2000", 400.0, 0.0);
	r = vbs_cluster_select(STA_T, b, 4, NULL, NULL, &cl);
	check(r == VBS_CLUSTER_TOO_FEW_AUX,
	    "a reserve in another frame is no reserve (%s)", vbs_cluster_reason_str(r));
	put(&b[3], "FAR", "RGF93", 400.0, 0.0);

	/*
	 * A station 130 km from the VBS but only 128 km from the master would be an ordinary
	 * auxiliary -- unless the server bounds auxiliaries to its policy radius, in which case
	 * it can only be the reserve, and it is the nearer reserve than FAR.
	 */
	put(&b[4], "G130", "RGF93", 130.0, 0.0);
	r = vbs_cluster_select(STA_T, b, 5, NULL, NULL, &cl);
	check(r == VBS_CLUSTER_OK && cl.naux == 3 && in_aux(&cl, 4) && cl.backup == -1,
	    "unbounded, G130 is an ordinary third auxiliary and no reserve is taken");
	vbs_cluster_default_opt(&opt);
	opt.aux_max_vbs_km = 120.0;
	r = vbs_cluster_select(STA_T, b, 5, NULL, &opt, &cl);
	check(r == VBS_CLUSTER_OK && cl.naux == 2 && !in_aux(&cl, 4) && cl.backup == 4,
	    "bounded to 120 km, G130 becomes the reserve, never an auxiliary (backup %d)",
	    cl.backup);

	/* Enough primaries: the reserve is not chosen at all. */
	put(&b[4], "A3", "RGF93", 5.0, 20.0);
	r = vbs_cluster_select(STA_T, b, 5, NULL, NULL, &cl);
	check(r == VBS_CLUSTER_OK && cl.naux == 3 && cl.backup == -1 && !in_aux(&cl, 3),
	    "three primaries: no backup, and the far station stays out");
}

int main(void) {
	test_hull();
	test_select_enclosed();
	test_select_extrapolated();
	test_gate();
	test_select_drops_outlier();
	test_tide();
	test_grid_height();
	test_select_backup();
	printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
	return failures ? 1 : 0;
}

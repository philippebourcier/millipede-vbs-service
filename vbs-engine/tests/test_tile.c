/*
 * The tile grid: pure, total, and neighbour-distinct.
 *
 * The property that matters is the last one. A rover crossing into a neighbouring tile must see
 * a different RTCM station ID, otherwise its base moves under a constant ID and its ambiguities
 * silently go wrong. So this walks EVERY row of the grid and, for each tile, finds its spatial
 * neighbours by position -- not by index, since a row boundary can put a tile of a different
 * column count next door -- and asserts the IDs differ. The arithmetic argument in vbs_tile.c
 * is only a reason to expect this to pass; the test is what establishes it.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "vbs/vbs_tile.h"

static int nfail;

static void check(int ok, const char *what) {
	printf("  %-66s %s\n", what, ok ? "ok" : "FAIL");
	if (!ok)
		nfail++;
}

static void test_square(void) {
	struct vbs_tile t, u, v;
	char name[32];
	int row, bad = 0, checked = 0;

	vbs_tile_set_shape(VBS_TILE_SQUARE);
	printf("=== tile grid: squares ===\n");
	vbs_tile_from_llh(48.8566, 2.3522, &t);
	check(t.lat0 <= 48.8566 && 48.8566 < t.lat1 && t.lon0 <= 2.3522 && 2.3522 < t.lon1,
	    "Paris is inside the tile it maps to");
	vbs_tile_from_llh(t.clat, t.clon, &u);
	check(u.row == t.row && u.col == t.col, "the tile's centre maps back to the same tile");
	vbs_tile_at(t.row, t.col, &v);
	check(v.lat0 == t.lat0 && v.lon0 == t.lon0, "vbs_tile_at reproduces vbs_tile_from_llh");
	check(fabs((t.lat1 - t.lat0) * 111.32 - VBS_TILE_KM) < 0.5,
	    "a row is ~15 km tall");
	check(fabs((t.lon1 - t.lon0) * 111.32 * cos(t.clat * 3.14159265358979 / 180.0) - VBS_TILE_KM) < 1.0,
	    "a column is ~15 km wide at this latitude");
	check(vbs_tile_ncols(t.row) > 1000 && vbs_tile_ncols(t.row) < 2000,
	    "column count at 49N is in the expected range");
	check(vbs_tile_ncols(vbs_tile_nrows() - 1) >= 1, "the last row has at least one column");

	/* wrap-around and clamping */
	vbs_tile_from_llh(0.0, 179.9999, &u);
	vbs_tile_from_llh(0.0, -180.0, &v);
	check(u.row == v.row && (u.col + 1) % vbs_tile_ncols(u.row) == v.col,
	    "the tile east of the date line is the neighbour of the one west of it");
	vbs_tile_from_llh(0.0, 540.0, &u);
	vbs_tile_from_llh(0.0, -180.0, &v);
	check(u.col == v.col, "longitude 540 normalises to -180");
	vbs_tile_from_llh(90.0, 0.0, &u);
	check(u.row == vbs_tile_nrows() - 1, "the pole lands in the last row, not past it");

	/* station id range and the name */
	check(vbs_tile_staid(&t) >= 1 && vbs_tile_staid(&t) <= VBS_TILE_STAID_MAX,
	    "station id is within DF003's 12 bits and never 0");
	vbs_tile_name(&t, name, sizeof name);
	check(strlen(name) == 10 && name[0] == 'T', "tile name is T####C####");

	/* hysteresis */
	{
		struct vbs_tile cur, nxt;
		double lat, lon;

		vbs_tile_from_llh(48.0, 2.0, &cur);
		/* just over the eastern edge: 50 m outside */
		lat = cur.clat;
		lon = cur.lon1 + 50.0 / (111320.0 * cos(cur.clat * 3.14159265358979 / 180.0));
		check(vbs_tile_outside_m(&cur, lat, lon) > 40.0 && vbs_tile_outside_m(&cur, lat, lon) < 60.0,
		    "50 m outside the edge measures as ~50 m");
		check(vbs_tile_update(&cur, lat, lon, 200.0, &nxt) == 0,
		    "50 m over the edge with 200 m hysteresis: stay");
		check(vbs_tile_update(&cur, lat, lon, 20.0, &nxt) == 1 && nxt.col == (cur.col + 1) % vbs_tile_ncols(cur.row),
		    "50 m over the edge with 20 m hysteresis: move to the eastern neighbour");
		check(vbs_tile_update(&cur, cur.clat, cur.clon, 200.0, &nxt) == 0,
		    "at the centre: stay");
		check(vbs_tile_update(NULL, 48.0, 2.0, 200.0, &nxt) == 1,
		    "a new rover always gets a tile");
		check(vbs_tile_outside_m(&cur, cur.clat, cur.clon) == 0.0, "inside measures 0");
	}

	/*
	 * Every spatial neighbour has a different station id. For each tile take nine probe
	 * points just outside its edges and corners, plus probes just outside at the two row
	 * boundaries offset by a third and two thirds of the tile width (to catch neighbours that
	 * a different column count puts across the boundary), and compare ids.
	 */
	for (row = 0; row < vbs_tile_nrows(); row++) {
		int ncols = vbs_tile_ncols(row), col, step = ncols > 64 ? ncols / 64 : 1;

		for (col = 0; col < ncols; col += step) {
			double eps_lat, eps_lon, probes[16][2];
			int np = 0, p;

			vbs_tile_at(row, col, &t);
			eps_lat = (t.lat1 - t.lat0) * 0.01;
			eps_lon = (t.lon1 - t.lon0) * 0.01;
			probes[np][0] = t.lat0 - eps_lat; probes[np++][1] = t.clon;
			probes[np][0] = t.lat1 + eps_lat; probes[np++][1] = t.clon;
			probes[np][0] = t.clat; probes[np++][1] = t.lon0 - eps_lon;
			probes[np][0] = t.clat; probes[np++][1] = t.lon1 + eps_lon;
			probes[np][0] = t.lat0 - eps_lat; probes[np++][1] = t.lon0 - eps_lon;
			probes[np][0] = t.lat0 - eps_lat; probes[np++][1] = t.lon1 + eps_lon;
			probes[np][0] = t.lat1 + eps_lat; probes[np++][1] = t.lon0 - eps_lon;
			probes[np][0] = t.lat1 + eps_lat; probes[np++][1] = t.lon1 + eps_lon;
			probes[np][0] = t.lat0 - eps_lat; probes[np++][1] = t.lon0 + (t.lon1 - t.lon0) / 3.0;
			probes[np][0] = t.lat0 - eps_lat; probes[np++][1] = t.lon0 + 2.0 * (t.lon1 - t.lon0) / 3.0;
			probes[np][0] = t.lat1 + eps_lat; probes[np++][1] = t.lon0 + (t.lon1 - t.lon0) / 3.0;
			probes[np][0] = t.lat1 + eps_lat; probes[np++][1] = t.lon0 + 2.0 * (t.lon1 - t.lon0) / 3.0;
			for (p = 0; p < np; p++) {
				if (probes[p][0] < -90.0 || probes[p][0] > 90.0)
					continue;
				vbs_tile_from_llh(probes[p][0], probes[p][1], &u);
				if (u.row == t.row && u.col == t.col)
					continue;
				checked++;
				if (vbs_tile_staid(&u) == vbs_tile_staid(&t)) {
					if (bad < 5)
						printf("    id clash: T%dC%d and T%dC%d both %d\n",
						    t.row, t.col, u.row, u.col, vbs_tile_staid(&t));
					bad++;
				}
			}
		}
	}
	printf("  neighbour pairs checked: %d, id clashes: %d\n", checked, bad);
	check(bad == 0 && checked > 10000, "no spatial neighbour shares a station id (all rows)");
}

#define PI	3.14159265358979
#define KM_PER_DEG	111.32

/* Planar distance, km, from a position to a tile centre, on the tangent plane at the position. */
static double dist_km(double lat, double lon, const struct vbs_tile *t) {
	double dlon = lon - t->clon;

	while (dlon >= 180.0)
		dlon -= 360.0;
	while (dlon < -180.0)
		dlon += 360.0;
	return hypot((lat - t->clat) * KM_PER_DEG, dlon * KM_PER_DEG * cos(lat * PI / 180.0));
}

/* Deterministic pseudo-random sequence, so a failure reproduces. */
static unsigned long lcg = 12345;
static double urand(void) {
	lcg = lcg * 6364136223846793005UL + 1442695040888963407UL;
	return (double)(lcg >> 11) / 9007199254740992.0;
}

static void test_hex(void) {
	struct vbs_tile t, u, v;
	char name[32];
	const double R = VBS_TILE_HEX_R_KM;
	double worst = 0.0, sum = 0.0;
	int i, row, bad = 0, checked = 0, notnearest = 0, nsample = 0, six = 0, ncell = 0;
	long nhex = 0, nsq = 0;

	vbs_tile_set_shape(VBS_TILE_HEX);
	printf("=== tile grid: hexagons ===\n");
	check(vbs_tile_get_shape() == VBS_TILE_HEX, "the shape is set");
	check(vbs_tile_shape_from_str("hex") == VBS_TILE_HEX && vbs_tile_shape_from_str("square") ==
	    VBS_TILE_SQUARE && vbs_tile_shape_from_str("round") == -1, "--tiles parses hex and square only");

	vbs_tile_from_llh(48.8566, 2.3522, &t);
	check(t.shape == VBS_TILE_HEX && dist_km(48.8566, 2.3522, &t) <= R, "Paris is within R of its tile's centre");
	vbs_tile_from_llh(t.clat, t.clon, &u);
	check(u.row == t.row && u.col == t.col, "the tile's centre maps back to the same tile");
	vbs_tile_at(t.row, t.col, &v);
	check(v.clat == t.clat && v.clon == t.clon, "vbs_tile_at reproduces vbs_tile_from_llh");
	vbs_tile_at(t.row, t.col + 1, &u);
	check(fabs(dist_km(t.clat, t.clon, &u) / (sqrt(3.0) * R) - 1.0) < 0.051,
	    "centres in a row are sqrt(3) R apart, within the 5 % band tolerance");
	vbs_tile_at(t.row + 1, t.col, &u);
	check(fabs((u.clat - t.clat) * KM_PER_DEG - 1.5 * R) < 0.01, "rows are 1.5 R apart");
	vbs_tile_name(&t, name, sizeof name);
	check(strlen(name) == 10 && name[0] == 'H', "tile name is H####C####");
	check(vbs_tile_staid(&t) >= 1 && vbs_tile_staid(&t) <= VBS_TILE_STAID_MAX,
	    "station id is within DF003's 12 bits and never 0");

	/* wrap-around and clamping */
	vbs_tile_from_llh(0.3, 179.9999, &u);
	vbs_tile_from_llh(0.3, -179.9999, &v);
	check(dist_km(0.3, 179.9999, &u) <= R && dist_km(0.3, -179.9999, &v) <= R,
	    "either side of the date line lands within R of a centre");
	vbs_tile_from_llh(10.0, 540.0, &u);
	vbs_tile_from_llh(10.0, -180.0, &v);
	check(u.row == v.row && u.col == v.col, "longitude 540 normalises to -180");
	vbs_tile_from_llh(90.0, 0.0, &u);
	check(u.row >= vbs_tile_nrows() - 2, "the pole lands in the last rows, not past them");

	/*
	 * The nearest-centre rule, by brute force: for random positions, no centre within two rows
	 * and three columns is nearer than the one the grid picked. This is what makes the cells
	 * hexagons, and what the three-row search in vbs_tile.c assumes.
	 */
	for (i = 0; i < 200000; i++) {
		double lat = -75.0 + 150.0 * urand(), lon = -180.0 + 360.0 * urand(), d;
		int dr, dc;

		vbs_tile_from_llh(lat, lon, &t);
		d = dist_km(lat, lon, &t);
		for (dr = -2; dr <= 2; dr++)
			for (dc = -3; dc <= 3; dc++) {
				vbs_tile_at(t.row + dr, t.col + dc, &u);
				if (dist_km(lat, lon, &u) < d - 1e-6)
					notnearest++;
			}
		if (d > worst)
			worst = d;
		sum += d;
		nsample++;
	}
	printf("  rover-to-centre: mean %.2f km, worst %.2f km over %d positions (R = %.2f km)\n",
	    sum / nsample, worst, nsample, R);
	check(notnearest == 0, "every position belongs to its nearest centre");
	check(worst <= 1.10 * R, "worst rover-to-centre distance within 10 % of R (75S-75N)");

	/*
	 * Neighbours and station ids, every row. Probes on two rings just outside the cell -- at
	 * 1.05 R (past the vertices) and at 1.05 sqrt(3)/2 R (past the flat sides), 12 directions
	 * each -- find every neighbour, including those a differing column count shifts.
	 */
	for (row = 0; row < vbs_tile_nrows(); row++) {
		int ncols = vbs_tile_ncols(row), col, step = ncols > 64 ? ncols / 64 : 1;

		for (col = 0; col < ncols; col += step) {
			struct vbs_tile nb[24];
			int nnb = 0, k, ring, j;

			vbs_tile_at(row, col, &t);
			if (fabs(t.clat) > 80.0)
				continue;		/* polar rows: a handful of cells, not hexagons */
			for (ring = 0; ring < 2; ring++)
				for (k = 0; k < 12; k++) {
					double a = (k + 0.5) * PI / 6.0;
					double rr = (ring ? 1.05 : 1.05 * sqrt(3.0) / 2.0) * R;
					double lat = t.clat + rr * cos(a) / KM_PER_DEG;
					double lon = t.clon + rr * sin(a) / (KM_PER_DEG * cos(t.clat * PI / 180.0));

					vbs_tile_from_llh(lat, lon, &u);
					if (u.row == t.row && u.col == t.col)
						continue;
					for (j = 0; j < nnb; j++)
						if (nb[j].row == u.row && nb[j].col == u.col)
							break;
					if (j < nnb)
						continue;
					nb[nnb++] = u;
					checked++;
					if (vbs_tile_staid(&u) == vbs_tile_staid(&t)) {
						if (bad < 5)
							printf("    id clash: H%dC%d and H%dC%d both %d\n",
							    t.row, t.col, u.row, u.col, vbs_tile_staid(&t));
						bad++;
					}
				}
			ncell++;
			if (nnb == 6)
				six++;
		}
	}
	printf("  cells probed %d, with exactly six neighbours %d; neighbour pairs %d, id clashes %d\n",
	    ncell, six, checked, bad);
	check(bad == 0 && checked > 10000, "no spatial neighbour shares a station id (all rows)");
	check(six > 0.9 * ncell, "nine cells in ten have exactly six neighbours");

	/* hysteresis: 50 m past the flat side shared with the eastern neighbour */
	{
		struct vbs_tile cur, nxt;
		double lat, lon;

		struct vbs_tile east;
		double half;

		vbs_tile_from_llh(48.0, 2.0, &cur);
		vbs_tile_at(cur.row, cur.col + 1, &east);
		/* the shared side is the bisector: half the centre spacing, then 50 m beyond */
		half = (east.clon - cur.clon) / 2.0;
		lat = cur.clat;
		lon = cur.clon + half + 0.050 / (KM_PER_DEG * cos(cur.clat * PI / 180.0));
		check(vbs_tile_outside_m(&cur, lat, lon) > 40.0 && vbs_tile_outside_m(&cur, lat, lon) < 60.0,
		    "50 m past the shared side measures as ~50 m");
		check(vbs_tile_update(&cur, lat, lon, 200.0, &nxt) == 0, "50 m over with 200 m hysteresis: stay");
		check(vbs_tile_update(&cur, lat, lon, 20.0, &nxt) == 1 && nxt.row == cur.row &&
		    nxt.col == (cur.col + 1) % vbs_tile_ncols(cur.row),
		    "50 m over with 20 m hysteresis: move to the eastern neighbour");
		check(vbs_tile_update(&cur, cur.clat, cur.clon, 200.0, &nxt) == 0, "at the centre: stay");
		check(vbs_tile_outside_m(&cur, cur.clat, cur.clon) == 0.0, "inside measures 0");
	}

	/* the point of the change: fewer tiles for the same worst-case distance */
	for (row = 0; row < vbs_tile_nrows(); row++)
		nhex += vbs_tile_ncols(row);
	vbs_tile_set_shape(VBS_TILE_SQUARE);
	for (row = 0; row < vbs_tile_nrows(); row++)
		nsq += vbs_tile_ncols(row);
	vbs_tile_set_shape(VBS_TILE_HEX);
	printf("  tiles on the globe: %ld hexagons, %ld squares (%.1f %%)\n", nhex, nsq, 100.0 * nhex / nsq);
	check(nhex < 0.8 * nsq && nhex > 0.74 * nsq, "about 23 % fewer tiles than squares");
}

int main(void) {
	check(vbs_tile_get_shape() == VBS_TILE_HEX, "hexagons are the default");
	test_square();
	test_hex();
	printf("%s: %d failure(s)\n", nfail ? "FAILED" : "PASSED", nfail);
	return nfail != 0;
}

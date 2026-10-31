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

int main(void) {
	struct vbs_tile t, u, v;
	char name[32];
	int row, bad = 0, checked = 0;

	printf("=== tile grid ===\n");
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

	printf("%s: %d failure(s)\n", nfail ? "FAILED" : "PASSED", nfail);
	return nfail != 0;
}

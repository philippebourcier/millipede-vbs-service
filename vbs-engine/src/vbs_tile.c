#include <math.h>
#include <stdio.h>

#include "rtklib.h"
#include "vbs/vbs_tile.h"

/*
 * Metres per degree of latitude, on the sphere that sizes the grid. The grid does not need the
 * ellipsoid: it only has to be the same function everywhere, and a cell that is 14.9 km rather
 * than 15.0 changes nothing downstream.
 */
#define M_PER_DEG	111320.0

static double row_height_deg(void) {
	return VBS_TILE_KM * 1000.0 / M_PER_DEG;
}

int vbs_tile_nrows(void) {
	return (int)ceil(180.0 / row_height_deg());
}

static double row_lat0(int row) {
	return -90.0 + row * row_height_deg();
}

int vbs_tile_ncols(int row) {
	double clat = row_lat0(row) + row_height_deg() / 2.0;
	double circ = 360.0 * cos(clat * D2R) * M_PER_DEG;	/* metres round this parallel */
	int n = (int)floor(circ / (VBS_TILE_KM * 1000.0) + 0.5);

	return n < 1 ? 1 : n;
}

static double norm_lon(double lon) {
	while (lon >= 180.0)
		lon -= 360.0;
	while (lon < -180.0)
		lon += 360.0;
	return lon;
}

void vbs_tile_at(int row, int col, struct vbs_tile *out) {
	int nrows = vbs_tile_nrows(), ncols;
	double dlat = row_height_deg(), dlon;

	if (row < 0)
		row = 0;
	if (row >= nrows)
		row = nrows - 1;
	ncols = vbs_tile_ncols(row);
	col = ((col % ncols) + ncols) % ncols;
	dlon = 360.0 / ncols;
	out->row = row;
	out->col = col;
	out->lat0 = row_lat0(row);
	out->lat1 = out->lat0 + dlat;
	if (out->lat1 > 90.0)
		out->lat1 = 90.0;
	out->lon0 = -180.0 + col * dlon;
	out->lon1 = out->lon0 + dlon;
	out->clat = (out->lat0 + out->lat1) / 2.0;
	out->clon = (out->lon0 + out->lon1) / 2.0;
}

void vbs_tile_from_llh(double lat_deg, double lon_deg, struct vbs_tile *out) {
	int row, col, ncols;

	if (lat_deg > 90.0)
		lat_deg = 90.0;
	if (lat_deg < -90.0)
		lat_deg = -90.0;
	lon_deg = norm_lon(lon_deg);
	row = (int)floor((lat_deg + 90.0) / row_height_deg());
	if (row >= vbs_tile_nrows())
		row = vbs_tile_nrows() - 1;
	ncols = vbs_tile_ncols(row);
	col = (int)floor((lon_deg + 180.0) / (360.0 / ncols));
	if (col >= ncols)
		col = ncols - 1;
	vbs_tile_at(row, col, out);
}

/*
 * (row, col) -> 1..4095, neighbour-distinct.
 *
 * The first attempt was a linear hash, row * 2671 + col mod 4095, and the test found ten
 * clashes: at the date line, column 0 of a row sits next to the LAST column of the row below,
 * and 2671 is exactly the column count there, so the two land on the same value. No linear
 * multiplier can avoid that for every row, because the set of column counts is dense enough
 * to hit any residue.
 *
 * So the ID is a 2-D colouring instead: 63 row classes times 65 column classes = 4095. Two
 * tiles share an ID only if their rows agree modulo 63 AND their columns agree modulo 65.
 * Adjacent rows never agree modulo 63, which settles every cross-row neighbour including the
 * ones a differing column count puts across the boundary. Within a row, neighbours differ by
 * one column -- distinct modulo 65 -- except across the wrap, where column 0 meets column
 * ncols - 1; when that last column happens to fall in class 0 it is moved to class 32, which
 * differs from both its west neighbour (class 64) and column 0 (class 0).
 *
 * test_tile probes every row and asserts the property rather than trusting this argument.
 */
#define STAID_ROWS	63
#define STAID_COLS	65

int vbs_tile_staid(const struct vbs_tile *t) {
	int rowpart = ((t->row % STAID_ROWS) + STAID_ROWS) % STAID_ROWS;
	int colpart = ((t->col % STAID_COLS) + STAID_COLS) % STAID_COLS;
	int ncols = vbs_tile_ncols(t->row);

	if (ncols > 1 && t->col == ncols - 1 && colpart == 0)
		colpart = STAID_COLS / 2;
	return 1 + rowpart * STAID_COLS + colpart;
}

const char *vbs_tile_name(const struct vbs_tile *t, char *buf, int len) {
	snprintf(buf, len, "T%04dC%04d", t->row, t->col);
	return buf;
}

double vbs_tile_outside_m(const struct vbs_tile *t, double lat_deg, double lon_deg) {
	double dlat = 0.0, dlon = 0.0, mlon = M_PER_DEG * cos(t->clat * D2R);

	lon_deg = norm_lon(lon_deg);
	if (lat_deg < t->lat0)
		dlat = t->lat0 - lat_deg;
	else if (lat_deg >= t->lat1)
		dlat = lat_deg - t->lat1;
	/*
	 * Longitude distance across the date line: a tile at lon0 = -180 has a neighbour at
	 * +179.9. Take the shorter way round.
	 */
	if (lon_deg < t->lon0)
		dlon = t->lon0 - lon_deg;
	else if (lon_deg >= t->lon1)
		dlon = lon_deg - t->lon1;
	if (dlon > 180.0)
		dlon = 360.0 - dlon;
	return hypot(dlat * M_PER_DEG, dlon * mlon);
}

int vbs_tile_update(const struct vbs_tile *cur, double lat_deg, double lon_deg,
    double hysteresis_m, struct vbs_tile *out) {
	struct vbs_tile t;

	vbs_tile_from_llh(lat_deg, lon_deg, &t);
	if (cur == NULL) {
		*out = t;
		return 1;
	}
	if (t.row == cur->row && t.col == cur->col)
		return 0;
	if (vbs_tile_outside_m(cur, lat_deg, lon_deg) <= hysteresis_m)
		return 0;
	*out = t;
	return 1;
}

void vbs_tile_xyz(const struct vbs_tile *t, double h_ellipsoidal, double *xyz) {
	double pos[3] = {t->clat * D2R, t->clon * D2R, h_ellipsoidal};

	pos2ecef(pos, xyz);
}

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "rtklib.h"
#include "vbs/vbs_tile.h"

/*
 * Metres per degree of latitude, on the sphere that sizes the grid. The grid does not need the
 * ellipsoid: it only has to be the same function everywhere, and a cell that is 14.9 km rather
 * than 15.0 changes nothing downstream.
 */
#define M_PER_DEG	111320.0
#define SQRT3		1.7320508075688772

/* One process, one grid: set at startup (vbsd --tiles), read by every function below. */
static int grid_shape = VBS_TILE_HEX;

int vbs_tile_set_shape(int shape) {
	if (shape != VBS_TILE_HEX && shape != VBS_TILE_SQUARE)
		return -1;
	grid_shape = shape;
	return 0;
}

int vbs_tile_get_shape(void) {
	return grid_shape;
}

int vbs_tile_shape_from_str(const char *s) {
	if (s == NULL)
		return -1;
	if (!strcmp(s, "hex") || !strcmp(s, "hexagon"))
		return VBS_TILE_HEX;
	if (!strcmp(s, "square"))
		return VBS_TILE_SQUARE;
	return -1;
}

/*
 * Row pitch and the target column spacing, metres. Square: both VBS_TILE_KM. Hexagon of
 * circumradius R ("pointy-top": vertices north and south): rows 1.5 R apart, centres sqrt(3) R
 * apart along a row, odd rows shifted by half that.
 */
static double row_pitch_m(int shape) {
	return shape == VBS_TILE_HEX ? 1.5 * VBS_TILE_HEX_R_KM * 1000.0 : VBS_TILE_KM * 1000.0;
}

static double col_pitch_m(int shape) {
	return shape == VBS_TILE_HEX ? SQRT3 * VBS_TILE_HEX_R_KM * 1000.0 : VBS_TILE_KM * 1000.0;
}

static double row_height_deg(int shape) {
	return row_pitch_m(shape) / M_PER_DEG;
}

static int nrows_of(int shape) {
	return (int)ceil(180.0 / row_height_deg(shape));
}

int vbs_tile_nrows(void) {
	return nrows_of(grid_shape);
}

static double row_lat0(int shape, int row) {
	return -90.0 + row * row_height_deg(shape);
}

/* The column count that keeps a cell at its nominal width on this row's central parallel. */
static int ideal_ncols(int shape, int row) {
	double clat = row_lat0(shape, row) + row_height_deg(shape) / 2.0;
	double circ = 360.0 * cos(clat * D2R) * M_PER_DEG;	/* metres round this parallel */
	int n = (int)floor(circ / col_pitch_m(shape) + 0.5);

	return n < 1 ? 1 : n;
}

/*
 * Hexagon rows share their column count in bands (D172). The half-column offset between
 * adjacent rows is what makes the cells hexagons; with a column count of its own per row --
 * about 4 columns more or less than the next row at 48 deg -- that offset drifts round the
 * globe, and wherever it drifts to zero the cells degrade to rectangles: measured, only 62 % of
 * cells had six neighbours and the worst rover-to-centre distance was 12.1 km. Holding the count
 * while the cell width stays within HEX_BAND_EPS of nominal keeps rows aligned inside a band;
 * only the seam between two bands is irregular. At 5 %: 98 % of cells have six neighbours, p99
 * 10.0 km, worst 11.4 km on seam rows only, and 21 % fewer tiles than squares.
 *
 * Built outwards from the equator, so the two hemispheres mirror; computed once.
 */
#define HEX_BAND_EPS	0.05
#define HEX_MAXROWS	2048

static int hex_tab[HEX_MAXROWS], hex_tab_ready;

static void hex_tab_build(void) {
	int nrows = nrows_of(VBS_TILE_HEX), eq = nrows / 2, r, n, cur;

	cur = 0;
	for (r = eq; r < nrows; r++) {
		n = ideal_ncols(VBS_TILE_HEX, r);
		if (cur == 0 || fabs((double)cur / n - 1.0) > HEX_BAND_EPS)
			cur = n;
		hex_tab[r] = cur;
	}
	cur = 0;
	for (r = eq - 1; r >= 0; r--) {
		n = ideal_ncols(VBS_TILE_HEX, r);
		if (cur == 0 || fabs((double)cur / n - 1.0) > HEX_BAND_EPS)
			cur = n;
		hex_tab[r] = cur;
	}
	hex_tab_ready = 1;
}

static int ncols_of(int shape, int row) {
	if (shape != VBS_TILE_HEX)
		return ideal_ncols(shape, row);
	if (!hex_tab_ready)
		hex_tab_build();
	return hex_tab[row];
}

int vbs_tile_ncols(int row) {
	return ncols_of(grid_shape, row);
}

static double norm_lon(double lon) {
	while (lon >= 180.0)
		lon -= 360.0;
	while (lon < -180.0)
		lon += 360.0;
	return lon;
}

/* Longitude of the first column's west edge: odd hexagon rows are shifted half a column east. */
static double row_lon_offset(int shape, int row, int ncols) {
	return shape == VBS_TILE_HEX && (row & 1) && ncols > 1 ? 0.5 * 360.0 / ncols : 0.0;
}

static void tile_at(int shape, int row, int col, struct vbs_tile *out) {
	int nrows = nrows_of(shape), ncols;
	double dlat = row_height_deg(shape), dlon, off;

	if (row < 0)
		row = 0;
	if (row >= nrows)
		row = nrows - 1;
	ncols = ncols_of(shape, row);
	col = ((col % ncols) + ncols) % ncols;
	dlon = 360.0 / ncols;
	off = row_lon_offset(shape, row, ncols);
	out->shape = shape;
	out->row = row;
	out->col = col;
	if (shape == VBS_TILE_HEX) {
		/*
		 * Centre at mid-row; the bounding box reaches R north and south (the vertices) and
		 * half a column east and west (the flat sides).
		 */
		double r_deg = VBS_TILE_HEX_R_KM * 1000.0 / M_PER_DEG;

		out->clat = row_lat0(shape, row) + dlat / 2.0;
		if (out->clat > 90.0)
			out->clat = (row_lat0(shape, row) + 90.0) / 2.0;
		out->clon = norm_lon(-180.0 + off + (col + 0.5) * dlon);
		out->lat0 = out->clat - r_deg;
		out->lat1 = out->clat + r_deg;
		if (out->lat0 < -90.0)
			out->lat0 = -90.0;
		if (out->lat1 > 90.0)
			out->lat1 = 90.0;
		out->lon0 = out->clon - dlon / 2.0;
		out->lon1 = out->clon + dlon / 2.0;
		return;
	}
	out->lat0 = row_lat0(shape, row);
	out->lat1 = out->lat0 + dlat;
	if (out->lat1 > 90.0)
		out->lat1 = 90.0;
	out->lon0 = -180.0 + col * dlon;
	out->lon1 = out->lon0 + dlon;
	out->clat = (out->lat0 + out->lat1) / 2.0;
	out->clon = (out->lon0 + out->lon1) / 2.0;
}

void vbs_tile_at(int row, int col, struct vbs_tile *out) {
	tile_at(grid_shape, row, col, out);
}

/*
 * Local east/north offset of (lat, lon) from a centre, metres, on the tangent plane at the
 * position's own latitude. Using the position's latitude (not the centre's) for every
 * candidate makes the comparison between candidates exact for a given position, which is all
 * the nearest-centre rule needs to be a function.
 */
static void local_en(double lat, double lon, double clat, double clon, double *e, double *n) {
	double dlon = norm_lon(lon - clon);

	*e = dlon * M_PER_DEG * cos(lat * D2R);
	*n = (lat - clat) * M_PER_DEG;
}

static void hex_from_llh(double lat, double lon, struct vbs_tile *out) {
	int nrows = nrows_of(VBS_TILE_HEX), row0, dr, dc, found = 0;
	double best = 0.0;

	row0 = (int)floor((lat + 90.0) / row_height_deg(VBS_TILE_HEX));
	if (row0 >= nrows)
		row0 = nrows - 1;
	/*
	 * A hexagon reaches R north and south of its centre, two thirds of the row pitch, so a
	 * position can only belong to its own row or an adjacent one; and within a row to the
	 * column under it or one of the two beside it.
	 */
	for (dr = -1; dr <= 1; dr++) {
		int row = row0 + dr, ncols, c0;
		double dlon, off;

		if (row < 0 || row >= nrows)
			continue;
		ncols = ncols_of(VBS_TILE_HEX, row);
		dlon = 360.0 / ncols;
		off = row_lon_offset(VBS_TILE_HEX, row, ncols);
		c0 = (int)floor(fmod(lon + 180.0 - off + 360.0, 360.0) / dlon);
		for (dc = -1; dc <= 1; dc++) {
			struct vbs_tile t;
			double e, n, d2;

			tile_at(VBS_TILE_HEX, row, c0 + dc, &t);
			local_en(lat, lon, t.clat, t.clon, &e, &n);
			d2 = e * e + n * n;
			/* ties go to the lower row, then the lower column: still a pure function */
			if (!found || d2 < best - 1e-6 ||
			    (fabs(d2 - best) <= 1e-6 && (t.row < out->row || (t.row == out->row && t.col < out->col)))) {
				*out = t;
				best = d2;
				found = 1;
			}
		}
	}
}

void vbs_tile_from_llh(double lat_deg, double lon_deg, struct vbs_tile *out) {
	int row, col, ncols;

	if (lat_deg > 90.0)
		lat_deg = 90.0;
	if (lat_deg < -90.0)
		lat_deg = -90.0;
	lon_deg = norm_lon(lon_deg);
	if (grid_shape == VBS_TILE_HEX) {
		hex_from_llh(lat_deg, lon_deg, out);
		return;
	}
	row = (int)floor((lat_deg + 90.0) / row_height_deg(VBS_TILE_SQUARE));
	if (row >= nrows_of(VBS_TILE_SQUARE))
		row = nrows_of(VBS_TILE_SQUARE) - 1;
	ncols = ncols_of(VBS_TILE_SQUARE, row);
	col = (int)floor((lon_deg + 180.0) / (360.0 / ncols));
	if (col >= ncols)
		col = ncols - 1;
	tile_at(VBS_TILE_SQUARE, row, col, out);
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
 * A hexagon's six neighbours are two in its row (one column away) and two in each adjacent
 * row, so the same argument covers it unchanged.
 *
 * test_tile probes every row of both grids and asserts the property rather than trusting this
 * argument.
 */
#define STAID_ROWS	63
#define STAID_COLS	65

int vbs_tile_staid(const struct vbs_tile *t) {
	int rowpart = ((t->row % STAID_ROWS) + STAID_ROWS) % STAID_ROWS;
	int colpart = ((t->col % STAID_COLS) + STAID_COLS) % STAID_COLS;
	int ncols = ncols_of(t->shape, t->row);

	if (ncols > 1 && t->col == ncols - 1 && colpart == 0)
		colpart = STAID_COLS / 2;
	return 1 + rowpart * STAID_COLS + colpart;
}

const char *vbs_tile_name(const struct vbs_tile *t, char *buf, int len) {
	snprintf(buf, len, "%c%04dC%04d", t->shape == VBS_TILE_HEX ? 'H' : 'T', t->row, t->col);
	return buf;
}

double vbs_tile_outside_m(const struct vbs_tile *t, double lat_deg, double lon_deg) {
	double dlat = 0.0, dlon = 0.0, mlon = M_PER_DEG * cos(t->clat * D2R);

	lon_deg = norm_lon(lon_deg);
	if (t->shape == VBS_TILE_HEX) {
		/*
		 * The cell is the set of positions nearer this centre than any other, so "how far
		 * outside" is the distance past the bisector between this centre C and the centre N
		 * the position now belongs to: (|pC|^2 - |pN|^2) / (2 |CN|), on the tangent plane at
		 * the position.
		 */
		struct vbs_tile n;
		double ec, nc, en, nn, dx, dy, cn;

		if (lat_deg > 90.0)
			lat_deg = 90.0;
		if (lat_deg < -90.0)
			lat_deg = -90.0;
		hex_from_llh(lat_deg, lon_deg, &n);
		if (n.row == t->row && n.col == t->col)
			return 0.0;
		local_en(lat_deg, lon_deg, t->clat, t->clon, &ec, &nc);
		local_en(lat_deg, lon_deg, n.clat, n.clon, &en, &nn);
		dx = ec - en;
		dy = nc - nn;
		cn = hypot(dx, dy);
		if (cn <= 0.0)
			return 0.0;
		return (ec * ec + nc * nc - en * en - nn * nn) / (2.0 * cn);
	}
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
	if (t.shape == cur->shape && t.row == cur->row && t.col == cur->col)
		return 0;
	if (t.shape == cur->shape && vbs_tile_outside_m(cur, lat_deg, lon_deg) <= hysteresis_m)
		return 0;
	*out = t;
	return 1;
}

void vbs_tile_xyz(const struct vbs_tile *t, double h_ellipsoidal, double *xyz) {
	double pos[3] = {t->clat * D2R, t->clon * D2R, h_ellipsoidal};

	pos2ecef(pos, xyz);
}

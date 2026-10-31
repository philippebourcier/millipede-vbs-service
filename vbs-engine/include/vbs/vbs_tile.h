#ifndef __VBS_TILE_H__
#define __VBS_TILE_H__

/*
 * The fixed global tile grid (spec §6b, decided 2026-09-12; hexagons by default since
 * 2026-10-31, D172).
 *
 * A tile is a pure function of position: the world is divided once into cells, so two rovers
 * at the same place get the same tile today and tomorrow, with no state and no negotiation.
 * Everything the VBS caster keys on -- the synthesis point, the RTCM station ID, which rovers
 * share a synthesis -- follows from the tile.
 *
 * Two shapes, chosen once per process with vbs_tile_set_shape():
 *
 *   VBS_TILE_HEX (default)  Rows equal in latitude, 1.5 R apart; columns sqrt(3) R apart, sized
 *                           per row; odd rows offset by half a column. A position belongs to the
 *                           nearest centre among its own row and the two adjacent ones, so the
 *                           cells are the Voronoi cells of that lattice: hexagons of circumradius
 *                           R = VBS_TILE_HEX_R_KM. R is the square's half-diagonal, so the worst
 *                           rover-to-VBS distance is unchanged (10.6 km) while a hexagon covers
 *                           2.60 R^2 = 292 km2 against the square's 225: about 23 % fewer tiles,
 *                           hence fewer syntheses, for the same area. Only three cells meet at a
 *                           vertex (four for squares), and the six neighbours are equidistant.
 *   VBS_TILE_SQUARE         Rows equal in latitude, VBS_TILE_KM tall; columns ~VBS_TILE_KM wide,
 *                           sized per row. The original grid.
 *
 * Either way the column count changes from row to row, so a cell's neighbour across a row
 * boundary is found by position, not by index.
 *
 * The station ID is derived from (row, col) so that every spatial neighbour of a tile gets a
 * different ID. That is the only property the ID needs: a rover leaving a tile must see its base
 * change, which RTCM signals by the station ID changing. Two tiles a thousand kilometres apart
 * sharing an ID is harmless, because no rover is ever on both. The same colouring serves both
 * shapes: every neighbour is in the same row (one column away) or in an adjacent row.
 */

#define VBS_TILE_KM		15.0
#define VBS_TILE_HEX_R_KM	(VBS_TILE_KM / 1.4142135623730951)	/* circumradius = square half-diagonal */
#define VBS_TILE_STAID_MAX	4095	/* DF003 is 12 bits; 0 is reserved */

#define VBS_TILE_HEX		0
#define VBS_TILE_SQUARE		1

struct vbs_tile {
	int shape;			/* VBS_TILE_HEX or VBS_TILE_SQUARE */
	int row, col;
	/*
	 * Square: the cell, [lat0, lat1) x [lon0, lon1). Hexagon: the cell's bounding box (it
	 * contains the cell; neighbours' boxes overlap it).
	 */
	double lat0, lat1;		/* degrees */
	double lon0, lon1;		/* degrees */
	double clat, clon;		/* centre, degrees */
};

/*
 * The grid shape every function below uses for new tiles. Set it once, before any tile is made
 * (vbsd does from --tiles); the default is VBS_TILE_HEX. Returns 0, or -1 for an unknown shape.
 */
int vbs_tile_set_shape(int shape);
int vbs_tile_get_shape(void);
/* "hex" / "square" -> VBS_TILE_*, -1 if neither. */
int vbs_tile_shape_from_str(const char *s);

/* The tile containing a position. Longitude is normalised to [-180, 180). */
void vbs_tile_from_llh(double lat_deg, double lon_deg, struct vbs_tile *out);

/* The tile at explicit grid coordinates; col is wrapped to the row's width. */
void vbs_tile_at(int row, int col, struct vbs_tile *out);

/* Columns in a row of the current shape. Rows near the poles have one. */
int vbs_tile_ncols(int row);
int vbs_tile_nrows(void);

/* RTCM reference-station ID for a tile, in 1..VBS_TILE_STAID_MAX. */
int vbs_tile_staid(const struct vbs_tile *t);

/* A stable, human-readable name: "H0987C1234" for a hexagon, "T0987C1234" for a square. */
const char *vbs_tile_name(const struct vbs_tile *t, char *buf, int len);

/*
 * How far a position lies OUTSIDE a tile, in metres -- 0 if inside. This is the quantity the
 * boundary hysteresis is applied to: a rover keeps its tile until it is more than the
 * hysteresis distance into another one, so a rover parked on an edge does not alternate.
 */
double vbs_tile_outside_m(const struct vbs_tile *t, double lat_deg, double lon_deg);

/*
 * Tile-change decision with hysteresis. `cur` is the rover's current tile (or NULL for a new
 * rover). Returns 1 and fills `out` if the rover should move to a different tile, 0 otherwise.
 */
int vbs_tile_update(const struct vbs_tile *cur, double lat_deg, double lon_deg,
	double hysteresis_m, struct vbs_tile *out);

/* ECEF of the tile centre at an ellipsoidal height, metres. */
void vbs_tile_xyz(const struct vbs_tile *t, double h_ellipsoidal, double *xyz);

#endif /* __VBS_TILE_H__ */

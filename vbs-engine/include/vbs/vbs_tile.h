#ifndef __VBS_TILE_H__
#define __VBS_TILE_H__

/*
 * The fixed global tile grid (spec §6b, decided 2026-09-12).
 *
 * A tile is a pure function of position: the world is divided once into cells of about
 * VBS_TILE_KM on a side, so two rovers at the same place get the same tile today and tomorrow,
 * with no state and no negotiation. Everything the VBS caster keys on -- the synthesis point,
 * the RTCM station ID, which rovers share a synthesis -- follows from the tile.
 *
 * Rows are equal in latitude. Columns are sized per row so a cell stays ~VBS_TILE_KM wide at
 * that latitude rather than shrinking toward the poles; the column count therefore changes from
 * row to row, and a cell's neighbour across a row boundary is found by position, not by index.
 *
 * The station ID is derived from (row, col) so that every spatial neighbour of a tile -- the
 * eight around it, and the two or three that a differing column count can put across a row
 * boundary -- gets a different ID. That is the only property the ID needs: a rover leaving a
 * tile must see its base change, which RTCM signals by the station ID changing. Two tiles a
 * thousand kilometres apart sharing an ID is harmless, because no rover is ever on both.
 */

#define VBS_TILE_KM		15.0
#define VBS_TILE_STAID_MAX	4095	/* DF003 is 12 bits; 0 is reserved */

struct vbs_tile {
	int row, col;
	double lat0, lat1;		/* degrees, [lat0, lat1) */
	double lon0, lon1;		/* degrees, [lon0, lon1) */
	double clat, clon;		/* centre, degrees */
};

/* The tile containing a position. Longitude is normalised to [-180, 180). */
void vbs_tile_from_llh(double lat_deg, double lon_deg, struct vbs_tile *out);

/* The tile at explicit grid coordinates; col is wrapped to the row's width. */
void vbs_tile_at(int row, int col, struct vbs_tile *out);

/* Columns in a row. Rows near the poles have one. */
int vbs_tile_ncols(int row);
int vbs_tile_nrows(void);

/* RTCM reference-station ID for a tile, in 1..VBS_TILE_STAID_MAX. */
int vbs_tile_staid(const struct vbs_tile *t);

/* A stable, human-readable name such as "T0987C1234", for logs and mountpoints. */
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

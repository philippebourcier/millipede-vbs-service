#include <math.h>
#include <float.h>
#include <string.h>

#include "vbs/vbs_basemeta.h"

/*
 * Base metadata, the quality gate, and D23's tide handling. See the header.
 */

/*
 * tidedisp()'s argument is a *bitmask*, not prcopt_t::tidecorr's enum -- D23 records the trap
 * in full: both MRTKLIB's and RTKLIB-EX's relative-mode zdres() pass tidecorr straight through,
 * so configuring tidecorr=2 there selects bit 2 alone (ocean loading only, solid and pole
 * *off*), which is the opposite of what its documentation says. The engine never goes through
 * tidecorr.
 *
 * Bit 1 alone: solid Earth tide, expressed against a conventional tide-free coordinate, which
 * is what RGF93/ETRS89 is. Bit 8 would express it against a mean-tide coordinate and is wrong
 * here; getting it backwards biases every rover by the permanent deformation, 2.4 cm at 43°N
 * to 4.5 cm at 50°N, nearly constant in space and therefore invisible to any leave-one-out
 * residual.
 */
#define TIDE_SOLID_ONLY		1

static const char *REASON[VBS_BASE_NREASON] = {
	"ok", "no coordinate", "frame unknown", "frame mismatch", "stale", "unhealthy",
	"coordinate outlier"
};

const char *vbs_base_reason_str(int reason) {
	return (reason >= 0 && reason < VBS_BASE_NREASON) ? REASON[reason] : "?";
}

void vbs_basemeta_default_opt(struct vbs_basemeta_opt *opt) {
	memset(opt, 0, sizeof *opt);
	/*
	 * 30 s: one epoch at the network's slowest streaming rate plus margin. A base that has
	 * gone quiet keeps its last coordinate, and using it would put a stale atmosphere into
	 * the corrections with nothing in the residuals to show the data was old.
	 */
	opt->max_age_s = 30.0;
	opt->apply_solid_tide = 1;
	opt->max_coord_resid_h_m = 0.010;	/* the network's published gate for E_N, E_E */
	opt->max_coord_resid_v_m = 0.020;	/* ... and for E_H */
}

int vbs_base_init(struct vbs_base *b, const char *name, const char *frame, const double *xyz) {
	memset(b, 0, sizeof *b);
	if (name == NULL || xyz == NULL)
		return -1;
	snprintf(b->name, sizeof b->name, "%s", name);
	if (frame != NULL)
		snprintf(b->frame, sizeof b->frame, "%s", frame);
	b->xyz[0] = xyz[0];
	b->xyz[1] = xyz[1];
	b->xyz[2] = xyz[2];
	if (norm(b->xyz, 3) > RE_WGS84 / 2.0)
		ecef2pos(b->xyz, b->pos);
	b->healthy = 1;
	b->reason = VBS_BASE_OK;
	return 0;
}

void vbs_base_set_coord_resid(struct vbs_base *b, const double *enu) {
	int i;

	if (b == NULL || enu == NULL)
		return;
	for (i = 0; i < 3; i++)
		b->coord_resid[i] = enu[i];
	b->coord_resid_known = 1;
}

int vbs_base_check(struct vbs_base *b, const char *cluster_frame, gtime_t now,
	const struct vbs_basemeta_opt *opt) {
	struct vbs_basemeta_opt defopt;

	if (opt == NULL) {
		vbs_basemeta_default_opt(&defopt);
		opt = &defopt;
	}
	/*
	 * Order matters: a base with no coordinate cannot be frame-checked, and a base with the
	 * wrong frame must be refused before its age is considered, so that the reported reason
	 * is the most fundamental one rather than whichever test happened to run first.
	 */
	if (norm(b->xyz, 3) <= RE_WGS84 / 2.0)
		return b->reason = VBS_BASE_NO_COORD;
	if (b->frame[0] == '\0')
		return b->reason = VBS_BASE_NO_FRAME;
	if (cluster_frame != NULL && cluster_frame[0] != '\0' &&
	    strcmp(b->frame, cluster_frame) != 0)
		return b->reason = VBS_BASE_FRAME_MISMATCH;
	if (!b->healthy)
		return b->reason = VBS_BASE_UNHEALTHY;
	if (b->coord_resid_known &&
	    ((opt->max_coord_resid_h_m > 0.0 &&
	      hypot(b->coord_resid[0], b->coord_resid[1]) > opt->max_coord_resid_h_m) ||
	     (opt->max_coord_resid_v_m > 0.0 &&
	      fabs(b->coord_resid[2]) > opt->max_coord_resid_v_m)))
		return b->reason = VBS_BASE_COORD_OUTLIER;
	if (opt->max_age_s > 0.0 && b->last_obs.time != 0 &&
	    fabs(timediff(now, b->last_obs)) > opt->max_age_s)
		return b->reason = VBS_BASE_STALE;
	return b->reason = VBS_BASE_OK;
}

void vbs_basemeta_tide_xyz(const double *xyz, gtime_t time,
	const struct vbs_basemeta_opt *opt, double *xyz_out) {
	struct vbs_basemeta_opt defopt;
	static erp_t erp;		/* zeroed: no pole tide is requested, so none is read */
	double dr[3];
	int i;

	if (opt == NULL) {
		vbs_basemeta_default_opt(&defopt);
		opt = &defopt;
	}
	for (i = 0; i < 3; i++)
		xyz_out[i] = xyz[i];
	if (!opt->apply_solid_tide || norm((double *)xyz, 3) <= RE_WGS84 / 2.0)
		return;
	tidedisp(gpst2utc(time), xyz, TIDE_SOLID_ONLY, &erp, NULL, dr);
	for (i = 0; i < 3; i++)
		xyz_out[i] += dr[i];
}

void vbs_base_instant_xyz(const struct vbs_base *b, gtime_t time,
	const struct vbs_basemeta_opt *opt, double *xyz_out) {
	vbs_basemeta_tide_xyz(b->xyz, time, opt, xyz_out);
}

double vbs_basemeta_geoid_sep(const double *pos) {
	double n = geoidh(pos);

	/*
	 * geoidh() returns 0.0 to mean "error", and 0.0 is also a real separation where the
	 * geoid crosses the ellipsoid. The two cannot be told apart from the return value, so
	 * exact zero is treated as unavailable: refusing a handful of legitimate points on a
	 * contour of measure zero costs nothing, while accepting an error as a value puts a
	 * ~50 m height blunder into a VBS coordinate with nothing downstream to catch it.
	 */
	if (n == 0.0 || !isfinite(n) || fabs(n) > VBS_GEOID_MAX_M)
		return NAN;
	return n;
}

int vbs_basemeta_grid_xyz(const double *pos_latlon, double dem_orthometric_h, int over_water,
	double *xyz_out, double *ellipsoidal_h_out) {
	double pos[3], n, h;

	pos[0] = pos_latlon[0];
	pos[1] = pos_latlon[1];
	pos[2] = 0.0;

	if (isnan(n = vbs_basemeta_geoid_sep(pos)))
		return -1;

	/*
	 * Clamp before converting, in orthometric space, because that is where "sea surface"
	 * means zero: the geoid *is* the zero of an orthometric height.
	 */
	h = dem_orthometric_h;
	if (over_water && h < 0.0)
		h = 0.0;

	pos[2] = h + n;
	/*
	 * The range tropmodel() will serve. A cell outside it must publish no VBS at all --
	 * D36 measured what happens when the model is asked anyway and its 0.0 is differenced
	 * against a real delay.
	 */
	if (pos[2] < -100.0 || pos[2] > 10000.0)
		return -1;

	if (ellipsoidal_h_out != NULL)
		*ellipsoidal_h_out = pos[2];
	if (xyz_out != NULL)
		pos2ecef(pos, xyz_out);
	return 0;
}

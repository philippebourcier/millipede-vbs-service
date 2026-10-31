#ifndef __VBS_BASEMETA_H__
#define __VBS_BASEMETA_H__

#include "rtklib.h"

/*
 * Physical base metadata and the quality gate -- spec §4.1's inputs, and D23's tide handling.
 *
 * A base enters the network only if the engine knows enough about it to use it safely: an
 * ECEF coordinate, the frame that coordinate is in, and recent data. `01-design-study.md` §3.4
 * found a station whose catalogue coordinate disagreed with its stream by enough to fail an
 * independent static baseline, so "we have a number for it" is not the same as "we can trust
 * it", and the gate is where that distinction lives.
 *
 * **Frames do not mix.** D13: clustering is same-frame only. Two bases whose coordinates are
 * in different realisations differ by the plate motion between their epochs -- centimetres to
 * decimetres across Europe -- and that difference would enter the network corrections as a
 * fictitious geometric error with no time signature to reveal it. The rule is global, not a
 * France-border caveat.
 *
 * **Tides (D23).** The instantaneous antenna position is the catalogue coordinate plus the
 * solid Earth tide, which reaches 11.7 mm of differential over a 160 km cluster. The catalogue
 * is RGF93 -- ETRS89, and therefore IERS *conventional tide-free* -- so `tidedisp()` is called
 * **without** bit 8, which would express the displacement against a mean-tide coordinate
 * instead. Pole tide and ocean loading are deliberately absent: measured at ≤0.34 mm and
 * removed by the §4.5 interpolation to ≤0.01 mm, they are below every threshold in this
 * project (`04-geophysical-corrections.md`).
 *
 * The same call must be applied to the VBS point, or the consistency that makes the tide
 * cancel at the rover is broken. That is why it lives here rather than inside the EKF.
 */

#define VBS_BASE_MAXNAME	32
#define VBS_BASE_MAXFRAME	16

/* Why a base was refused. Reported rather than silently dropped. */
#define VBS_BASE_OK		0
#define VBS_BASE_NO_COORD	1	/* no usable ECEF position */
#define VBS_BASE_NO_FRAME	2	/* frame unknown: cannot be checked against the cluster */
#define VBS_BASE_FRAME_MISMATCH	3	/* frame differs from the cluster's (D13) */
#define VBS_BASE_STALE		4	/* no data recently enough to be trusted */
#define VBS_BASE_UNHEALTHY	5	/* flagged by the operator or by residual screening */
#define VBS_BASE_COORD_OUTLIER	6	/* its declared position disagrees with its neighbours (D56) */
#define VBS_BASE_NREASON	7

const char *vbs_base_reason_str(int reason);

struct vbs_base {
	char name[VBS_BASE_MAXNAME];		/* mountpoint */
	char frame[VBS_BASE_MAXFRAME];		/* "RGF93", "ETRF2000", ... */
	double xyz[3];				/* catalogue ECEF, metres */
	double pos[3];				/* geodetic {lat, lon, h}, rad/m */
	gtime_t last_obs;			/* last epoch seen, for the staleness test */
	int healthy;				/* cleared by residual screening (§4.3 FDE) */
	int reason;				/* VBS_BASE_* as of the last gate evaluation */
	/*
	 * How far this station's declared position sits from where its neighbours put it, local
	 * ENU in metres, and whether that has been measured at all (D56).
	 *
	 * A declared coordinate is what the base broadcasts in RTCM 1005/1006, not what the base
	 * is: on the live network one station of five was found seven to thirteen centimetres
	 * high, six times the network's own published acceptance gate. A wrong coordinate does not
	 * cancel at the rover -- it enters as a satellite-dependent geometric error along each
	 * line of sight, lands in that baseline's ambiguities, and is then interpolated into
	 * every VBS the cluster serves.
	 *
	 * The engine cannot measure this itself; it comes from the leave-one-out network
	 * integrity monitor (D59), which compares each station against the rest continuously,
	 * or from station metadata. Unmeasured is not the same as good, so `known` is separate:
	 * a station nobody has checked is admitted, and the monitor is what changes that.
	 */
	double coord_resid[3];
	int coord_resid_known;
};

struct vbs_basemeta_opt {
	double max_age_s;		/* older than this and the base is stale */
	int apply_solid_tide;		/* D23; on by default */
	/*
	 * Coordinate acceptance, horizontal and vertical, metres. Defaulted to the thresholds
	 * a typical network's build-base documentation sets for a station's published position --
	 * "below 10 mm for E_N and E_E, and below 20 mm for E_H" -- so a station this engine
	 * refuses is one the network's own rules already refuse. 0 disables either test.
	 */
	double max_coord_resid_h_m;
	double max_coord_resid_v_m;
};

void vbs_basemeta_default_opt(struct vbs_basemeta_opt *opt);

/*
 * Register a base. `frame` may be NULL or empty, in which case the base will fail the gate
 * with VBS_BASE_NO_FRAME rather than being silently assumed to match. Returns 0 on success.
 */
int vbs_base_init(struct vbs_base *b, const char *name, const char *frame, const double *xyz);

/*
 * Record how far this base's declared position was measured to be from its neighbours', in
 * local ENU metres. Until this is called the base counts as unchecked, not as good.
 */
void vbs_base_set_coord_resid(struct vbs_base *b, const double *enu);

/*
 * Evaluate the gate for one base against a reference frame (NULL: accept any single frame).
 * Sets and returns b->reason; VBS_BASE_OK means usable.
 */
int vbs_base_check(struct vbs_base *b, const char *cluster_frame, gtime_t now,
	const struct vbs_basemeta_opt *opt);

/*
 * Instantaneous antenna position: the catalogue coordinate plus the solid Earth tide (D23).
 *
 * Apply this to every physical base *and* to the VBS point, from the same model at the same
 * epoch. Applying it to only one side is worse than applying it to neither: the tide is
 * common-mode across a cluster and cancels when treated consistently, so an inconsistent
 * application converts a term that would have vanished into a real error.
 */
void vbs_base_instant_xyz(const struct vbs_base *b, gtime_t time,
	const struct vbs_basemeta_opt *opt, double *xyz_out);

/* Same, for an arbitrary point (the VBS), so the two cannot drift apart. */
void vbs_basemeta_tide_xyz(const double *xyz, gtime_t time,
	const struct vbs_basemeta_opt *opt, double *xyz_out);

/*
 * --- Grid-VBS heights (D15, D37) --------------------------------------------------------
 *
 * A Grid-VBS point sits at a cell centre with its height from a DEM lookup, and DEM values are
 * **orthometric** -- heights above the geoid. Everything downstream wants **ellipsoidal**: the
 * 1005 coordinate, the geometric range, and the troposphere model. In France the two differ by
 * about 50 m, which is ~1.5 cm of zenith troposphere delivered straight into the rover's
 * vertical -- the component D17 gives the loosest target, and therefore the one least likely
 * to be noticed going wrong.
 */

#define VBS_GEOID_MAX_M		110.0	/* |N| never exceeds this anywhere on Earth */

/*
 * Geoid separation N at a position, in metres: ellipsoidal = orthometric + N.
 *
 * Returns NAN, not 0.0, when the model cannot answer. `geoidh()` signals failure by returning
 * 0.0, which is a legitimate separation in the places where the geoid crosses the ellipsoid --
 * the third in-band sentinel found in this library (after `phase_m()` in D21 and `tropmodel()`
 * in D36), and guarded rather than trusted. An implausible magnitude is rejected the same way.
 */
double vbs_basemeta_geoid_sep(const double *pos);

/*
 * ECEF position of a Grid-VBS point from its cell centre and DEM height.
 *
 * `dem_orthometric_h` is the DEM value as read. `over_water` clamps the point to the sea
 * surface: some DEMs return **bathymetry**, so an offshore cell can yield a height hundreds of
 * metres below the geoid -- which lands inside the range where `tropmodel()` returns its 0.0
 * sentinel, and D36 measured what that produces (a -4.79 m "troposphere correction"). A VBS
 * over water sits at the sea surface, never beneath it.
 *
 * Returns 0 on success, -1 if no geoid separation was available or the result is outside the
 * range the troposphere model can serve -- in which case no VBS may be published for that
 * cell, rather than one with a silently wrong height.
 */
int vbs_basemeta_grid_xyz(const double *pos_latlon, double dem_orthometric_h, int over_water,
	double *xyz_out, double *ellipsoidal_h_out);

#endif /* __VBS_BASEMETA_H__ */

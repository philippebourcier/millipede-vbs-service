#ifndef __VBS_POLICY_H__
#define __VBS_POLICY_H__

/*
 * The serving decision (spec §6c), ported from data/scripts/vbs_policy.py.
 *
 * The Python module is the reference: its docstrings carry the evidence for every threshold,
 * and this file must produce the same rule id for the same input. tests/test_policy.c checks
 * that on the same fixtures the Python self-test uses, and tools/vbs_policy_dump lets the two
 * be compared over any base cache and any set of positions. Change a constant here without
 * changing it there and the comparison fails, which is the point of having two.
 *
 * Rules are evaluated in order and the first that fires decides.
 *
 *   R0  density      fewer than MIN_BASES_IN_RADIUS bases within CLUSTER_RADIUS_KM   -> 404
 *   R1  covered      nearest base < COVERED_KM, or MAX_REDUNDANT+ bases within
 *                    REDUNDANCY_KM of the nearest (below the Mosaic floor)         -> 404
 *   R2  homogeneous  a single-family cluster whose nearest station is inside that
 *                    family's envelope (UM980 13-30 km, D133; Mosaic 28-60 km, D135);
 *                    nearest qualifying family wins, Mosaic breaks an exact tie     -> serve
 *   R2s capability   a station of that cluster delivers fewer than MIN_GPS_BANDS      -> 404
 *   R3  all-F9P      an F9P cluster and no Mosaic or UM980 cluster within reach       -> 404
 *   R4  mixed        the nearest five span more than one family                       -> 404
 *   R5  out of range nearest base outside every envelope                              -> 404
 */

#define VBS_FAM_MOSAIC	0
#define VBS_FAM_UM980	1
#define VBS_FAM_F9P	2
#define VBS_FAM_OTHER	3
#define VBS_NFAM	4

/* One base as the policy sees it. */
struct vbs_policy_base {
	const char *mp;
	int family;			/* VBS_FAM_* */
	double xyz[3];
	int gps_bands;			/* distinct GPS bands among L1/L2/L5 delivered; 0 = unknown */
	const char *fw;			/* receiver firmware from RTCM 1033; NULL/empty = unknown */
	int live;			/* 0: not currently streaming, ignored entirely */
};

/*
 * Firmware builds known to produce unusable observations. Checked ALONGSIDE the band count, not
 * instead of it (D159): the band check catches an unknown bad build by its symptom after enough
 * epochs, and this catches a known one at connect time. Dropping either would lose a case the
 * other cannot see.
 *
 * Maintenance: this is a blocklist, so it only ever refuses what is already known bad. A build
 * absent from it is not thereby endorsed.
 */
int vbs_policy_fw_blocked(const char *fw);

/*
 * The build number in a firmware string ("2.7.0 R4.10Build17548" -> 17548), or -1 when it has
 * none. The one parser for it: the blocklist above and the per-build grouping of vbs_c2wcal both
 * read it here, so the two cannot disagree on which build a base runs. A RINEX header truncates
 * the string by one digit (D159); the result is then wrong as a version but still distinct per
 * build, which is all a grouping needs.
 */
long vbs_policy_fw_build(const char *fw);

#define VBS_POLICY_MAXCLUSTER	5

struct vbs_decision {
	int serve;
	char rule[4];			/* "R0".."R5"; R2s reports as "R2" with serve = 0 */
	int family;			/* the serving family, when serve */
	int dispersive;			/* apply the dispersive correction */
	int prefer_l2c;			/* GPS L2 from L2C: Mosaic yes (D138), UM980 no (D140) */
	int iono_sys;			/* dispersive on these constellations only (SYS_*); 0 = all */
	int glonass;
	double nearest_km;		/* to the nearest live base of any family */
	int cluster[VBS_POLICY_MAXCLUSTER];	/* indices into the caller's base array, nearest first */
	int ncluster;
	char reason[160];
};

double vbs_policy_cluster_radius_km(void);
double vbs_policy_backup_radius_km(void);	/* "2 + 1 backup": how far a reserve may be */

/*
 * The distance window R2 applies to a family's nearest station: 1 and the bounds when the
 * family is served at all, 0 otherwise. Exported so that tools measuring what the server
 * would do read the window from here rather than restating it.
 */
int vbs_policy_envelope(int family, double *lo_km, double *hi_km);

/* Receiver-type string -> family, the same mapping as the Python. */
int vbs_policy_family(const char *rectype);
const char *vbs_policy_family_str(int family);

/*
 * vbs_policy_decide() flags. VBS_POLICY_CHECK_BANDS enables R2s; a base whose gps_bands is 0
 * (unknown) is not refused on that ground, only one measured below the minimum.
 * VBS_POLICY_NO_BACKUP withholds the "2 + 1 backup" reserve: no density relief from it, no
 * short family admitted on it, no reserve candidate offered to the cluster module.
 * VBS_POLICY_FORCE is for field tests only (vbsd --force): R1 is skipped and the R2 envelopes
 * become the whole 120 km reach, so a homogeneous Mosaic or UM980 cluster is served wherever one
 * can form. R0, R2s, R3 and R4 still refuse, and the reason starts with "forced:". The Python
 * reference has no such flag: production never sets it.
 */
#define VBS_POLICY_CHECK_BANDS	0x1
#define VBS_POLICY_NO_BACKUP	0x2
#define VBS_POLICY_FORCE	0x4

/* Decide for a position. */
void vbs_policy_decide(const double *rover_xyz, const struct vbs_policy_base *bases, int nbase,
	int flags, struct vbs_decision *out);

#endif /* __VBS_POLICY_H__ */

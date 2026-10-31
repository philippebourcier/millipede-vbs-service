#ifndef __VBS_CLUSTER_H__
#define __VBS_CLUSTER_H__

#include "vbs/vbs_basemeta.h"

/*
 * Cluster and master selection -- spec §4.1.
 *
 * Given a VBS point and the available bases, choose a master and its auxiliaries, and say
 * honestly whether the VBS is inside the cluster or outside it and by how far.
 *
 * **The extrapolation metric is the deliverable, not a footnote.** `01-design-study.md` §4
 * found that 81 % of the area where VBS is needed cannot form an enclosing cluster, so
 * extrapolation is the ordinary case on this network rather than a degraded corner. D17 item 4
 * requires `vbs.max_extrapolation_km` to be *calibrated* as the largest out-of-hull distance
 * still meeting the targets, and beyond it the caster falls back to nearest-base pass-through
 * rather than serving degraded corrections. Both numbers this module publishes -- distance
 * outside the hull, and the interpolation leverage that distance implies -- exist to be
 * compared against that calibration at run time.
 *
 * **Master selection.** The master is the station whose observations are *translated* (§4.6),
 * so every VBS observable is its observable plus corrections. Two consequences drive the
 * choice: the master must be near the VBS, because Δρ and the correction magnitudes both grow
 * with that distance; and the master's signal set is the VBS's signal set, since nothing can
 * be synthesized that the master did not observe. A rover with more capable hardware than the
 * master therefore gains nothing from it, which is why the signal count enters the score.
 *
 * **Same frame only (D13).** Bases whose frame differs from the master's are excluded, not
 * transformed. The transform between realisations is small, epoch-dependent, and would enter
 * the corrections as a fictitious geometric error with no time signature to reveal it.
 */

#define VBS_CLUSTER_MAXAUX	12

struct vbs_cluster_opt {
	int max_aux;			/* auxiliaries to keep, at most VBS_CLUSTER_MAXAUX */
	int min_aux;			/* below this, no VBS can be formed */
	double max_baseline_km;		/* refuse auxiliaries further than this from the master */
	double max_extrapolation_km;	/* 0 disables the cap (calibration runs need it off) */
	double min_separation_km;	/* two bases closer than this add no geometry */
	/*
	 * "2 + 1 backup". When fewer than min_aux auxiliaries qualify but at least
	 * backup_min_aux do, the cluster is still admitted with one more station of the master's
	 * frame held in reserve -- the nearest remaining one within backup_max_km of the master,
	 * which may be far, because it is meant to be used almost never. It is processed
	 * continuously (an ambiguity cannot be fixed the moment it is needed) but the
	 * interpolator keeps it out of the fit unless the primaries for a satellite fall below
	 * vbs_interp_opt.min_aux. 0 disables it.
	 */
	int backup_min_aux;
	double backup_max_km;
	/*
	 * A regular auxiliary must lie within this distance of the VBS; 0 = no limit. The server
	 * sets it to the policy's cluster radius, because the policy may now offer one station
	 * from beyond that radius as a backup candidate, and without this bound the auxiliary loop
	 * would take it as an ordinary auxiliary -- a distant station permanently in the fit,
	 * which is exactly what the reserve exists to avoid.
	 */
	double aux_max_vbs_km;
};

void vbs_cluster_default_opt(struct vbs_cluster_opt *opt);

struct vbs_cluster {
	int master;			/* index into the caller's base array */
	int aux[VBS_CLUSTER_MAXAUX];
	int naux;

	int inside_hull;		/* VBS within the convex hull of master + auxiliaries */
	double extrapolation_km;	/* 0 when inside; distance to the hull otherwise */
	double leverage;		/* interpolation variance inflation at the VBS */
	double max_baseline_km;		/* longest master->auxiliary baseline */
	double master_dist_km;

	char frame[VBS_BASE_MAXFRAME];
	int nsignal_master;		/* signals the master carries, i.e. the VBS's ceiling */

	/*
	 * The "2 + 1 backup" station, or -1. Never in aux[]: the geometry above (hull,
	 * extrapolation, leverage) is that of the primaries, which is what the fit uses whenever
	 * they are all present.
	 */
	int backup;
	double backup_dist_km;		/* VBS to backup */
	double backup_baseline_km;	/* master to backup */
};

/* Why no cluster could be formed. */
#define VBS_CLUSTER_OK			0
#define VBS_CLUSTER_NO_MASTER		1
#define VBS_CLUSTER_TOO_FEW_AUX		2
#define VBS_CLUSTER_TOO_FAR		3	/* beyond max_extrapolation_km */
#define VBS_CLUSTER_NREASON		4

const char *vbs_cluster_reason_str(int reason);

/*
 * Select a cluster for `vbs_xyz` from `bases` (only those whose gate reason is VBS_BASE_OK are
 * considered). `nsignal` may be NULL; when given, nsignal[i] is the number of canonical
 * signals base i carries, and is used to prefer a richer master.
 *
 * Returns VBS_CLUSTER_OK and fills `out`, or a VBS_CLUSTER_* reason.
 */
int vbs_cluster_select(const double *vbs_xyz, const struct vbs_base *bases, int nbase,
	const int *nsignal, const struct vbs_cluster_opt *opt, struct vbs_cluster *out);

/*
 * Distance from a point to the convex hull of a set of points, all in a local ENU plane, in
 * the same units as the input. Returns 0 when the point is inside. Exposed because the
 * extrapolation metric is worth testing directly rather than only through cluster selection.
 */
double vbs_cluster_hull_distance(const double *e, const double *n, int m,
	double pe, double pn, int *inside);

#endif /* __VBS_CLUSTER_H__ */

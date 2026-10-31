#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_cluster.h"

/*
 * Cluster and master selection -- spec §4.1. See the header.
 */

#define MAXPT	(VBS_CLUSTER_MAXAUX + 1)

static const char *REASON[VBS_CLUSTER_NREASON] = {
	"ok", "no usable master", "too few auxiliaries", "beyond max_extrapolation_km"
};

const char *vbs_cluster_reason_str(int reason) {
	return (reason >= 0 && reason < VBS_CLUSTER_NREASON) ? REASON[reason] : "?";
}

void vbs_cluster_default_opt(struct vbs_cluster_opt *opt) {
	memset(opt, 0, sizeof *opt);
	opt->max_aux = 4;
	opt->min_aux = 3;		/* a plane has three parameters (§4.5) */
	/*
	 * Two primaries and a reserve (decided 2026-10-02). Two auxiliaries plus
	 * the master anchor still determine the plane exactly, and the reserve keeps a satellite
	 * corrected when one primary loses it. Far is acceptable for a station used this rarely.
	 */
	opt->backup_min_aux = 2;
	/*
	 * This caps the *network* baseline, master->auxiliary, which for an enclosed cluster is
	 * roughly twice the nominal spacing: a "150 km" case means the target is 150 km from
	 * each station, so two stations on opposite sides are ~300 km apart. Setting this to the
	 * spacing bin refused every long enclosed case in the dataset outright, which is how the
	 * confusion surfaced. 320 km admits the dataset's longest geometry; §4.3 warns that
	 * fixing degrades over this range, and measuring that degradation is the point.
	 */
	opt->max_baseline_km = 320.0;
	/*
	 * 0 by design: the cap exists to be *calibrated* (D17 item 4), and a calibration run
	 * that silently refused the cases beyond the current guess would only ever confirm it.
	 * The caster applies the cap at run time; the engine reports the distance.
	 */
	opt->max_extrapolation_km = 0.0;
	opt->min_separation_km = 1.0;
	/*
	 * Twice the auxiliary cap. In the server the effective limit is the policy's own backup
	 * radius (vbs_policy.c), which decides which stations are offered at all; this bound only
	 * keeps a caller that hands the module a whole network from reaching across a continent.
	 */
	opt->backup_max_km = 2.0 * opt->max_baseline_km;
}

static double km(const double *a, const double *b) {
	double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};

	return norm(d, 3) / 1000.0;
}

/* Monotone-chain convex hull. Writes indices into `hull`, returns the count. */
static int convex_hull(const double *e, const double *n, int m, int *hull) {
	int idx[MAXPT], i, k = 0, t;

	for (i = 0; i < m; i++)
		idx[i] = i;
	/* Sort by (e, n) -- m is at most 13, so insertion sort is the honest choice. */
	for (i = 1; i < m; i++) {
		int cur = idx[i], j = i - 1;

		while (j >= 0 && (e[idx[j]] > e[cur] ||
		    (e[idx[j]] == e[cur] && n[idx[j]] > n[cur]))) {
			idx[j + 1] = idx[j];
			j--;
		}
		idx[j + 1] = cur;
	}
#define CROSS(o, a, b)	((e[a] - e[o]) * (n[b] - n[o]) - (n[a] - n[o]) * (e[b] - e[o]))
	for (i = 0; i < m; i++) {
		while (k >= 2 && CROSS(hull[k - 2], hull[k - 1], idx[i]) <= 0.0)
			k--;
		hull[k++] = idx[i];
	}
	for (i = m - 2, t = k + 1; i >= 0; i--) {
		while (k >= t && CROSS(hull[k - 2], hull[k - 1], idx[i]) <= 0.0)
			k--;
		hull[k++] = idx[i];
	}
#undef CROSS
	return k > 1 ? k - 1 : k;	/* the last point repeats the first */
}

/* Distance from p to the segment ab. */
static double seg_dist(double ae, double an, double be, double bn, double pe, double pn) {
	double de = be - ae, dn = bn - an, len2 = de * de + dn * dn, t;

	if (len2 <= 0.0)
		return hypot(pe - ae, pn - an);
	t = ((pe - ae) * de + (pn - an) * dn) / len2;
	t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
	return hypot(pe - (ae + t * de), pn - (an + t * dn));
}

double vbs_cluster_hull_distance(const double *e, const double *n, int m,
	double pe, double pn, int *inside) {
	int hull[MAXPT * 2], nh, i, in = 1;
	double best = 1e18;

	if (inside != NULL)
		*inside = 0;
	if (m < 3) {
		/*
		 * Fewer than three points span no area, so nothing is "inside" -- report the
		 * distance to the nearest of them rather than pretending containment.
		 */
		for (i = 0; i < m; i++) {
			double d = hypot(pe - e[i], pn - n[i]);

			if (d < best)
				best = d;
		}
		return m > 0 ? best : 0.0;
	}
	nh = convex_hull(e, n, m, hull);
	if (nh < 3) {			/* collinear: no area, same treatment */
		for (i = 0; i < m; i++) {
			double d = hypot(pe - e[i], pn - n[i]);

			if (d < best)
				best = d;
		}
		return best;
	}
	/*
	 * The hull comes out counter-clockwise, so an interior point is left of every edge.
	 * The same loop measures the distance to each edge, which is what the point needs when
	 * it turns out to be outside.
	 */
	for (i = 0; i < nh; i++) {
		int a = hull[i], b = hull[(i + 1) % nh];
		double cross = (e[b] - e[a]) * (pn - n[a]) - (n[b] - n[a]) * (pe - e[a]);
		double d = seg_dist(e[a], n[a], e[b], n[b], pe, pn);

		if (cross < 0.0)
			in = 0;
		if (d < best)
			best = d;
	}
	if (inside != NULL)
		*inside = in;
	return in ? 0.0 : best;
}

/* Interpolation leverage at the VBS: A^-1[0,0] of the unit-weight plane fit (see vbs_interp). */
static double leverage(const double *e, const double *n, int m, double pe, double pn) {
	double A[9] = {0};
	int i, r, c;

	for (i = 0; i < m; i++) {
		double g[3] = {1.0, e[i] - pe, n[i] - pn};

		for (r = 0; r < 3; r++)
			for (c = 0; c < 3; c++)
				A[r + c * 3] += g[r] * g[c];
	}
	if (m < 3 || matinv(A, 3) != 0 || A[0] <= 0.0)
		return m > 0 ? 1.0 / m : 0.0;
	return A[0];
}

int vbs_cluster_select(const double *vbs_xyz, const struct vbs_base *bases, int nbase,
	const int *nsignal, const struct vbs_cluster_opt *opt, struct vbs_cluster *out) {
	struct vbs_cluster_opt defopt;
	double e[MAXPT], n[MAXPT], vpos[3], rel[3], denu[3];
	int cand[64], ncand = 0, i, j, best = -1, maxaux;
	double bestscore = -1e18;

	if (opt == NULL) {
		vbs_cluster_default_opt(&defopt);
		opt = &defopt;
	}
	maxaux = opt->max_aux > VBS_CLUSTER_MAXAUX ? VBS_CLUSTER_MAXAUX : opt->max_aux;
	memset(out, 0, sizeof *out);
	out->master = -1;
	out->backup = -1;		/* memset made it 0, which is a valid index */

	for (i = 0; i < nbase && ncand < (int)(sizeof cand / sizeof cand[0]); i++)
		if (bases[i].reason == VBS_BASE_OK)
			cand[ncand++] = i;
	if (ncand == 0)
		return VBS_CLUSTER_NO_MASTER;

	/*
	 * Master: nearest to the VBS, with a bonus for carrying more signals. The distance term
	 * dominates -- 1 km of extra separation costs more than one extra signal is worth --
	 * because Δρ and the interpolated corrections both grow with it, while an extra signal
	 * only helps rovers that can use it. The bonus breaks ties between comparable stations
	 * in favour of the one that does not cap what the VBS can offer.
	 */
	for (j = 0; j < ncand; j++) {
		double d = km(vbs_xyz, bases[cand[j]].xyz);
		double score = -d + (nsignal != NULL ? 0.5 * nsignal[cand[j]] : 0.0);

		if (score > bestscore) {
			bestscore = score;
			best = cand[j];
		}
	}
	if (best < 0)
		return VBS_CLUSTER_NO_MASTER;
	out->master = best;
	out->master_dist_km = km(vbs_xyz, bases[best].xyz);
	snprintf(out->frame, sizeof out->frame, "%s", bases[best].frame);
	out->nsignal_master = nsignal != NULL ? nsignal[best] : 0;

	/*
	 * Auxiliaries: nearest to the *VBS*, same frame, not co-located with an already chosen
	 * station, and within max_baseline_km of the master.
	 *
	 * Nearest to the VBS rather than to the master, because what the auxiliaries are for is
	 * §4.5's interpolation *at the VBS* -- choosing them around the master optimises the
	 * wrong point and, for a VBS off to one side of the master, would happily pick a set
	 * that leaves it outside the hull. The master-baseline cap stays as a separate guard,
	 * since that length is what drives ambiguity resolution difficulty.
	 *
	 * Co-located pairs (the network has several, a metre or two apart) add a row to the
	 * plane fit without adding geometry, inflating apparent redundancy while leaving the fit
	 * as badly determined as before.
	 */
	for (;;) {
		int pick = -1;
		double pickd = 1e18;

		for (j = 0; j < ncand; j++) {
			int b = cand[j], taken = 0;
			double d;

			if (b == out->master)
				continue;
			if (strcmp(bases[b].frame, out->frame) != 0)
				continue;
			for (i = 0; i < out->naux; i++)
				if (out->aux[i] == b)
					taken = 1;
			if (taken)
				continue;
			if (km(bases[out->master].xyz, bases[b].xyz) > opt->max_baseline_km)
				continue;
			if (km(bases[out->master].xyz, bases[b].xyz) < opt->min_separation_km)
				continue;
			d = km(vbs_xyz, bases[b].xyz);
			if (opt->aux_max_vbs_km > 0.0 && d > opt->aux_max_vbs_km)
				continue;	/* beyond reach: a backup at most, never an auxiliary */
			for (i = 0; i < out->naux; i++)
				if (km(bases[out->aux[i]].xyz, bases[b].xyz) <
				    opt->min_separation_km)
					taken = 1;
			if (taken)
				continue;
			if (d < pickd) {
				pickd = d;
				pick = b;
			}
		}
		if (pick < 0 || out->naux >= maxaux)
			break;
		out->aux[out->naux++] = pick;
		{
			double bl = km(bases[out->master].xyz, bases[pick].xyz);

			if (bl > out->max_baseline_km)
				out->max_baseline_km = bl;
		}
	}
	if (out->naux < opt->min_aux) {
		int pick = -1;
		double pickd = 1e18;

		if (opt->backup_min_aux <= 0 || out->naux < opt->backup_min_aux)
			return VBS_CLUSTER_TOO_FEW_AUX;
		/*
		 * "2 + 1 backup": the loop above ran out of stations within max_baseline_km, so look
		 * once more, out to backup_max_km, for the nearest remaining station of the same
		 * frame that adds geometry. Same rules as an auxiliary otherwise -- nearest to the
		 * VBS, not co-located with anything already chosen.
		 */
		for (j = 0; j < ncand; j++) {
			int b = cand[j], taken = 0;
			double d, bl;

			if (b == out->master || strcmp(bases[b].frame, out->frame) != 0)
				continue;
			for (i = 0; i < out->naux; i++)
				if (out->aux[i] == b ||
				    km(bases[out->aux[i]].xyz, bases[b].xyz) < opt->min_separation_km)
					taken = 1;
			if (taken)
				continue;
			bl = km(bases[out->master].xyz, bases[b].xyz);
			if (bl > opt->backup_max_km || bl < opt->min_separation_km)
				continue;
			d = km(vbs_xyz, bases[b].xyz);
			if (d < pickd) {
				pickd = d;
				pick = b;
			}
		}
		if (pick < 0)
			return VBS_CLUSTER_TOO_FEW_AUX;
		out->backup = pick;
		out->backup_dist_km = pickd;
		out->backup_baseline_km = km(bases[out->master].xyz, bases[pick].xyz);
	}

	/*
	 * Containment is evaluated over master *and* auxiliaries: the master is a corner of the
	 * cluster like any other, and leaving it out would report a VBS sitting right on the
	 * master as an extrapolation whenever the auxiliaries happened to lie to one side.
	 */
	ecef2pos(vbs_xyz, vpos);
	for (i = 0; i < 3; i++)
		rel[i] = bases[out->master].xyz[i] - vbs_xyz[i];
	ecef2enu(vpos, rel, denu);
	e[0] = denu[0] / 1000.0;
	n[0] = denu[1] / 1000.0;
	for (j = 0; j < out->naux; j++) {
		for (i = 0; i < 3; i++)
			rel[i] = bases[out->aux[j]].xyz[i] - vbs_xyz[i];
		ecef2enu(vpos, rel, denu);
		e[j + 1] = denu[0] / 1000.0;
		n[j + 1] = denu[1] / 1000.0;
	}
	out->extrapolation_km = vbs_cluster_hull_distance(e, n, out->naux + 1, 0.0, 0.0,
	    &out->inside_hull);
	out->leverage = leverage(e, n, out->naux + 1, 0.0, 0.0);

	if (opt->max_extrapolation_km > 0.0 &&
	    out->extrapolation_km > opt->max_extrapolation_km)
		return VBS_CLUSTER_TOO_FAR;
	return VBS_CLUSTER_OK;
}

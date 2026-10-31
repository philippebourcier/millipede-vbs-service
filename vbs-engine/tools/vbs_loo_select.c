/*
 * vbs_loo_select -- the server's own cluster choice, made once per station with that station
 * left out.
 *
 * Calibrating `vbs.max_extrapolation_km` (D17 item 3: "the largest distance outside the hull
 * where fix rate >= 90 % and false fix <= 1 %") needs many leave-one-out cases whose targets
 * sit at different distances outside their cluster's hull. Choosing those clusters by hand, or
 * in Python, would calibrate the cap against clusters the server never forms: on 2026-10-02 a
 * Python re-implementation of the serving rules disagreed with vbsd on 283 of 5849 tiles. So
 * every step that decides membership is the engine's own:
 *
 *   family     vbs_policy_family() on the receiver string, the mapping vbsd uses;
 *   candidates the nearest VBS_POLICY_MAXCLUSTER stations of that family within
 *              vbs_policy_cluster_radius_km() -- what R2 hands to the cluster module;
 *   cluster    vbs_cluster_select() over those candidates, default options.
 *
 * Two departures, both deliberate. The target is removed from the network, which is the point
 * of leave-one-out. And the extrapolation cap is off, because it is what is being calibrated.
 * The policy's commercial rules (R1 "already covered", the family distance envelope) are NOT
 * applied as filters: the calibration needs the inside-hull cases R1 would discard as its
 * reference, and the envelope is reported per case (`in_envelope`) so the analysis can
 * restrict to what is served without the tool deciding it.
 *
 * Input, one station per line ('#' starts a comment):
 *
 *     MP MARKER X Y Z FRAME RECTYPE [SUSPECT]
 *
 * RECTYPE is the receiver type with spaces written as underscores. SUSPECT marks a station
 * whose declared receiver its observations contradict (an "F9P" carrying band 5); it is never
 * part of a network, since its family is unknown, but it can still be a target. FRAME: a
 * cluster only takes stations in its master's frame (D13), so a file mixing frames calibrates
 * nothing about frames and the caller should not mix them.
 *
 * Output: one JSON object per target on stdout -- the cluster, or why none could be formed.
 *
 * --family MIX takes UM980 and Mosaic stations together (no F9P): the cluster a mixed-family
 * server would form. vbsd refuses every mix (R4), so this is a study of what such a rule would
 * serve, not of what is served; it reports no envelope.
 *
 * Usage:
 *   vbs_loo_select --family MOSAIC|UM980|F9P|MIX [--whole-network] [--min-aux N]
 *                  [--max-aux N] [--target MP] FILE
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_basemeta.h"
#include "vbs/vbs_cluster.h"
#include "vbs/vbs_policy.h"

#define MAXST		4096

struct st {
	char mp[VBS_BASE_MAXNAME];
	char marker[16];
	char frame[VBS_BASE_MAXFRAME];
	char rec[64];
	double xyz[3];
	int family;			/* VBS_FAM_*, VBS_FAM_OTHER when suspect */
	int suspect;
};

static struct st S[MAXST];

struct ranked {
	int i;
	double d;
};

static double km(const double *a, const double *b) {
	return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
	    (a[2] - b[2]) * (a[2] - b[2])) / 1000.0;
}

static int cmp_ranked(const void *a, const void *b) {
	double x = ((const struct ranked *)a)->d, y = ((const struct ranked *)b)->d;

	return x < y ? -1 : x > y ? 1 : 0;
}

static int load(const char *path) {
	char line[512], flag[16];
	FILE *f = fopen(path, "r");
	int n = 0, k;

	if (f == NULL) {
		perror(path);
		return -1;
	}
	while (fgets(line, sizeof line, f) != NULL) {
		char *c = strchr(line, '#');
		struct st *s = &S[n];

		if (c != NULL)
			*c = '\0';
		if (n >= MAXST) {
			fprintf(stderr, "%s: more than %d stations\n", path, MAXST);
			break;
		}
		flag[0] = '\0';
		k = sscanf(line, "%31s %15s %lf %lf %lf %15s %63s %15s", s->mp, s->marker,
		    &s->xyz[0], &s->xyz[1], &s->xyz[2], s->frame, s->rec, flag);
		if (k < 7)
			continue;
		for (c = s->rec; *c; c++)
			if (*c == '_')
				*c = ' ';
		s->suspect = !strcmp(flag, "SUSPECT");
		s->family = s->suspect ? VBS_FAM_OTHER : vbs_policy_family(s->rec);
		n++;
	}
	fclose(f);
	return n;
}

/* MIX: UM980 and Mosaic together; otherwise the one family */
static int in_network(int fam, int family, int mix) {
	return mix ? fam == VBS_FAM_UM980 || fam == VBS_FAM_MOSAIC : fam == family;
}

static void emit_member(const struct st *s, const double *vbs, int first) {
	printf("%s{\"mp\":\"%s\",\"marker\":\"%s\",\"xyz\":[%.4f,%.4f,%.4f],\"family\":\"%s\","
	    "\"dist_km\":%.3f}", first ? "" : ",", s->mp, s->marker, s->xyz[0], s->xyz[1],
	    s->xyz[2], vbs_policy_family_str(s->family), km(s->xyz, vbs));
}

int main(int argc, char **argv) {
	static struct ranked rk[MAXST];
	struct vbs_cluster_opt opt;
	struct vbs_base bases[MAXST];
	int idx[MAXST];
	const char *only = NULL, *path = NULL, *famname = NULL;
	int family = -1, mix = 0, whole = 0, n, t, i, f;
	double radius = vbs_policy_cluster_radius_km(), lo = 0.0, hi = 0.0;

	vbs_cluster_default_opt(&opt);
	opt.max_extrapolation_km = 0.0;		/* the quantity under calibration: never cap it */
	opt.aux_max_vbs_km = radius;		/* as vbsd: a reserve from beyond reach is never an auxiliary */
	for (i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--family") && i + 1 < argc)
			famname = argv[++i];
		else if (!strcmp(argv[i], "--whole-network"))
			whole = 1;
		else if (!strcmp(argv[i], "--min-aux") && i + 1 < argc)
			opt.min_aux = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--max-aux") && i + 1 < argc)
			opt.max_aux = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--target") && i + 1 < argc)
			only = argv[++i];
		else if (argv[i][0] != '-')
			path = argv[i];
		else {
			fprintf(stderr, "unknown option %s\n", argv[i]);
			return 2;
		}
	}
	for (f = 0; famname != NULL && f < VBS_NFAM; f++)
		if (!strcmp(famname, vbs_policy_family_str(f)))
			family = f;
	mix = famname != NULL && !strcmp(famname, "MIX");
	if (path == NULL || (family < 0 && !mix)) {
		fprintf(stderr, "usage: vbs_loo_select --family MOSAIC|UM980|F9P|MIX "
		    "[--whole-network] [--min-aux N] [--max-aux N] [--target MP] FILE\n");
		return 2;
	}
	if ((n = load(path)) <= 0)
		return 1;
	if (!mix)
		vbs_policy_envelope(family, &lo, &hi);

	for (t = 0; t < n; t++) {
		const struct st *T = &S[t];
		struct vbs_cluster cl;
		int nr = 0, nb = 0, in_radius = 0, r;
		double nearest = -1.0;

		if (only != NULL && strcmp(only, T->mp) != 0)
			continue;
		/* the family's stations other than the target, nearest first */
		for (i = 0; i < n; i++)
			if (i != t && in_network(S[i].family, family, mix)) {
				rk[nr].i = i;
				rk[nr].d = km(S[i].xyz, T->xyz);
				nr++;
			}
		qsort(rk, nr, sizeof rk[0], cmp_ranked);
		if (nr > 0)
			nearest = rk[0].d;
		for (i = 0; i < nr; i++)
			if (rk[i].d <= radius)
				in_radius++;
		/* what R2 hands the cluster module: the nearest few within reach */
		for (i = 0; i < nr; i++) {
			if (!whole && (rk[i].d > radius || nb >= VBS_POLICY_MAXCLUSTER))
				break;
			vbs_base_init(&bases[nb], S[rk[i].i].mp, S[rk[i].i].frame, S[rk[i].i].xyz);
			idx[nb++] = rk[i].i;
		}
		/* ... and, when there is room, the "2 + 1 backup" candidate from beyond reach */
		if (!whole && nb < VBS_POLICY_MAXCLUSTER)
			for (; i < nr; i++)
				if (rk[i].d > radius && rk[i].d <= vbs_policy_backup_radius_km()) {
					vbs_base_init(&bases[nb], S[rk[i].i].mp, S[rk[i].i].frame,
					    S[rk[i].i].xyz);
					idx[nb++] = rk[i].i;
					break;
				}
		printf("{\"target\":\"%s\",\"target_marker\":\"%s\",\"target_xyz\":[%.4f,%.4f,%.4f],"
		    "\"target_family\":\"%s\",\"target_rec\":\"%s\",\"network_family\":\"%s\","
		    "\"family_in_radius\":%d,\"candidates\":%d,\"nearest_family_km\":%.3f,"
		    "\"in_envelope\":%d,",
		    T->mp, T->marker, T->xyz[0], T->xyz[1], T->xyz[2],
		    T->suspect ? "SUSPECT" : vbs_policy_family_str(T->family), T->rec,
		    mix ? "MIX" : vbs_policy_family_str(family), in_radius, nb, nearest,
		    hi > 0.0 && nearest >= lo && nearest <= hi);
		r = vbs_cluster_select(T->xyz, bases, nb, NULL, &opt, &cl);
		if (r != VBS_CLUSTER_OK) {
			printf("\"ok\":false,\"reason\":\"%s\"}\n", vbs_cluster_reason_str(r));
			continue;
		}
		printf("\"ok\":true,\"frame\":\"%s\",\"naux\":%d,\"inside_hull\":%d,"
		    "\"extrapolation_km\":%.3f,\"leverage\":%.4f,\"max_baseline_km\":%.3f,"
		    "\"master_dist_km\":%.3f,\"master\":", cl.frame, cl.naux, cl.inside_hull,
		    cl.extrapolation_km, cl.leverage, cl.max_baseline_km, cl.master_dist_km);
		emit_member(&S[idx[cl.master]], T->xyz, 1);
		printf(",\"aux\":[");
		for (i = 0; i < cl.naux; i++)
			emit_member(&S[idx[cl.aux[i]]], T->xyz, i == 0);
		printf("],\"backup\":");
		if (cl.backup >= 0)
			emit_member(&S[idx[cl.backup]], T->xyz, 1);
		else
			printf("null");
		printf("}\n");
	}
	return 0;
}

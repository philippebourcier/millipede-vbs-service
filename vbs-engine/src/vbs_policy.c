#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <strings.h>

#include "rtklib.h"
#include "vbs/vbs_policy.h"

/*
 * Constants, verbatim from data/scripts/vbs_policy.py. The evidence for each is there; here
 * they are only values, and the cross-check test is what keeps them equal.
 */
#define COVERED_KM		13.0
#define CLUSTER_RADIUS_KM	120.0
#define MIN_BASES_IN_RADIUS	4
#define REDUNDANCY_KM		12.0
#define MAX_REDUNDANT		2
/*
 * R1's redundancy clause stops here. It used to read the Mosaic floor, which was the same 20 km;
 * D135 raised that floor to 28 km for an unrelated reason, and R1 must not follow it.
 */
#define REDUNDANT_MAX_KM	20.0
#define MIN_GPS_BANDS		3

/*
 * Builds whose observations are unusable. D159: on these, SIGNALGROUP 2 declares a GPS L2C column
 * and never fills it, while tracking L2 P(Y) on about a fifth of satellite-records. Build 17548
 * runs the identical configuration on 77 stations with no defect, so this is the firmware and not
 * the configuration.
 *
 * Matched on the BUILD NUMBER, extracted and compared as an integer rather than by substring:
 * "Build1183" is a prefix of "Build11833", so a strstr blocklist would refuse a build nobody has
 * shown to be bad. The RINEX header truncates this string at 20 characters and drops the last
 * digit, which is how the number was first misreported -- another reason not to match prefixes.
 */
static const long fw_bad_builds[] = { 11833, 13504 };

long vbs_policy_fw_build(const char *fw) {
	const char *p;

	if (fw == NULL || (p = strstr(fw, "Build")) == NULL)
		return -1;
	p += 5;
	if (*p < '0' || *p > '9')
		return -1;
	return strtol(p, NULL, 10);
}

int vbs_policy_fw_blocked(const char *fw) {
	long n = vbs_policy_fw_build(fw);
	size_t i;

	if (n < 0)
		return 0;		/* unknown firmware is not a finding */
	for (i = 0; i < sizeof fw_bad_builds / sizeof fw_bad_builds[0]; i++)
		if (n == fw_bad_builds[i])
			return 1;
	return 0;
}
#define MIN_AUXILIARIES		3
/*
 * "2 + 1 backup" -- decided 2026-10-02. A homogeneous family with only its
 * master and two auxiliaries within CLUSTER_RADIUS_KM still forms a cluster when one more of
 * its stations lies within BACKUP_RADIUS_KM: that station is processed but held in reserve,
 * and fitted only for a satellite whose primaries fall short (vbs_interp). The radius is
 * generous on purpose -- the reserve should almost never be used, and when it is, a far
 * station still beats losing the satellite. Mixed and all-F9P clusters are refused before
 * this matters, so only the served families (Mosaic, UM980) can use it.
 */
#define BACKUP_MIN_AUX		2
#define BACKUP_RADIUS_KM	250.0

/*
 * UM980 13-40 km with the dispersive correction on Galileo + BeiDou only (D167; was 13-30 km
 * without it, D133, D142), Mosaic 28-60 km (D135); prefer_l2c Mosaic only (D138, D140)
 */
static const struct {
	double lo, hi;
	int dispersive, prefer_l2c;
	int iono_sys;		/* with dispersive: the constellations carrying it, SYS_*; 0 = all */
} ENVELOPE[VBS_NFAM] = {
	[VBS_FAM_MOSAIC] = {28.0, 60.0, 1, 1, 0},
	[VBS_FAM_UM980] = {13.0, 40.0, 1, 0, SYS_GAL | SYS_CMP},
	[VBS_FAM_F9P] = {0.0, 0.0, 0, 0, 0},
	[VBS_FAM_OTHER] = {0.0, 0.0, 0, 0, 0},
};

double vbs_policy_cluster_radius_km(void) { return CLUSTER_RADIUS_KM; }
double vbs_policy_backup_radius_km(void) { return BACKUP_RADIUS_KM; }

int vbs_policy_envelope(int family, double *lo_km, double *hi_km) {
	if (family < 0 || family >= VBS_NFAM || ENVELOPE[family].hi <= 0.0)
		return 0;
	if (lo_km != NULL)
		*lo_km = ENVELOPE[family].lo;
	if (hi_km != NULL)
		*hi_km = ENVELOPE[family].hi;
	return 1;
}

int vbs_policy_family(const char *rectype) {
	char u[64];
	int i;

	if (rectype == NULL)
		return VBS_FAM_OTHER;
	for (i = 0; i < (int)sizeof u - 1 && rectype[i]; i++)
		u[i] = (char)toupper((unsigned char)rectype[i]);
	u[i] = '\0';
	if (strstr(u, "MOSAIC"))
		return VBS_FAM_MOSAIC;
	if (strstr(u, "UM98") || strstr(u, "UNICORE"))
		return VBS_FAM_UM980;
	if (strstr(u, "F9P"))
		return VBS_FAM_F9P;
	return VBS_FAM_OTHER;
}

const char *vbs_policy_family_str(int family) {
	static const char *const names[VBS_NFAM] = {"MOSAIC", "UM980", "F9P", "OTHER"};

	return family >= 0 && family < VBS_NFAM ? names[family] : "?";
}

static double km(const double *a, const double *b) {
	return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
	    (a[2] - b[2]) * (a[2] - b[2])) / 1000.0;
}

struct ranked {
	int idx;
	double d;
};

static int cmp_ranked(const void *a, const void *b) {
	double da = ((const struct ranked *)a)->d, db = ((const struct ranked *)b)->d;

	return da < db ? -1 : da > db ? 1 : 0;
}

static void refuse(struct vbs_decision *out, const char *rule, const char *reason) {
	memset(out, 0, sizeof *out);
	out->serve = 0;
	snprintf(out->rule, sizeof out->rule, "%s", rule);
	snprintf(out->reason, sizeof out->reason, "%s", reason);
}

void vbs_policy_decide(const double *rover_xyz, const struct vbs_policy_base *bases, int nbase,
    int flags, struct vbs_decision *out) {
	struct ranked *all;
	int nall = 0, i, k;
	double nearest;
	int total_in_radius = 0, redundant = 0;
	/* per family: the members within reach, nearest first */
	struct ranked *fam[VBS_NFAM];
	int nfam[VBS_NFAM] = {0, 0, 0, 0};
	int has[VBS_NFAM] = {0, 0, 0, 0};
	int n_reach = 0, has_mix, backup_ok = 0;
	int bk[VBS_NFAM] = {-1, -1, -1, -1};	/* nearest beyond reach, within the backup radius */
	char buf[160];

	memset(out, 0, sizeof *out);
	if (nbase <= 0) {
		refuse(out, "R0", "no bases in range");
		return;
	}
	all = (struct ranked *)malloc(sizeof *all * nbase);
	for (k = 0; k < VBS_NFAM; k++)
		fam[k] = (struct ranked *)malloc(sizeof *all * nbase);
	for (i = 0; i < nbase; i++) {
		if (!bases[i].live)
			continue;
		all[nall].idx = i;
		all[nall].d = km(rover_xyz, bases[i].xyz);
		nall++;
	}
	if (nall == 0) {
		refuse(out, "R0", "no bases in range");
		goto done;
	}
	qsort(all, nall, sizeof *all, cmp_ranked);
	nearest = all[0].d;
	out->nearest_km = nearest;

	/*
	 * density(): the radius count is over the MEASURED families only (a receiver the ladders
	 * never tested cannot form a cluster we know anything about); the redundancy count is over
	 * every live base, because a base of any make is an alternative for direct RTK.
	 */
	for (i = 0; i < nall; i++) {
		int f = bases[all[i].idx].family;

		if (all[i].d <= nearest + REDUNDANCY_KM)
			redundant++;
		if (all[i].d <= CLUSTER_RADIUS_KM) {
			n_reach++;
			if (f != VBS_FAM_OTHER) {
				total_in_radius++;
				fam[f][nfam[f]++] = all[i];
			}
		}
	}

	/*
	 * The backup candidates: all[] is sorted, so the first beyond reach is the nearest. A
	 * station R2s would refuse is passed over for the next one out: offered as the reserve, it
	 * would veto a cluster that serves without it (T0990C0990: CSEL, two GPS bands, refused a
	 * Mosaic cluster of four).
	 */
	for (i = 0; i < nall && !(flags & VBS_POLICY_NO_BACKUP); i++) {
		const struct vbs_policy_base *b = &bases[all[i].idx];
		int f = b->family;

		if (all[i].d <= CLUSTER_RADIUS_KM || all[i].d > BACKUP_RADIUS_KM ||
		    f == VBS_FAM_OTHER || bk[f] >= 0)
			continue;
		if ((flags & VBS_POLICY_CHECK_BANDS) && ((b->gps_bands > 0 &&
		    b->gps_bands < MIN_GPS_BANDS) || vbs_policy_fw_blocked(b->fw)))
			continue;
		bk[f] = i;
	}
	for (k = 0; k < 2; k++)			/* MOSAIC, UM980: the families that can be served */
		if (nfam[k] >= BACKUP_MIN_AUX + 1 && bk[k] >= 0)
			backup_ok = 1;

	/* R0 -- density */
	if (total_in_radius < MIN_BASES_IN_RADIUS && !backup_ok) {
		snprintf(buf, sizeof buf, "only %d base(s) within %.0f km: no cluster is possible here",
		    total_in_radius, CLUSTER_RADIUS_KM);
		refuse(out, "R0", buf);
		out->nearest_km = nearest;
		goto done;
	}

	/* R1 -- covered (a field test forces past it, to measure the VBS where bases are dense) */
	if (nearest < COVERED_KM && !(flags & VBS_POLICY_FORCE)) {
		snprintf(buf, sizeof buf, "nearest base %.1f km < %.0f km: direct RTK wins outright",
		    nearest, COVERED_KM);
		refuse(out, "R1", buf);
		out->nearest_km = nearest;
		goto done;
	}
	if (redundant >= MAX_REDUNDANT && nearest < REDUNDANT_MAX_KM && !(flags & VBS_POLICY_FORCE)) {
		snprintf(buf, sizeof buf, "%d bases within %.0f km of the nearest (%.1f km): not a coverage gap",
		    redundant, REDUNDANCY_KM, nearest);
		refuse(out, "R1", buf);
		out->nearest_km = nearest;
		goto done;
	}

	/* candidate clusters: a family needs MIN_AUXILIARIES + 1 of its own within reach */
	for (k = 0; k < VBS_NFAM; k++)
		has[k] = nfam[k] >= MIN_AUXILIARIES + 1 ||
		    (k < 2 && nfam[k] >= BACKUP_MIN_AUX + 1 && bk[k] >= 0);
	has_mix = n_reach >= MIN_AUXILIARIES + 1;
	if (!has[VBS_FAM_MOSAIC] && !has[VBS_FAM_UM980] && !has[VBS_FAM_F9P] && !has_mix) {
		snprintf(buf, sizeof buf, "fewer than %d auxiliaries in any composition", MIN_AUXILIARIES);
		refuse(out, "R0", buf);
		out->nearest_km = nearest;
		goto done;
	}

	/*
	 * R2 -- the nearest qualifying homogeneous family, Mosaic breaking an exact tie. The
	 * Python sorts on (round(d, 1), 0 for Mosaic else 1), so ties are decided at 100 m.
	 */
	{
		int best = -1;
		double bestd = 0.0;

		for (k = 0; k < 2; k++) {		/* MOSAIC = 0, UM980 = 1: in tiebreak order */
			double d;

			if (!has[k])
				continue;
			d = fam[k][0].d;
			if (!(flags & VBS_POLICY_FORCE) && (d < ENVELOPE[k].lo || d > ENVELOPE[k].hi))
				continue;
			if (best < 0 || floor(d * 10.0 + 0.5) < floor(bestd * 10.0 + 0.5)) {
				best = k;
				bestd = d;
			}
		}
		if (best >= 0) {
			int n = nfam[best] < VBS_POLICY_MAXCLUSTER ? nfam[best] : VBS_POLICY_MAXCLUSTER;
			int with_backup = 0;

			/*
			 * Offer the backup candidate whenever there is room for it, not only when
			 * the family is short: a member refused later -- an outage, a gate -- then
			 * costs a satellite's reserve rather than the whole cluster. The cluster
			 * module takes it only when the primaries fall below min_aux, and the server
			 * bounds regular auxiliaries to CLUSTER_RADIUS_KM so it cannot become one.
			 */
			if (n < VBS_POLICY_MAXCLUSTER && bk[best] >= 0) {
				fam[best][n] = all[bk[best]];
				n++;
				with_backup = 1;
			}

			if (flags & VBS_POLICY_CHECK_BANDS)
				for (i = 0; i < n; i++) {
					const struct vbs_policy_base *b = &bases[fam[best][i].idx];

					if (b->gps_bands > 0 && b->gps_bands < MIN_GPS_BANDS) {
						snprintf(buf, sizeof buf, "%s cluster at %.1f km but %s emits "
						    "%d GPS band(s): the cascade cannot use it",
						    vbs_policy_family_str(best), bestd, b->mp, b->gps_bands);
						refuse(out, "R2", buf);
						out->nearest_km = nearest;
						out->family = best;
						goto done;
					}
					if (vbs_policy_fw_blocked(b->fw)) {
						snprintf(buf, sizeof buf, "%s cluster at %.1f km but %s runs "
						    "firmware \"%s\", a build that declares GPS L2C and "
						    "never fills it (D159)",
						    vbs_policy_family_str(best), bestd, b->mp, b->fw);
						refuse(out, "R2", buf);
						out->nearest_km = nearest;
						out->family = best;
						goto done;
					}
				}
			out->serve = 1;
			snprintf(out->rule, sizeof out->rule, "R2");
			out->family = best;
			out->dispersive = ENVELOPE[best].dispersive;
			out->prefer_l2c = ENVELOPE[best].prefer_l2c;
			out->iono_sys = ENVELOPE[best].iono_sys;
			out->glonass = 1;
			out->nearest_km = nearest;
			for (i = 0; i < n; i++)
				out->cluster[i] = fam[best][i].idx;
			out->ncluster = n;
			snprintf(out->reason, sizeof out->reason, "%shomogeneous %s at %.1f km%s",
			    (flags & VBS_POLICY_FORCE) ? "forced: " : "", vbs_policy_family_str(best), bestd,
			    with_backup && nfam[best] < MIN_AUXILIARIES + 1 ? " (2 + 1 backup)" : "");
			goto done;
		}
	}

	/* R3 -- all-F9P */
	if (has[VBS_FAM_F9P] && !has[VBS_FAM_MOSAIC] && !has[VBS_FAM_UM980]) {
		refuse(out, "R3", "all-F9P cluster: 0 wins in 12 scorable pairs across 13-28 km");
		out->nearest_km = nearest;
		goto done;
	}

	/* R4 -- the nearest five within reach span more than one family */
	if (has_mix) {
		int seen[VBS_NFAM] = {0, 0, 0, 0}, nfams = 0, taken = 0;

		for (i = 0; i < nall && taken < MIN_AUXILIARIES + 2; i++) {
			if (all[i].d > CLUSTER_RADIUS_KM)
				continue;
			if (!seen[bases[all[i].idx].family]) {
				seen[bases[all[i].idx].family] = 1;
				nfams++;
			}
			taken++;
		}
		if (nfams > 1) {
			refuse(out, "R4", "mixed cluster: refused on three independent tests, and mixing "
			    "forfeits network GLONASS");
			out->nearest_km = nearest;
			goto done;
		}
	}

	/* R5 */
	snprintf(buf, sizeof buf, "nearest base %.1f km is outside every family's envelope "
	    "(UM980 %.0f-%.0f, Mosaic %.0f-%.0f)", nearest,
	    ENVELOPE[VBS_FAM_UM980].lo, ENVELOPE[VBS_FAM_UM980].hi,
	    ENVELOPE[VBS_FAM_MOSAIC].lo, ENVELOPE[VBS_FAM_MOSAIC].hi);
	refuse(out, "R5", buf);
	out->nearest_km = nearest;
done:
	free(all);
	for (k = 0; k < VBS_NFAM; k++)
		free(fam[k]);
}

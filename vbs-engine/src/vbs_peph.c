#include <math.h>
#include <string.h>

#include "vbs/vbs_peph.h"

static int peph_sysidx(int sat) {
	switch (satsys(sat, NULL)) {
	case SYS_GPS: return 0;
	case SYS_GLO: return 1;
	case SYS_GAL: return 2;
	case SYS_CMP: return 3;
	}
	return -1;
}

int vbs_peph_load(nav_t *nav, const char *const *sp3, int nsp3, const char *atx) {
	int i;

	if (nsp3 <= 0)
		return 0;
	if (atx == NULL) {
		fprintf(stderr, "SP3 orbits need the satellite antennas (--sat-atx): SP3 is the "
		    "centre of mass, broadcast the phase centre\n");
		return -1;
	}
	for (i = 0; i < nsp3; i++) {
		int before = nav->ne;

		readsp3(sp3[i], nav, 0);
		if (nav->ne <= before) {
			fprintf(stderr, "no SP3 epochs read from %s\n", sp3[i]);
			return -1;
		}
	}
	if (!readsap(atx, nav->peph[0].time, nav)) {
		fprintf(stderr, "cannot read satellite antennas from %s\n", atx);
		return -1;
	}
	return nav->ne;
}

void vbs_peph_substitute(const nav_t *nav, const obsd_t *obs, int n, double *rs,
	const double *a, const double *b, struct vbs_peph_stats *st) {
	static double rp[MAXOBS * 6], dtp[MAXOBS * 2], varp[MAXOBS];
	static int svhp[MAXOBS];
	int i, j, s;

	if (n > MAXOBS)
		n = MAXOBS;
	satposs(obs[0].time, obs, n, nav, EPHOPT_PREC, rp, dtp, varp, svhp);
	for (i = 0; i < n; i++) {
		double *r = rs + i * 6, *p = rp + i * 6, d[3], ea[3], eb[3], u[3], proj = 0.0;

		if (norm(r, 3) <= 0.0)
			continue;	/* no broadcast position: the broadcast run had none either */
		s = peph_sysidx(obs[i].sat);
		if (norm(p, 3) <= 0.0) {
			if (st != NULL && s >= 0)
				st->kept[s]++;
			continue;
		}
		for (j = 0; j < 3; j++)
			d[j] = p[j] - r[j];
		if (st != NULL && s >= 0 && geodist(r, a, ea) > 0.0 && geodist(r, b, eb) > 0.0 &&
		    normv3(r, u)) {
			for (j = 0; j < 3; j++)
				proj += d[j] * (ea[j] - eb[j]);
			st->replaced[s]++;
			st->sum_d2[s] += dot(d, d, 3);
			st->sum_radial[s] += dot(d, u, 3);
			st->sum_proj2[s] += proj * proj;
			if (fabs(proj) > st->max_proj[s])
				st->max_proj[s] = fabs(proj);
		}
		memcpy(r, p, 6 * sizeof *r);
	}
}

void vbs_peph_report(FILE *fp, const struct vbs_peph_stats *st, const char *geometry) {
	static const char *name[VBS_PEPH_NSYS] = {"GPS", "GLO", "GAL", "BDS"};
	int s;

	for (s = 0; s < VBS_PEPH_NSYS; s++) {
		long n = st->replaced[s];

		if (n == 0 && st->kept[s] == 0)
			continue;
		fprintf(fp, "# sp3 %s: %ld positions replaced, %ld kept broadcast (not in the SP3); "
		    "SP3 - broadcast %.3f m RMS, radial mean %+.3f m; on the %s geometry "
		    "%.2f mm RMS, %.2f mm max\n", name[s], n, st->kept[s],
		    n ? sqrt(st->sum_d2[s] / n) : 0.0, n ? st->sum_radial[s] / n : 0.0, geometry,
		    n ? 1e3 * sqrt(st->sum_proj2[s] / n) : 0.0, 1e3 * st->max_proj[s]);
	}
}

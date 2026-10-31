/*
 * vbs_loo -- leave-one-out residual measurement (spec S8 step 1, data/README.md S5 step 1).
 *
 * Takes a master base and a target base, synthesizes a VBS at the *target's own coordinates*
 * from the master's observations, and double-differences the result against the target's real
 * observations. Because the VBS sits exactly on the target, the geometric range cancels and
 * what is left is the error the network corrections have to supply, plus an integer ambiguity
 * per satellite pair (removed as an arc mean) and noise.
 *
 * With corrections forced to zero -- which is all this tool can do today, since nothing
 * computes them yet -- the numbers are the *uncorrected* error budget as a function of
 * baseline length. That is the thing the network EKF has to remove, and it is worth measuring
 * before designing the filter that removes it.
 *
 * Deliberately takes coordinates and file paths as arguments rather than reading
 * MANIFEST.json: that file is validation-dataset metadata, while the engine's real `basemeta`
 * module reads the network catalogue and the caster's API (D6/D18).
 * Teaching the engine to parse a test fixture would be building the wrong module. A Python
 * driver iterates the cases and aggregates.
 *
 * Usage:
 *   vbs_loo <master.rnx> <mx> <my> <mz> <target.rnx> <tx> <ty> <tz> <nav.rnx> [label]
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_ddres.h"
#include "vbs/vbs_peph.h"
#include "vbs/vbs_signal.h"
#include "vbs/vbs_synth.h"

#define ELEV_MASK_DEG	15.0

static const char *sys_name(int sys) {
	switch (sys) {
	case SYS_GPS: return "GPS";
	case SYS_GAL: return "GAL";
	case SYS_CMP: return "BDS";
	case SYS_GLO: return "GLO";
	}
	return "???";
}

static const char *comb_name(const struct vbs_dd_stat *s) {
	static char buf[16];

	switch (s->comb) {
	case VBS_DD_GEOMFREE: return "geom-free";
	case VBS_DD_IONOFREE: return "iono-free";
	case VBS_DD_BAND:
	default:
		snprintf(buf, sizeof buf, "L%d", s->band);
		return buf;
	}
}

/* Records belonging to the epoch starting at index *i; advances *i past them. */
static int epoch_slice(const obs_t *obs, int *i, const obsd_t **out) {
	int start = *i, n = 0;

	if (start >= obs->n)
		return 0;
	*out = &obs->data[start];
	while (start + n < obs->n &&
	    fabs(timediff(obs->data[start + n].time, obs->data[start].time)) <= 1e-3)
		n++;
	*i = start + n;
	return n;
}

int main(int argc, char **argv) {
	static obs_t mobs, tobs;
	static nav_t nav;
	static obsd_t vbs[MAXOBS];
	static double rs[6 * MAXOBS], dts[2 * MAXOBS], var[MAXOBS];
	static int svh[MAXOBS];
	struct vbs_ddres *dd;
	struct vbs_dd_result res;
	gtime_t t0 = {0, 0.0};
	double mxyz[3], txyz[3], baseline;
	const char *label;
	int mi = 0, ti = 0, nvbs, g, nmatched = 0, i;
	/* --sp3 / --sat-atx: the same measurement on precise orbits (D163, vbs_peph.h) */
	const char *sp3path[VBS_PEPH_MAXFILE], *sat_atx = NULL;
	int nsp3 = 0;
	static struct vbs_peph_stats sp3st;

	if (argc < 10) {
		fprintf(stderr, "usage: %s <master.rnx> <mx> <my> <mz> <target.rnx> "
		    "<tx> <ty> <tz> <nav.rnx> [label] [--sp3 <orbits.sp3> ... --sat-atx <igs20.atx>]\n",
		    argv[0]);
		return 2;
	}
	mxyz[0] = atof(argv[2]); mxyz[1] = atof(argv[3]); mxyz[2] = atof(argv[4]);
	txyz[0] = atof(argv[6]); txyz[1] = atof(argv[7]); txyz[2] = atof(argv[8]);
	label = "";
	i = 10;
	if (argc > 10 && strncmp(argv[10], "--", 2) != 0)
		label = argv[i++];
	for (; i < argc; i++) {
		if (!strcmp(argv[i], "--sp3") && i + 1 < argc && nsp3 < VBS_PEPH_MAXFILE)
			sp3path[nsp3++] = argv[++i];
		else if (!strcmp(argv[i], "--sat-atx") && i + 1 < argc)
			sat_atx = argv[++i];
		else {
			fprintf(stderr, "unknown or malformed argument: %s\n", argv[i]);
			return 2;
		}
	}

	baseline = sqrt((mxyz[0] - txyz[0]) * (mxyz[0] - txyz[0]) +
	                (mxyz[1] - txyz[1]) * (mxyz[1] - txyz[1]) +
	                (mxyz[2] - txyz[2]) * (mxyz[2] - txyz[2]));

	memset(&mobs, 0, sizeof mobs);
	memset(&tobs, 0, sizeof tobs);
	memset(&nav, 0, sizeof nav);

	if (readrnxt(argv[1], 1, t0, t0, 0.0, "", &mobs, &nav, NULL) < 0 || mobs.n <= 0) {
		fprintf(stderr, "cannot read master %s\n", argv[1]);
		return 1;
	}
	if (readrnxt(argv[5], 2, t0, t0, 0.0, "", &tobs, &nav, NULL) < 0 || tobs.n <= 0) {
		fprintf(stderr, "cannot read target %s\n", argv[5]);
		return 1;
	}
	if (readrnxt(argv[9], 1, t0, t0, 0.0, "", NULL, &nav, NULL) < 0) {
		fprintf(stderr, "cannot read nav %s\n", argv[9]);
		return 1;
	}
	if (vbs_peph_load(&nav, sp3path, nsp3, sat_atx) < 0)
		return 1;
	if ((dd = vbs_ddres_new(ELEV_MASK_DEG)) == NULL) {
		fprintf(stderr, "out of memory\n");
		return 1;
	}

	/*
	 * Walk both files in step. The two stations do not always have the same epochs (the
	 * dataset has real outages, data/README.md S7b), so epochs are matched by time and
	 * anything unmatched is simply skipped.
	 */
	while (mi < mobs.n && ti < tobs.n) {
		const obsd_t *mep, *tep;
		int nm, nt;
		double dt;

		/* Peek both epoch times; only advance the one that is behind. */
		dt = timediff(mobs.data[mi].time, tobs.data[ti].time);
		if (dt < -1e-3) {
			epoch_slice(&mobs, &mi, &mep);
			continue;
		}
		if (dt > 1e-3) {
			epoch_slice(&tobs, &ti, &tep);
			continue;
		}
		nm = epoch_slice(&mobs, &mi, &mep);
		nt = epoch_slice(&tobs, &ti, &tep);
		if (nm == 0 || nt == 0)
			break;
		if (nm > MAXOBS)
			nm = MAXOBS;
		nmatched++;

		/* Satellite positions from the master's own observations. */
		satposs(mep[0].time, mep, nm, &nav, EPHOPT_BRDC, rs, dts, var, svh);
		if (nsp3 > 0)
			vbs_peph_substitute(&nav, mep, nm, rs, mxyz, txyz, &sp3st);

		nvbs = vbs_synth_epoch(mep, nm, rs, mxyz, txyz, NULL, 0, &nav, vbs);
		if (nvbs <= 0)
			continue;

		/*
		 * The synthesized set is compacted relative to the master's (GLONASS dropped),
		 * so the satellite positions no longer line up by index. Rebuild them for the
		 * synthesized set, which is what the DD code needs elevations for.
		 */
		satposs(vbs[0].time, vbs, nvbs, &nav, EPHOPT_BRDC, rs, dts, var, svh);

		vbs_ddres_add_epoch(dd, vbs, nvbs, tep, nt, rs, txyz, &nav);
	}

	if (vbs_ddres_finish(dd, &res) < 0) {
		fprintf(stderr, "cannot reduce residuals\n");
		vbs_ddres_free(dd);
		return 1;
	}

	printf("# leave-one-out residual, corrections forced to zero\n");
	printf("# label=%s baseline_km=%.2f epochs_matched=%d epochs_used=%d elev_mask=%.0f\n",
	    label, baseline / 1000.0, nmatched, res.nepoch, ELEV_MASK_DEG);
	printf("# %-4s %-10s %8s %6s %10s %10s %10s\n",
	    "sys", "comb", "samples", "arcs", "rms_m", "p95_m", "max_m");
	for (g = 0; g < res.nstat; g++) {
		const struct vbs_dd_stat *s = &res.stat[g];

		printf("%-6s %-10s %8d %6d %10.4f %10.4f %10.4f\n",
		    sys_name(s->sys), comb_name(s), s->nsample, s->narc, s->rms, s->p95, s->max);
	}
	if (nsp3 > 0)
		vbs_peph_report(stdout, &sp3st, "master-target");
	if (getenv("VBS_LOO_OUTLIERS") != NULL)
		vbs_ddres_dump_outliers(dd, atof(getenv("VBS_LOO_OUTLIERS")), 25, stdout);

	vbs_ddres_free(dd);
	free(mobs.data);
	free(tobs.data);
	freenav(&nav, 0xFF);
	return 0;
}

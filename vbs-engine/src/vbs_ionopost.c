#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_ionopost.h"
#include "vbs/vbs_signal.h"

/* See the header for what each option is and why it exists. */

#define RE_KM		6371.0
#define SHELL_KM	350.0
#define NDATUM		5		/* G R E C J, the datum columns */
#define MAXPAR		(NDATUM + 6)
#define IDW_D0_KM	50.0		/* inverse-distance softening */
#define IDW_MAX_KM	400.0		/* residuals farther than this do not contribute */
#define IDW_PRIOR	(1.0 / (200.0 * 200.0))	/* shrinks a lone far residual toward the fit */

int vbs_ionopost_parse_sys(const char *s) {
	int m = 0;

	for (; *s; s++)
		switch (*s) {
		case 'G': m |= SYS_GPS; break;
		case 'R': m |= SYS_GLO; break;
		case 'E': m |= SYS_GAL; break;
		case 'C': m |= SYS_CMP; break;
		case 'J': m |= SYS_QZS; break;
		default: return -1;
		}
	return m;
}

int vbs_ionopost_active(const struct vbs_ionopost_opt *opt) {
	return opt->need_band5_sys || (opt->scale_gps > 0.0 && opt->scale_gps != 1.0) ||
	    opt->model_sys || opt->drop_sys || opt->shrink_m > 0.0;
}

static int sysi(int sys) {
	switch (sys) {
	case SYS_GPS: return 0;
	case SYS_GLO: return 1;
	case SYS_GAL: return 2;
	case SYS_CMP: return 3;
	case SYS_QZS: return 4;
	}
	return -1;
}

static int master_has_band5(int sat, const obsd_t *mep, int nm) {
	int k, j;

	for (k = 0; k < nm; k++) {
		if (mep[k].sat != sat)
			continue;
		for (j = 0; j < NFREQ + NEXOBS; j++)
			if (mep[k].L[j] != 0.0 && vbs_signal_band(mep[k].code[j]) == 5)
				return 1;
		return 0;
	}
	return 0;
}

/* Pierce point offset (north, east, km) from the VBS and the shell mapping function. */
static void pierce(double az, double el, double *pn, double *pe, double *map) {
	double z = PI / 2.0 - el, sz = RE_KM / (RE_KM + SHELL_KM) * sin(z), zp = asin(sz);
	double psi = z - zp;

	*pn = psi * RE_KM * cos(az);
	*pe = psi * RE_KM * sin(az);
	*map = 1.0 / cos(zp);
}

/* dI on its constellation's datum frequency <-> the same slant delay on L1 */
static double to_l1(int sys, double dI) {
	double f = vbs_synth_ref_freq(sys);

	return f > 0.0 ? dI * (f / FREQ1) * (f / FREQ1) : 0.0;
}

static double from_l1(int sys, double y) {
	double f = vbs_synth_ref_freq(sys);

	return f > 0.0 ? y * (FREQ1 / f) * (FREQ1 / f) : 0.0;
}

/* Design row of the shell model (without the datum columns) for one pierce point. */
static int shell_row(int mode, double pn, double pe, double map, double *h) {
	double n = pn / 100.0, e = pe / 100.0;	/* 100 km units keep the normal matrix sane */
	int k = 0;

	h[k++] = map;
	h[k++] = map * n;
	h[k++] = map * e;
	if (mode == VBS_IONOPOST_QUAD) {
		h[k++] = map * n * n;
		h[k++] = map * e * e;
		h[k++] = map * n * e;
	}
	return k;
}

static int cmp_dbl(const void *a, const void *b) {
	double x = *(const double *)a, y = *(const double *)b;

	return x < y ? -1 : x > y;
}

/* An own, kept correction of constellation index si that the shrink applies to. */
static int shrinkable(const struct vbs_corr *c, int withheld, const struct vbs_ionopost_opt *opt,
	int si) {
	int sys = satsys(c->sat, NULL);

	return c->valid && (withheld == 0 || withheld == 3) && sysi(sys) == si &&
	    (!opt->sys_mask || (sys & opt->sys_mask));
}

int vbs_ionopost_apply(struct vbs_corr *corr, int n, const double *azel, const double *vpos,
	const obsd_t *mep, int nm, const nav_t *nav, const struct vbs_ionopost_opt *opt,
	struct vbs_ionopost_stats *stats) {
	static int withheld[MAXOBS];
	static double pn[MAXOBS], pe[MAXOBS], map[MAXOBS], res[MAXOBS];
	int k, nchg = 0;

	(void)vpos;
	(void)nav;
	if (n > MAXOBS)
		n = MAXOBS;
	for (k = 0; k < n; k++) {
		int sys = satsys(corr[k].sat, NULL);

		withheld[k] = 0;
		if (!corr[k].valid)
			continue;
		if (opt->sys_mask && !(sys & opt->sys_mask)) {
			if (corr[k].dI != 0.0)
				nchg++;
			corr[k].dI = 0.0;
			withheld[k] = 1;
		} else if ((sys & opt->need_band5_sys) && !master_has_band5(corr[k].sat, mep, nm)) {
			corr[k].dI = 0.0;
			withheld[k] = 1;
			nchg++;
			if (stats)
				stats->band5_withheld++;
		} else if (sys == SYS_GPS && opt->scale_gps > 0.0 && opt->scale_gps != 1.0) {
			corr[k].dI *= opt->scale_gps;
			nchg++;
		}
	}

	if (opt->model_sys && azel != NULL) {
		static double A[MAXPAR * MAXOBS], y[MAXOBS], x[MAXPAR], Q[MAXPAR * MAXPAR];
		static int fit_idx[MAXOBS];
		int col_of[NDATUM], ncol = 0, nfit = 0, npar, i, nshell;
		double h[6], ss = 0.0;

		for (i = 0; i < NDATUM; i++)
			col_of[i] = -1;
		/* the datum columns: one per constellation in the fit set */
		for (k = 0; k < n; k++) {
			int sys = satsys(corr[k].sat, NULL), si = sysi(sys);
			double el;

			if (!corr[k].valid || withheld[k] || (sys & opt->model_sys) || si < 0)
				continue;
			if ((el = azel[(corr[k].sat - 1) * 2 + 1]) <= 0.0)
				continue;
			if (col_of[si] < 0)
				col_of[si] = ncol++;
			fit_idx[nfit++] = k;
		}
		nshell = opt->model_mode == VBS_IONOPOST_QUAD ? 6 : 3;
		npar = ncol + nshell;
		if (nfit < npar + 2) {
			if (stats)
				stats->fit_refused++;
			/* no model: withheld model-system corrections stay withheld */
		} else {
			for (i = 0; i < nfit; i++) {
				struct vbs_corr *c = &corr[fit_idx[i]];
				int sys = satsys(c->sat, NULL), j;
				double w;

				pierce(azel[(c->sat - 1) * 2], azel[(c->sat - 1) * 2 + 1], &pn[i], &pe[i],
				    &map[i]);
				w = 1.0 / map[i];	/* sigma grows with the mapping function */
				for (j = 0; j < npar; j++)
					A[j + i * npar] = 0.0;
				A[col_of[sysi(sys)] + i * npar] = w;
				shell_row(opt->model_mode, pn[i], pe[i], map[i], h);
				for (j = 0; j < nshell; j++)
					A[ncol + j + i * npar] = w * h[j];
				y[i] = w * to_l1(sys, c->dI);
			}
			if (lsq(A, y, npar, nfit, x, Q) != 0) {
				if (stats)
					stats->fit_refused++;
			} else {
				static double off[MAXOBS];
				int noff[NDATUM], si;
				double med[NDATUM];

				for (i = 0; i < nfit; i++) {
					struct vbs_corr *c = &corr[fit_idx[i]];
					double f = x[col_of[sysi(satsys(c->sat, NULL))]];
					int j;

					shell_row(opt->model_mode, pn[i], pe[i], map[i], h);
					for (j = 0; j < nshell; j++)
						f += x[ncol + j] * h[j];
					res[i] = to_l1(satsys(c->sat, NULL), c->dI) - f;
					ss += res[i] * res[i];
				}
				if (stats) {
					stats->fits++;
					stats->fit_rms_sum += sqrt(ss / nfit);
				}
				/* the prediction, and each kept own correction's offset to it */
				for (si = 0; si < NDATUM; si++)
					noff[si] = 0;
				for (k = 0; k < n; k++) {
					int sys = satsys(corr[k].sat, NULL), j;
					double tn, te, tm, p = 0.0;

					off[k] = 0.0;
					if (!corr[k].valid || !(sys & opt->model_sys) || sysi(sys) < 0 ||
					    azel[(corr[k].sat - 1) * 2 + 1] <= 0.0)
						continue;
					pierce(azel[(corr[k].sat - 1) * 2], azel[(corr[k].sat - 1) * 2 + 1],
					    &tn, &te, &tm);
					shell_row(opt->model_mode, tn, te, tm, h);
					for (j = 0; j < nshell; j++)
						p += x[ncol + j] * h[j];
					if (opt->model_mode == VBS_IONOPOST_IDW) {
						double sw = IDW_PRIOR, swr = 0.0;

						for (i = 0; i < nfit; i++) {
							double d2 = (tn - pn[i]) * (tn - pn[i]) +
							    (te - pe[i]) * (te - pe[i]), w;

							if (d2 > IDW_MAX_KM * IDW_MAX_KM)
								continue;
							w = 1.0 / (d2 + IDW_D0_KM * IDW_D0_KM);
							sw += w;
							swr += w * res[i] / map[i];
						}
						p += tm * swr / sw;
					}
					p = from_l1(sys, p);
					if (withheld[k]) {
						corr[k].dI = p;
						nchg++;
						if (stats)
							stats->filled++;
						withheld[k] = 2;	/* filled */
					} else {
						off[k] = corr[k].dI - p;
						si = sysi(sys);
						noff[si]++;
						withheld[k] = 3;	/* own kept, prediction in off */
					}
				}
				if (opt->model_check_m > 0.0)
					for (si = 0; si < NDATUM; si++) {
						static double tmp[MAXOBS];
						int m = 0;

						if (noff[si] == 0)
							continue;
						for (k = 0; k < n; k++)
							if (withheld[k] == 3 && sysi(satsys(corr[k].sat, NULL)) == si)
								tmp[m++] = off[k];
						qsort(tmp, m, sizeof *tmp, cmp_dbl);
						med[si] = m % 2 ? tmp[m / 2] : 0.5 * (tmp[m / 2 - 1] + tmp[m / 2]);
						for (k = 0; k < n; k++) {
							if (withheld[k] != 3 || sysi(satsys(corr[k].sat, NULL)) != si)
								continue;
							if (fabs(off[k] - med[si]) > opt->model_check_m) {
								/* the prediction, on this constellation's own datum */
								corr[k].dI = corr[k].dI - off[k] + med[si];
								nchg++;
								if (stats)
									stats->replaced++;
							} else if (stats)
								stats->kept++;
						}
					}
			}
		}
	}

	if (opt->shrink_m > 0.0) {
		static double tmp[MAXOBS];
		int si;

		/*
		 * On each satellite's departure from its constellation's median: dI carries the
		 * constellation datum (receiver code biases, metres), which the rover's double
		 * difference cancels and which would otherwise swamp the weight (measured: mean weight
		 * 1.000 on day 251 before this). Own corrections only: a filled one is a model already.
		 */
		for (si = 0; si < NDATUM; si++) {
			int m = 0;
			double med;

			for (k = 0; k < n; k++)
				if (shrinkable(&corr[k], withheld[k], opt, si))
					tmp[m++] = to_l1(satsys(corr[k].sat, NULL), corr[k].dI);
			if (m < 2)
				continue;
			qsort(tmp, m, sizeof *tmp, cmp_dbl);
			med = m % 2 ? tmp[m / 2] : 0.5 * (tmp[m / 2 - 1] + tmp[m / 2]);
			for (k = 0; k < n; k++) {
				int sys = satsys(corr[k].sat, NULL);
				double d, w;

				if (!shrinkable(&corr[k], withheld[k], opt, si))
					continue;
				d = to_l1(sys, corr[k].dI) - med;
				w = d * d / (d * d + opt->shrink_m * opt->shrink_m);
				corr[k].dI = from_l1(sys, med + d * w);
				nchg++;
				if (stats) {
					stats->shrunk++;
					stats->shrink_w_sum += w;
				}
			}
		}
	}

	if (opt->drop_sys)
		for (k = 0; k < n; k++)
			if (corr[k].valid && (satsys(corr[k].sat, NULL) & opt->drop_sys)) {
				corr[k].valid = 0;
				nchg++;
			}
	return nchg;
}

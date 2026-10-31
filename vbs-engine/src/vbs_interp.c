#include <math.h>
#include <string.h>

#include "vbs/vbs_interp.h"

/*
 * Weighted plane fit at the VBS point -- spec §4.5. See the header for why each piece is
 * there.
 */

void vbs_interp_default_opt(struct vbs_interp_opt *opt) {
	memset(opt, 0, sizeof *opt);
	/*
	 * 0.05 m: below this the correction is worth using, above it the state has not
	 * converged enough to help. D25 measured the double-differenced ionosphere at
	 * 5.5-9.4 mm once ambiguities are fixed and 42-76 mm float, so this threshold is
	 * roughly "fixed, or nearly so".
	 */
	opt->max_sigma = 0.05;
	/*
	 * The per-satellite fit floor, deliberately separate from the cluster module's admission
	 * floor (vbs_cluster_opt.min_aux, 3). Two auxiliaries plus the master anchor still give
	 * three points for the plane's three parameters, so a satellite one auxiliary has lost
	 * keeps its correction instead of leaving the VBS and re-entering a few epochs later --
	 * which costs the rover that satellite's ambiguity each time. Decided
	 * 2026-10-02, together with the "2 + 1 backup" admission.
	 */
	opt->min_aux = 2;
	opt->humi = 0.7;
	opt->iono_scale = 1.0;
	/*
	 * On by default (D69). The master's own correction is zero by definition, and spending
	 * the fit's redundancy on the one sample that carries no error is free accuracy: measured,
	 * it took M30's false fixes at 5 km offset from 1.14 % to 0.73 % -- through the 1 % target
	 * -- and at 46-75 km it removed a 93 cm wrong-integer offset outright, lifting that case
	 * from 2.1 % fixed to 24.7 %.
	 */
	opt->anchor_master = 1;
	/*
	 * Collocation by default (D70). Measured against the plane fit at M30's 10 km tile
	 * corner: 90.5 % fixed with 5.89 % false becomes 98.0 % with 0.60 % -- the plane's
	 * curvature error was most of the horizontal budget, and the 10 km tile
	 * only works with this. At the synthesis point it is also better (H95 4.8 -> 3.9 cm).
	 */
	opt->method = VBS_INTERP_METHOD_COLLOCATION;
	/*
	 * 60 km and 5 cm: the correlation length is a little longer than the clusters being
	 * served, so nearby auxiliaries dominate without distant ones being discarded, and the
	 * field sigma is the scale of the corrections themselves. Both are starting values to be
	 * calibrated, not asserted.
	 */
	opt->corr_len_km = 60.0;
	opt->sig_field_m = 0.05;
}

/*
 * Fit d = a + b·E + c·N weighted by w, evaluate at E = N = 0 (the VBS), and report the fitted
 * value together with `var`, the variance of that value.
 *
 * Note what `var` is and is not. Because the weights are 1/σ², A⁻¹[0,0] comes out in m² -- it
 * is the propagated variance of the interpolated correction, which is what §4.4/§4.5 want
 * carried forward. It is *not* the dimensionless extrapolation metric: that has to be
 * computed from geometry alone (see geom_leverage below), or a cluster of precise stations
 * would look like a safe extrapolation purely because its corrections were tight.
 *
 * Returns 1 on success, 0 if the geometry is degenerate. On degeneracy it falls back to the
 * weighted mean, which is the right answer when the auxiliaries are collinear: a plane is not
 * determined, but a level is.
 */
static int wplane(const double *e, const double *n, const double *v, const double *w, int m,
	double *value, double *var) {
	double A[9] = {0}, b[3] = {0}, x[3], sw = 0.0, sv = 0.0;
	int i, r, c;

	for (i = 0; i < m; i++) {
		double g[3] = {1.0, e[i], n[i]};

		for (r = 0; r < 3; r++) {
			for (c = 0; c < 3; c++)
				A[r + c * 3] += w[i] * g[r] * g[c];
			b[r] += w[i] * g[r] * v[i];
		}
		sw += w[i];
		sv += w[i] * v[i];
	}
	if (m >= 3 && matinv(A, 3) == 0 && A[0] > 0.0) {
		matmul("NN", 3, 1, 3, 1.0, A, b, 0.0, x);
		*value = x[0];
		*var = A[0];
		return 1;
	}
	if (sw <= 0.0)
		return 0;
	*value = sv / sw;
	*var = 1.0 / sw;
	return 1;
}

/*
 * The consistency gate's counters and its leave-one-out test (vbs_interp_opt.iono_gate_m).
 */
static long gate_checked, gate_dropped, gate_withheld;

void vbs_interp_gate_stats(long *checked, long *dropped, long *withheld) {
	if (checked != NULL)
		*checked = gate_checked;
	if (dropped != NULL)
		*dropped = gate_dropped;
	if (withheld != NULL)
		*withheld = gate_withheld;
}

void vbs_interp_gate_reset(void) {
	gate_checked = gate_dropped = gate_withheld = 0;
}

/*
 * Misfit of sample `skip`: fit a weighted plane through every other sample plus the master
 * anchor (value 0 at E = N = 0, weight the mean of the others, as the fit itself anchors), and
 * return v[skip] minus the plane's value there. Positions are master-relative, in km. Returns 0
 * when the remaining points do not determine a plane.
 */
static int loo_misfit(const double *e, const double *n, const double *v, const double *w,
	int m, int skip, double *misfit) {
	double A[9] = {0}, b[3] = {0}, x[3], wa = 0.0;
	int i, r, c, cnt = 0;

	for (i = 0; i < m; i++)
		if (i != skip) {
			wa += w[i];
			cnt++;
		}
	if (cnt < 2)
		return 0;
	wa /= cnt;
	for (i = 0; i <= m; i++) {
		double g[3], wi, vi;

		if (i == skip)
			continue;
		g[0] = 1.0;
		g[1] = i < m ? e[i] : 0.0;
		g[2] = i < m ? n[i] : 0.0;
		wi = i < m ? w[i] : wa;
		vi = i < m ? v[i] : 0.0;
		for (r = 0; r < 3; r++) {
			for (c = 0; c < 3; c++)
				A[r + c * 3] += wi * g[r] * g[c];
			b[r] += wi * g[r] * vi;
		}
	}
	if (matinv(A, 3) != 0)
		return 0;
	matmul("NN", 3, 1, 3, 1.0, A, b, 0.0, x);
	*misfit = v[skip] - (x[0] + x[1] * e[skip] + x[2] * n[skip]);
	return 1;
}

/*
 * Least-squares collocation at the VBS (D70).
 *
 * `e`/`n` are sample positions in km relative to the VBS (so the prediction point is the
 * origin), `v` the sampled corrections, `sig` each sample's own standard deviation. The field
 * covariance is sig_f^2 * exp(-d / L): isotropic, exponential, zero mean -- the zero mean being
 * correct here because the corrections are differences against the master and the master itself
 * is one of the samples, with value zero.
 *
 * Returns 1 and writes the predicted value and its error variance, or 0 if the system is
 * singular. The caller falls back to the plane fit in that case rather than inventing a value.
 */
static int collocate(const double *e, const double *n, const double *v, const double *sig,
	int m, double sig_f, double corr_len, double *value, double *var) {
	double C[VBS_INTERP_MAXAUX + 1][VBS_INTERP_MAXAUX + 1], c0[VBS_INTERP_MAXAUX + 1];
	double A[(VBS_INTERP_MAXAUX + 1) * (VBS_INTERP_MAXAUX + 1)], x[VBS_INTERP_MAXAUX + 1];
	double sf2 = sig_f * sig_f;
	int i, j;

	if (m < 1 || m > VBS_INTERP_MAXAUX + 1 || corr_len <= 0.0 || sig_f <= 0.0)
		return 0;
	for (i = 0; i < m; i++) {
		double d0 = sqrt(e[i] * e[i] + n[i] * n[i]);

		c0[i] = sf2 * exp(-d0 / corr_len);
		for (j = 0; j < m; j++) {
			double de = e[i] - e[j], dn = n[i] - n[j];

			C[i][j] = sf2 * exp(-sqrt(de * de + dn * dn) / corr_len);
			if (i == j)
				C[i][j] += sig[i] * sig[i];	/* the sample's own noise */
		}
	}
	for (i = 0; i < m; i++)
		for (j = 0; j < m; j++)
			A[i + j * m] = C[i][j];
	if (matinv(A, m) != 0)
		return 0;
	/* x = (C+N)^-1 v, then value = c0' x */
	for (i = 0; i < m; i++) {
		x[i] = 0.0;
		for (j = 0; j < m; j++)
			x[i] += A[i + j * m] * v[j];
	}
	*value = 0.0;
	for (i = 0; i < m; i++)
		*value += c0[i] * x[i];
	/* Error variance: sig_f^2 - c0' (C+N)^-1 c0. */
	*var = sf2;
	for (i = 0; i < m; i++) {
		double t = 0.0;

		for (j = 0; j < m; j++)
			t += A[i + j * m] * c0[j];
		*var -= c0[i] * t;
	}
	if (*var < 0.0)
		*var = 0.0;
	return 1;
}

/*
 * The extrapolation metric §4.5 asks to be published: the variance inflation of the fitted
 * value at the VBS from *geometry alone*, i.e. A⁻¹[0,0] of the same plane fit with unit
 * weights. It is ~1/m inside a well-spread cluster and grows without bound as the fit is
 * asked to extrapolate or as the auxiliaries approach collinearity.
 *
 * Dimensionless on purpose. Deriving it from the weighted fit instead would make a cluster of
 * precise stations look like a safe extrapolation, and make a well-enclosed cluster of noisy
 * ones look dangerous -- both backwards, because what `vbs.max_extrapolation_km` has to be
 * calibrated against is where the VBS sits relative to the cluster, not how good the
 * corrections happen to be this epoch.
 *
 * `04-geophysical-corrections.md` measured it on the real LOO geometries: 0.25-0.35 for the
 * enclosed cases against 6-12 for the extrapolated ones.
 */
static double geom_leverage(const double *e, const double *n, int m) {
	double A[9] = {0};
	int i, r, c;

	for (i = 0; i < m; i++) {
		double g[3] = {1.0, e[i], n[i]};

		for (r = 0; r < 3; r++)
			for (c = 0; c < 3; c++)
				A[r + c * 3] += g[r] * g[c];
	}
	if (m < 3 || matinv(A, 3) != 0 || A[0] <= 0.0)
		return m > 0 ? 1.0 / m : 0.0;
	return A[0];
}

/*
 * Constellation index for the datum, 0..NDATUM-1. Grouping is per constellation because the
 * bias being removed is the receivers' differential code bias, which is a property of the
 * signal pair each constellation is tracked on, not of the individual satellite.
 */
#define NDATUM 8

static int datum_idx(int sat) {
	switch (satsys(sat, NULL)) {
	case SYS_GPS: return 0;
	case SYS_GLO: return 1;
	case SYS_GAL: return 2;
	case SYS_CMP: return 3;
	case SYS_QZS: return 4;
	case SYS_IRN: return 5;
	case SYS_SBS: return 6;
	default: return 7;
	}
}

/* Is this auxiliary's correction for `sat` one the fit would accept? */
static int usable(const struct vbs_interp_aux *a, int c, double max_sigma) {
	double sI, sT;

	if (c < 0 || !a->corr[c].valid)
		return 0;
	sI = a->sd_I != NULL ? a->sd_I[c] : 0.0;
	sT = a->sd_T != NULL ? a->sd_T[c] : 0.0;
	return !(max_sigma > 0.0 && (sI > max_sigma || sT > max_sigma));
}

/* Find satellite `sat` in an auxiliary's correction list, or -1. */
static int find_corr(const struct vbs_interp_aux *a, int sat) {
	int i;

	for (i = 0; i < a->ncorr; i++)
		if (a->corr[i].sat == sat)
			return i;
	return -1;
}

int vbs_interp_epoch(gtime_t time, const double *master_xyz, const double *vbs_xyz,
	const struct vbs_interp_aux *aux, int naux, const struct vbs_interp_opt *opt,
	const double *azel, struct vbs_corr *out, struct vbs_interp_quality *qual,
	int maxout) {
	struct vbs_interp_opt defopt;
	double mpos[3], vpos[3], rel[3], denu[3];
	double ve, vn;
	double level[VBS_INTERP_MAXAUX][NDATUM], lw[VBS_INTERP_MAXAUX][NDATUM];
	int nout = 0, k, i;

	if (opt == NULL) {
		vbs_interp_default_opt(&defopt);
		opt = &defopt;
	}
	if (naux > VBS_INTERP_MAXAUX)
		naux = VBS_INTERP_MAXAUX;
	if (naux < 1)
		return 0;

	ecef2pos(master_xyz, mpos);
	ecef2pos(vbs_xyz, vpos);
	for (i = 0; i < 3; i++)
		rel[i] = vbs_xyz[i] - master_xyz[i];
	ecef2enu(mpos, rel, denu);
	ve = denu[0] / 1000.0;			/* km: keeps the normal matrix conditioned */
	vn = denu[1] / 1000.0;

	/*
	 * Each auxiliary's ionosphere datum per constellation (D55). Estimated only from the
	 * satellites *every* participating auxiliary carries this epoch: a datum taken over
	 * each auxiliary's own satellite set would differ between auxiliaries simply because
	 * the sets differ, putting back a version of the bias it exists to remove.
	 *
	 * A weighted mean over that common set rather than one reference satellite, so losing
	 * a satellite moves the datum by (its residual)/n -- about a centimetre on eight
	 * satellites -- instead of stepping every correction of the constellation by the
	 * difference between two reference satellites.
	 */
	if (opt->datum == VBS_INTERP_DATUM_COMMON) {
		int d, s;

		memset(level, 0, sizeof level);
		memset(lw, 0, sizeof lw);
		for (s = 0; s < MAXSAT; s++) {
			int sat = s + 1, d2 = datum_idx(sat), all = 1;
			double w, v;

			/*
			 * The common set is what the PRIMARIES all hold: a backup that is absent
			 * most of the time must not shrink it. A backup still gets its own level, over
			 * those of the same satellites it does hold, so that on the rare epoch it is
			 * used its corrections sit on the same datum as the primaries'.
			 */
			for (k = 0; k < naux; k++)
				if (!aux[k].backup &&
				    !usable(&aux[k], find_corr(&aux[k], sat), opt->max_sigma)) {
					all = 0;
					break;
				}
			if (!all)
				continue;
			for (k = 0; k < naux; k++) {
				int c = find_corr(&aux[k], sat);
				double sI;

				if (aux[k].backup && !usable(&aux[k], c, opt->max_sigma))
					continue;
				sI = aux[k].sd_I != NULL ? aux[k].sd_I[c] : 0.0;

				w = sI > 0.0 ? 1.0 / (sI * sI) : 1.0;
				v = aux[k].corr[c].dI;
				level[k][d2] += w * v;
				lw[k][d2] += w;
			}
		}
		for (k = 0; k < naux; k++)
			for (d = 0; d < NDATUM; d++)
				level[k][d] = lw[k][d] > 0.0 ? level[k][d] / lw[k][d] : 0.0;
	}

	/* One plane fit per satellite, over whichever auxiliaries carry it. */
	for (k = 0; k < naux && nout < maxout; k++) {
		for (i = 0; i < aux[k].ncorr && nout < maxout; i++) {
			int sat = aux[k].corr[i].sat, j, m = 0, seen = 0;
			double e[VBS_INTERP_MAXAUX], n[VBS_INTERP_MAXAUX];
			double vI[VBS_INTERP_MAXAUX], vT[VBS_INTERP_MAXAUX];
			double wI[VBS_INTERP_MAXAUX], wT[VBS_INTERP_MAXAUX];
			double dI, dT, varI, varT, lev = 0.0;
			int used_backup = 0, withhold_I = 0;

			/* Emit each satellite once: skip if an earlier auxiliary owned it. */
			for (j = 0; j < k; j++)
				if (find_corr(&aux[j], sat) >= 0) {
					seen = 1;
					break;
				}
			if (seen)
				continue;

			/*
			 * Two passes: the primaries, then -- only if they fall below the fit floor for
			 * this satellite -- the reserve. Held back otherwise, because a distant station
			 * in a fit that does not need it only adds its noise: M30 widened from four
			 * auxiliaries to six, with the two extra at 45 and 51 km, went from 98.2 % fixed
			 * and 0.35 % false to 65.3 % and 16.7 % (D70).
			 */
			for (int pass = 0; pass < 2; pass++) {
				if (pass == 1 && m >= opt->min_aux)
					break;
				for (j = 0; j < naux; j++) {
					int c = find_corr(&aux[j], sat);
					double sI, sT, arel[3], aenu[3];

					if ((aux[j].backup != 0) != pass)
						continue;
					if (c < 0 || !aux[j].corr[c].valid)
						continue;
					sI = aux[j].sd_I != NULL ? aux[j].sd_I[c] : 0.0;
					sT = aux[j].sd_T != NULL ? aux[j].sd_T[c] : 0.0;
					if (opt->max_sigma > 0.0 &&
					    (sI > opt->max_sigma || sT > opt->max_sigma))
						continue;

					for (int q = 0; q < 3; q++)
						arel[q] = aux[j].xyz[q] - master_xyz[q];
					ecef2enu(mpos, arel, aenu);
					e[m] = aenu[0] / 1000.0;
					n[m] = aenu[1] / 1000.0;
					vI[m] = aux[j].corr[c].dI -
					    (opt->datum == VBS_INTERP_DATUM_COMMON ?
					    level[j][datum_idx(sat)] : 0.0);
					vT[m] = aux[j].corr[c].dT;
					/*
					 * 1/σ², so an auxiliary whose ambiguities never fixed fades out
					 * of the fit instead of having to be cut by a threshold. A
					 * caller that supplies no variances gets an unweighted fit.
					 */
					wI[m] = sI > 0.0 ? 1.0 / (sI * sI) : 1.0;
					wT[m] = sT > 0.0 ? 1.0 / (sT * sT) : 1.0;
					if (pass == 1)
						used_backup = 1;	/* it actually entered the fit */
					m++;
				}
			}
			if (m < opt->min_aux)
				continue;

			/* The consistency gate (vbs_interp_opt.iono_gate_m); off by default. */
			if (opt->iono_gate_m > 0.0 && m >= 3) {
				gate_checked++;
				for (;;) {
					double worst = 0.0, r;
					int wj = -1;

					for (j = 0; j < m; j++)
						if (loo_misfit(e, n, vI, wI, m, j, &r) && fabs(r) > worst) {
							worst = fabs(r);
							wj = j;
						}
					if (wj < 0 || worst <= opt->iono_gate_m)
						break;
					if (m <= 3) {
						withhold_I = 1;
						gate_withheld++;
						break;
					}
					for (j = wj; j < m - 1; j++) {
						e[j] = e[j + 1];
						n[j] = n[j + 1];
						vI[j] = vI[j + 1];
						vT[j] = vT[j + 1];
						wI[j] = wI[j + 1];
						wT[j] = wT[j + 1];
					}
					m--;
					gate_dropped++;
				}
			}

			/*
			 * The plane is fitted in the master-relative frame and evaluated at the
			 * VBS offset, not at the origin: the corrections are differences against
			 * the master, so the field is zero at the master by construction and the
			 * fit must be allowed to say so rather than be forced through it.
			 */
			{
				double eo[VBS_INTERP_MAXAUX + 1], no[VBS_INTERP_MAXAUX + 1];
				double vI2[VBS_INTERP_MAXAUX + 1], vT2[VBS_INTERP_MAXAUX + 1];
				double wI2[VBS_INTERP_MAXAUX + 1], wT2[VBS_INTERP_MAXAUX + 1];
				int mm = m;

				for (j = 0; j < m; j++) {
					eo[j] = e[j] - ve;
					no[j] = n[j] - vn;
					vI2[j] = vI[j];
					vT2[j] = vT[j];
					wI2[j] = wI[j];
					wT2[j] = wT[j];
				}
				/*
				 * Anchor the fit at the master (D69). These corrections are
				 * master-relative, so the field is *exactly* zero at the master --
				 * that is a definition, not an observation. With three auxiliaries
				 * and three plane parameters the free fit passes through the three
				 * of them exactly and throws that away, leaving the value at the VBS
				 * determined by three points and nothing else. Adding the master as a
				 * point of value zero spends the redundancy on the one sample that
				 * carries no error at all.
				 */
				if (opt->anchor_master && mm < VBS_INTERP_MAXAUX + 1) {
					double wa = 0.0;

					for (j = 0; j < m; j++)
						wa += wI[j];
					eo[mm] = -ve;
					no[mm] = -vn;
					vI2[mm] = 0.0;
					vT2[mm] = 0.0;
					wI2[mm] = m > 0 ? wa / m : 1.0;
					wa = 0.0;
					for (j = 0; j < m; j++)
						wa += wT[j];
					wT2[mm] = m > 0 ? wa / m : 1.0;
					mm++;
				}
				if (opt->method == VBS_INTERP_METHOD_COLLOCATION) {
					double sI2[VBS_INTERP_MAXAUX + 1];
					double sT2[VBS_INTERP_MAXAUX + 1];
					int q;

					for (q = 0; q < mm; q++) {
						sI2[q] = wI2[q] > 0.0 ? sqrt(1.0 / wI2[q]) : 1.0;
						sT2[q] = wT2[q] > 0.0 ? sqrt(1.0 / wT2[q]) : 1.0;
					}
					if (!collocate(eo, no, vI2, sI2, mm, opt->sig_field_m,
					    opt->corr_len_km, &dI, &varI) ||
					    !collocate(eo, no, vT2, sT2, mm, opt->sig_field_m,
					    opt->corr_len_km, &dT, &varT)) {
						if (!wplane(eo, no, vI2, wI2, mm, &dI, &varI))
							continue;
						if (!wplane(eo, no, vT2, wT2, mm, &dT, &varT))
							continue;
					}
				} else {
					if (!wplane(eo, no, vI2, wI2, mm, &dI, &varI))
						continue;
					if (!wplane(eo, no, vT2, wT2, mm, &dT, &varT))
						continue;
				}
				lev = geom_leverage(eo, no, mm);
			}

			/*
			 * The troposphere model term, computed rather than interpolated (§4.5),
			 * per satellite along its own line of sight. `dT` from the plane fit is
			 * the interpolated *residual* only -- vbs_netekf returns the residual
			 * alone for exactly this reason.
			 */
			{
				double dT_model = 0.0;

				if (azel != NULL && sat > 0 && sat <= MAXSAT) {
					const double *ae = azel + (sat - 1) * 2;

					/*
					 * tropmodel() returns 0.0 for an out-of-range input --
					 * height below -100 m or above 10 km, or a non-positive
					 * elevation -- and 0.0 is a *value* here, not a flag. A
					 * difference taken against it silently yields minus the
					 * other station's entire slant delay: measured at -4.79 m
					 * for a VBS placed below sea level, which would be applied
					 * to every observation of that satellite as if it were a
					 * troposphere correction.
					 *
					 * Same class as D21's `phase_m()` returning 0.0 for "no
					 * observation", and the rule from that case stands: treat
					 * any in-band sentinel in this engine as a defect, and
					 * check both terms before differencing them.
					 */
					double tv = ae[1] > 0.0 ?
					    tropmodel(time, vpos, ae, opt->humi) : 0.0;
					double tm = ae[1] > 0.0 ?
					    tropmodel(time, mpos, ae, opt->humi) : 0.0;

					if (tv > 0.0 && tm > 0.0)
						dT_model = tv - tm;
				}
				out[nout].dT = dT_model +
				    (opt->apply_nondispersive ? dT : 0.0);
				if (qual != NULL)
					qual[nout].dT_height = dT_model;
			}
			out[nout].sat = sat;
			out[nout].dI = withhold_I ? 0.0 :
			    dI * (opt->iono_scale > 0.0 ? opt->iono_scale : 1.0);
			out[nout].valid = 1;
			if (qual != NULL) {
				qual[nout].sat = sat;
				qual[nout].naux = m;
				qual[nout].backup_used = used_backup;
				qual[nout].leverage = lev;
				qual[nout].sd_I = sqrt(varI > 0.0 ? varI : 0.0);
				qual[nout].sd_T = sqrt(varT > 0.0 ? varT : 0.0);
			}
			nout++;
		}
	}
	return nout;
}

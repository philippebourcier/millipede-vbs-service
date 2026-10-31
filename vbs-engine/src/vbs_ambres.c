#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_ambres.h"

/*
 * The §4.3 cascade. See the header for why it exists and what each step is worth.
 */

void vbs_ambres_default_opt(struct vbs_ambres_opt *opt) {
	memset(opt, 0, sizeof *opt);
	opt->min_ratio = 3.0;		/* the preflight's working value (data/README.md §6) */
	opt->min_elev_deg = 20.0;	/* above the filter's own 15° mask: fixing is stricter */
	opt->max_float_sigma = 0.5;
	opt->q_inflate = 1.0;	/* cycles; beyond this a step is not worth attempting */
	opt->fix_var = 1e-6;		/* a held integer, not a measurement */
	opt->min_fixed = 2;
	opt->partial = 1;
	opt->max_step = VBS_AMB_STEP_NL;	/* the whole cascade, unless a caller stops it short */
	opt->if_narrow_lane = 0;		/* D67; off until measured, then it should be the default */
	opt->if_min_epochs = 20;
	opt->if_max_sigma = 0.15;		/* cycles of the narrow lane: 1.6 cm on GPS L1/L2 */
	opt->if_sys_m = 0.005;			/* troposphere residual that averaging cannot remove */
	opt->if_max_frac = 0.25;		/* how near an integer the mean must sit */
	opt->cascade_by_lambda = 1;		/* D74; 0 restores the pre-D74 table order */
	opt->use_band6 = 1;			/* D76; 0 drops Galileo E6 / BeiDou B3I again */
	opt->skip_wl12_sys = 0;			/* D129 diagnostic; see the header */
	opt->wl12_derive_sys = 0;		/* D131; see the header */
	opt->one_ref_sys = 0;			/* D132; see the header */
}

/*
 * The cascade steps, longest wavelength first.
 *
 * Bands are RINEX digits, so the same table covers every constellation that happens to carry
 * them; a constellation that does not simply yields no DD set for that step. BeiDou's digits
 * are the ones worth reading carefully -- 2 is B1I (1561 MHz) and 1 is B1C (1575 MHz), and
 * 7 is B2b, 6 is B3I (see vbs_signal_band_freq).
 */
static const struct vbs_amb_comb CASCADE[] = {
	/* Extra-wide lane: needs three bands. GPS L2-L5, Galileo E5a-E5b, BeiDou B2a-B2b. */
	{"EWL 2-5", 2, {2, 5}, {+1.0, -1.0}},
	{"EWL 5-7", 2, {5, 7}, {+1.0, -1.0}},
	/*
	 * Band 6 -- Galileo E6 and BeiDou B3I (D76). Neither appeared here before, so their
	 * ambiguity states were estimated by the filter and never fixed by anything: carried as
	 * unconstrained states, and their information thrown away. Both are well populated on
	 * this network -- E6 on 93-94 % of Galileo satellite-epochs at a UM980 (the Mosaic-X5
	 * does not stream it at all) and B3I on 99-100 % of BeiDou at *both* receivers.
	 *
	 * They are worth having for their wavelengths, not merely for the extra observable:
	 * E5a-E6 and B2a-B3I are 2.93 and 3.26 m, E6-E5b and B3I-B2b are 4.19 and 4.88 m, and
	 * the band-6 wide lanes at 0.98-1.02 m are *longer* than every wide lane already in the
	 * table. No new states are needed -- VBS_NET_MAXBAND is 5, and band 6 was already
	 * occupying one of those slots.
	 */
	{"EWL 6-7", 2, {6, 7}, {+1.0, -1.0}},
	{"EWL 5-6", 2, {5, 6}, {+1.0, -1.0}},
	/* Wide lane. */
	{"WL 1-2", 2, {1, 2}, {+1.0, -1.0}},
	{"WL 2-6", 2, {2, 6}, {+1.0, -1.0}},
	{"WL 1-6", 2, {1, 6}, {+1.0, -1.0}},
	{"WL 1-7", 2, {1, 7}, {+1.0, -1.0}},
	{"WL 1-5", 2, {1, 5}, {+1.0, -1.0}},
	/* Narrow lane: the individual carriers. */
	{"NL 1", 1, {1}, {+1.0}},
	{"NL 2", 1, {2}, {+1.0}},
	{"NL 5", 1, {5}, {+1.0}},
	{"NL 6", 1, {6}, {+1.0}},
	{"NL 7", 1, {7}, {+1.0}},
};

/* Is this the L1/L2 wide lane? Used only by the `skip_wl12_sys` measurement toggle. */
static int is_wl12(const struct vbs_amb_comb *comb) {
	return comb->nband == 2 && comb->band[0] == 1 && comb->band[1] == 2;
}

/* Does this combination need band 6? Used only by the `use_band6` measurement toggle. */
static int needs_band6(const struct vbs_amb_comb *comb) {
	int k;

	for (k = 0; k < comb->nband; k++)
		if (comb->band[k] == 6)
			return 1;
	return 0;
}
#define NCASCADE	((int)(sizeof CASCADE / sizeof CASCADE[0]))

/* Row index of a cascade combination by name, or -1. Looked up rather than hardcoded so that
 * inserting a row cannot silently repoint the derivation at the wrong combination. */
static int comb_row(const char *name) {
	int c;

	for (c = 0; c < NCASCADE; c++)
		if (!strcmp(CASCADE[c].name, name))
			return c;
	return -1;
}

/*
 * Constellations the cascade runs over. GLONASS is included (D90); whether it is actually
 * processed is decided by `sys_mask`, which the caller gates on cluster receiver homogeneity --
 * its FDMA inter-channel biases only cancel in a double difference when both ends have matching
 * hardware delays.
 */
static const int SYSLIST[] = {SYS_GPS, SYS_GLO, SYS_GAL, SYS_CMP};
#define NSYSLIST	((int)(sizeof SYSLIST / sizeof SYSLIST[0]))

/*
 * Wavelength of a combination on one constellation, or 0 if that constellation does not carry
 * every band the combination needs.
 *
 * The sign follows sum(coef·f), so a combination whose bands happen to be listed
 * low-frequency-first comes out negative -- Galileo and BeiDou "EWL 5-7" both do. That is
 * harmless, because the integer search works in cycles and the same lambda converts the result
 * back, but it does mean every comparison here has to be on the magnitude.
 */
static double comb_lambda(int sys, const struct vbs_amb_comb *comb) {
	double sumf = 0.0;
	int k;

	for (k = 0; k < comb->nband; k++) {
		double f = vbs_signal_band_freq(sys, comb->band[k]);

		if (f == 0.0)
			return 0.0;
		sumf += comb->coef[k] * f;
	}
	return fabs(sumf) < 1.0 ? 0.0 : CLIGHT / sumf;
}

/*
 * Which cascade bucket a combination belongs to -- from its wavelength, not from its row in
 * the table (D74).
 *
 * The table is keyed on RINEX band digits so that one set of rows covers every constellation,
 * but the same digit pair is a different physical combination on each: {2,5} is GPS L2-L5 at
 * 5.86 m and BeiDou B1I-B2a at 0.78 m. Classifying by row put that 0.78 m BeiDou combination
 * in the extra-wide-lane bucket, where `max_step` admits it as an essentially free integer.
 *
 * 2 m is well clear of both groups: the genuine extra-wide lanes present on this network are
 * 5.86, 9.77 and 20.93 m, and every true wide lane is 0.75-0.86 m.
 */
#define VBS_AMB_EWL_MIN_M	2.0

static int step_of_lambda(const struct vbs_amb_comb *comb, double lam) {
	if (comb->nband == 1)
		return VBS_AMB_STEP_NL;
	return fabs(lam) >= VBS_AMB_EWL_MIN_M ? VBS_AMB_STEP_EWL : VBS_AMB_STEP_WL;
}

/*
 * The pre-D74 bucket classification: by row in the table. Kept so the change can be measured
 * against the behaviour it replaced from a single binary (`--cascade-order table`), rather than
 * against numbers recorded by an older build on a different day.
 */
static int step_of_row(int i) {
	/*
	 * Row indices of the *original* table, so that `--cascade-order table` still reproduces
	 * pre-D74 behaviour after D76 inserted rows. Band-6 rows had no pre-D74 behaviour to
	 * reproduce; they are classified by wavelength like everything else.
	 */
	if (needs_band6(&CASCADE[i]))
		/*
		 * Galileo stands in for "a constellation that has band 6" here: this arm exists
		 * only to reproduce pre-D74 row classification, and band 6 had no pre-D74 rows.
		 */
		return step_of_lambda(&CASCADE[i], comb_lambda(SYS_GAL, &CASCADE[i]));
	if (CASCADE[i].nband == 1)
		return VBS_AMB_STEP_NL;
	return CASCADE[i].band[0] == 1 ? VBS_AMB_STEP_WL : VBS_AMB_STEP_EWL;
}

/* One (constellation, combination) to attempt, with the wavelength that orders it. */
struct cascade_step {
	int c;				/* row in CASCADE, i.e. the comb_id */
	int s;				/* index into SYSLIST */
	int step;			/* VBS_AMB_STEP_*, derived from lam */
	double lam;			/* magnitude, metres */
};

/*
 * The cascade's running order for one set of constellations: every (constellation,
 * combination) the constellations can actually form, longest wavelength first.
 *
 * `out` must hold NCASCADE * NSYSLIST entries. Returns how many were written.
 */
static int build_order(int sys_mask, int by_lambda, int use_band6, int skip_wl12_sys,
	struct cascade_step *out) {
	int n = 0, c, s, i, j;

	for (c = 0; c < NCASCADE; c++)
		for (s = 0; s < NSYSLIST; s++) {
			double lam;

			if (!(SYSLIST[s] & sys_mask))
				continue;
			if (!use_band6 && needs_band6(&CASCADE[c]))
				continue;
			if ((SYSLIST[s] & skip_wl12_sys) && is_wl12(&CASCADE[c]))
				continue;
			if ((lam = comb_lambda(SYSLIST[s], &CASCADE[c])) == 0.0)
				continue;		/* constellation lacks a band */
			out[n].c = c;
			out[n].s = s;
			out[n].lam = fabs(lam);
			out[n].step = by_lambda ? step_of_lambda(&CASCADE[c], lam)
						: step_of_row(c);
			n++;
		}
	if (!by_lambda)
		return n;			/* table order, as built */
	/*
	 * Insertion sort, descending. At most 27 entries, and stable on ties, so
	 * equal-wavelength rows (every constellation's NL 1, for one) keep table order and the
	 * sequence stays deterministic.
	 */
	for (i = 1; i < n; i++) {
		struct cascade_step tmp = out[i];

		for (j = i; j > 0 && out[j - 1].lam < tmp.lam; j--)
			out[j] = out[j - 1];
		out[j] = tmp;
	}
	return n;
}

/*
 * Integer least squares on the DD float vector.
 *
 * MRTKLIB's lambda() does the decorrelation and search (spec §0: do not re-implement
 * well-tested GNSS primitives; §4.3 names LAMBDA/MLAMBDA explicitly). It returns the two best
 * candidates and their residual norms, which is exactly what the ratio test needs.
 *
 * Writes the best integer vector to `fix` and the ratio to *ratio. Returns 1 on success.
 */
static int ils(const double *a, const double *Q, int n, double *fix, double *ratio) {
	double *F = NULL, s[2];
	int ok = 0;

	if (n < 1)
		return 0;
	if ((F = (double *)malloc(sizeof(double) * n * 2)) == NULL)
		return 0;
	if (lambda(n, 2, a, Q, F, s) == 0 && s[0] > 0.0) {
		int i;

		for (i = 0; i < n; i++)
			fix[i] = F[i];
		/*
		 * s[] are squared residual norms of the best and second-best candidates, so the
		 * conventional ratio is s[1]/s[0]. A perfect fit gives s[0] == 0; clamp so the
		 * test reports a large ratio rather than dividing by zero.
		 */
		*ratio = s[0] > 1e-12 ? s[1] / s[0] : 999.9;
		ok = 1;
	}
	free(F);
	return ok;
}

/*
 * Index (into the compressed set) of the ambiguity with the largest float variance -- the one
 * partial resolution drops first, since it is the one least likely to be carrying a correct
 * integer and most likely to be dragging the ratio down.
 */
static int worst_index(const double *Q, int n) {
	int i, worst = -1;
	double vmax = -1.0;

	for (i = 0; i < n; i++)
		if (Q[i + i * n] > vmax) {
			vmax = Q[i + i * n];
			worst = i;
		}
	return worst;
}

/* Extract the sub-vector and sub-matrix selected by keep[] into a/Qs. Returns its size. */
static int compress(const double *a0, const double *Q0, int n, const int *keep,
	double *a, double *Qs, int *map) {
	int i, j, ni = 0, nj;

	for (i = 0; i < n; i++) {
		if (!keep[i])
			continue;
		a[ni] = a0[i];
		map[ni] = i;
		ni++;
	}
	for (i = 0; i < ni; i++)
		for (j = 0, nj = 0; j < n; j++) {
			if (!keep[j])
				continue;
			Qs[i + nj * ni] = Q0[map[i] + j * n];
			nj++;
		}
	return ni;
}

/*
 * Resolve one (constellation, combination). Returns the number of ambiguities constrained.
 */
static int resolve_one(struct vbs_netekf *ekf, const struct vbs_ambres_opt *opt, int sys,
	const struct vbs_amb_comb *comb, int comb_id, int step, int prefer_ref,
	struct vbs_ambres_stats *stats, struct vbs_amb_fix *fixed, int maxfixed, int *nfixed) {
	struct vbs_dd_amb dd[VBS_AMB_MAXDD];
	static double Q[VBS_AMB_MAXDD * VBS_AMB_MAXDD];
	static double Qs[VBS_AMB_MAXDD * VBS_AMB_MAXDD];
	double a0[VBS_AMB_MAXDD], a[VBS_AMB_MAXDD], fix[VBS_AMB_MAXDD], ratio = 0.0;
	int keep[VBS_AMB_MAXDD], map[VBS_AMB_MAXDD];
	int n, i, ni, nfix = 0, dropped = 0;

	if ((n = vbs_netekf_dd_set_ref(ekf, sys, comb, prefer_ref, dd, VBS_AMB_MAXDD, Q)) < 2)
		return 0;
	if (opt->q_inflate != 1.0)
		for (i = 0; i < n * n; i++)
			Q[i] *= opt->q_inflate;

	/*
	 * Screen before searching. A low satellite or a wildly uncertain float ambiguity does
	 * not just fail to fix -- included in the search it distorts the decorrelation and can
	 * drag a whole otherwise-good subset below the ratio threshold.
	 */
	for (i = 0, ni = 0; i < n; i++) {
		a0[i] = dd[i].value;
		keep[i] = dd[i].el >= opt->min_elev_deg * D2R &&
		    Q[i + i * n] > 0.0 && sqrt(Q[i + i * n]) <= opt->max_float_sigma;
		ni += keep[i];
	}
	if (ni < opt->min_fixed)
		return 0;

	stats->attempts[step]++;

	for (;;) {
		int m = compress(a0, Q, n, keep, a, Qs, map);

		if (m < opt->min_fixed)
			break;
		if (!ils(a, Qs, m, fix, &ratio))
			break;
		if (ratio >= opt->min_ratio) {
			for (i = 0; i < m; i++)
				if (vbs_netekf_constrain_dd(ekf, &dd[map[i]], comb, fix[i],
				    opt->fix_var)) {
					nfix++;
					if (fixed != NULL && nfixed != NULL &&
					    *nfixed < maxfixed) {
						struct vbs_amb_fix *f = &fixed[(*nfixed)++];

						f->sat = dd[map[i]].sat;
						f->refsat = dd[map[i]].refsat;
						f->sys = sys;
						f->comb_id = comb_id;
						f->step = step;
						f->integer = fix[i];
						f->value = a0[map[i]];
						f->sigma = sqrt(Q[map[i] + map[i] * n]);
					}
				}
			stats->last_ratio[step] = ratio;
			stats->last_ndd[step] = n;
			stats->last_nfix[step] = nfix;
			stats->fixed[step] += nfix;
			if (dropped)
				stats->partial_used[step]++;
			return nfix;
		}
		stats->ratio_fail[step]++;
		if (!opt->partial)
			break;
		/* Drop the least certain and try the remaining subset (§4.3 partial AR). */
		if ((i = worst_index(Qs, m)) < 0)
			break;
		keep[map[i]] = 0;
		dropped++;
	}
	stats->last_ratio[step] = ratio;
	stats->last_ndd[step] = n;
	stats->last_nfix[step] = 0;
	return 0;
}

/* Defined below, next to the explanation of why the raw carrier cannot serve at network range. */
static int resolve_narrow_if(struct vbs_netekf *ekf, const struct vbs_ambres_opt *opt, int sys,
	const struct vbs_amb_fix *fix, int nfix, struct vbs_ambres_stats *stats,
	struct vbs_amb_fix *fixed, int maxfixed, int *nfixed);

/*
 * Derive the L1/L2 wide lane from the two clean ones instead of resolving it (D131).
 *
 * N(1-2) = N1 - N2 = (N1 - N5) - (N2 - N5) = N(1-5) - N(2-5), exactly, as integers. D129
 * measured WL 1-2's float sitting 0.055-0.111 cycles from integer on a UM980 when the cascade
 * rounds, with 10-27 % of decisions within 0.05 cycles of the 0.25 tolerance, because it is
 * the only cascade row built from the L1/L2 Melbourne-Wubbena combination and a UM980 reports
 * L2C under the L2W label so that code carries two different satellite biases. WL 1-5 and
 * EWL 2-5 on the same baselines and epochs measure 0.000.
 *
 * The two inputs are read back from the FILTER, not from the caller's `fixed` array, which is
 * a diagnostic the caller sizes (the pipeline passes 256). A first version read that array and
 * produced a byte-identical emit on one cluster while changing 34 % of the integers on
 * another; the buffer bound is NOT the explanation -- `vbs_ambres_resolve` resets it every
 * epoch and a busy epoch holds about ninety entries -- so why that version was inert on some
 * clusters is not established, and this one does not depend on the answer. A constrained
 * ambiguity's float sits on its integer with the variance of the constraint, so the filter
 * itself carries both the value and the evidence that it was fixed; `FIXED_SIGMA_CYC` is that
 * test, far below the 0.25-cycle rounding tolerance and far above `fix_var`. Reading the
 * filter is also the honest source: it is the thing the constraint actually acted on.
 */
#define FIXED_SIGMA_CYC 0.05

static int derived_input(const struct vbs_netekf *ekf, int sys, int row, int sat, int refsat,
	double *integer) {
	struct vbs_dd_amb dd[VBS_AMB_MAXDD];
	static double Q[VBS_AMB_MAXDD * VBS_AMB_MAXDD];
	int n, i;

	if (row < 0 || (n = vbs_netekf_dd_set(ekf, sys, &CASCADE[row], dd, VBS_AMB_MAXDD, Q)) < 1)
		return 0;
	for (i = 0; i < n; i++) {
		if (dd[i].sat != sat || dd[i].refsat != refsat)
			continue;
		if (!(Q[i + i * n] > 0.0) || sqrt(Q[i + i * n]) > FIXED_SIGMA_CYC)
			return 0;		/* not fixed: rounding it would be the very error
						 * this function exists to avoid */
		*integer = floor(dd[i].value + 0.5);
		return 1;
	}
	return 0;
}

static int resolve_wl12_derived(struct vbs_netekf *ekf, const struct vbs_ambres_opt *opt,
	int sys, int prefer_ref, struct vbs_ambres_stats *stats, struct vbs_amb_fix *fixed,
	int maxfixed, int *nfixed) {
	struct vbs_dd_amb dd[VBS_AMB_MAXDD];
	static double Q[VBS_AMB_MAXDD * VBS_AMB_MAXDD];
	int c12 = comb_row("WL 1-2"), c15 = comb_row("WL 1-5"), c25 = comb_row("EWL 2-5");
	int n, i, nfix = 0;

	if (c12 < 0 || c15 < 0 || c25 < 0)
		return 0;
	if ((n = vbs_netekf_dd_set_ref(ekf, sys, &CASCADE[c12], prefer_ref, dd, VBS_AMB_MAXDD, Q)) < 1)
		return 0;
	for (i = 0; i < n; i++) {
		double n15, n25;

		/*
		 * The reference satellite must be the same in all three sets. dd_set picks the
		 * highest-elevation satellite carrying every band the combination needs, and WL 1-2
		 * needs bands 1 and 2 where WL 1-5 needs 1 and 5, so the two can differ. A
		 * difference taken across two references is not the wide lane, so it is refused
		 * rather than turned into a plausible wrong integer.
		 */
		if (!derived_input(ekf, sys, c15, dd[i].sat, dd[i].refsat, &n15) ||
		    !derived_input(ekf, sys, c25, dd[i].sat, dd[i].refsat, &n25))
			continue;
		stats->attempts[VBS_AMB_STEP_WL]++;
		if (!vbs_netekf_constrain_dd(ekf, &dd[i], &CASCADE[c12], n15 - n25, opt->fix_var))
			continue;
		nfix++;
		stats->fixed[VBS_AMB_STEP_WL]++;
		if (fixed != NULL && nfixed != NULL && *nfixed < maxfixed) {
			struct vbs_amb_fix *f = &fixed[(*nfixed)++];

			f->sat = dd[i].sat;
			f->refsat = dd[i].refsat;
			f->sys = sys;
			f->comb_id = c12;
			f->step = VBS_AMB_STEP_WL;
			f->integer = n15 - n25;
			f->value = dd[i].value;	/* the float it did NOT round, kept for the dump */
			f->sigma = sqrt(Q[i + i * n]);
		}
	}
	return nfix;
}

int vbs_ambres_resolve(struct vbs_netekf *ekf, const struct vbs_ambres_opt *opt,
	int sys_mask, struct vbs_ambres_stats *stats,
	struct vbs_amb_fix *fixed, int maxfixed, int *nfixed) {
	struct vbs_ambres_opt defopt;
	int ref_of[NSYSLIST];
	int total = 0, s;

	if (nfixed != NULL)
		*nfixed = 0;

	if (opt == NULL) {
		vbs_ambres_default_opt(&defopt);
		opt = &defopt;
	}
	/*
	 * Longest wavelength first, and each fix is applied before the next step runs -- that
	 * ordering is the entire mechanism. Fixing the extra-wide lane tightens the individual
	 * ambiguities through the covariance, which is what makes the narrow lane reachable at
	 * baselines where it would otherwise not be.
	 *
	 * The order is therefore built here from the combinations' actual wavelengths on the
	 * constellations in play, rather than being taken from the table's row order (D74). The
	 * rows are keyed on band digits, and the same digits are a different combination on each
	 * constellation, so row order delivered that promise only for GPS. On BeiDou it ran
	 * 0.78 m -> 9.77 m -> 20.93 m: integers committed on the *shortest* wavelength in the set
	 * before the longest was attempted, with the 20.93 m B1C-B1I wide lane -- the easiest
	 * integer available anywhere in GNSS -- reached third, by which point the cascade had
	 * already constrained the ambiguities it feeds.
	 *
	 * The receiver-dependent part is what makes this matter in practice. A Mosaic-X5 on this
	 * network reports BeiDou bands 1, 2, 6 and 7, so the 0.78 m row (which needs band 5) is
	 * absent for it and its cascade started, correctly, at 20.93 m. A UM980 also reports
	 * B2a, so the row existed and its cascade started at 0.78 m: the receiver supplying
	 * *more* signals was penalised for supplying them. Deriving the order from the physics
	 * of whatever bands a station actually provides removes the whole class of problem
	 * instead of special-casing a receiver model.
	 */
	/*
	 * D132: one reference satellite per constellation for the whole cascade, when asked.
	 * Chosen from the union of every band the table uses, so the satellite picked is the one
	 * that can serve the most combinations; the rest fall back per combination as before.
	 */
	for (s = 0; s < NSYSLIST; s++) {
		int bands[NCASCADE * 2], nb = 0, c, k, q;

		ref_of[s] = 0;
		if (!(SYSLIST[s] & sys_mask & opt->one_ref_sys))
			continue;
		for (c = 0; c < NCASCADE; c++)
			for (k = 0; k < CASCADE[c].nband; k++) {
				for (q = 0; q < nb; q++)
					if (bands[q] == CASCADE[c].band[k])
						break;
				if (q == nb && nb < (int)(sizeof bands / sizeof bands[0]))
					bands[nb++] = CASCADE[c].band[k];
			}
		ref_of[s] = vbs_netekf_best_ref(ekf, SYSLIST[s], bands, nb);
	}

	{
		struct cascade_step ord[NCASCADE * NSYSLIST];
		int nord = build_order(sys_mask, opt->cascade_by_lambda, opt->use_band6,
		    opt->skip_wl12_sys | opt->wl12_derive_sys, ord), i;

		for (i = 0; i < nord; i++) {
			if (ord[i].step > opt->max_step)
				continue;
			/*
			 * With the ionosphere-free narrow lane in use the raw-carrier steps are
			 * skipped entirely: they are the ones that cannot work past ~45 km, and
			 * leaving them in would apply their wrong integers as hard constraints
			 * before the IF step ran.
			 */
			if (opt->if_narrow_lane && ord[i].step == VBS_AMB_STEP_NL)
				continue;
			total += resolve_one(ekf, opt, SYSLIST[ord[i].s], &CASCADE[ord[i].c],
			    ord[i].c, ord[i].step, ref_of[ord[i].s], stats, fixed, maxfixed,
			    nfixed);
		}
	}
	/* D131: the L1/L2 wide lane, taken as N(1-5) - N(2-5) now that both are fixed. */
	if (opt->wl12_derive_sys)
		for (s = 0; s < NSYSLIST; s++)
			if (SYSLIST[s] & sys_mask & opt->wl12_derive_sys)
				total += resolve_wl12_derived(ekf, opt, SYSLIST[s], ref_of[s],
				    stats, fixed, maxfixed, nfixed);

	if (opt->if_narrow_lane && fixed != NULL && nfixed != NULL)
		for (s = 0; s < NSYSLIST; s++) {
			if (!(SYSLIST[s] & sys_mask))
				continue;
			total += resolve_narrow_if(ekf, opt, SYSLIST[s], fixed, *nfixed, stats,
			    fixed, maxfixed, nfixed);
		}
	return total;
}

/* Exact entry for (sat, comb), whatever reference it used, or NULL. */
static const struct vbs_amb_fix *find_raw(const struct vbs_amb_fix *list, int n, int sat,
	int comb_id) {
	int i;

	for (i = 0; i < n; i++)
		if (list[i].sat == sat && list[i].comb_id == comb_id)
			return &list[i];
	return NULL;
}

/*
 * The DD ambiguity of (sat, ref) on one baseline, re-referenced if that baseline chose a
 * different reference satellite.
 *
 * A DD ambiguity is a difference of per-station single differences, so for any three
 * satellites N(s,r) = N(s,k) − N(r,k): changing the reference is exact arithmetic on integers,
 * not an approximation. Skipping a triangle whenever the three baselines disagreed on a
 * reference therefore threw away checks that were perfectly well defined — and reference
 * changes are frequent, since each filter picks the highest satellite independently.
 *
 * Returns 1 and writes the value, or 0 if this baseline lacks a fix for either satellite
 * against its own reference.
 */
static int dd_against(const struct vbs_amb_fix *list, int n, int sat, int ref, int comb_id,
	double *out) {
	const struct vbs_amb_fix *fs = find_raw(list, n, sat, comb_id), *fr;

	if (fs == NULL)
		return 0;
	if (fs->refsat == ref) {
		*out = fs->integer;
		return 1;
	}
	/* This baseline used k = fs->refsat; N(s,r) = N(s,k) − N(r,k). */
	if ((fr = find_raw(list, n, ref, comb_id)) == NULL || fr->refsat != fs->refsat)
		return 0;
	*out = fs->integer - fr->integer;
	return 1;
}

int vbs_ambres_loop_check(const struct vbs_amb_fix *ma, int nma,
	const struct vbs_amb_fix *mb, int nmb,
	const struct vbs_amb_fix *ab, int nab,
	double tol_cycles, struct vbs_amb_loop_stats *stats,
	struct vbs_amb_fix *bad, int maxbad) {
	int i, nbad = 0;
	/* Per-satellite: did a wide- / extra-wide-lane combination misclose here? First pass. */
	static unsigned char wl_failed[MAXSAT], ewl_failed[MAXSAT];

	if (tol_cycles <= 0.0)
		tol_cycles = 0.25;	/* half the distance to the next integer, with margin */

	memset(wl_failed, 0, sizeof wl_failed);
	memset(ewl_failed, 0, sizeof ewl_failed);
	for (i = 0; i < nma; i++) {
		double a, b, c;

		if (ma[i].step != VBS_AMB_STEP_WL && ma[i].step != VBS_AMB_STEP_EWL)
			continue;
		if (ma[i].sat <= 0 || ma[i].sat > MAXSAT)
			continue;
		if (!dd_against(ma, nma, ma[i].sat, ma[i].refsat, ma[i].comb_id, &a) ||
		    !dd_against(mb, nmb, ma[i].sat, ma[i].refsat, ma[i].comb_id, &b) ||
		    !dd_against(ab, nab, ma[i].sat, ma[i].refsat, ma[i].comb_id, &c))
			continue;
		if (fabs(c - (b - a)) > tol_cycles) {
			if (ma[i].step == VBS_AMB_STEP_EWL)
				ewl_failed[ma[i].sat - 1] = 1;
			else
				wl_failed[ma[i].sat - 1] = 1;
		}
	}

	for (i = 0; i < nma; i++) {
		int ref = ma[i].refsat, sat = ma[i].sat, cid = ma[i].comb_id;
		double nma_v, nmb_v, nab_v, mis;

		/*
		 * All three are expressed against MA's reference, re-referencing where a baseline
		 * chose its own. The identity is exact in any common reference; what it cannot
		 * tolerate is the three being expressed against *different* ones.
		 */
		if (!dd_against(ma, nma, sat, ref, cid, &nma_v) ||
		    !dd_against(mb, nmb, sat, ref, cid, &nmb_v) ||
		    !dd_against(ab, nab, sat, ref, cid, &nab_v)) {
			if (stats != NULL)
				stats->skipped_ref++;
			continue;	/* not fixed on all three, or not re-referenceable */
		}
		mis = nab_v - (nmb_v - nma_v);
		if (stats != NULL) {
			double a = fabs(mis), frac = fabs(a - floor(a + 0.5));
			int bin, st = ma[i].step;

			stats->checked++;
			if (a > stats->worst_misclosure)
				stats->worst_misclosure = a;
			if (a <= tol_cycles)
				bin = VBS_LOOP_BIN_CLOSED;
			else if (frac > 0.1)
				bin = VBS_LOOP_BIN_NONINT;
			else if (a < 1.5)
				bin = VBS_LOOP_BIN_1CYC;
			else if (a < 2.5)
				bin = VBS_LOOP_BIN_2CYC;
			else
				bin = VBS_LOOP_BIN_LARGE;
			if (st >= 0 && st < VBS_AMB_NSTEP)
				stats->hist[st][bin]++;
			if (st == VBS_AMB_STEP_NL && sat > 0 && sat <= MAXSAT) {
				int wl = wl_failed[sat - 1], ewl = ewl_failed[sat - 1];

				if (bin == VBS_LOOP_BIN_LARGE) {
					if (wl) stats->nl_tail_with_wl++;
					else    stats->nl_tail_without_wl++;
					if (ewl) stats->nl_tail_with_ewl++;
					else     stats->nl_tail_without_ewl++;
				} else if (bin == VBS_LOOP_BIN_1CYC) {
					if (wl) stats->nl_one_with_wl++;
					else    stats->nl_one_without_wl++;
				}
			}
			if (st == VBS_AMB_STEP_WL && sat > 0 && sat <= MAXSAT &&
			    bin == VBS_LOOP_BIN_LARGE) {
				if (ewl_failed[sat - 1]) stats->wl_tail_with_ewl++;
				else                     stats->wl_tail_without_ewl++;
			}
		}
		if (fabs(mis) > tol_cycles) {
			if (stats != NULL)
				stats->misclosed++;
			if (bad != NULL && nbad < maxbad)
				bad[nbad++] = ma[i];
		}
	}
	return nbad;
}

/*
 * --- D67: the narrow lane resolved on the ionosphere-free combination ---------------------
 *
 * The cascade's last steps fix the individual carriers, which needs the slant ionosphere known
 * to a fraction of a 19 cm wavelength. Past about 45 km it never is, and the measured result is
 * unambiguous: at 46-75 km the wide lanes misclose *zero* times beyond two cycles while the raw
 * carriers miss 6 755 times. Server-side network RTK does not attempt this. With the wide-lane
 * integer already known it resolves the remaining integer on the ionosphere-free combination,
 * which contains no ionosphere at all.
 *
 * For bands a and b of a satellite pair, writing the de-geometried single differences as zL and
 * the combination weights w = 1/lambda^2:
 *
 *   IF = (w_a*zL_a - w_b*zL_b) / (w_a - w_b)
 *
 * whose ambiguity content works out to lambda_NL*N_a + lambda_NL*(lam_a/(lam_b-lam_a))*N_WL,
 * with lambda_NL = lam_a*lam_b/(lam_a+lam_b) -- 10.7 cm for GPS L1/L2. Double-differenced, the
 * clock and the receiver biases are gone and only the troposphere remains, which is already a
 * state. So with N_WL substituted, N_a follows from the measurements with no ionosphere in the
 * model and no dependence on baseline length.
 *
 * One epoch of it is too noisy to round: the IF combination amplifies phase noise about
 * sixfold, which on a 10.7 cm wavelength is a third of a cycle. The float is constant over an
 * arc, so it is averaged down instead -- the standard practice, and what makes the difference
 * between a step that works at 30 km and one that works at 100.
 */
static int resolve_narrow_if(struct vbs_netekf *ekf, const struct vbs_ambres_opt *opt, int sys,
	const struct vbs_amb_fix *fix, int nfix, struct vbs_ambres_stats *stats,
	struct vbs_amb_fix *fixed, int maxfixed, int *nfixed) {
	int k, nfixed_here = 0;

	for (k = 0; k < nfix; k++) {
		const struct vbs_amb_fix *w = &fix[k];
		const struct vbs_amb_comb *cw, *cn = NULL;
		double za, zb, ra, rb, lama, lamb, mua, mub, wa, wb;
		double mfs, mfr, lnl, ifs, ifr, ddif, flt, mean, integer;
		double sig1, signl;
		int c2, cn_id = -1, cnt = 0;
		struct vbs_dd_amb dd;

		if (w->sys != sys || w->comb_id < 0 || w->comb_id >= NCASCADE)
			continue;
		cw = &CASCADE[w->comb_id];
		if (!(cw->nband == 2 && cw->coef[0] == 1.0 && cw->coef[1] == -1.0))
			continue;		/* wide lanes only: this step consumes their integers */
		/* The plain-carrier entry for band a, whose slot in the fix list we reuse. */
		for (c2 = 0; c2 < NCASCADE; c2++)
			if (CASCADE[c2].nband == 1 && CASCADE[c2].band[0] == cw->band[0] &&
			    CASCADE[c2].coef[0] == 1.0) {
				cn = &CASCADE[c2];
				cn_id = c2;
			}
		if (cn == NULL)
			continue;
		if (!vbs_netekf_phase_sd(ekf, w->sat, cw->band[0], &za, &mua, &lama) ||
		    !vbs_netekf_phase_sd(ekf, w->sat, cw->band[1], &zb, &mub, &lamb) ||
		    !vbs_netekf_phase_sd(ekf, w->refsat, cw->band[0], &ra, NULL, NULL) ||
		    !vbs_netekf_phase_sd(ekf, w->refsat, cw->band[1], &rb, NULL, NULL))
			continue;
		if (!vbs_netekf_slot_mapfw(ekf, w->sat, &mfs, NULL) ||
		    !vbs_netekf_slot_mapfw(ekf, w->refsat, &mfr, NULL))
			continue;
		if (lamb <= lama)
			continue;		/* band b must be the lower frequency */
		wa = 1.0 / (lama * lama);
		wb = 1.0 / (lamb * lamb);
		if (fabs(wa - wb) < 1e-9)
			continue;
		lnl = lama * lamb / (lama + lamb);
		ifs = (wa * za - wb * zb) / (wa - wb);
		ifr = (wa * ra - wb * rb) / (wa - wb);
		ddif = (ifs - ifr) - (mfs - mfr) * vbs_netekf_ztd(ekf);
		flt = (ddif - lnl * (lama / (lamb - lama)) * w->integer) / lnl;
		if (!vbs_netekf_if_average(ekf, w->sat, w->refsat, cw->band[0], flt, &mean, &cnt))
			continue;
		/*
		 * Noise of the double-differenced IF float in cycles of the narrow lane, divided
		 * down by the number of epochs averaged. The IF coefficients are w/(wa-wb), and
		 * eight phase measurements enter a double difference.
		 */
		sig1 = sqrt(wa * wa + wb * wb) / fabs(wa - wb) * 2.0 * vbs_netekf_sig_phase(ekf);
		/*
		 * Averaging beats down phase noise and nothing else. What remains is the residual
		 * troposphere the ZTD state did not absorb, which is correlated over minutes, so
		 * dividing the whole sigma by sqrt(n) claims a precision the float does not have:
		 * measured against the raw-carrier step, that optimism put 1.66 % of the integers
		 * a cycle out, and a wrong integer held at 1e-6 cycles^2 is worse than no integer.
		 * The floor is what averaging cannot remove.
		 */
		signl = sqrt(sig1 * sig1 / cnt + opt->if_sys_m * opt->if_sys_m) / lnl;
		integer = floor(mean + 0.5);
		/*
		 * Round only when the mean is both precise and unambiguous. 0.25 cycles of the
		 * narrow lane is 2.7 cm, comfortably inside what the averaged IF resolves, and a
		 * float sitting near a half integer is refused rather than guessed.
		 */
		if (cnt < opt->if_min_epochs || signl > opt->if_max_sigma ||
		    fabs(mean - integer) > opt->if_max_frac)
			continue;
		stats->attempts[VBS_AMB_STEP_NL]++;
		dd.sat = w->sat;
		dd.refsat = w->refsat;
		dd.sys = sys;
		dd.value = mean;
		dd.lam = lnl;
		dd.el = 0.0;
		if (!vbs_netekf_constrain_dd(ekf, &dd, cn, integer, opt->fix_var))
			continue;
		nfixed_here++;
		stats->fixed[VBS_AMB_STEP_NL]++;
		if (fixed != NULL && nfixed != NULL && *nfixed < maxfixed) {
			struct vbs_amb_fix *f = &fixed[(*nfixed)++];

			f->sat = w->sat;
			f->refsat = w->refsat;
			f->sys = sys;
			f->comb_id = cn_id;
			f->step = VBS_AMB_STEP_NL;
			f->integer = integer;
			f->value = mean;
			f->sigma = signl;
		}
	}
	return nfixed_here;
}

/*
 * --- D66: the ionosphere as an output of ambiguity resolution -----------------------------
 *
 * The filter estimates a slant ionosphere per satellite as a random-walk state, and
 * `vbs_netekf_corrections()` refuses to emit any correction whose variance has not come down
 * to 0.05 m. Both are wrong at network range. Server-side network RTK fixes the wide lane with
 * Melbourne-Wübbena (geometry-free and first-order ionosphere-free, so it holds over hundreds
 * of kilometres), then resolves the remaining integer on the ionosphere-free combination, and
 * only *then* reads the ionosphere off the geometry-free combination with those integers
 * substituted. Estimating the ionosphere in order to fix, as this engine did, inverts the
 * dependency and is why every baseline past ~45 km failed: the state can never converge to
 * 5 cm at that range, so either the gate rejects the satellite or -- worse -- a too-tight
 * random walk lets the filter believe an ionosphere it has not measured and the error lands in
 * the ambiguities.
 *
 * With N1 and N2 known integers in the double difference, between bands a and b of satellites
 * s and r:
 *
 *   DD_GF = (zL_a - zL_b)_s - (zL_a - zL_b)_r = (mu_b - mu_a)·DD_I + lam_a·N_a - lam_b·N_b
 *
 * because the clock and troposphere cancel between bands and the receiver phase biases cancel
 * between satellites. So DD_I follows directly, to phase noise, with no state at all.
 *
 * What cannot be recovered this way is the reference satellite's own level: a double difference
 * is blind to it. That is exactly the datum D55 measured and showed to be constellation-common
 * to within 0.02 cm, and therefore to cancel in the double differences a rover forms -- so the
 * filter's estimate is used for the reference and the result is correct up to a term no rover
 * can see.
 */

/* The DD integers of one satellite pair on bands a and b, from the wide lane plus one carrier. */
static int dd_integers(const struct vbs_amb_fix *fix, int nfix, int sat, int ref,
	int band_a, int band_b, double *Na, double *Nb) {
	double wl = 0.0, na = 0.0;
	int have_wl = 0, have_na = 0, k;

	for (k = 0; k < nfix; k++) {
		const struct vbs_amb_fix *f = &fix[k];
		const struct vbs_amb_comb *c;

		if (f->sat != sat || f->refsat != ref || f->comb_id < 0 || f->comb_id >= NCASCADE)
			continue;
		c = &CASCADE[f->comb_id];
		if (c->nband == 2 && c->band[0] == band_a && c->band[1] == band_b &&
		    c->coef[0] == 1.0 && c->coef[1] == -1.0) {
			wl = f->integer;		/* N_a - N_b */
			have_wl = 1;
		} else if (c->nband == 1 && c->band[0] == band_a && c->coef[0] == 1.0) {
			na = f->integer;
			have_na = 1;
		}
	}
	if (!have_wl || !have_na)
		return 0;
	*Na = na;
	*Nb = na - wl;
	return 1;
}

/*
 * Why ion_from_fixes replaces nothing, counted rather than guessed.
 *
 * The feature is inert in the shipped configuration (0 of 1 750 839 corrections) and reached
 * only when --if-narrow-lane is removed, which independently ruins the solution. Isolating it
 * needs to know WHICH of its two preconditions fails: no usable reference satellite for the
 * constellation, or a reference with no wide-lane/narrow-lane pair against it. These counters
 * answer that in one run; `vbs_ambres_ion_diag()` reads them.
 */
static long ion_diag_norev, ion_diag_nopair, ion_diag_ok;

void vbs_ambres_ion_diag(long *no_ref, long *no_pair, long *replaced) {
	if (no_ref != NULL)
		*no_ref = ion_diag_norev;
	if (no_pair != NULL)
		*no_pair = ion_diag_nopair;
	if (replaced != NULL)
		*replaced = ion_diag_ok;
}

int vbs_ambres_ion_from_fixes(const struct vbs_netekf *ekf, const struct vbs_amb_fix *fix,
	int nfix, struct vbs_corr *corr, int ncorr, double *sd_I) {
	int k, n = 0, si;
	int sysref[NSYSLIST];

	if (ekf == NULL || corr == NULL)
		return 0;
	/*
	 * One reference satellite per constellation, for the whole epoch.
	 *
	 * A double difference cannot see its reference's own ionosphere, so every satellite
	 * referred to the same reference carries the same unknown level -- which is what D55
	 * measured to be invisible to a rover. Let each satellite keep the reference it happened
	 * to be fixed against and that guarantee is gone: each carries a *different* unknown, and
	 * the rover sees per-satellite bias. Measured, that error is not subtle -- it took M30
	 * from 98.4 % fixed to 13.2 % with 86 % of those fixes wrong and a 12.6 cm vertical bias,
	 * while the loop closure, which is blind to a per-satellite constant, still read 0.01 %.
	 *
	 * The most-used reference of each constellation is chosen, and satellites fixed against
	 * any other are left to the filter rather than mixed in.
	 */
	for (si = 0; si < NSYSLIST; si++) {
		int best = 0, bestn = 0, j;

		for (j = 0; j < nfix; j++) {
			int cand = fix[j].refsat, cnt = 0, q;

			if (fix[j].sys != SYSLIST[si] || cand == 0)
				continue;
			for (q = 0; q < nfix; q++)
				if (fix[q].sys == SYSLIST[si] && fix[q].refsat == cand)
					cnt++;
			if (cnt > bestn) {
				bestn = cnt;
				best = cand;
			}
		}
		sysref[si] = best;
	}
	for (k = 0; k < ncorr; k++) {
		double gfs, gfr, mua, mub, lama, lamb, za, zb, ra, rb, ddi, sig;
		double Na = 0.0, Nb = 0.0;
		int ref = 0, band_a = 0, band_b = 0, j, refidx = -1;

		/*
		 * A reference this satellite has *both* a wide lane and a carrier fixed against.
		 * Each cascade step picks its own highest-elevation reference per constellation,
		 * so taking the first fix's reference finds a pair that does not exist and the
		 * satellite is silently dropped -- which is what collapsed coverage to 2.3
		 * satellites per epoch the first time this ran.
		 */
		{
			int sysk = satsys(corr[k].sat, NULL), q;

			ref = 0;
			for (q = 0; q < NSYSLIST; q++)
				if (SYSLIST[q] == sysk)
					ref = sysref[q];
			if (ref == 0 || ref == corr[k].sat) {
				ion_diag_norev++;
				continue;
			}
		}
		for (j = 0; j < nfix && band_a == 0; j++) {
			int cand = fix[j].refsat, c2;

			if (fix[j].sat != corr[k].sat || cand != ref)
				continue;
			for (c2 = 0; c2 < NCASCADE; c2++) {
				double ma, mb;

				if (!(CASCADE[c2].nband == 2 && CASCADE[c2].coef[0] == 1.0 &&
				    CASCADE[c2].coef[1] == -1.0))
					continue;
				if (!dd_integers(fix, nfix, corr[k].sat, cand, CASCADE[c2].band[0],
				    CASCADE[c2].band[1], &Na, &Nb))
					continue;
				if (!vbs_netekf_phase_sd(ekf, corr[k].sat, CASCADE[c2].band[0], NULL, &ma, NULL) ||
				    !vbs_netekf_phase_sd(ekf, corr[k].sat, CASCADE[c2].band[1], NULL, &mb, NULL) ||
				    !vbs_netekf_phase_sd(ekf, cand, CASCADE[c2].band[0], NULL, NULL, NULL) ||
				    !vbs_netekf_phase_sd(ekf, cand, CASCADE[c2].band[1], NULL, NULL, NULL))
					continue;
				/*
				 * Reject a pair too close in frequency to separate the ionosphere
				 * *here*, while other pairs can still be tried. Checking it after
				 * the pair was committed instead abandoned the satellite: BeiDou's
				 * B1C-B1I sit 14 MHz apart, giving mu_b - mu_a of about 0.018, so
				 * whenever that pair was picked first the satellite's ionosphere
				 * was dropped rather than taken from one of the wide pairs the
				 * same satellite also carries.
				 */
				if (fabs(mb - ma) < 0.05)
					continue;
				ref = cand;
				band_a = CASCADE[c2].band[0];
				band_b = CASCADE[c2].band[1];
				break;
			}
		}
		if (band_a == 0) {
			ion_diag_nopair++;	/* a reference exists, but this satellite has no
						 * wide-lane/narrow-lane pair fixed against it */
			continue;
		}
		if (ref == 0)
			continue;
		for (j = 0; j < ncorr; j++)
			if (corr[j].sat == ref)
				refidx = j;
		if (refidx < 0)
			continue;
		/*
		 * The lowest two bands both satellites carry. Lowest because the wide lane the
		 * cascade fixes first is built on them, and because their ionosphere separation
		 * (mu_b - mu_a) is the largest, which is what divides the phase noise down.
		 */
		if (!vbs_netekf_phase_sd(ekf, corr[k].sat, band_a, &za, &mua, &lama) ||
		    !vbs_netekf_phase_sd(ekf, corr[k].sat, band_b, &zb, &mub, &lamb) ||
		    !vbs_netekf_phase_sd(ekf, ref, band_a, &ra, NULL, NULL) ||
		    !vbs_netekf_phase_sd(ekf, ref, band_b, &rb, NULL, NULL))
			continue;
		if (fabs(mub - mua) < 0.05)
			continue;		/* bands too close to separate the ionosphere */
		gfs = za - zb;
		gfr = ra - rb;
		ddi = ((gfs - gfr) - lama * Na + lamb * Nb) / (mub - mua);
		/*
		 * Eight phase measurements enter a double-differenced geometry-free combination,
		 * so its noise is sqrt(8)·sigma_phase, divided by the ionosphere separation.
		 */
		sig = 2.828427 * vbs_netekf_sig_phase(ekf) / fabs(mub - mua);
		corr[k].dI = corr[refidx].dI + ddi;
		corr[k].valid = 1;
		if (sd_I != NULL)
			sd_I[k] = sig;
		ion_diag_ok++;
		n++;
	}
	return n;
}

/*
 * The cascade order the resolver will use, for tests and diagnostics. Reports the same array
 * `vbs_ambres_resolve` runs, built by the same function -- a test against a reimplementation
 * of the ordering would pass while the engine did something else.
 */
int vbs_ambres_cascade_order(int sys_mask, int by_lambda, int use_band6,
	struct vbs_amb_step_info *out, int maxout) {
	struct cascade_step ord[NCASCADE * NSYSLIST];
	int n = build_order(sys_mask, by_lambda, use_band6, 0, ord), i;

	for (i = 0; i < n && i < maxout; i++) {
		out[i].sys = SYSLIST[ord[i].s];
		out[i].step = ord[i].step;
		out[i].lam = ord[i].lam;
		out[i].name = CASCADE[ord[i].c].name;
	}
	return n < maxout ? n : maxout;
}

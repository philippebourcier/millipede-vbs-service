#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_netekf.h"

/*
 * Network EKF for one baseline -- spec §4.3/§4.4. See the header for the state vector and
 * why it is single-difference.
 */

/* Constellations that get their own receiver-clock-difference state. */
#define NSYSIDX		5

/*
 * State layout. Ionosphere and ambiguity states are indexed by *slot*, not by satellite
 * number: MAXSAT is ~250 with every constellation enabled, and carrying a 250*(1+3) state
 * covariance would be 1000x1000 doubles per baseline for the ~35 satellites actually in view.
 * Slots are allocated when a satellite rises and released when it sets.
 */
#define IX_ZTD		0
/*
 * Horizontal troposphere gradients of the *difference* between the two stations, north and
 * east (D68). One zenith scalar per baseline cannot represent a tilted wet atmosphere: the
 * unmodelled part is elevation- and azimuth-dependent, reaches several centimetres on low
 * satellites, and -- because it biases every baseline of a cluster alike -- passes loop closure
 * untouched while corrupting the corrections. It is what stopped the ionosphere-free narrow
 * lane from working once the ionosphere itself was handled.
 */
#define IX_GRN		1
#define IX_GRE		2
#define IX_CLK		3
#define IX_ION		(IX_CLK + NSYSIDX)
#define IX_AMB		(IX_ION + VBS_NET_MAXSAT)
#define NX		(IX_AMB + VBS_NET_MAXSAT * VBS_NET_MAXBAND)

#define ION(slot)		(IX_ION + (slot))
#define AMB(slot, b)		(IX_AMB + (slot) * VBS_NET_MAXBAND + (b))

/*
 * One observation with its geometry already computed. The epoch is processed in two passes --
 * prepare, then update -- so that the receiver clock difference can be seeded from the whole
 * epoch's code observations before the first Kalman update runs. See the comment at the seed.
 */
struct prep {
	int slot;			/* satellite slot */
	int si;				/* constellation index, for the clock state */
	int bi;				/* band index within the slot, for the ambiguity state */
	double mu;			/* (f_ref/f)², the ionosphere coefficient */
	double lam;			/* carrier wavelength, m */
	double mapfw;			/* wet mapping function, the ΔZTD partial */
	double mgn, mge;		/* gradient partials: m_grad·cos(az), m_grad·sin(az) */
	double hclk;			/* receiver-clock partial, 1 - rdot_aux/c (see below) */
	double zP, zL;			/* code and phase single differences, m, de-geometried */
	double var_p, var_l;
};

struct slot {
	int sat;			/* 0 when free */
	int band[VBS_NET_MAXBAND];	/* RINEX band digits, ascending; 0 when unused */
	/*
	 * This satellite's carrier frequency on each of those bands, Hz.
	 *
	 * Held per satellite because GLONASS is FDMA: its satellites are separated by frequency,
	 * so a combination's wavelength is a property of the satellite and not of the
	 * constellation. The observation model already worked this way -- it takes each
	 * frequency from sat2freq() -- but the double-difference set computed one wavelength for
	 * a whole constellation, which is what confined the engine to CDMA systems.
	 */
	double freq[VBS_NET_MAXBAND];
	int nband;
	gtime_t last;
	double el;			/* master-side elevation, rad, last epoch */
	double az;			/* master-side azimuth, rad, last epoch */
	double mapg;			/* gradient mapping function, last epoch */
	double mapfw;			/* wet mapping function, last epoch */
	double dT_model;		/* a priori slant tropo difference, m, last epoch */
	long epochs_tracked;		/* consecutive epochs since this slot was allocated */
	long arc_master, arc_aux;	/* station arc ids the ambiguities rest on (D45) */
	int seen_this_epoch;
	int phase_rejects;		/* consecutive epochs with every phase of this satellite gated out */
	int phase_accepted;		/* any phase accepted this epoch */
};

struct vbs_netekf {
	double master_xyz[3], aux_xyz[3];
	double master_pos[3], aux_pos[3];	/* geodetic, for the tropo model */
	struct vbs_netekf_opt opt;

	double x[NX];
	double P[NX * NX];

	struct slot slot[VBS_NET_MAXSAT];
	int slot_of_sat[MAXSAT];		/* 1-based slot index, 0 = none */

	struct vbs_sigset sigset;
	struct prep prep[VBS_NET_MAXSAT * VBS_NET_MAXBAND];
	int nprep;			/* entries filled by the last epoch, for D66's geometry-free readout */
	/*
	 * Running mean of the ionosphere-free narrow-lane float, per satellite (D67). The float
	 * is constant over an arc, so averaging it down is free precision -- which is how network
	 * AR reaches long baselines: one epoch of an IF observable is too noisy to round, a few
	 * minutes of it is not. Reset when either station's arc for that satellite changes.
	 */
	struct if_acc {
		int sat, refsat, band;
		long arc_m, arc_a;
		double sum;
		int n;
	} ifacc[VBS_NET_MAXSAT];

	double sq_code, sq_phase, norm_code, norm_phase;
	long n_code, n_phase;

	gtime_t last_time;
	int have_last;
	struct vbs_netekf_stats stats;
};

void vbs_netekf_default_opt(struct vbs_netekf_opt *opt) {
	memset(opt, 0, sizeof *opt);
	opt->sys_mask = SYS_GPS | SYS_GAL | SYS_CMP;
	opt->elev_mask_deg = 15.0;
	opt->humi = 0.7;
	/*
	 * Zenith sigmas. The phase figure is deliberately loose for carrier phase: it is not
	 * the receiver's noise but the modelling error this filter does not have states for
	 * (per-satellite orbit and asymmetric troposphere, see the header), and setting it to
	 * true carrier noise would make the filter trust a model it should not.
	 */
	opt->sig_code = 0.6;
	opt->sig_phase = 0.006;
	opt->q_ztd = 1e-8;		/* ~6 mm/hour of ZTD difference wander */
	opt->q_ion = 1e-6;		/* ~5.5 mm per 30 s of slant ionosphere wander */
	opt->var_ztd0 = 0.3 * 0.3;
	/*
	 * Gradients are a few millimetres in normal weather and a couple of centimetres in a
	 * front; 30 mm of prior admits both without letting the state absorb anything else, and
	 * the random walk matches the ZTD's own timescale.
	 */
	opt->var_grad0 = 0.002 * 0.002;
	opt->q_grad = 1e-12;
	opt->max_mapg = 10.0;
	/*
	 * On by default (D69). Measured on both working cases: M30 false fixes 1.20 % -> 0.40 %
	 * at an unchanged 98 % fix rate, E1-0 92.6 % -> 93.3 % fixed with one fewer never-fixed
	 * window and time-to-first-fix p95 from beyond the window to 257 s. The states stay small
	 * -- the gain is that the filter stops forcing an azimuthally tilted atmosphere into the
	 * ambiguities.
	 */
	opt->est_gradients = 1;
	opt->var_ion0 = 0.5 * 0.5;
	/*
	 * Uncertainty around the per-epoch median code estimate of the clock difference, not
	 * the clock difference itself -- the state is seeded from the data every epoch, so this
	 * only has to cover the median's own error against code noise and multipath.
	 */
	opt->var_clk = 10.0 * 10.0;
	opt->var_amb0 = 1e8;
	opt->max_innov_sigma = 4.0;
	opt->gate_after_epochs = 10;
	opt->reject_streak_reset = 3;
	opt->min_snr_dbhz = 0.0;
	opt->prefer_l2c = 0;
	opt->drop_band5 = 0;
	opt->drop_band2 = 0;
}

static int cmp_double(const void *a, const void *b) {
	double x = *(const double *)a, y = *(const double *)b;

	return x < y ? -1 : x > y ? 1 : 0;
}

static int sysidx(int sys) {
	switch (sys) {
	case SYS_GPS: return 0;
	case SYS_GAL: return 1;
	case SYS_CMP: return 2;
	case SYS_GLO: return 3;
	case SYS_QZS: return 4;
	}
	return -1;
}

struct vbs_netekf *vbs_netekf_new(const double *master_xyz, const double *aux_xyz,
	const struct vbs_netekf_opt *opt) {
	struct vbs_netekf *this = (struct vbs_netekf *)calloc(1, sizeof *this);
	int i;

	if (this == NULL)
		return NULL;
	if (opt != NULL)
		this->opt = *opt;
	else
		vbs_netekf_default_opt(&this->opt);

	for (i = 0; i < 3; i++) {
		this->master_xyz[i] = master_xyz[i];
		this->aux_xyz[i] = aux_xyz[i];
	}
	ecef2pos(this->master_xyz, this->master_pos);
	ecef2pos(this->aux_xyz, this->aux_pos);

	this->x[IX_ZTD] = 0.0;
	this->P[IX_ZTD + IX_ZTD * NX] = this->opt.var_ztd0;
	this->x[IX_GRN] = this->x[IX_GRE] = 0.0;
	this->P[IX_GRN + IX_GRN * NX] = this->opt.var_grad0;
	this->P[IX_GRE + IX_GRE * NX] = this->opt.var_grad0;
	return this;
}

void vbs_netekf_free(struct vbs_netekf *this) {
	if (this == NULL)
		return;
	free(this);
}

/* Zero row and column `i` and set its variance -- how a state is (re)initialised. */
static void reset_state(struct vbs_netekf *this, int i, double value, double var) {
	int j;

	for (j = 0; j < NX; j++)
		this->P[i + j * NX] = this->P[j + i * NX] = 0.0;
	this->x[i] = value;
	this->P[i + i * NX] = var;
}

static int slot_alloc(struct vbs_netekf *this, int sat) {
	int s;

	if (sat <= 0 || sat > MAXSAT)
		return -1;
	if (this->slot_of_sat[sat - 1] > 0)
		return this->slot_of_sat[sat - 1] - 1;
	for (s = 0; s < VBS_NET_MAXSAT; s++)
		if (this->slot[s].sat == 0) {
			int b;

			memset(&this->slot[s], 0, sizeof this->slot[s]);
			this->slot[s].sat = sat;
			this->slot_of_sat[sat - 1] = s + 1;
			reset_state(this, ION(s), 0.0, this->opt.var_ion0);
			for (b = 0; b < VBS_NET_MAXBAND; b++)
				reset_state(this, AMB(s, b), 0.0, this->opt.var_amb0);
			return s;
		}
	return -1;
}

static void slot_release(struct vbs_netekf *this, int s) {
	int b;

	if (this->slot[s].sat > 0)
		this->slot_of_sat[this->slot[s].sat - 1] = 0;
	reset_state(this, ION(s), 0.0, 0.0);
	for (b = 0; b < VBS_NET_MAXBAND; b++)
		reset_state(this, AMB(s, b), 0.0, 0.0);
	memset(&this->slot[s], 0, sizeof this->slot[s]);
}

/*
 * Wavelength of a combination for one satellite, from the frequencies its slot recorded, or 0
 * if the satellite lacks a band or the sum degenerates.
 */
static double slot_comb_lam(const struct slot *sl, const struct vbs_amb_comb *comb,
	const int *bi) {
	double sumf = 0.0;
	int k;

	for (k = 0; k < comb->nband; k++) {
		double f = sl->freq[bi[k]];

		if (f == 0.0)
			return 0.0;
		sumf += comb->coef[k] * f;
	}
	return fabs(sumf) < 1.0 ? 0.0 : CLIGHT / sumf;
}

/* Index within the slot for a band digit, allocating one if there is room. -1 if full. */
static int band_index(struct slot *sl, int band) {
	int b;

	for (b = 0; b < sl->nband; b++)
		if (sl->band[b] == band)
			return b;
	if (sl->nband >= VBS_NET_MAXBAND)
		return -1;
	sl->band[sl->nband] = band;
	return sl->nband++;
}

/*
 * Sequential scalar Kalman update. h is sparse -- at most four non-zero partials -- so it is
 * passed as (index, value) pairs rather than a dense NX vector, which is the difference
 * between an O(NX) and an O(NX²) inner loop per observation.
 *
 * Returns 0 if the observation was rejected by the FDE gate, 1 if it was applied.
 */
static int kf_update(struct vbs_netekf *this, const int *ix, const double *hv, int nh,
	double z, double r, int gate, double *acc_sq, double *acc_norm, long *acc_n) {
	double Ph[NX], hPh = 0.0, innov = z, s, K;
	int i, j, k;

	for (k = 0; k < nh; k++)
		innov -= hv[k] * this->x[ix[k]];

	for (i = 0; i < NX; i++) {
		double sum = 0.0;

		for (k = 0; k < nh; k++)
			sum += this->P[i + ix[k] * NX] * hv[k];
		Ph[i] = sum;
	}
	for (k = 0; k < nh; k++)
		hPh += hv[k] * Ph[ix[k]];

	s = hPh + r;
	if (s <= 0.0 || !isfinite(s) || !isfinite(innov))
		return 0;
	/*
	 * FDE (§4.3). The gate is applied against the *predicted* innovation variance, before
	 * the observation touches the state, so a bad observation is judged by a state it has
	 * not yet corrupted. It is held off for the first few epochs because a state still at
	 * its prior has a huge predicted variance and would let anything through, and then a
	 * tight variance afterwards would reject everything.
	 */
	if (gate && fabs(innov) > this->opt.max_innov_sigma * sqrt(s))
		return 0;

	if (acc_sq != NULL) {
		*acc_sq += innov * innov;
		*acc_norm += innov * innov / s;
		(*acc_n)++;
	}

	K = 1.0 / s;
	for (i = 0; i < NX; i++)
		this->x[i] += Ph[i] * K * innov;
	/*
	 * Rank-1 update P -= K·Ph·Phᵀ, one column at a time. Updating only the upper triangle
	 * and mirroring costs half the multiplies but makes one of the two writes stride NX
	 * doubles -- 2.3 KB -- so every element is a cache miss; perf put 82 % of the 1 Hz emit
	 * in that loop. Walking full columns doubles the arithmetic and removes the strided
	 * access; the per-element expression is kept verbatim, so the result is bit-identical
	 * to the triangular form (Ph[i]·Ph[j] commutes exactly, and P is symmetric on entry).
	 * A column whose Ph[j] is zero is untouched, exactly as a zero Ph[i] skipped a row.
	 */
	for (j = 0; j < NX; j++) {
		double *col = this->P + j * NX, phj = Ph[j];

		if (phj == 0.0)
			continue;
		for (i = 0; i < NX; i++)
			col[i] -= Ph[i] * phj * K;
	}
	return 1;
}

/* Slot carrying a given band with usable phase and code, or -1. */
static int is_l2c(uint8_t code) {
	return code == CODE_L2S || code == CODE_L2L || code == CODE_L2X;
}

static int obs_slot_of_band(const obsd_t *o, int band, int prefer_l2c) {
	int j, first = -1;

	for (j = 0; j < NFREQ + NEXOBS; j++) {
		if (o->code[j] == CODE_NONE || vbs_signal_band(o->code[j]) != band ||
		    o->L[j] == 0.0 || o->P[j] == 0.0)
			continue;
		if (!prefer_l2c || band != 2 || satsys(o->sat, NULL) != SYS_GPS)
			return j;
		if (is_l2c(o->code[j]))
			return j;
		if (first < 0)
			first = j;
	}
	return first;
}

uint8_t vbs_netekf_band_code(const obsd_t *o, int band, int prefer_l2c) {
	int j = obs_slot_of_band(o, band, prefer_l2c);

	return j < 0 ? CODE_NONE : o->code[j];
}

static const obsd_t *find_sat(const obsd_t *obs, int n, int sat) {
	int i;

	for (i = 0; i < n; i++)
		if (obs[i].sat == sat)
			return &obs[i];
	return NULL;
}

static void time_update(struct vbs_netekf *this, double dt) {
	int i, s;

	/* Receiver clock difference is white noise: no receiver clock is predictable. */
	for (i = 0; i < NSYSIDX; i++)
		reset_state(this, IX_CLK + i, 0.0, this->opt.var_clk);

	if (dt <= 0.0)
		return;
	this->P[IX_ZTD + IX_ZTD * NX] += this->opt.q_ztd * dt;
	this->P[IX_GRN + IX_GRN * NX] += this->opt.q_grad * dt;
	this->P[IX_GRE + IX_GRE * NX] += this->opt.q_grad * dt;
	for (s = 0; s < VBS_NET_MAXSAT; s++) {
		if (this->slot[s].sat == 0)
			continue;
		/*
		 * Ionospheric process noise grows toward the horizon: a low satellite's slant
		 * path crosses far more of the ionosphere, so the same vertical change moves
		 * the slant delay by the obliquity factor. Using the wet mapping function as
		 * that factor is an approximation, but it is the right shape and it is already
		 * computed.
		 */
		double f = this->slot[s].mapfw > 1.0 ? this->slot[s].mapfw : 1.0;
		double q = this->opt.q_ion_gps > 0.0 && satsys(this->slot[s].sat, NULL) == SYS_GPS ?
		    this->opt.q_ion_gps : this->opt.q_ion;

		this->P[ION(s) + ION(s) * NX] += q * f * f * dt;
	}
}

int vbs_netekf_epoch(struct vbs_netekf *this, gtime_t time,
	const obsd_t *master, int nm, const double *rs, const double *dts,
	const obsd_t *aux, int na, const nav_t *nav,
	const long *arc_master, const long *arc_aux) {
	int i, s, b, k, nprep, used = 0, gate;
	double dt;

	if (this == NULL || master == NULL || aux == NULL ||
	    arc_master == NULL || arc_aux == NULL)
		return -1;

	/* Log the baseline's common signal set once it is known (§4.3). */
	if (this->sigset.n == 0)
		vbs_sigset_common(master, nm, aux, na, this->opt.sys_mask, &this->sigset);

	dt = this->have_last ? fabs(timediff(time, this->last_time)) : 0.0;
	time_update(this, dt);
	this->last_time = time;
	this->have_last = 1;
	this->stats.epochs++;
	gate = this->stats.epochs > this->opt.gate_after_epochs;

	for (s = 0; s < VBS_NET_MAXSAT; s++)
		this->slot[s].seen_this_epoch = 0;

	nprep = 0;
	for (i = 0; i < nm; i++) {
		const obsd_t *om = &master[i], *oa;
		double rm, ra, e[3], azel_m[2], azel_a[2], drho, mapfw_m, mapfw_a, mapfw;
		double dT_model, fref, rdot, hclk;
		struct slot *sl;
		int sat = om->sat, sys = satsys(sat, NULL), si = sysidx(sys), islot;

		if (si < 0 || !(sys & this->opt.sys_mask))
			continue;
		if ((oa = find_sat(aux, na, sat)) == NULL)
			continue;

		/* Geometry at both ends, from the same satellite position (§4.6 reuse). */
		if ((rm = geodist(rs + (sat - 1) * 6, this->master_xyz, e)) <= 0.0) {
			this->stats.obs_rejected[VBS_NET_REJ_GEOM]++;
			continue;
		}
		satazel(this->master_pos, e, azel_m);
		if (azel_m[1] < this->opt.elev_mask_deg * D2R) {
			this->stats.obs_rejected[VBS_NET_REJ_ELEV]++;
			continue;
		}
		if ((ra = geodist(rs + (sat - 1) * 6, this->aux_xyz, e)) <= 0.0) {
			this->stats.obs_rejected[VBS_NET_REJ_GEOM]++;
			continue;
		}
		satazel(this->aux_pos, e, azel_a);

		/*
		 * Receiver clock partial, and the reason it is not simply 1.
		 *
		 * An observation tagged at the nominal epoch was actually made when the
		 * receiver's own clock read that epoch, i.e. at true time t - dt_rcv. The
		 * geometry therefore has to be evaluated at t - dt_rcv, not at t, and getting
		 * that wrong costs rdot * dt_rcv. Two independent receivers sit at whatever
		 * offset their steering leaves, and -- measured on a 29 km F9P baseline -- they drift
		 * apart at ~7.4 us/s until the receiver resets itself by exactly one
		 * millisecond. At that jump the SD geometry error steps by rdot * 1 ms, up to
		 * +-0.8 m, the same on every band and different on every satellite. Before this
		 * term existed the filter rejected 39% of all phase observations from the jump
		 * onward, because the gate kept out the very observations that would have
		 * corrected the state.
		 *
		 * Rather than detect and repair the jump as a special case, fold it into the
		 * model: the clock state already carries c*dt in the observation, so the
		 * geometry effect is exactly the same state with an extra partial of -rdot/c.
		 * Drift and jumps are then both handled by the ordinary update, with no
		 * threshold to tune and nothing to detect.
		 *
		 * The master is taken as the time reference (dt_master := 0). That is not an
		 * approximation being swept away: the corrections are defined relative to the
		 * master, the VBS is synthesized from the master's observations, and so the VBS
		 * inherits the master's clock exactly as a physical base would. The convention
		 * only has to be self-consistent, and it is.
		 */
		/*
		 * d/dt of geodist() itself, term by term, rather than a range-rate formula
		 * borrowed from elsewhere: geodist returns |Xs - Xr| + OMGE(Xs_x Xr_y -
		 * Xs_y Xr_x)/c with `e` pointing receiver->satellite, so with a static receiver
		 * the derivative is vs.e plus the derivative of that Sagnac term. (The Sagnac
		 * rate is ~5 mm/s against a range rate of up to 800 m/s, so it changes nothing
		 * numerically -- it is here because a term that is easy to get backwards is
		 * worth deriving rather than recalling.)
		 */
		rdot = dot(rs + (sat - 1) * 6 + 3, e, 3)
		     + OMGE / CLIGHT * (rs[(sat - 1) * 6 + 3] * this->aux_xyz[1]
				      - rs[(sat - 1) * 6 + 4] * this->aux_xyz[0]);
		hclk = 1.0 - rdot / CLIGHT;

		/*
		 * The single-difference geometry. Not ra - rm: that would use the satellite
		 * position at the *master's* transmission time for both ends, and over a 150 km
		 * baseline the ~500 us light-time difference moves the satellite ~2 m.
		 * vbs_synth_drho() re-converges the transmission time for the second station,
		 * and using it here is what keeps the filter and the synthesis consistent.
		 */
		drho = vbs_synth_drho(rs + (sat - 1) * 6, rs + (sat - 1) * 6 + 3, this->master_xyz,
		    this->aux_xyz);
		if (!isfinite(drho)) {
			this->stats.obs_rejected[VBS_NET_REJ_GEOM]++;
			continue;
		}

		/*
		 * A priori troposphere: the model difference is removed, only the residual
		 * estimated.
		 *
		 * Both terms are checked before differencing. tropmodel() returns 0.0 for an
		 * out-of-range input rather than signalling failure, so a difference taken
		 * against it silently becomes minus the other station's whole slant delay --
		 * metres, applied as though it were an atmosphere model. A surveyed base is
		 * unlikely to sit below -100 m, which is exactly why this would have sat here
		 * undetected; the same construct in vbs_interp was reached in practice and
		 * measured at -4.79 m (D36).
		 */
		{
			double t_aux = tropmodel(time, this->aux_pos, azel_a, this->opt.humi);
			double t_mst = tropmodel(time, this->master_pos, azel_m, this->opt.humi);

			if (t_aux <= 0.0 || t_mst <= 0.0) {
				this->stats.obs_rejected[VBS_NET_REJ_GEOM]++;
				continue;
			}
			dT_model = t_aux - t_mst;
		}
		tropmapf(time, this->master_pos, azel_m, &mapfw_m);
		tropmapf(time, this->aux_pos, azel_a, &mapfw_a);
		mapfw = 0.5 * (mapfw_m + mapfw_a);

		if ((islot = slot_alloc(this, sat)) < 0) {
			this->stats.obs_rejected[VBS_NET_REJ_NOSLOT]++;
			continue;
		}
		sl = &this->slot[islot];
		sl->seen_this_epoch = 1;
		sl->epochs_tracked++;
		sl->last = time;
		sl->el = azel_m[1];
		sl->az = azel_m[0];
		/*
		 * Chen and Herring's gradient mapping, m_g(e) = 1/(sin e · tan e + 0.0032): the
		 * standard partial for a horizontal gradient, and the same one the IGS troposphere
		 * products are estimated with, so the states mean the same thing as theirs.
		 */
		sl->mapg = 1.0 / (sin(azel_m[1]) * tan(azel_m[1]) + 0.0032);
		/*
		 * The mapping reaches 30 at 10 degrees, so an unconstrained gradient state can
		 * absorb a metre of low-elevation signal that belongs to the ionosphere or the
		 * ambiguities -- measured, that took M30 from 0.78 % misclosure to 20.8 %. Capping
		 * it keeps the gradient a troposphere term rather than a sink for everything that
		 * happens to be low in the sky.
		 */
		if (sl->mapg > this->opt.max_mapg)
			sl->mapg = this->opt.max_mapg;
		sl->mapfw = mapfw;
		sl->dT_model = dT_model;

		/*
		 * Arc boundaries come from the stations, not from here. A change in either
		 * station's arc id since the previous epoch means the phase is no longer
		 * continuous at that station, and the ambiguities on this baseline are reset --
		 * in the same epoch as every other filter that includes the station, because
		 * they all read the same id.
		 */
		if (sl->arc_master >= 0 &&
		    (arc_master[sat - 1] != sl->arc_master || arc_aux[sat - 1] != sl->arc_aux)) {
			if (arc_master[sat - 1] != sl->arc_master)
				this->stats.slips_master++;
			else
				this->stats.slips_aux++;
			for (b = 0; b < VBS_NET_MAXBAND; b++)
				reset_state(this, AMB(islot, b), 0.0, this->opt.var_amb0);
		}
		sl->arc_master = arc_master[sat - 1];
		sl->arc_aux = arc_aux[sat - 1];

		if ((fref = vbs_synth_ref_freq(sys)) == 0.0)
			continue;

		for (b = 0; b < this->sigset.n; b++) {
			const struct vbs_signal *sig = &this->sigset.sig[b];
			struct prep *p;
			int jm, ja, bi;
			double f, sinel;

			if (sig->sys != sys)
				continue;
			if (this->opt.drop_band5 && sig->band == 5)
				continue;
			if (this->opt.drop_band2 && sig->band == 2 && sys == SYS_GPS)
				continue;
			if (nprep >= (int)(sizeof this->prep / sizeof this->prep[0]))
				break;
			if ((jm = obs_slot_of_band(om, sig->band, this->opt.prefer_l2c)) < 0)
				continue;
			if ((ja = obs_slot_of_band(oa, sig->band, this->opt.prefer_l2c)) < 0)
				continue;
			if ((f = sat2freq(sat, om->code[jm], nav)) == 0.0)
				continue;
			if (this->opt.min_snr_dbhz > 0.0 &&
			    ((om->SNR[jm] > 0 && om->SNR[jm] * SNR_UNIT < this->opt.min_snr_dbhz) ||
			     (oa->SNR[ja] > 0 && oa->SNR[ja] * SNR_UNIT < this->opt.min_snr_dbhz))) {
				this->stats.obs_rejected[VBS_NET_REJ_SNR]++;
				continue;
			}
			if ((bi = band_index(sl, sig->band)) < 0)
				continue;

			sinel = sin(azel_m[1]);
			if (sinel < 0.05)
				sinel = 0.05;

			sl->freq[bi] = f;
			p = &this->prep[nprep++];
			p->slot = islot;
			p->si = si;
			p->bi = bi;
			p->lam = CLIGHT / f;
			p->mu = (fref / f) * (fref / f);
			p->mapfw = mapfw;
		p->mgn = sl->mapg * cos(sl->az);
		p->mge = sl->mapg * sin(sl->az);
			p->hclk = hclk;
			p->zP = (oa->P[ja] - om->P[jm]) - drho - dT_model;
			p->zL = (oa->L[ja] - om->L[jm]) * p->lam - drho - dT_model;
			{
				double sc = sys == SYS_GPS && this->opt.sig_code_gps > 0.0 ?
				    this->opt.sig_code_gps : this->opt.sig_code;

				p->var_p = sc * sc * (1.0 + 1.0 / (sinel * sinel));
			}
			p->var_l = this->opt.sig_phase * this->opt.sig_phase *
			    (1.0 + 1.0 / (sinel * sinel));
		}
	}

	/*
	 * Seed the receiver clock difference from this epoch's own code observations, before
	 * any update and before the FDE gate can act.
	 *
	 * Without this the filter cannot start. The clock difference between two independent
	 * receivers is whatever their steering happens to leave -- commonly up to a millisecond,
	 * which is 300 km -- while every other state in the vector is a centimetre-scale
	 * quantity. Resetting the clock to zero with any sane prior variance makes the first
	 * innovation astronomically large, the gate rejects it, and because the gate rejected it
	 * the state never learns the clock: the filter starves itself. Measured on a 29 km F9P baseline
	 * before this was added: 177 386 observations rejected against 2 228 used.
	 *
	 * The median is the right estimator here rather than the mean: it costs one sort and it
	 * is immune to the handful of gross code outliers that a real base always has, which is
	 * the same robustness §4.3 asks for from the gate itself.
	 */
	for (i = 0; i < NSYSIDX; i++) {
		double v[sizeof this->prep / sizeof this->prep[0]];
		int nv = 0;

		for (k = 0; k < nprep; k++)
			if (this->prep[k].si == i)
				v[nv++] = this->prep[k].zP;
		if (nv == 0)
			continue;
		qsort(v, nv, sizeof *v, cmp_double);
		reset_state(this, IX_CLK + i, v[nv / 2], this->opt.var_clk);
	}

	/*
	 * Code first, then phase. The clock seed above is only as good as the median of noisy
	 * code, so its variance is still metres when the first update runs; running every code
	 * observation first brings the clock down to centimetres, and only then is the gate
	 * tight enough to mean anything on phase. Interleaving them would give the first few
	 * phase observations a metre-wide gate -- which is exactly where an undetected slip
	 * would slip through.
	 */
	for (k = 0; k < nprep; k++) {
		const struct prep *p = &this->prep[k];
		int ix[5] = {IX_CLK + p->si, IX_ZTD, ION(p->slot), IX_GRN, IX_GRE};
		double hv[5] = {p->hclk, p->mapfw, +p->mu, p->mgn, p->mge};
		int nh = this->opt.est_gradients ? 5 : 3;

		if (kf_update(this, ix, hv, nh, p->zP, p->var_p, gate,
		    &this->sq_code, &this->norm_code, &this->n_code))
			used++;
		else
			this->stats.obs_rejected[VBS_NET_REJ_CODE]++;
	}
	for (k = 0; k < nprep; k++) {
		const struct prep *p = &this->prep[k];
		int ix[6] = {IX_CLK + p->si, IX_ZTD, ION(p->slot), AMB(p->slot, p->bi),
		    IX_GRN, IX_GRE};
		double hv[6] = {p->hclk, p->mapfw, -p->mu, p->lam, p->mgn, p->mge};
		int nh = this->opt.est_gradients ? 6 : 4;

		if (kf_update(this, ix, hv, nh, p->zL, p->var_l, gate,
		    &this->sq_phase, &this->norm_phase, &this->n_phase)) {
			used++;
			this->slot[p->slot].phase_accepted = 1;
		} else
			this->stats.obs_rejected[VBS_NET_REJ_PHASE]++;
	}
	/*
	 * A satellite whose phase is gated out epoch after epoch is not being filtered: its
	 * ambiguities sit where the last accepted epoch left them, and any fix on them is an
	 * assertion nobody is checking. That is what an undetected phase jump looks like from
	 * in here (D49: T30 held such fixes for hours). Treat a streak like a slip: reset the
	 * ambiguities so the next epoch re-initialises from the phase as it now is.
	 */
	for (s = 0; s < VBS_NET_MAXSAT; s++) {
		struct slot *sl = &this->slot[s];

		if (sl->sat == 0 || !sl->seen_this_epoch)
			continue;
		if (sl->phase_accepted) {
			sl->phase_rejects = 0;
		} else if (++sl->phase_rejects >= this->opt.reject_streak_reset) {
			for (b = 0; b < VBS_NET_MAXBAND; b++)
				reset_state(this, AMB(s, b), 0.0, this->opt.var_amb0);
			sl->phase_rejects = 0;
			this->stats.reject_resets++;
		}
		sl->phase_accepted = 0;
	}

	/* Release slots for satellites that have set, so the state stays compact. */
	for (s = 0; s < VBS_NET_MAXSAT; s++) {
		if (this->slot[s].sat == 0 || this->slot[s].seen_this_epoch)
			continue;
		if (fabs(timediff(time, this->slot[s].last)) > 120.0)
			slot_release(this, s);
	}

	this->stats.obs_used += used;
	this->stats.nsat_active = 0;
	this->stats.nband_max = 0;
	for (s = 0; s < VBS_NET_MAXSAT; s++)
		if (this->slot[s].sat != 0) {
			this->stats.nsat_active++;
			if (this->slot[s].nband > this->stats.nband_max)
				this->stats.nband_max = this->slot[s].nband;
		}
	this->nprep = nprep;		/* kept for vbs_netekf_phase_sd(), D66 */
	this->stats.ztd = this->x[IX_ZTD];
	this->stats.ztd_sigma = sqrt(this->P[IX_ZTD + IX_ZTD * NX]);
	this->stats.innov_rms_code = this->n_code ? sqrt(this->sq_code / this->n_code) : 0.0;
	this->stats.innov_rms_phase = this->n_phase ? sqrt(this->sq_phase / this->n_phase) : 0.0;
	this->stats.innov_norm_code = this->n_code ? sqrt(this->norm_code / this->n_code) : 0.0;
	this->stats.innov_norm_phase = this->n_phase ? sqrt(this->norm_phase / this->n_phase) : 0.0;
	return used;
}

int vbs_netekf_if_average(struct vbs_netekf *this, int sat, int refsat, int band, double value,
	double *mean, int *count) {
	int s, free_i = -1, i;
	const struct slot *sl;

	if (this == NULL || (s = this->slot_of_sat[sat - 1] - 1) < 0)
		return 0;
	sl = &this->slot[s];
	for (i = 0; i < VBS_NET_MAXSAT; i++) {
		struct if_acc *a = &this->ifacc[i];

		if (a->sat == 0) {
			if (free_i < 0)
				free_i = i;
			continue;
		}
		if (a->sat != sat)
			continue;
		/* A different reference, band or arc is a different quantity: start again. */
		if (a->refsat != refsat || a->band != band ||
		    a->arc_m != sl->arc_master || a->arc_a != sl->arc_aux) {
			a->refsat = refsat;
			a->band = band;
			a->arc_m = sl->arc_master;
			a->arc_a = sl->arc_aux;
			a->sum = 0.0;
			a->n = 0;
		}
		a->sum += value;
		a->n++;
		*mean = a->sum / a->n;
		*count = a->n;
		return 1;
	}
	if (free_i < 0)
		return 0;
	{
		struct if_acc *a = &this->ifacc[free_i];

		a->sat = sat;
		a->refsat = refsat;
		a->band = band;
		a->arc_m = sl->arc_master;
		a->arc_a = sl->arc_aux;
		a->sum = value;
		a->n = 1;
		*mean = value;
		*count = 1;
	}
	return 1;
}

int vbs_netekf_slot_mapfw(const struct vbs_netekf *this, int sat, double *mapfw, double *el) {
	int s;

	if (this == NULL || sat <= 0 || sat > MAXSAT ||
	    (s = this->slot_of_sat[sat - 1] - 1) < 0 || this->slot[s].sat == 0)
		return 0;
	if (mapfw != NULL)
		*mapfw = this->slot[s].mapfw;
	if (el != NULL)
		*el = this->slot[s].el;
	return 1;
}

double vbs_netekf_ztd(const struct vbs_netekf *this) {
	return this != NULL ? this->x[IX_ZTD] : 0.0;
}

int vbs_netekf_phase_sd(const struct vbs_netekf *this, int sat, int band, double *zL,
	double *mu, double *lam) {
	int k;

	if (this == NULL)
		return 0;
	for (k = 0; k < this->nprep; k++) {
		const struct prep *p = &this->prep[k];
		const struct slot *sl = &this->slot[p->slot];

		if (sl->sat != sat || p->bi < 0 || p->bi >= sl->nband || sl->band[p->bi] != band)
			continue;
		if (zL != NULL)
			*zL = p->zL;
		if (mu != NULL)
			*mu = p->mu;
		if (lam != NULL)
			*lam = p->lam;
		return 1;
	}
	return 0;
}

double vbs_netekf_sig_phase(const struct vbs_netekf *this) {
	return this != NULL ? this->opt.sig_phase : 0.0;
}

int vbs_netekf_corrections(const struct vbs_netekf *this, struct vbs_corr *out, int maxout,
	double max_sigma) {
	int s, n = 0;

	if (this == NULL || out == NULL)
		return 0;
	/*
	 * A threshold at or above the prior cannot distinguish a state the data has informed
	 * from one that is still sitting at its initial variance, so it is not a gate at all --
	 * it would wave through a satellite that has just risen and mark its untouched prior as
	 * a valid correction. Measured on the closed-loop test: at max_sigma = 1.0 m against a
	 * 0.5 m prior the worst recovered ΔI was 61.8 mm out; at 0.05 m it was 6.8 mm; at
	 * 0.02 m it was exact. Refuse the meaningless case outright rather than silently
	 * tightening it, so a caller that asks for nonsense is told, not quietly corrected.
	 */
	if (!(max_sigma < sqrt(this->opt.var_ion0)))
		return 0;
	for (s = 0; s < VBS_NET_MAXSAT && n < maxout; s++) {
		const struct slot *sl = &this->slot[s];
		double var_i, var_t;

		if (sl->sat == 0 || !sl->seen_this_epoch)
			continue;
		var_i = this->P[ION(s) + ION(s) * NX];
		var_t = sl->mapfw * sl->mapfw * this->P[IX_ZTD + IX_ZTD * NX];

		out[n].sat = sl->sat;
		out[n].dI = this->x[ION(s)];
		/*
		 * Non-dispersive: the estimated *residual* only, not the a priori model
		 * difference.
		 *
		 * §4.5 requires the model term between master and VBS to be computed from the
		 * model and only the residual interpolated. Returning model+residual here made
		 * that impossible downstream: the plane fit interpolated the auxiliaries' model
		 * differences to the VBS, and `vbs_interp` then added its own computed height
		 * term on top, counting the same physics twice. The paired decomposition caught
		 * it as a *negative* alpha on the non-dispersive term (-1.45 GAL, -2.15 GPS at
		 * 29 km) -- the correction anti-correlated with the error it was meant to
		 * remove, which noise alone cannot produce.
		 */
		/*
		 * The gradient belongs in the correction, not only in the filter. Estimating it
		 * sharpens the ambiguities -- measured, it cut the rover's false fixes from 1.20 %
		 * to 0.40 % -- but if it is left out of what is emitted, the synthesized
		 * observations still carry the tilt, and the rover sees it as a horizontal bias
		 * that grows with the master-to-VBS displacement.
		 */
		out[n].dT = sl->mapfw * this->x[IX_ZTD];
		if (this->opt.est_gradients)
			out[n].dT += sl->mapg * (cos(sl->az) * this->x[IX_GRN] +
			    sin(sl->az) * this->x[IX_GRE]);
		out[n].valid = (var_i >= 0.0 && var_t >= 0.0 &&
		    sqrt(var_i) <= max_sigma && sqrt(var_t) <= max_sigma);
		n++;
	}
	return n;
}

int vbs_netekf_correction_sigmas(const struct vbs_netekf *this, const struct vbs_corr *corr,
	int ncorr, double *sd_I, double *sd_T) {
	int k, n = 0;

	for (k = 0; k < ncorr; k++) {
		int sat = corr[k].sat, s;

		sd_I[k] = sd_T[k] = 0.0;
		if (sat <= 0 || sat > MAXSAT)
			continue;
		if ((s = this->slot_of_sat[sat - 1] - 1) < 0)
			continue;
		sd_I[k] = sqrt(this->P[ION(s) + ION(s) * NX]);
		sd_T[k] = this->slot[s].mapfw * sqrt(this->P[IX_ZTD + IX_ZTD * NX]);
		n++;
	}
	return n;
}

void vbs_netekf_get_stats(const struct vbs_netekf *this, struct vbs_netekf_stats *out) {
	*out = this->stats;
}

const struct vbs_sigset *vbs_netekf_sigset(const struct vbs_netekf *this) {
	return &this->sigset;
}

/* --- Ambiguity resolution support (§4.3) ---------------------------------------------- */

/* Slot for a satellite, or -1. */
static int slot_of(const struct vbs_netekf *this, int sat) {
	if (sat <= 0 || sat > MAXSAT)
		return -1;
	return this->slot_of_sat[sat - 1] - 1;
}

/*
 * Band indices within a slot for every band the combination needs, or 0 if the satellite
 * does not carry them all. A combination is only formable where both ends tracked every band
 * it uses, which is exactly the per-baseline signal-set question §4.3 asks to be answered
 * from the data.
 */
static int comb_indices(const struct slot *sl, const struct vbs_amb_comb *comb, int *bi) {
	int k, b;

	for (k = 0; k < comb->nband; k++) {
		for (b = 0; b < sl->nband; b++)
			if (sl->band[b] == comb->band[k])
				break;
		if (b >= sl->nband)
			return 0;
		bi[k] = b;
	}
	return 1;
}


/*
 * The satellite of `sys` that carries the most of `bands`, highest elevation breaking ties, or
 * 0. Used to hold ONE reference satellite across every step of a constellation's cascade.
 *
 * Why that matters. `vbs_netekf_dd_set()` picks the highest-elevation satellite carrying the
 * bands of the combination in front of it, so each step gets its own reference. The wide lane
 * and the narrow lane derived from it do share one (resolve_narrow_if inherits w->refsat), but
 * different wide lanes do not, and `vbs_ambres_ion_from_fixes()` can only use satellites fixed
 * against a single common reference -- mixing them gives each satellite a different unknown
 * ionosphere level, which D66 measured taking M30 from 98.4 % fixed to 13.2 % with a 12.6 cm
 * vertical bias. The result was that 64 % of corrections had no usable wide-lane/narrow-lane
 * pair and the feature replaced nothing at all (D131).
 *
 * Carrying the *most* of the bands rather than all of them is deliberate: requiring all would
 * often return nothing (only about eighteen GPS satellites transmit L5), and a combination
 * whose bands the held reference lacks simply falls back to its own choice, which is no worse
 * than the behaviour this replaces.
 */
int vbs_netekf_best_ref(const struct vbs_netekf *this, int sys, const int *bands, int nband) {
	int s, k, b, best = 0, bestn = -1;
	double bestel = -1.0;

	if (this == NULL || bands == NULL || nband <= 0)
		return 0;
	for (s = 0; s < VBS_NET_MAXSAT; s++) {
		const struct slot *sl = &this->slot[s];
		int have = 0;

		if (sl->sat == 0 || !sl->seen_this_epoch || satsys(sl->sat, NULL) != sys)
			continue;
		for (k = 0; k < nband; k++)
			for (b = 0; b < sl->nband; b++)
				if (sl->band[b] == bands[k]) {
					have++;
					break;
				}
		if (have > bestn || (have == bestn && sl->el > bestel)) {
			bestn = have;
			bestel = sl->el;
			best = sl->sat;
		}
	}
	return bestn > 0 ? best : 0;
}

int vbs_netekf_dd_set(const struct vbs_netekf *this, int sys, const struct vbs_amb_comb *comb,
	struct vbs_dd_amb *out, int maxout, double *Q) {
	return vbs_netekf_dd_set_ref(this, sys, comb, 0, out, maxout, Q);
}

/*
 * As above, but `prefer_sat` (0 for none) is used as the double-difference reference when it is
 * a valid candidate for this combination -- seen this epoch, carrying every band, with a
 * frequency for each. Otherwise the highest-elevation candidate is used exactly as before, so
 * passing 0 reproduces the previous behaviour byte for byte.
 */
int vbs_netekf_dd_set_ref(const struct vbs_netekf *this, int sys, const struct vbs_amb_comb *comb,
	int prefer_sat, struct vbs_dd_amb *out, int maxout, double *Q) {
	int cand[VBS_NET_MAXSAT], bi[VBS_NET_MAXBAND], refbi[VBS_NET_MAXBAND];
	int ncand = 0, refslot = -1, preferslot = -1, s, i, j, k, l, n = 0;
	double bestel = -1.0, preferel = -1.0;

	if (this == NULL || comb == NULL || comb->nband <= 0)
		return 0;

	/*
	 * The combination's wavelength is 1 / sum(coef_i / lambda_i), i.e. c / sum(coef_i f_i),
	 * and it is a property of the *satellite*: for GLONASS the frequencies differ per
	 * satellite, so this is computed per candidate below from the frequencies the slot
	 * recorded. It used to be computed once from the constellation's nominal frequencies,
	 * which is what confined the filter to CDMA systems.
	 *
	 * A CDMA constellation still gets one value for all of its satellites, arrived at the
	 * same way, so nothing changes for GPS, Galileo or BeiDou.
	 */
	if (vbs_synth_ref_freq(sys) == 0.0 && sys != SYS_GLO)
		return 0;

	/* Candidates: satellites of this system carrying every band, seen this epoch. */
	for (s = 0; s < VBS_NET_MAXSAT; s++) {
		const struct slot *sl = &this->slot[s];

		if (sl->sat == 0 || !sl->seen_this_epoch)
			continue;
		if (satsys(sl->sat, NULL) != sys)
			continue;
		if (!comb_indices(sl, comb, bi))
			continue;
		if (slot_comb_lam(sl, comb, bi) == 0.0)
			continue;	/* no frequency for one of its bands */
		cand[ncand++] = s;
		if (prefer_sat != 0 && sl->sat == prefer_sat) {
			preferslot = s;
			preferel = sl->el;
		}
		if (sl->el > bestel) {
			bestel = sl->el;
			refslot = s;
		}
	}
	if (preferslot >= 0) {
		refslot = preferslot;
		bestel = preferel;	/* `el` is the pair's lower elevation, so it must follow
					 * the reference actually used, not the discarded one */
	}
	if (refslot < 0 || ncand < 2)
		return 0;
	comb_indices(&this->slot[refslot], comb, refbi);

	/*
	 * DD value and covariance. Each DD is c·(N_sat − N_ref) over the combination, so with
	 * G the (n × NX) map from state to DD, Q = G P Gᵀ. G has at most 2·nband non-zeros per
	 * row, so the product is formed directly rather than by building G.
	 */
	for (i = 0; i < ncand && n < maxout; i++) {
		const struct slot *sl = &this->slot[cand[i]];
		double v = 0.0;

		if (cand[i] == refslot)
			continue;
		comb_indices(sl, comb, bi);
		for (k = 0; k < comb->nband; k++)
			v += comb->coef[k] * (this->x[AMB(cand[i], bi[k])]
					    - this->x[AMB(refslot, refbi[k])]);
		out[n].sat = sl->sat;
		out[n].refsat = this->slot[refslot].sat;
		out[n].sys = sys;
		out[n].value = v;
		out[n].lam = slot_comb_lam(sl, comb, bi);
		out[n].el = sl->el < bestel ? sl->el : bestel;
		cand[n] = cand[i];		/* compact, keeping cand[] parallel to out[] */
		n++;
	}

	for (i = 0; i < n; i++) {
		int bi_i[VBS_NET_MAXBAND];

		comb_indices(&this->slot[cand[i]], comb, bi_i);
		for (j = 0; j < n; j++) {
			int bi_j[VBS_NET_MAXBAND];
			double q = 0.0;

			comb_indices(&this->slot[cand[j]], comb, bi_j);
			for (k = 0; k < comb->nband; k++)
				for (l = 0; l < comb->nband; l++) {
					double ci = comb->coef[k], cj = comb->coef[l];
					int ai = AMB(cand[i], bi_i[k]), ar = AMB(refslot, refbi[k]);
					int aj = AMB(cand[j], bi_j[l]), br = AMB(refslot, refbi[l]);

					q += ci * cj * (this->P[ai + aj * NX] - this->P[ai + br * NX]
						      - this->P[ar + aj * NX] + this->P[ar + br * NX]);
				}
			Q[i + j * n] = q;
		}
	}
	return n;
}

int vbs_netekf_constrain_dd(struct vbs_netekf *this, const struct vbs_dd_amb *a,
	const struct vbs_amb_comb *comb, double integer, double var) {
	int ix[2 * VBS_NET_MAXBAND], bi[VBS_NET_MAXBAND], refbi[VBS_NET_MAXBAND];
	double hv[2 * VBS_NET_MAXBAND];
	int s = slot_of(this, a->sat), r = slot_of(this, a->refsat), k, nh = 0;

	if (s < 0 || r < 0)
		return 0;
	if (!comb_indices(&this->slot[s], comb, bi) ||
	    !comb_indices(&this->slot[r], comb, refbi))
		return 0;
	for (k = 0; k < comb->nband; k++) {
		ix[nh] = AMB(s, bi[k]);
		hv[nh++] = comb->coef[k];
		ix[nh] = AMB(r, refbi[k]);
		hv[nh++] = -comb->coef[k];
	}
	/* Never gated: a held integer is an assertion about the world, not an observation. */
	return kf_update(this, ix, hv, nh, integer, var, 0, NULL, NULL, NULL);
}

int vbs_netekf_dd_ion_sigma(const struct vbs_netekf *this, int sys, double *out, int maxout) {
	int s, ref = -1, n = 0;
	double bestel = -1.0;

	for (s = 0; s < VBS_NET_MAXSAT; s++) {
		const struct slot *sl = &this->slot[s];

		if (sl->sat == 0 || !sl->seen_this_epoch || satsys(sl->sat, NULL) != sys)
			continue;
		if (sl->nband < 2)	/* single-band: ΔI is not separable at all (D22) */
			continue;
		if (sl->el > bestel) {
			bestel = sl->el;
			ref = s;
		}
	}
	if (ref < 0)
		return 0;
	for (s = 0; s < VBS_NET_MAXSAT && n < maxout; s++) {
		const struct slot *sl = &this->slot[s];
		double v;

		if (s == ref || sl->sat == 0 || !sl->seen_this_epoch ||
		    satsys(sl->sat, NULL) != sys || sl->nband < 2)
			continue;
		v = this->P[ION(s) + ION(s) * NX] + this->P[ION(ref) + ION(ref) * NX]
		  - 2.0 * this->P[ION(s) + ION(ref) * NX];
		out[n++] = v > 0.0 ? sqrt(v) : 0.0;
	}
	return n;
}

int vbs_netekf_satinfo(const struct vbs_netekf *this, struct vbs_netekf_satinfo *out,
	int maxout) {
	int s, b, n = 0;

	for (s = 0; s < VBS_NET_MAXSAT && n < maxout; s++) {
		const struct slot *sl = &this->slot[s];
		double worst = 0.0;

		if (sl->sat == 0 || !sl->seen_this_epoch)
			continue;
		for (b = 0; b < sl->nband; b++) {
			double v = this->P[AMB(s, b) + AMB(s, b) * NX];

			if (v > worst)
				worst = v;
		}
		out[n].sat = sl->sat;
		out[n].sys = satsys(sl->sat, NULL);
		out[n].nband = sl->nband;
		out[n].epochs_tracked = sl->epochs_tracked;
		out[n].el_deg = sl->el * R2D;
		out[n].sd_ion = sqrt(this->P[ION(s) + ION(s) * NX]);
		out[n].amb_sigma_max = sqrt(worst);
		n++;
	}
	return n;
}

void vbs_netekf_unfix_sat(struct vbs_netekf *this, int sat) {
	int s = slot_of(this, sat), b;

	if (s < 0)
		return;
	for (b = 0; b < VBS_NET_MAXBAND; b++)
		reset_state(this, AMB(s, b), 0.0, this->opt.var_amb0);
}

void vbs_satpos_by_sat(const obsd_t *obs, int n, const double *rs_rec, double *rs_sat) {
	int i;

	memset(rs_sat, 0, sizeof(double) * MAXSAT * 6);
	for (i = 0; i < n; i++)
		if (obs[i].sat >= 1 && obs[i].sat <= MAXSAT)
			memcpy(rs_sat + (obs[i].sat - 1) * 6, rs_rec + i * 6, sizeof(double) * 6);
}

int vbs_netekf_station_arcs(const struct vbs_netekf *this, int sat, long *arc_master,
	long *arc_aux) {
	int s = slot_of(this, sat);

	if (s < 0)
		return 0;
	*arc_master = this->slot[s].arc_master;
	*arc_aux = this->slot[s].arc_aux;
	return 1;
}

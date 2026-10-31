#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_ambres.h"
#include "vbs/vbs_basemeta.h"
#include "vbs/vbs_interp.h"
#include "vbs/vbs_lock.h"
#include "vbs/vbs_netekf.h"
#include "vbs/vbs_preproc.h"
#include "vbs/vbs_rtcm_out.h"
#include "vbs/vbs_session.h"
#include "vbs/vbs_signal.h"
#include "vbs/vbs_synth.h"

/*
 * Everything below tracks tools/vbs_pipeline.c's epoch loop, minus its diagnostics. Where the
 * pipeline offers a flag, the value used here is the one the shipped configuration passes; the
 * comment names the flag so the two can be compared line by line. tests/test_session.c is the
 * proof that they still agree.
 */

#define NCOMB		16			/* as tools/vbs_pipeline.c */

static int sysidx(int sat) {
	switch (satsys(sat, NULL)) {
	case SYS_GPS: return 0;
	case SYS_GLO: return 1;
	case SYS_GAL: return 2;
	case SYS_CMP: return 3;
	default: return -1;
	}
}
#define CORR_GATE_M	0.05			/* vbs_netekf_corrections() sigma gate, as the pipeline */

struct slot {
	gtime_t t;
	obsd_t obs[MAXOBS];
	int n;
	double arrived;
	int used;			/* consumed, or superseded: free for reuse */
	int filled;
};

struct station {
	struct vbs_session_station def;
	struct slot ring[VBS_SESSION_RING];
	struct vbs_preproc *pp;
	long arc[MAXSAT];
	/* auxiliaries only */
	struct vbs_netekf *ekf;
	struct vbs_amb_fix fix[256];
	int nfix;
	struct vbs_corr corr[MAXOBS];
	double sd_I[MAXOBS], sd_T[MAXOBS];
	int ncorr;
	double last_fix[MAXSAT][NCOMB];
};

struct vbs_session {
	struct vbs_session_opt opt;
	struct station sta[VBS_SESSION_MAXSTA];
	int nsta, master;
	double vbs_xyz[3];
	nav_t *nav;

	struct vbs_netekf_opt nopt;
	struct vbs_ambres_opt aopt;
	struct vbs_interp_opt iopt;
	struct vbs_basemeta_opt bopt;
	struct vbs_lock *lock;
	struct vbs_rtcm_out *rtcm;
	int glo_fcn_installed;

	long last_arc_m[MAXSAT];
	gtime_t last_processed;		/* master epoch time last run (or attempted) */
	int have_processed;
	gtime_t last_pos_sent;
	int pos_sent_once;
	double last_keepalive_unix;
	long epochs_out_at_keepalive;
	int converged;
	double last_vbs_epoch_unix;

	unsigned char out[VBS_SESSION_OUTBUF];
	int nout;

	obsd_t last_synth[MAXOBS];
	int nlast_synth;
	gtime_t last_synth_time;

	struct vbs_session_stats stats;
	struct vbs_ambres_stats ar;
	long slip_by_sat[MAXSAT];
	double last_el[MAXSAT];
};

void vbs_session_default_opt(struct vbs_session_opt *opt) {
	memset(opt, 0, sizeof *opt);
	opt->sys_mask = SYS_GPS | SYS_GAL | SYS_CMP;
	opt->dispersive = 1;
	opt->prefer_l2c = 0;
	opt->interval_s = 1.0;
	/* the shipped flag set: --out-elmask 25 --if-narrow-lane --if-sys 0.020 --if-max-frac 0.12
	 * --if-min-epochs 120 --if-max-sigma 0.08 --q-ion 1e-4, with --ar-corr-time 30 from the
	 * harness (--net-corr-time 30) */
	opt->out_elmask_deg = 25.0;
	opt->ar_corr_time_s = 30.0;
	opt->net_ratio = 0.0;
	opt->if_narrow_lane = 1;
	opt->if_min_epochs = 120;
	opt->if_max_sigma = 0.08;
	opt->if_sys_m = 0.020;
	opt->if_max_frac = 0.12;
	opt->q_ion = 1e-4;
	opt->staid = 1;
	opt->pos_period_s = 10.0;
	opt->max_wait_s = 0.8;
	opt->reconverge_s = 30.0;
}

struct vbs_session *vbs_session_new(const double *vbs_xyz, const struct vbs_session_station *sta,
    int nsta, int master, const struct vbs_session_opt *opt, nav_t *nav) {
	struct vbs_session *this;
	struct vbs_lock_opt lopt;
	int i, a, c;

	if (nsta < 2 || nsta > VBS_SESSION_MAXSTA || master < 0 || master >= nsta || nav == NULL)
		return NULL;
	this = (struct vbs_session *)calloc(1, sizeof *this);
	if (this == NULL)
		return NULL;
	this->opt = *opt;
	this->nsta = nsta;
	this->master = master;
	this->nav = nav;
	memcpy(this->vbs_xyz, vbs_xyz, sizeof this->vbs_xyz);

	/* options, mapped exactly as tools/vbs_pipeline.c maps its flags */
	vbs_interp_default_opt(&this->iopt);
	if (opt->interp_min_aux > 0)
		this->iopt.min_aux = opt->interp_min_aux;
	vbs_netekf_default_opt(&this->nopt);
	this->nopt.min_snr_dbhz = 0.0;
	this->nopt.prefer_l2c = this->opt.prefer_l2c;
	if (opt->q_ion > 0.0)
		this->nopt.q_ion = opt->q_ion;
	this->nopt.drop_band5 = 0;
	this->nopt.drop_band2 = 0;
	this->nopt.sys_mask = opt->sys_mask;
	vbs_ambres_default_opt(&this->aopt);
	if (opt->ar_corr_time_s > 0.0 && opt->interval_s > 0.0 && opt->interval_s < opt->ar_corr_time_s)
		this->aopt.q_inflate = opt->ar_corr_time_s / opt->interval_s;
	if (opt->net_ratio > 0.0)
		this->aopt.min_ratio = opt->net_ratio;
	this->aopt.if_narrow_lane = opt->if_narrow_lane;
	if (opt->if_min_epochs > 0)
		this->aopt.if_min_epochs = opt->if_min_epochs;
	if (opt->if_max_sigma > 0.0)
		this->aopt.if_max_sigma = opt->if_max_sigma;
	if (opt->if_sys_m >= 0.0)
		this->aopt.if_sys_m = opt->if_sys_m;
	if (opt->if_max_frac > 0.0)
		this->aopt.if_max_frac = opt->if_max_frac;
	vbs_basemeta_default_opt(&this->bopt);
	vbs_lock_default_opt(&lopt);
	if (3.0 * opt->interval_s > lopt.max_gap_s)
		lopt.max_gap_s = 3.0 * opt->interval_s;
	if ((this->lock = vbs_lock_new(&lopt)) == NULL) {
		vbs_session_free(this);
		return NULL;
	}
	if ((this->rtcm = vbs_rtcm_out_new(opt->staid, vbs_xyz)) == NULL) {
		vbs_session_free(this);
		return NULL;
	}

	for (i = 0; i < nsta; i++) {
		struct station *s = &this->sta[i];

		s->def = sta[i];
		if ((s->pp = vbs_preproc_new(opt->interval_s, 3.0)) == NULL) {
			vbs_session_free(this);
			return NULL;
		}
		for (a = 0; a < VBS_SESSION_RING; a++)
			s->ring[a].used = 1;
		if (i == master)
			continue;
		if ((s->ekf = vbs_netekf_new(sta[master].xyz, sta[i].xyz, &this->nopt)) == NULL) {
			vbs_session_free(this);
			return NULL;
		}
		for (a = 0; a < MAXSAT; a++)
			for (c = 0; c < NCOMB; c++)
				s->last_fix[a][c] = NAN;
	}
	return this;
}

void vbs_session_free(struct vbs_session *this) {
	int i;

	if (this == NULL)
		return;
	for (i = 0; i < this->nsta; i++) {
		vbs_preproc_free(this->sta[i].pp);
		vbs_netekf_free(this->sta[i].ekf);
	}
	vbs_lock_free(this->lock);
	vbs_rtcm_out_free(this->rtcm);
	free(this);
}

int vbs_session_refresh_glo(struct vbs_session *this) {
	int n = vbs_rtcm_out_set_glo_fcn(this->rtcm, this->nav);

	if (n > 0)
		this->glo_fcn_installed = 1;
	return n;
}

/* ---- the epoch ring ---------------------------------------------------------------------- */

static struct slot *find_slot(struct station *s, gtime_t t) {
	int i;

	for (i = 0; i < VBS_SESSION_RING; i++)
		if (!s->ring[i].used && s->ring[i].filled &&
		    fabs(timediff(s->ring[i].t, t)) <= 1e-3)
			return &s->ring[i];
	return NULL;
}

int vbs_session_push(struct vbs_session *this, int station, const obsd_t *obs, int n,
    double now_unix) {
	struct station *s;
	struct slot *dst = NULL;
	int i;
	gtime_t t;

	if (station < 0 || station >= this->nsta || n <= 0)
		return -1;
	s = &this->sta[station];
	t = obs[0].time;
	/*
	 * An epoch older than the one last processed can never be used; the filter has moved on.
	 * Drop it here rather than let it sit in the ring pretending to be pending.
	 */
	if (this->have_processed && timediff(t, this->last_processed) <= 1e-3)
		return -1;
	if ((dst = find_slot(s, t)) == NULL) {
		/* a free slot, else the oldest live one -- if it is older than this epoch */
		struct slot *oldest = NULL;

		for (i = 0; i < VBS_SESSION_RING; i++) {
			if (s->ring[i].used || !s->ring[i].filled) {
				dst = &s->ring[i];
				break;
			}
			if (oldest == NULL || timediff(s->ring[i].t, oldest->t) < 0.0)
				oldest = &s->ring[i];
		}
		if (dst == NULL) {
			if (oldest == NULL || timediff(oldest->t, t) >= 0.0)
				return -1;
			dst = oldest;
		}
	}
	if (n > MAXOBS)
		n = MAXOBS;
	memcpy(dst->obs, obs, n * sizeof *obs);
	dst->n = n;
	dst->t = t;
	dst->arrived = now_unix;
	dst->used = 0;
	dst->filled = 1;
	if (station == this->master)
		this->stats.epochs_in++;
	return 0;
}

/* The earliest pending master epoch, or NULL. */
static struct slot *next_master_epoch(struct vbs_session *this) {
	struct station *m = &this->sta[this->master];
	struct slot *best = NULL;
	int i;

	for (i = 0; i < VBS_SESSION_RING; i++)
		if (!m->ring[i].used && m->ring[i].filled &&
		    (best == NULL || timediff(m->ring[i].t, best->t) < 0.0))
			best = &m->ring[i];
	return best;
}

static int has_later_master_epoch(struct vbs_session *this, gtime_t t) {
	struct station *m = &this->sta[this->master];
	int i;

	for (i = 0; i < VBS_SESSION_RING; i++)
		if (!m->ring[i].used && m->ring[i].filled && timediff(m->ring[i].t, t) > 1e-3)
			return 1;
	return 0;
}

/* Retire every slot of every station at or before t: nothing can use them now. */
static void retire_through(struct vbs_session *this, gtime_t t) {
	int i, k;

	for (i = 0; i < this->nsta; i++)
		for (k = 0; k < VBS_SESSION_RING; k++)
			if (this->sta[i].ring[k].filled && !this->sta[i].ring[k].used &&
			    timediff(this->sta[i].ring[k].t, t) <= 1e-3)
				this->sta[i].ring[k].used = 1;
}

/* ---- output ------------------------------------------------------------------------------ */

static int emit(struct vbs_session *this, const unsigned char *b, int n) {
	if (n <= 0)
		return 0;
	if (this->nout + n > (int)sizeof this->out) {
		/* nobody is draining: keep the newest, the old bytes are stale anyway */
		this->nout = 0;
		if (n > (int)sizeof this->out)
			return -1;
	}
	memcpy(this->out + this->nout, b, n);
	this->nout += n;
	this->stats.bytes_out += n;
	return n;
}

int vbs_session_output(struct vbs_session *this, unsigned char *buf, int len) {
	int n = this->nout < len ? this->nout : len;

	if (n <= 0)
		return 0;
	memcpy(buf, this->out, n);
	memmove(this->out, this->out + n, this->nout - n);
	this->nout -= n;
	return n;
}

/* ---- one epoch --------------------------------------------------------------------------- */

static int process_epoch(struct vbs_session *this, struct slot *ms, double now_unix) {
	static double rs[MAXOBS * 6], dts[MAXOBS * 2], var[MAXOBS], rs_sat[MAXSAT * 6];
	static int svh[MAXOBS];
	static struct vbs_slip mslip[MAXOBS];
	static int mslip_sat[MAXOBS], refixed[MAXSAT];
	static struct vbs_interp_aux iaux[VBS_SESSION_MAXSTA];
	static struct vbs_corr vcorr[MAXOBS];
	static struct vbs_interp_quality vqual[MAXOBS];
	static double azel[MAXSAT * 2];
	static obsd_t synth[MAXOBS];
	static unsigned char buf[8192];
	const obsd_t *mep = ms->obs;
	int nm = ms->n, nmslip, nmslip_sat = 0, nrefixed = 0, ncl = 0, nvalid, nsynth, i, k, len;
	double mxyz[3], vxyz[3], vpos[3];
	gtime_t t = ms->t;
	int incomplete = 0;

	if (nm > MAXOBS)
		nm = MAXOBS;

	satposs(t, mep, nm, this->nav, EPHOPT_BRDC, rs, dts, var, svh);
	vbs_satpos_by_sat(mep, nm, rs, rs_sat);
	for (k = 0; k < nm; k++)
		if (sysidx(mep[k].sat) >= 0) {
			this->stats.in_sys[sysidx(mep[k].sat)]++;
			if (rs[k * 6] == 0.0 && rs[k * 6 + 1] == 0.0 && rs[k * 6 + 2] == 0.0)
				this->stats.nopos_sys[sysidx(mep[k].sat)]++;
		}

	/* the master's station-level arcs, read once and shared by every filter (D45) */
	nmslip = vbs_preproc_epoch(this->sta[this->master].pp, mep, nm, this->nav, mslip, MAXOBS);
	vbs_preproc_arc_ids(this->sta[this->master].pp, this->sta[this->master].arc);
	for (k = 0; k < nmslip && k < MAXOBS; k++) {
		int si = sysidx(mslip[k].sat), b;

		mslip_sat[nmslip_sat++] = mslip[k].sat;
		if (si >= 0)
			for (b = 0; b < 5; b++)
				if (mslip[k].reason & (1 << b))
					this->stats.slip_sys[si][b]++;
		/*
		 * Only constellations this session actually uses. The master's stream carries SBAS
		 * too (EGNOS PRN 123 and 136 over France), the preprocessor tracks them like any
		 * other satellite, and they slip often -- but no filter ever sees them, so counting
		 * their slips here would report a problem the engine does not have. Left in the data
		 * path deliberately: excluding them there would change the measured pipeline for no
		 * gain beyond a few microseconds.
		 */
		if (mslip[k].sat >= 1 && mslip[k].sat <= MAXSAT && !(mslip[k].reason & VBS_SLIP_NEW) &&
		    si >= 0)
			this->slip_by_sat[mslip[k].sat - 1]++;
	}
	for (k = 0; k < nm; k++) {
		int si = sysidx(mep[k].sat), f;

		if (si < 0)
			continue;
		for (f = 0; f < NFREQ + NEXOBS; f++) {
			if (mep[k].L[f] == 0.0)
				continue;
			this->stats.phase_sys[si]++;
			if (mep[k].LLI[f] & 1)
				this->stats.lli1_sys[si]++;
			if (mep[k].LLI[f] & 2)
				this->stats.lli2_sys[si]++;
		}
	}

	/* tide-corrected instantaneous positions, same model, same epoch, both ends (D23) */
	vbs_basemeta_tide_xyz(this->sta[this->master].def.xyz, t, &this->bopt, mxyz);
	vbs_basemeta_tide_xyz(this->vbs_xyz, t, &this->bopt, vxyz);
	ecef2pos(vxyz, vpos);

	for (i = 0; i < this->nsta; i++) {
		struct station *s = &this->sta[i];
		struct slot *as;

		if (i == this->master)
			continue;
		s->ncorr = 0;
		if ((as = find_slot(s, t)) == NULL) {
			incomplete = 1;
			continue;
		}
		vbs_preproc_epoch(s->pp, as->obs, as->n, this->nav, NULL, 0);
		vbs_preproc_arc_ids(s->pp, s->arc);
		vbs_netekf_epoch(s->ekf, t, mep, nm, rs_sat, dts, as->obs, as->n, this->nav,
		    this->sta[this->master].arc, s->arc);
		vbs_ambres_resolve(s->ekf, &this->aopt, this->opt.sys_mask, &this->ar,
		    s->fix, 256, &s->nfix);
		s->ncorr = vbs_netekf_corrections(s->ekf, s->corr, MAXOBS, CORR_GATE_M);
		vbs_netekf_correction_sigmas(s->ekf, s->corr, s->ncorr, s->sd_I, s->sd_T);
		for (k = 0; k < s->ncorr; k++)
			if (s->corr[k].valid && sysidx(s->corr[k].sat) >= 0)
				this->stats.aux_corr_sys[sysidx(s->corr[k].sat)]++;

		memcpy(iaux[ncl].xyz, s->def.xyz, sizeof iaux[ncl].xyz);
		iaux[ncl].backup = s->def.backup;	/* iaux is static: set it every time */
		iaux[ncl].corr = s->corr;
		iaux[ncl].sd_I = s->sd_I;
		iaux[ncl].sd_T = s->sd_T;
		iaux[ncl].ncorr = s->ncorr;
		ncl++;
	}
	if (incomplete)
		this->stats.epochs_incomplete++;

	/* §5.3 case 2: a network ambiguity re-fixed to a different integer, or lost */
	for (i = 0; i < this->nsta; i++) {
		struct station *s = &this->sta[i];
		int q;

		if (i == this->master || s->ekf == NULL)
			continue;
		for (q = 0; q < s->nfix; q++) {
			const struct vbs_amb_fix *f = &s->fix[q];
			double *slot;

			if (f->sat <= 0 || f->sat > MAXSAT || f->comb_id < 0 || f->comb_id >= NCOMB)
				continue;
			slot = &s->last_fix[f->sat - 1][f->comb_id];
			if (!isnan(*slot) && *slot != f->integer) {
				int r, dup = 0;

				for (r = 0; r < nrefixed; r++)
					if (refixed[r] == f->sat)
						dup = 1;
				if (!dup && nrefixed < MAXSAT)
					refixed[nrefixed++] = f->sat;
			}
			*slot = f->integer;
		}
	}
	/* an arc change at the master resets every baseline: the ambiguity is lost */
	for (k = 0; k < MAXSAT; k++) {
		long arc = this->sta[this->master].arc[k];

		if (arc != this->last_arc_m[k]) {
			int r, dup = 0;

			for (r = 0; r < nrefixed; r++)
				if (refixed[r] == k + 1)
					dup = 1;
			if (!dup && nrefixed < MAXSAT && this->last_arc_m[k] != 0)
				refixed[nrefixed++] = k + 1;
			for (i = 0; i < this->nsta; i++) {
				int c;

				for (c = 0; c < NCOMB; c++)
					this->sta[i].last_fix[k][c] = NAN;
			}
		}
		this->last_arc_m[k] = arc;
	}

	/* elevations at the VBS for the troposphere model term, from the same geometry */
	memset(azel, 0, sizeof azel);
	for (k = 0; k < nm; k++) {
		double ev[3], ae[2];

		if (mep[k].sat <= 0 || mep[k].sat > MAXSAT)
			continue;
		if (geodist(rs + k * 6, vxyz, ev) <= 0.0)
			continue;
		satazel(vpos, ev, ae);
		azel[(mep[k].sat - 1) * 2] = ae[0];
		azel[(mep[k].sat - 1) * 2 + 1] = ae[1];
		this->last_el[mep[k].sat - 1] = ae[1] * R2D;
	}

	nvalid = vbs_interp_epoch(t, mxyz, vxyz, iaux, ncl, &this->iopt, azel, vcorr, vqual, MAXOBS);
	if (nvalid <= 0) {
		this->stats.epochs_no_corr++;
		return 0;
	}
	/* how often the reserve is actually needed: meant to be almost never */
	for (k = 0; k < nvalid; k++)
		if (vcorr[k].valid && vqual[k].backup_used)
			this->stats.sats_backup++;
	/* --out-elmask (D52) */
	if (this->opt.out_elmask_deg > 0.0)
		for (k = 0; k < nvalid; k++) {
			int sat = vcorr[k].sat;

			if (sat >= 1 && sat <= MAXSAT && vcorr[k].valid &&
			    azel[(sat - 1) * 2 + 1] < this->opt.out_elmask_deg * D2R)
				vcorr[k].valid = 0;
		}
	/* --no-iono: the dispersive component withheld (UM980) */
	if (!this->opt.dispersive)
		for (k = 0; k < nvalid; k++)
			vcorr[k].dI = 0.0;
	else if (vbs_ionopost_active(&this->opt.iono_post)) {
		struct vbs_ionopost_opt ipo = this->opt.iono_post;

		ipo.sys_mask = this->opt.iono_sys_mask;
		vbs_ionopost_apply(vcorr, nvalid, azel, vpos, mep, nm, this->nav, &ipo, NULL);
	} else if (this->opt.iono_sys_mask)
		for (k = 0; k < nvalid; k++)
			if (!(satsys(vcorr[k].sat, NULL) & this->opt.iono_sys_mask))
				vcorr[k].dI = 0.0;
	if (!this->opt.dispersive && this->opt.iono_post.drop_sys)
		for (k = 0; k < nvalid; k++)
			if (satsys(vcorr[k].sat, NULL) & this->opt.iono_post.drop_sys)
				vcorr[k].valid = 0;

	for (k = 0; k < nvalid; k++)
		if (vcorr[k].valid && sysidx(vcorr[k].sat) >= 0)
			this->stats.corr_sys[sysidx(vcorr[k].sat)]++;
	nsynth = vbs_synth_epoch(mep, nm, rs, mxyz, vxyz, vcorr, nvalid, this->nav, synth);
	if (nsynth <= 0) {
		this->stats.epochs_no_corr++;
		return 0;
	}
	for (k = 0; k < nsynth; k++)
		if (sysidx(synth[k].sat) >= 0)
			this->stats.out_sys[sysidx(synth[k].sat)]++;

	/* §5.3: stamp the output's lock-time signalling, in place on `synth` */
	vbs_lock_update(this->lock, t, mep, nm, synth, nsynth, vcorr, nvalid,
	    mslip_sat, nmslip_sat, refixed, nrefixed);

	memcpy(this->last_synth, synth, nsynth * sizeof *synth);
	this->nlast_synth = nsynth;
	this->last_synth_time = t;

	if (!this->converged) {
		/* the switch: the tile's own station messages lead, then its observations */
		this->converged = 1;
		this->stats.converged = 1;
		this->stats.mode_switches++;
		this->pos_sent_once = 0;
	}
	this->last_vbs_epoch_unix = now_unix;

	/* station messages every pos_period_s, and before the first observations */
	if (!this->pos_sent_once || timediff(t, this->last_pos_sent) >= this->opt.pos_period_s) {
		if ((len = vbs_rtcm_out_1005(this->rtcm, buf, sizeof buf)) > 0)
			emit(this, buf, len);
		if ((len = vbs_rtcm_out_station(this->rtcm, (this->opt.sys_mask & SYS_GLO) != 0,
		    buf, sizeof buf)) > 0)
			emit(this, buf, len);
		this->last_pos_sent = t;
		this->pos_sent_once = 1;
		this->last_keepalive_unix = now_unix;
		this->epochs_out_at_keepalive = this->stats.epochs_out + 1;
	}
	if ((this->opt.sys_mask & SYS_GLO) && !this->glo_fcn_installed)
		vbs_session_refresh_glo(this);
	len = vbs_rtcm_out_msm(this->rtcm, t, synth, nsynth, buf, sizeof buf);
	if (len < 0) {
		/* GLONASS without channel numbers: refuse the epoch rather than emit it wrong (D90) */
		this->stats.epochs_no_corr++;
		return 0;
	}
	emit(this, buf, len);
	this->stats.epochs_out++;
	this->stats.sats_out += nsynth;
	this->stats.glo_dropped = vbs_rtcm_out_glo_dropped(this->rtcm);
	return 1;
}

int vbs_session_run(struct vbs_session *this, double now_unix) {
	int emitted = 0;
	static unsigned char kbuf[1024];

	for (;;) {
		struct slot *ms = next_master_epoch(this);
		int complete = 1, i;

		if (ms == NULL)
			break;
		for (i = 0; i < this->nsta; i++)
			if (i != this->master && find_slot(&this->sta[i], ms->t) == NULL)
				complete = 0;
		if (!complete && !has_later_master_epoch(this, ms->t) &&
		    now_unix - ms->arrived < this->opt.max_wait_s)
			break;		/* give the auxiliaries a moment longer */
		emitted += process_epoch(this, ms, now_unix);
		this->last_processed = ms->t;
		this->have_processed = 1;
		retire_through(this, ms->t);
	}
	/*
	 * Lost convergence: no VBS epoch for reconverge_s after having served some. Back to
	 * pass-through rather than silence; the base's own station messages announce the change.
	 */
	if (this->converged && this->opt.reconverge_s > 0.0 &&
	    now_unix - this->last_vbs_epoch_unix > this->opt.reconverge_s) {
		this->converged = 0;
		this->stats.converged = 0;
		this->stats.mode_switches++;
	}
	/*
	 * Keepalive in VBS mode: the station messages on their period even when no epoch passed
	 * the gate, so the caster pulling this stream never takes a quiet upstream for a dead one.
	 */
	if (this->converged && now_unix - this->last_keepalive_unix >= this->opt.pos_period_s) {
		if (this->stats.epochs_out == this->epochs_out_at_keepalive) {
			int n = vbs_session_station_msgs(this, kbuf, sizeof kbuf);

			if (n > 0)
				emit(this, kbuf, n);
		}
		this->last_keepalive_unix = now_unix;
		this->epochs_out_at_keepalive = this->stats.epochs_out;
	}
	return emitted;
}

int vbs_session_passthrough(struct vbs_session *this, const unsigned char *buf, int n) {
	if (this->converged || n <= 0)
		return 0;
	this->stats.passthrough_bytes += n;
	return emit(this, buf, n);
}

int vbs_session_converged(const struct vbs_session *this) { return this->converged; }

int vbs_session_last_synth(const struct vbs_session *this, const obsd_t **obs, gtime_t *time) {
	*obs = this->last_synth;
	*time = this->last_synth_time;
	return this->nlast_synth;
}

void vbs_session_get_stats(const struct vbs_session *this, struct vbs_session_stats *out) {
	static const int SYS[4] = {SYS_GPS, SYS_GLO, SYS_GAL, SYS_CMP};
	int i, k;

	*out = this->stats;
	memset(out->rej, 0, sizeof out->rej);
	for (k = 0; k < 3; k++) {
		int best = -1;

		for (i = 0; i < MAXSAT; i++) {
			int taken = 0, q;

			for (q = 0; q < k; q++)
				if (out->worst_sat[q] == i + 1)
					taken = 1;
			if (!taken && this->slip_by_sat[i] > 0 && (best < 0 || this->slip_by_sat[i] > this->slip_by_sat[best]))
				best = i;
		}
		out->worst_sat[k] = best >= 0 ? best + 1 : 0;
		out->worst_n[k] = best >= 0 ? this->slip_by_sat[best] : 0;
		out->worst_el[k] = best >= 0 ? this->last_el[best] : -1.0;
	}
	for (i = 0; i < this->nsta; i++) {
		struct vbs_netekf_stats ns;
		const struct vbs_sigset *ss;

		if (this->sta[i].ekf == NULL)
			continue;
		vbs_netekf_get_stats(this->sta[i].ekf, &ns);
		for (k = 0; k < VBS_NET_NREJ && k < 8; k++)
			out->rej[k] += ns.obs_rejected[k];
		if ((ss = vbs_netekf_sigset(this->sta[i].ekf)) != NULL)
			for (k = 0; k < 4; k++)
				if (out->bands_sys[k] == 0)
					out->bands_sys[k] = vbs_sigset_nbands(ss, SYS[k]);
		if (out->fsat_sys[0] == 0 && out->fsat_sys[1] == 0 && out->fsat_sys[2] == 0 && out->fsat_sys[3] == 0) {
			static struct vbs_netekf_satinfo si[VBS_NET_MAXSAT];
			int n = vbs_netekf_satinfo(this->sta[i].ekf, si, VBS_NET_MAXSAT), q;

			for (q = 0; q < n; q++) {
				int sx = sysidx(si[q].sat);

				if (sx < 0)
					continue;
				out->fsat_sys[sx]++;
				out->fion_sys[sx] += si[q].sd_ion;
				out->famb_sys[sx] += si[q].amb_sigma_max;
				if (si[q].epochs_tracked > out->ftrack_sys[sx])
					out->ftrack_sys[sx] = si[q].epochs_tracked;
			}
			for (k = 0; k < 4; k++)
				if (out->fsat_sys[k] > 0) {
					out->fion_sys[k] /= out->fsat_sys[k];
					out->famb_sys[k] /= out->fsat_sys[k];
				}
		}
	}
}

int vbs_session_station_msgs(struct vbs_session *this, unsigned char *buf, int len) {
	int used = 0, n;

	if ((n = vbs_rtcm_out_1005(this->rtcm, (uint8_t *)buf, len)) > 0)
		used += n;
	if ((n = vbs_rtcm_out_station(this->rtcm, (this->opt.sys_mask & SYS_GLO) != 0,
	    (uint8_t *)buf + used, len - used)) > 0)
		used += n;
	return used;
}

void vbs_session_set_master_station(struct vbs_session *this, const sta_t *sta) {
	vbs_rtcm_out_set_station(this->rtcm, sta);
}


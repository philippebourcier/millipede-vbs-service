#ifndef __VBS_SESSION_H__
#define __VBS_SESSION_H__

/*
 * One tile's synthesis: the Phase 1 pipeline core, lifted out of tools/vbs_pipeline.c and
 * driven epoch by epoch instead of file by file (spec §6b "one shared synthesis per tile").
 *
 * The contract that matters is that this is the SAME computation the measurements were made
 * with. tools/vbs_pipeline.c produced every number in the Phase 1 report; if the live path
 * diverged from it, the report would describe a different engine from the one deployed.
 * tests/test_session.c therefore drives both over the same RINEX and requires the synthesized
 * observations to agree epoch for epoch -- the only proof that lifting the code did not change
 * it.
 *
 * Observations arrive by push, one epoch per station, from whatever source: live RTCM decoded by
 * vbs_feed, or RINEX in the test. vbs_session_run() then processes every master epoch that is
 * complete (all auxiliaries present) or has waited long enough, in time order, and appends the
 * resulting RTCM to an output buffer the caller drains to its clients.
 */

#include "rtklib.h"
#include "vbs/vbs_ionopost.h"

#define VBS_SESSION_MAXSTA	8		/* master plus auxiliaries */
#define VBS_SESSION_RING	8		/* epochs buffered per station */
#define VBS_SESSION_OUTBUF	16384

struct vbs_session_station {
	char name[32];
	double xyz[3];			/* catalogue ECEF, metres */
	char rectype[32];
	/*
	 * The "2 + 1 backup" reserve (vbs_cluster.h). Processed every epoch like any auxiliary,
	 * because its ambiguities cannot be fixed the moment it is needed, but kept out of the
	 * fit unless a satellite's primaries fall below the interpolator's floor.
	 */
	int backup;
};

struct vbs_session_opt {
	int sys_mask;			/* constellations; SYS_GLO only on a homogeneous cluster */
	int dispersive;			/* apply the ionospheric correction (Mosaic yes, UM980 no) */
	/*
	 * With `dispersive` on, the constellations that carry it (SYS_* mask); 0 = every one. The
	 * others are served as --no-iono serves them. vbs_pipeline --iono-sys, D167: a UM980 cluster's
	 * code-levelled GPS L2 and GLONASS ionosphere is worse than none, its Galileo and BeiDou is not.
	 */
	int iono_sys_mask;
	/* D168 candidates for the constellations iono_sys_mask leaves out (vbs_ionopost); all off */
	struct vbs_ionopost_opt iono_post;
	/*
	 * Take GPS L2 from L2C rather than the first L2 slot (Mosaic yes, UM980 not measured).
	 *
	 * D138: on the Mosaic ladder this dominates the shipped recipe on all four rungs inside
	 * the envelope -- more on-spec time AND less harmful distance on every one, 11 wins and
	 * 0 losses over 20 rung-days, total harmful distance 27 108 -> 11 661 m. It has not been
	 * run on the UM980 ladder when this was written. It has been since (D140): on the three
	 * clean UM980 rungs it is INERT -- 1 win, 0 losses, 14 trades over 15 rung-days, every
	 * difference in the third significant figure -- because a UM980's `2W` column already
	 * carries L2C, so preferring L2C changes which column is read and not which signal. On a
	 * Mosaic `2W` is genuine semi-codeless, so there it changes the observable and pays.
	 * `vbs_policy.ENVELOPE` carries the same per-family split.
	 */
	int prefer_l2c;
	double interval_s;		/* nominal observation interval */
	double out_elmask_deg;		/* satellites below this are not emitted (D52) */
	double ar_corr_time_s;		/* --ar-corr-time: search-side covariance inflation (D48) */
	double net_ratio;		/* network ratio-test threshold; 0 = engine default */
	int if_narrow_lane;		/* --if-narrow-lane */
	int if_min_epochs;
	double if_max_sigma;
	double if_sys_m;
	double if_max_frac;
	double q_ion;
	/*
	 * Fewest auxiliaries a satellite's correction is fitted from; 0 = the interpolator's
	 * default (2, which lets the "2 + 1 backup" reserve stand in). 3 is the floor every
	 * result before 2026-10-02 was measured with.
	 */
	int interp_min_aux;
	int staid;			/* RTCM station id (from the tile) */
	double pos_period_s;		/* how often to send 1005 (+1033, 1230) */
	double max_wait_s;		/* run an incomplete epoch after this long */
	double reconverge_s;		/* no VBS epoch for this long after converging: back to pass-through */
};

/* The shipped configuration: the flag set every Phase 1 result was measured with. */
void vbs_session_default_opt(struct vbs_session_opt *opt);

struct vbs_session_stats {
	int converged;			/* serving the synthesis (1) or passing a base through (0) */
	long mode_switches;
	long passthrough_bytes;
	long epochs_in;			/* master epochs received */
	long epochs_out;		/* VBS epochs emitted */
	long epochs_incomplete;		/* run without every auxiliary */
	long epochs_no_corr;		/* skipped: no converged correction */
	long sats_out;			/* satellites emitted, summed */
	long sats_backup;		/* ... of which the reserve was needed to correct */
	long bytes_out;
	long glo_dropped;
	/* per constellation (GPS, GLO, GAL, BDS): master observations in, corrections valid at the
	 * tile, satellites emitted -- to see which system is being lost where */
	long in_sys[4], corr_sys[4], out_sys[4];
	long aux_corr_sys[4];		/* auxiliary corrections passing the gate, before interpolation */
	int bands_sys[4];		/* common bands per constellation on the first baseline; 0 = none */
	long nopos_sys[4];		/* master observations with no satellite position (no ephemeris) */
	long phase_sys[4], lli1_sys[4], lli2_sys[4];	/* master phases, with LLI lock-loss, with half-cycle flag */
	long slip_sys[4][5];		/* master slips by system and reason bit (LLI, GF, MW, GAP, NEW) */
	int worst_sat[3];		/* the three satellites with most master slips, and their counts */
	long worst_n[3];
	double worst_el[3];		/* their elevation at the last epoch, degrees (-1 unknown) */
	/* the first auxiliary filter's view, per system: satellites in state, mean ionosphere sigma (m),
	 * mean worst ambiguity sigma (cycles), longest tracked run (epochs) */
	int fsat_sys[4];
	double fion_sys[4], famb_sys[4];
	long ftrack_sys[4];
	long rej[8];			/* filter rejections by VBS_NET_REJ_* reason, summed over auxiliaries */
};

struct vbs_session;

/*
 * Create a session for a VBS at `vbs_xyz`, from `sta[0..nsta)` with `master` as the master.
 * `nav` is shared and owned by the caller; it must stay valid and may be updated between calls
 * (ephemerides arriving on the base streams).
 */
struct vbs_session *vbs_session_new(const double *vbs_xyz, const struct vbs_session_station *sta,
	int nsta, int master, const struct vbs_session_opt *opt, nav_t *nav);
void vbs_session_free(struct vbs_session *this);

/*
 * Deliver one epoch of one station. `now_unix` is the wall clock of arrival, used for the
 * incomplete-epoch timeout. Returns 0, or -1 if the epoch was refused (ring full of newer data).
 */
int vbs_session_push(struct vbs_session *this, int station, const obsd_t *obs, int n,
	double now_unix);

/*
 * Process every runnable master epoch. Returns the number of VBS epochs emitted; their RTCM is
 * appended to the output buffer.
 */
int vbs_session_run(struct vbs_session *this, double now_unix);

/*
 * Pass-through (2026-09-13: waiting 3 minutes for a first fix is unacceptable).
 *
 * The network filters need about three minutes of a new cluster before the first correction
 * passes its gate, and a rover must not sit on an empty stream for that. Until the synthesis
 * has converged, the session hands its clients a real base's own RTCM verbatim -- the nearest
 * live base, which is exactly what NEAR would have given -- and the moment the first VBS epoch
 * is ready it switches: the tile's 1005/1033/1230, then the synthesized observations. The
 * rover sees one base change, the event it already handles when NEAR switches base.
 *
 * Verbatim, not re-encoded and not translated: a geometric translation of the master's
 * observations to the tile would carry the master-tile atmosphere unannounced under a
 * short-baseline ARP, and a rover weighting its ionosphere by baseline length would fix wrong
 * integers on it. The base's own stream under its own ARP tells the rover the truth.
 *
 * Bytes given here are emitted only while not converged; once converged they are dropped, and
 * if the synthesis then produces nothing for `reconverge_s` the session goes back to pass-through
 * (one more base change, announced by the base's own station messages).
 */
int vbs_session_passthrough(struct vbs_session *this, const unsigned char *buf, int n);
int vbs_session_converged(const struct vbs_session *this);

/* Drain the output buffer: copies up to `len` bytes to `buf`, returns how many. */
int vbs_session_output(struct vbs_session *this, unsigned char *buf, int len);

/*
 * Debug hook for the identity test: the synthesized observations of the LAST emitted epoch,
 * before encoding. Returns the count, and sets *time.
 */
int vbs_session_last_synth(const struct vbs_session *this, const obsd_t **obs, gtime_t *time);

void vbs_session_get_stats(const struct vbs_session *this, struct vbs_session_stats *out);

/*
 * Install the GLONASS channel numbers from the shared nav into the encoder (D90). Call when
 * ephemerides arrive; cheap, idempotent. Returns the number installed.
 */
int vbs_session_refresh_glo(struct vbs_session *this);

/*
 * Encode the station messages (1005, 1033, 1230 when GLONASS is served) into `buf` now, for a
 * client that has just attached and cannot wait for the next periodic set. Returns bytes.
 */
int vbs_session_station_msgs(struct vbs_session *this, unsigned char *buf, int len);

/* Give the encoder the master's station record (1033/1230 content); see vbs_rtcm_out_set_station. */
void vbs_session_set_master_station(struct vbs_session *this, const sta_t *sta);

#endif /* __VBS_SESSION_H__ */

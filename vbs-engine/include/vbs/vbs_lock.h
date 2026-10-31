#ifndef __VBS_LOCK_H__
#define __VBS_LOCK_H__

#include "rtklib.h"

#include "vbs/vbs_synth.h"

/*
 * Lock-time tracking for the VBS output -- spec §5.3.
 *
 * MSM has no RINEX-style loss-of-lock flag. The only way to tell a rover that phase continuity
 * is broken is the **Lock Time Indicator**: resetting it, or letting it decrease, is the
 * signal. A rover that does not see the reset will carry its ambiguity across the
 * discontinuity and keep a fix that is now wrong by whole cycles -- which is worse than losing
 * the fix, because the position stays confident.
 *
 * §5.3 lists four events that must reset it:
 *
 *   1. the master reports a slip on that signal -- the VBS is the master's observation
 *      translated, so the master's discontinuity is the VBS's discontinuity;
 *   2. the network ambiguity for that satellite is re-fixed or lost -- the correction then
 *      steps by an integer number of cycles, which is exactly a slip as far as the rover is
 *      concerned, even though nothing happened at the master;
 *   3. the interpolated correction jumps by more than a threshold, e.g. because an auxiliary
 *      dropped out of the cluster and the plane fit changed underneath;
 *   4. the satellite reappears after exclusion by FDE.
 *
 * Cases 2-4 are the ones that matter here and the reason this is not simply a passthrough of
 * the master's lock time: the VBS carries a *computed* phase, so it can break continuity when
 * the master did not. **This must be in place before §8 step 4**, because a missed reset makes
 * the rover hold a wrong fix -- inflating the fix rate and deflating the false-fix rate, both
 * in the flattering direction, which is the worst way for an acceptance test to be wrong.
 *
 * **How the signal actually reaches the wire.** MRTKLIB's encoder derives the MSM Lock Time
 * Indicator from `obsd_t::LLI`: `locktime()` restarts the lock epoch whenever `LLI & 1` is
 * set, and otherwise reports the elapsed seconds. So the job here is to decide *when* to set
 * that bit on the synthesized observation; the encoder does the rest. The seconds counted here
 * are for diagnosis and for the tests -- they are what the rover will be told.
 *
 * **The master's slips do not come from the master's LLI.** D21/D22 measured that flag on this
 * network: it is derived from MSM lock-time indicators, and on the storm day it reported *zero*
 * slips while the geometry-free and Melbourne-Wübbena combinations found 1 991. Case 1 is
 * therefore driven by the slips `vbs_preproc` actually detects, with the master's LLI kept only
 * as an additional trigger. Trusting the flag alone would reproduce, on the VBS output, exactly
 * the failure that D21 found in the input.
 */

/* Why continuity was broken, for diagnosis. A reset usually has one cause, but not always. */
#define VBS_LOCK_MASTER_SLIP	0x01	/* §5.3 case 1 */
#define VBS_LOCK_AMB_REFIX	0x02	/* §5.3 case 2 */
#define VBS_LOCK_CORR_JUMP	0x04	/* §5.3 case 3 */
#define VBS_LOCK_FDE_RETURN	0x08	/* §5.3 case 4 */
#define VBS_LOCK_NEW		0x10	/* first observation: no continuity to break */

struct vbs_lock_opt {
	/*
	 * `vbs.correction_jump_threshold` (§5.3 case 3), metres. A correction that moves by
	 * more than this between epochs is treated as a discontinuity rather than as the
	 * atmosphere. It has to sit above the real epoch-to-epoch atmospheric change --
	 * `03-error-budget.md` measured the ionospheric term at 1.4-12.5 cm over a *baseline*,
	 * and its change between two epochs 1 s apart is far smaller -- and below the smallest
	 * step a re-fix can produce, which is one cycle of the longest wavelength involved.
	 */
	double jump_threshold_m;
	double max_gap_s;		/* longer than this and continuity is not assumed */
};

void vbs_lock_default_opt(struct vbs_lock_opt *opt);

struct vbs_lock_stats {
	long resets[5];			/* indexed by bit position of VBS_LOCK_* */
	long epochs;
	long tracked;			/* satellite-signal slots carried */
};

struct vbs_lock;

struct vbs_lock *vbs_lock_new(const struct vbs_lock_opt *opt);
void vbs_lock_free(struct vbs_lock *this);

/*
 * Update lock times for one epoch and write them into the synthesized observations.
 *
 * `master` is the master's epoch (its own lock times are the baseline, §5.3's last line),
 * `synth`/`nsynth` the synthesized VBS observations to stamp, `corr`/`ncorr` the interpolated
 * corrections for this epoch, and `refixed`/`nrefixed` the satellites whose network ambiguity
 * changed since the previous epoch (§5.3 case 2) -- pass NULL/0 if none.
 *
 * Returns the number of satellite-signals whose lock was reset.
 */
int vbs_lock_update(struct vbs_lock *this, gtime_t time,
	const obsd_t *master, int nmaster, obsd_t *synth, int nsynth,
	const struct vbs_corr *corr, int ncorr,
	const int *master_slip, int nmaster_slip,
	const int *refixed, int nrefixed);

void vbs_lock_get_stats(const struct vbs_lock *this, struct vbs_lock_stats *out);

/*
 * Lock time in seconds for one satellite/slot as of the last update, or -1 if untracked.
 * Exposed so the encoder and the tests can read the same number the rover will see.
 */
double vbs_lock_seconds(const struct vbs_lock *this, int sat, int slot);

#endif /* __VBS_LOCK_H__ */

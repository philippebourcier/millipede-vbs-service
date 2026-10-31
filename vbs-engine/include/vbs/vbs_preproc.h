#ifndef __VBS_PREPROC_H__
#define __VBS_PREPROC_H__

#include "rtklib.h"

#include "vbs/vbs_signal.h"

/*
 * Per-station preprocessing: cycle-slip detection -- spec S4.2.
 *
 * The spec asks for slip detection from "geometry-free and Melbourne-Wubbena combinations,
 * plus lock time". All three, deliberately, because each is blind to something the others
 * catch:
 *
 *   - Geometry-free (L_a - L_b) cancels geometry, clocks and troposphere, leaving ionosphere
 *     plus the ambiguity difference. It sees almost any slip -- except a slip whose cycle
 *     counts on the two bands happen to move it by about the ionosphere's own step. For GPS
 *     L1/L2 the classic blind spot is a 9:7 slip (9 cycles L1, 7 cycles L2), which shifts GF
 *     by under a millimetre.
 *   - Melbourne-Wubbena (wide-lane phase minus narrow-lane code) cancels geometry, clocks,
 *     troposphere *and* ionosphere, leaving the wide-lane ambiguity plus code noise. It sees
 *     exactly the slips GF is blind to, because those change the wide-lane count. It is noisy
 *     -- code noise is tens of centimetres against a 0.86 m wide-lane wavelength -- so it is
 *     used against a running mean rather than epoch to epoch.
 *   - Lock time / LLI is the receiver's own opinion. On this network it is not trustworthy on
 *     its own: docs/vbs/03-error-budget.md measured jumps of 50-120 m between consecutive
 *     epochs with LLI clear, because the dataset's flags are derived from MSM lock-time
 *     indicators (data/README.md S4). It is kept as a trigger, never as the only one.
 *
 * A satellite with only one usable band (BeiDou on this network, S3.1 of the design study)
 * has neither GF nor MW available. That is reported rather than papered over: such a
 * satellite can only be slip-checked by lock time, and the caller should know its slip
 * detection is weaker.
 */

/* Why a slip was declared. A slip usually trips more than one of these. */
#define VBS_SLIP_LLI	0x01	/* receiver's loss-of-lock flag */
#define VBS_SLIP_GF	0x02	/* geometry-free jump */
#define VBS_SLIP_MW	0x04	/* Melbourne-Wubbena jump */
#define VBS_SLIP_GAP	0x08	/* too long since this satellite was last seen */
#define VBS_SLIP_NEW	0x10	/* first observation of an arc: no continuity to defend */

struct vbs_slip {
	int sat;
	int band;		/* VBS_BAND_*, or 0 when the cause is not band-specific */
	unsigned reason;	/* OR of VBS_SLIP_* */
};

/* What the detectors could actually be applied to, for an honest coverage report. */
struct vbs_preproc_stats {
	long epochs;
	long sat_epochs;		/* satellite-epochs seen */
	long dual_band;			/* of those, with two usable bands (GF and MW possible) */
	long single_band;		/* only lock time available */
	long slips;			/* satellite-epochs declared slipped */
	long slip_lli, slip_gf, slip_mw, slip_gap, slip_new;
	long slip_missed_by_lli;	/* GF or MW fired, LLI did not -- the point of S4.2 */
	long slip_lli_only;		/* LLI fired, neither GF nor MW agreed */
};

struct vbs_preproc;

/*
 * sampling_s is the nominal epoch interval; a gap longer than gap_factor * sampling_s ends
 * the arc. Pass 0 for gap_factor to take the default (3).
 */
struct vbs_preproc *vbs_preproc_new(double sampling_s, double gap_factor);
void vbs_preproc_free(struct vbs_preproc *this);

/*
 * Process one epoch of one station. Writes at most maxout slip records and returns how many
 * were written, or -1 on error. Satellites with no slip produce no record.
 */
int vbs_preproc_epoch(struct vbs_preproc *this, const obsd_t *obs, int n, const nav_t *nav,
	struct vbs_slip *out, int maxout);

void vbs_preproc_get_stats(const struct vbs_preproc *this, struct vbs_preproc_stats *out);

/*
 * Station-level arc identity for one satellite: a counter that advances on every slip declared
 * for it, of any cause. Two epochs with the same arc id are on the same continuous phase arc
 * *at this station*.
 *
 * This is the single source of truth for arcs (D44/D45). There is one preprocessor per
 * station, run over that station's whole stream, and every consumer -- each network filter
 * that includes the station, the loop-closure check, and the §5.3 lock signalling on the VBS
 * output -- takes its arc boundaries from here. The alternative, each baseline filter running
 * its own detector over its own epoch subset, made the same slip land in different filters at
 * different epochs or not at all: closure then differenced across arc boundaries and produced
 * integer "misclosures" equal to the slip size (D44), and there was no single answer to which
 * detector drove the output LLI.
 *
 * Copied for every satellite into `out` (MAXSAT entries, indexed by satellite number - 1);
 * a satellite never seen reads 0.
 */
void vbs_preproc_arc_ids(const struct vbs_preproc *this, long *out);
/*
 * Melbourne-Wübbena scatter of a satellite over its current arc, in wide-lane cycles, or
 * NAN while the arc is too short to say (D53: a data-quality gate for stations whose
 * measurements wander; the C/N0 mask was a proxy for this).
 */
double vbs_preproc_mw_sigma(const struct vbs_preproc *this, int sat);

#endif /* __VBS_PREPROC_H__ */

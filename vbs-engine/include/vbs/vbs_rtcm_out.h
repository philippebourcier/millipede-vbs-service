#ifndef __VBS_RTCM_OUT_H__
#define __VBS_RTCM_OUT_H__

#include <stdint.h>
#include <stdio.h>

#include "rtklib.h"

/*
 * RTCM output for a VBS point -- spec S5.1.
 *
 * Wraps MRTKLIB's RTCM3 encoder rather than packing bits here: spec S0 forbids
 * re-implementing well-tested GNSS primitives, and MRTKLIB's encode_rtcm3() already covers
 * every MSM type we need (1071-1077, 1091-1097, 1121-1127) plus 1005/1006/1033/1230. The
 * bit-level MSM layout the caster carries in caster/rtcm.c (for its MSM7->MSM4 downconversion)
 * was the fallback if that had not been the case; it is not needed.
 *
 * What this module owns is the *policy* around the encoder, which is where the spec's
 * requirements actually live:
 *   - one superset stream per VBS point (S5.1): every message a rover might need, in one
 *     buffer, because rovers ignore what they do not use;
 *   - the virtual station ID comes from a reserved range and is frozen for the session
 *     (S3.3/S5.1), never reused by a physical base;
 *   - 1005 carries exactly the VBS point coordinates and never changes mid-session (S3.3).
 */

/* Message types emitted per constellation. MSM7 per spec S5.1. */
#define VBS_MSM7_GPS	1077
#define VBS_MSM7_GLO	1087
#define VBS_MSM7_GAL	1097
#define VBS_MSM7_BDS	1127

struct vbs_rtcm_out;

/*
 * Create an encoder for one VBS point.
 *
 * staid must come from the configured reserved range (vbs.station_id_range) -- this module
 * does not allocate it, because allocation is a property of the VBS point's lifetime, not of
 * the encoder.
 * xyz is the VBS point's ECEF position in metres; it is what 1005 will carry, verbatim.
 */
struct vbs_rtcm_out *vbs_rtcm_out_new(int staid, const double *xyz);
void vbs_rtcm_out_free(struct vbs_rtcm_out *this);

/*
 * Encode one 1005 (station coordinates) message.
 * Returns the number of bytes appended to buf, or -1 on error/overflow.
 */
int vbs_rtcm_out_1005(struct vbs_rtcm_out *this, uint8_t *buf, int buflen);

/*
 * Encode one epoch of observations as MSM7, one message per constellation present.
 *
 * Satellites of systems with no configured message type are skipped. The sync flag is set on
 * every message but the last of the epoch, as RTCM requires for a multi-message epoch.
 * Returns the number of bytes appended to buf, or -1 on error/overflow.
 */
int vbs_rtcm_out_msm(struct vbs_rtcm_out *this, gtime_t time,
	const obsd_t *obs, int n, uint8_t *buf, int buflen);

/*
 * Give the encoder the GLONASS frequency channel numbers (D90).
 *
 * MSM messages for GLONASS carry each satellite's channel number so the rover can work out its
 * frequency, and the encoder takes it from its own `nav`. Without it `fcn_glo()` returns -1 and
 * the encoder **silently drops that satellite from the message** -- the observations are
 * synthesized, handed to the encoder, and never reach the wire. The zero-baseline round-trip
 * caught exactly that: 16 observables lost between encode and decode.
 *
 * Call once after vbs_rtcm_out_new() with a nav that carries GLONASS ephemerides. Returns the
 * number of channel numbers installed.
 */
int vbs_rtcm_out_set_glo_fcn(struct vbs_rtcm_out *this, const nav_t *nav);

/*
 * How many GLONASS satellites have been skipped because their channel number was not known.
 *
 * Zero is the expected value. A non-zero, slowly growing count means some satellite's ephemeris
 * is missing from the nav handed to vbs_rtcm_out_set_glo_fcn() -- worth reporting, not worth
 * dropping the stream for. Forgetting set_glo_fcn() entirely is a different matter and makes
 * vbs_rtcm_out_msm() return -1 the first time GLONASS appears.
 */
long vbs_rtcm_out_glo_dropped(const struct vbs_rtcm_out *this);

/*
 * Give the encoder the master's station record: receiver and antenna descriptors for 1033, the
 * GLONASS code-phase alignment and biases for 1230.
 *
 * The VBS advertises the MASTER's hardware, not its own, because a synthesized observation is
 * the master's observation plus a correction: every receiver-dependent bias in it -- the
 * GLONASS inter-channel bias above all -- is the master's. A rover that matches biases by
 * receiver type (which is what 1033 exists for) must be told the master's type to match them
 * correctly. Copies the record; the caller's copy may go away.
 */
void vbs_rtcm_out_set_station(struct vbs_rtcm_out *this, const sta_t *sta);

/*
 * Encode the periodic station messages: 1033, plus 1230 when `glonass` is set. Sent beside the
 * 1005 every few seconds; a rover joining mid-stream needs all three before it can use the
 * observations. Returns bytes appended, or -1.
 */
int vbs_rtcm_out_station(struct vbs_rtcm_out *this, int glonass, uint8_t *buf, int buflen);

#endif /* __VBS_RTCM_OUT_H__ */

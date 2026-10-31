#ifndef __VBS_FEED_H__
#define __VBS_FEED_H__

/*
 * A physical base's live stream, decoded (spec §6b, base-quality monitor).
 *
 * One feed per mountpoint: an NTRIP client session (vbs_ntrip) whose bytes go straight into
 * MRTKLIB's RTCM 3 decoder. Complete observation epochs are handed to the owner; ephemerides
 * are merged into the shared nav; the station record (1005/1006 position, 1033 receiver and
 * antenna, 1230 GLONASS biases) is kept for the sessions that will advertise it.
 *
 * Health is the only judgement made here: a feed that has delivered no epoch for `stale_s` is
 * unhealthy, and the daemon elects a replacement (the closest remaining station) and re-forms
 * the cluster. Nothing about correction quality is decided at this level -- D112 found no
 * stream statistic that predicts it.
 */

#include <time.h>

#include "rtklib.h"
#include "vbs/vbs_ntrip.h"

struct vbs_feed_opt {
	const char *host;
	int port;
	const char *user;
	const char *password;
	double stale_s;			/* no epoch for this long: unhealthy */
};

struct vbs_feed_stats {
	long epochs;
	long ephemerides;
	long station_msgs;
	long bytes;
	long decode_errors;
	double last_epoch_unix;		/* 0 until the first */
	gtime_t last_epoch_time;	/* GNSS time of the last epoch */
};

struct vbs_feed;

/*
 * Open a feed. `nav` is the shared ephemeris store, allocated index-style (see
 * vbs_feed_nav_init); every feed writes its ephemerides into it. Non-blocking.
 */
struct vbs_feed *vbs_feed_open(const struct vbs_feed_opt *opt, const char *mountpoint, nav_t *nav);
void vbs_feed_close(struct vbs_feed *this);

int vbs_feed_fd(const struct vbs_feed *this);
short vbs_feed_events(const struct vbs_feed *this);

/*
 * What the owner hears from a feed. `on_epoch` gets each complete observation epoch (valid only
 * during the call). `on_message` gets every complete, CRC-checked RTCM frame as it came off the
 * wire -- the bytes a rover would have received from this base directly -- for the pass-through
 * a tile serves until its own synthesis has converged. Either may be NULL.
 */
struct vbs_feed_sink {
	void (*on_epoch)(struct vbs_feed *feed, const obsd_t *obs, int n, void *arg);
	void (*on_message)(struct vbs_feed *feed, int type, const unsigned char *frame, int len, void *arg);
	void *arg;
};

/*
 * Advance the feed: read bytes, decode, deliver. Returns the number of epochs delivered, or -1
 * on a connection error (the feed reconnects by itself).
 */
int vbs_feed_pump(struct vbs_feed *this, double now_unix, const struct vbs_feed_sink *sink);

/*
 * The base's most recent station messages -- 1005 or 1006, 1033, 1230 -- as received, so a
 * rover attaching mid-stream gets them at once instead of waiting for the base's next period.
 * Copies into `buf`, returns the bytes, 0 if none seen yet.
 */
int vbs_feed_station_msgs(const struct vbs_feed *this, unsigned char *buf, int len);
int vbs_feed_staid(const struct vbs_feed *this);		/* the base's own DF003 */

/* 1: delivering; 0: stale, never delivered, or disconnected. */
int vbs_feed_healthy(const struct vbs_feed *this, double now_unix);

const sta_t *vbs_feed_station(const struct vbs_feed *this);	/* from 1005/1006/1033/1230 */
int vbs_feed_has_position(const struct vbs_feed *this);		/* a 1005/1006 was seen */
void vbs_feed_get_stats(const struct vbs_feed *this, struct vbs_feed_stats *out);

/*
 * Prepare a shared nav for index-style ephemeris storage (eph[sat-1 + MAXSAT*set],
 * geph[prn-1]), the layout MRTKLIB's decoder fills and its satposs() reads. Returns 0, or -1.
 */
int vbs_feed_nav_init(nav_t *nav);
void vbs_feed_nav_free(nav_t *nav);

/*
 * Feed one buffer of RTCM bytes without a socket -- what replay and the tests use, and the
 * same decode path the live pump takes.
 */
int vbs_feed_input(struct vbs_feed *this, const unsigned char *buf, int n, double now_unix,
	const struct vbs_feed_sink *sink);

/* A feed with no socket, for replay: `mountpoint` names it, `nav` as above. */
struct vbs_feed *vbs_feed_open_offline(const char *mountpoint, nav_t *nav, double stale_s);

#endif /* __VBS_FEED_H__ */

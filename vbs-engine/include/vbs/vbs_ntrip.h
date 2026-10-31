#ifndef __VBS_NTRIP_H__
#define __VBS_NTRIP_H__

#include <stdio.h>
#include <sys/types.h>

/*
 * NTRIP client sessions to the physical bases, with verbatim recording and replay --
 * `02-architecture.md` §6 and spec §8's replay requirement.
 *
 * The contract that matters is that **recording and replay are the same code path as live**.
 * A recorded stream is the bytes exactly as they arrived, plus an index of when each read
 * returned; replay hands those same byte groupings back to the same decoder at the same
 * relative times. Anything the engine does differently between live and replay is a bug the
 * validation cannot see, so there is deliberately no separate "file mode" decoder.
 *
 * Why raw bytes rather than decoded observations: RTCM carries state the decoder builds up
 * (ephemerides, reference station records, MSM signal masks). Recording decoded output would
 * freeze one decoder's interpretation into the archive, and D19/D20's open item -- whether
 * MRTKLIB reads Galileo E5b correctly from a live stream -- is precisely a question about the
 * decoder. That question can only be asked of an archive that kept the original bytes.
 *
 * The caster cannot record streams (`00-codebase-map.md` §10) and does not need to learn how:
 * this is an ordinary NTRIP client, so it works against any caster.
 */

struct vbs_ntrip_opt {
	const char *host;
	int port;
	const char *mountpoint;
	const char *user;		/* NULL: anonymous (an open caster needs no authentication) */
	const char *password;
	const char *record_dir;		/* NULL: do not record */
	int reconnect_delay_s;		/* 0: take the default */
};

/* Per-stream counters, for a recording run that has to be trusted hours later. */
struct vbs_ntrip_stats {
	long long bytes;
	long connects;
	long disconnects;
	long errors;
	time_t first_byte;		/* 0 until the stream actually delivers */
	time_t last_byte;
};

struct vbs_ntrip;

/*
 * Create a session. Does not block: the connection is established as vbs_ntrip_pump() runs,
 * so a caller can drive many streams from one poll loop.
 */
struct vbs_ntrip *vbs_ntrip_open(const struct vbs_ntrip_opt *opt);
void vbs_ntrip_close(struct vbs_ntrip *this);

/*
 * The socket to poll, or -1 while the session is waiting out a reconnect delay. Re-read it
 * after every vbs_ntrip_pump(): a reconnect replaces it.
 */
int vbs_ntrip_fd(const struct vbs_ntrip *this);

/* What to poll for -- POLLOUT while a non-blocking connect is in flight, POLLIN otherwise. */
short vbs_ntrip_events(const struct vbs_ntrip *this);

/*
 * Advance the session: finish a connect, send the request, parse the response header, read
 * stream bytes. Received payload is appended to the recording (if any) and passed to `sink`
 * when one is given -- the same bytes, in the same grouping, that replay will reproduce.
 *
 * Returns the number of payload bytes handled, 0 if nothing happened, -1 on an error that
 * closed the connection (the session schedules its own reconnect; this is not fatal).
 */
int vbs_ntrip_pump(struct vbs_ntrip *this, time_t now,
	void (*sink)(const char *mountpoint, const void *buf, size_t n, void *arg), void *arg);

void vbs_ntrip_get_stats(const struct vbs_ntrip *this, struct vbs_ntrip_stats *out);
const char *vbs_ntrip_mountpoint(const struct vbs_ntrip *this);

/*
 * Replay a recording written by vbs_ntrip_pump(). Calls `sink` with the same byte groupings
 * the live session saw. `speed` scales the inter-chunk delay: 0 replays as fast as possible
 * (what a test wants), 1.0 in real time (what a latency measurement wants).
 *
 * Returns the number of chunks replayed, or -1 if the recording could not be opened.
 */
long vbs_ntrip_replay(const char *record_dir, const char *mountpoint, double speed,
	void (*sink)(const char *mountpoint, const void *buf, size_t n, void *arg), void *arg);

#endif /* __VBS_NTRIP_H__ */

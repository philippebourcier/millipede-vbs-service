/*
 * The RTCM feed, on a real recorded stream.
 *
 * testdata/rtcm holds seven minutes of three live bases, recorded byte for
 * byte by tools/vbs_record. Replaying it through vbs_feed exercises exactly the decode path the
 * daemon uses live -- vbs_ntrip_replay hands the bytes over in the same groupings the socket
 * did -- so what this checks is that the live path turns a real stream into epochs, ephemerides
 * and a station record, at the rate the base actually sends them.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "vbs/vbs_feed.h"

#define RECORD_DIR	"testdata/rtcm"

static int nfail;

static void check(int ok, const char *what) {
	printf("  %-66s %s\n", what, ok ? "ok" : "FAIL");
	if (!ok)
		nfail++;
}

struct ctx {
	struct vbs_feed *feed;
	long epochs, sats, glo, gaps;
	long frames, frame_bytes, msm, pos, t1033, t1230, badtype;
	gtime_t last;
	double clock;
};

static void on_message(struct vbs_feed *feed, int type, const unsigned char *frame, int len, void *arg) {
	struct ctx *c = (struct ctx *)arg;

	(void)feed;
	c->frames++;
	c->frame_bytes += len;
	if (frame[0] != 0xd3 || len < 6)
		c->badtype++;
	if (type >= 1071 && type <= 1137)
		c->msm++;
	else if (type == 1005 || type == 1006)
		c->pos++;
	else if (type == 1033)
		c->t1033++;
	else if (type == 1230)
		c->t1230++;
}

static void on_epoch(struct vbs_feed *feed, const obsd_t *obs, int n, void *arg) {
	struct ctx *c = (struct ctx *)arg;
	int i;

	(void)feed;
	c->epochs++;
	c->sats += n;
	for (i = 0; i < n; i++)
		if (satsys(obs[i].sat, NULL) == SYS_GLO)
			c->glo++;
	if (c->epochs > 1 && fabs(timediff(obs[0].time, c->last) - 1.0) > 0.01)
		c->gaps++;
	c->last = obs[0].time;
}

static void sink(const char *mp, const void *buf, size_t n, void *arg) {
	struct ctx *c = (struct ctx *)arg;
	struct vbs_feed_sink fs = {on_epoch, on_message, c};

	(void)mp;
	c->clock += 0.001;
	vbs_feed_input(c->feed, (const unsigned char *)buf, (int)n, c->clock, &fs);
}

int main(void) {
	static nav_t nav;
	struct ctx c;
	struct vbs_feed_stats st;
	long chunks;
	const sta_t *sta;
	int i, neph = 0, ngeph = 0;

	printf("=== RTCM feed on the recorded T30 capture ===\n");
	if (access(RECORD_DIR "/STA_R.rtcm3", R_OK) != 0) {
		printf("  capture not present under %s\n", RECORD_DIR);
		printf("SKIPPED: recorded capture not available\n");
		return 0;
	}
	check(vbs_feed_nav_init(&nav) == 0, "shared nav allocated");
	memset(&c, 0, sizeof c);
	c.feed = vbs_feed_open_offline("STA_R", &nav, 30.0);
	check(c.feed != NULL, "offline feed created");
	chunks = vbs_ntrip_replay(RECORD_DIR, "STA_R", 0.0, sink, &c);
	vbs_feed_get_stats(c.feed, &st);
	printf("  chunks %ld, bytes %ld, epochs %ld (%.1f sats/epoch, %ld GLONASS obs), ephemerides %ld, "
	    "station msgs %ld, decode errors %ld, gaps %ld\n",
	    chunks, st.bytes, st.epochs, c.epochs ? (double)c.sats / c.epochs : 0.0, c.glo,
	    st.ephemerides, st.station_msgs, st.decode_errors, c.gaps);
	check(chunks > 1000, "the recording replayed");
	check(st.epochs > 300 && st.epochs < 500, "about seven minutes of 1 Hz epochs decoded");
	check(c.sats > st.epochs * 8, "a realistic satellite count per epoch");
	check(st.ephemerides > 0, "ephemerides arrived on the stream");
	check(st.station_msgs > 0 && vbs_feed_has_position(c.feed), "1005/1006 seen: the base declared its position");
	check(c.gaps < st.epochs / 20, "epochs are consecutive at 1 Hz (few gaps)");
	check(st.decode_errors == 0, "no decode errors on a clean recording");
	sta = vbs_feed_station(c.feed);
	printf("  station: pos %.1f %.1f %.1f rec '%s' ant '%s'\n", sta->pos[0], sta->pos[1], sta->pos[2],
	    sta->rectype, sta->antdes);
	check(fabs(sta->pos[0]) > 1e6, "station ECEF is plausible");
	for (i = 0; i < nav.n; i++)
		if (nav.eph[i].sat > 0)
			neph++;
	for (i = 0; i < nav.ng; i++)
		if (nav.geph[i].sat > 0)
			ngeph++;
	printf("  shared nav: %d ephemerides, %d GLONASS\n", neph, ngeph);
	check(neph > 10, "ephemerides merged into the shared nav, index-style");
	{
		unsigned char sm[256];
		int n = vbs_feed_station_msgs(c.feed, sm, sizeof sm);

		printf("  raw frames %ld (%ld bytes: %.1f%% of the stream), MSM %ld, 1005/6 %ld, 1033 %ld, 1230 %ld, staid %d\n",
		    c.frames, c.frame_bytes, st.bytes ? 100.0 * c.frame_bytes / st.bytes : 0.0,
		    c.msm, c.pos, c.t1033, c.t1230, vbs_feed_staid(c.feed));
		check(c.frames > 1000 && c.badtype == 0, "every complete frame delivered raw, preamble intact");
		check(c.frame_bytes > st.bytes * 0.98, "the raw frames account for the whole stream (nothing lost)");
		check(c.msm >= st.epochs * 3, "MSM frames: several per epoch");
		check(c.pos > 0 && n > 0 && sm[0] == 0xd3, "the base's station messages are cached for late joiners");
	}
	check(vbs_feed_healthy(c.feed, c.clock) == 1, "healthy right after the last epoch");
	check(vbs_feed_healthy(c.feed, c.clock + 60.0) == 0, "stale a minute later");
	vbs_feed_close(c.feed);
	vbs_feed_nav_free(&nav);

	printf("%s: %d failure(s)\n", nfail ? "FAILED" : "PASSED", nfail);
	return nfail != 0;
}

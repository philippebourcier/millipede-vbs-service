#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_feed.h"

struct vbs_feed {
	char mountpoint[64];
	struct vbs_ntrip *ntrip;	/* NULL for an offline feed */
	rtcm_t rtcm;
	nav_t *nav;
	double stale_s;
	int have_pos;
	struct vbs_feed_stats stats;
	/* the pump's sink, stashed so the byte reader can reach it */
	struct vbs_feed_sink sink;
	double now;
	int epochs_this_pump;
	/* the latest station messages, raw, for a late-joining rover */
	unsigned char msg_pos[64], msg_1033[128], msg_1230[64];
	int len_pos, len_1033, len_1230;
};

#define RTCM_MAXFRAME	1200

/* Message type of a complete frame: 12 bits after the 24-bit header. */
static int frame_type(const unsigned char *b) {
	return ((b[3] << 4) | (b[4] >> 4)) & 0xfff;
}

static void keep_station_msg(struct vbs_feed *this, int type, const unsigned char *f, int len) {
	if ((type == 1005 || type == 1006) && len <= (int)sizeof this->msg_pos) {
		memcpy(this->msg_pos, f, len);
		this->len_pos = len;
	} else if (type == 1033 && len <= (int)sizeof this->msg_1033) {
		memcpy(this->msg_1033, f, len);
		this->len_1033 = len;
	} else if (type == 1230 && len <= (int)sizeof this->msg_1230) {
		memcpy(this->msg_1230, f, len);
		this->len_1230 = len;
	}
}

int vbs_feed_nav_init(nav_t *nav) {
	int i;
	eph_t eph0;
	geph_t geph0;

	/* an empty entry: sat 0, so satposs() never selects it; iode/iodc -1 as MRTKLIB's own */
	memset(&eph0, 0, sizeof eph0);
	eph0.iode = eph0.iodc = -1;
	memset(&geph0, 0, sizeof geph0);
	geph0.iode = -1;
	memset(nav, 0, sizeof *nav);
	nav->eph = (eph_t *)calloc(MAXSAT * 2, sizeof(eph_t));
	nav->geph = (geph_t *)calloc(MAXPRNGLO > 0 ? MAXPRNGLO : 1, sizeof(geph_t));
	if (nav->eph == NULL || nav->geph == NULL) {
		vbs_feed_nav_free(nav);
		return -1;
	}
	for (i = 0; i < MAXSAT * 2; i++)
		nav->eph[i] = eph0;
	for (i = 0; i < MAXPRNGLO; i++)
		nav->geph[i] = geph0;
	nav->n = nav->nmax = MAXSAT * 2;
	nav->ng = nav->ngmax = MAXPRNGLO;
	return 0;
}

void vbs_feed_nav_free(nav_t *nav) {
	free(nav->eph);
	free(nav->geph);
	nav->eph = NULL;
	nav->geph = NULL;
	nav->n = nav->ng = 0;
}

static struct vbs_feed *feed_new(const char *mountpoint, nav_t *nav, double stale_s) {
	struct vbs_feed *this = (struct vbs_feed *)calloc(1, sizeof *this);

	if (this == NULL)
		return NULL;
	snprintf(this->mountpoint, sizeof this->mountpoint, "%s", mountpoint);
	if (!init_rtcm(&this->rtcm)) {
		free(this);
		return NULL;
	}
	this->nav = nav;
	this->stale_s = stale_s > 0.0 ? stale_s : 30.0;
	return this;
}

struct vbs_feed *vbs_feed_open(const struct vbs_feed_opt *opt, const char *mountpoint, nav_t *nav) {
	struct vbs_feed *this = feed_new(mountpoint, nav, opt->stale_s);
	struct vbs_ntrip_opt no;

	if (this == NULL)
		return NULL;
	memset(&no, 0, sizeof no);
	no.host = opt->host;
	no.port = opt->port;
	no.mountpoint = this->mountpoint;
	no.user = opt->user;
	no.password = opt->password;
	no.record_dir = NULL;
	no.reconnect_delay_s = 0;
	if ((this->ntrip = vbs_ntrip_open(&no)) == NULL) {
		free_rtcm(&this->rtcm);
		free(this);
		return NULL;
	}
	return this;
}

struct vbs_feed *vbs_feed_open_offline(const char *mountpoint, nav_t *nav, double stale_s) {
	return feed_new(mountpoint, nav, stale_s);
}

void vbs_feed_close(struct vbs_feed *this) {
	if (this == NULL)
		return;
	if (this->ntrip)
		vbs_ntrip_close(this->ntrip);
	free_rtcm(&this->rtcm);
	free(this);
}

int vbs_feed_fd(const struct vbs_feed *this) {
	return this->ntrip ? vbs_ntrip_fd(this->ntrip) : -1;
}

short vbs_feed_events(const struct vbs_feed *this) {
	return this->ntrip ? vbs_ntrip_events(this->ntrip) : 0;
}

/*
 * Merge one ephemeris into the shared store. The decoder has written it into its own nav at
 * the index-style position; copy that one entry across. Other feeds may carry the same
 * satellite -- the latest arrival wins, which is what a receiver does too.
 */
static void merge_eph(struct vbs_feed *this) {
	int sat = this->rtcm.ephsat, prn, set = this->rtcm.ephset;

	if (sat <= 0 || sat > MAXSAT)
		return;
	if (satsys(sat, &prn) == SYS_GLO) {
		if (prn >= 1 && prn <= MAXPRNGLO && this->rtcm.nav.geph != NULL)
			this->nav->geph[prn - 1] = this->rtcm.nav.geph[prn - 1];
		return;
	}
	if (set < 0 || set > 1)
		set = 0;
	if (this->rtcm.nav.eph != NULL)
		this->nav->eph[sat - 1 + MAXSAT * set] = this->rtcm.nav.eph[sat - 1 + MAXSAT * set];
}

static void handle_ret(struct vbs_feed *this, int ret) {
	switch (ret) {
	case 1:				/* observation epoch complete */
		this->stats.epochs++;
		this->stats.last_epoch_unix = this->now;
		this->stats.last_epoch_time = this->rtcm.time;
		this->epochs_this_pump++;
		if (this->sink.on_epoch && this->rtcm.obs.n > 0)
			this->sink.on_epoch(this, this->rtcm.obs.data, this->rtcm.obs.n, this->sink.arg);
		break;
	case 2:				/* ephemeris */
		this->stats.ephemerides++;
		merge_eph(this);
		break;
	case 5:				/* station: 1005/1006/1033/1230 and friends */
		this->stats.station_msgs++;
		if (this->rtcm.sta.pos[0] != 0.0 || this->rtcm.sta.pos[1] != 0.0)
			this->have_pos = 1;
		break;
	case -1:
		this->stats.decode_errors++;
		break;
	default:
		break;
	}
}

/*
 * One byte through the decoder, watching for the frame boundary. input_rtcm3() zeroes nbyte
 * when a frame is complete (parity checked or not); a frame that decoded to nothing -- a type
 * the decoder does not handle -- still returns 0, so completeness is read from the byte count,
 * not from the return value, and the parity is re-checked here before the frame is passed on
 * as the base's own bytes.
 */
static void feed_byte(struct vbs_feed *this, unsigned char b) {
	int had = this->rtcm.nbyte, ret;

	ret = input_rtcm3(&this->rtcm, b);
	if (had > 0 && this->rtcm.nbyte == 0) {
		/* rtcm.len is header + payload (the decoder adds the 3 header bytes); parity follows */
		int len = this->rtcm.len;

		if (len >= 6 && len + 3 <= RTCM_MAXFRAME &&
		    rtk_crc24q(this->rtcm.buff, len) == getbitu(this->rtcm.buff, len * 8, 24)) {
			int type = frame_type(this->rtcm.buff);

			keep_station_msg(this, type, this->rtcm.buff, len + 3);
			if (this->sink.on_message)
				this->sink.on_message(this, type, this->rtcm.buff, len + 3, this->sink.arg);
		}
	}
	handle_ret(this, ret);
}

int vbs_feed_input(struct vbs_feed *this, const unsigned char *buf, int n, double now_unix,
    const struct vbs_feed_sink *sink) {
	int i;

	this->sink = *sink;
	this->now = now_unix;
	this->epochs_this_pump = 0;
	this->stats.bytes += n;
	for (i = 0; i < n; i++)
		feed_byte(this, buf[i]);
	return this->epochs_this_pump;
}

static void byte_sink(const char *mountpoint, const void *buf, size_t n, void *arg) {
	struct vbs_feed *this = (struct vbs_feed *)arg;
	const unsigned char *b = (const unsigned char *)buf;
	size_t i;

	(void)mountpoint;
	this->stats.bytes += (long)n;
	for (i = 0; i < n; i++)
		feed_byte(this, b[i]);
}

int vbs_feed_pump(struct vbs_feed *this, double now_unix, const struct vbs_feed_sink *sink) {
	int r;

	if (this->ntrip == NULL)
		return 0;
	this->sink = *sink;
	this->now = now_unix;
	this->epochs_this_pump = 0;
	r = vbs_ntrip_pump(this->ntrip, (time_t)now_unix, byte_sink, this);
	if (r < 0)
		return -1;
	return this->epochs_this_pump;
}

int vbs_feed_station_msgs(const struct vbs_feed *this, unsigned char *buf, int len) {
	int used = 0;

	if (this->len_pos && used + this->len_pos <= len) {
		memcpy(buf + used, this->msg_pos, this->len_pos);
		used += this->len_pos;
	}
	if (this->len_1033 && used + this->len_1033 <= len) {
		memcpy(buf + used, this->msg_1033, this->len_1033);
		used += this->len_1033;
	}
	if (this->len_1230 && used + this->len_1230 <= len) {
		memcpy(buf + used, this->msg_1230, this->len_1230);
		used += this->len_1230;
	}
	return used;
}

int vbs_feed_staid(const struct vbs_feed *this) { return this->rtcm.staid; }

int vbs_feed_healthy(const struct vbs_feed *this, double now_unix) {
	if (this->stats.last_epoch_unix <= 0.0)
		return 0;
	return now_unix - this->stats.last_epoch_unix <= this->stale_s;
}

const sta_t *vbs_feed_station(const struct vbs_feed *this) { return &this->rtcm.sta; }
int vbs_feed_has_position(const struct vbs_feed *this) { return this->have_pos; }

void vbs_feed_get_stats(const struct vbs_feed *this, struct vbs_feed_stats *out) {
	*out = this->stats;
}

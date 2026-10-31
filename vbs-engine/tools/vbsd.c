/*
 * vbsd -- the VBS caster (spec §6b).
 *
 * An NTRIP server addressed by position. A client (in practice the millipede caster, acting for
 * a rover) asks for
 *
 *	GET /<LAT>/<LON>/<ALT>
 *
 * and then sends the rover's GGA on the same connection. This daemon answers 404 when the
 * serving rules refuse the position -- no cluster, area already covered, all-F9P, mixed,
 * out of range -- and otherwise streams RTCM synthesized for the rover's TILE: a cell of the
 * fixed global grid, with a station ID derived from the cell, shared by every rover on it.
 * Tile changes are read from the GGA and handled here; the client never reconnects for movement.
 *
 * Structure, all in one poll loop:
 *
 *   feeds     one per physical base in use, shared by every tile that uses that base
 *   sessions  one per occupied tile: the Phase 1 pipeline, driven by the feeds
 *   clients   attached to a session; receive its RTCM as it is produced
 *
 * Base health: a feed that stops delivering is dropped from the cache's live view, and every
 * session using it is re-formed from the rules with the closest remaining stations. Nothing
 * here judges correction quality (D112 found no statistic that predicts it); the rules are a
 * prior and the report says so.
 *
 * Usage:
 *   vbsd --cache base_cache.json --coords COORDINATES.json --caster HOST[:PORT]
 *        [--user U --password P] [--listen PORT] [--hysteresis-m 300] [--stale-s 30]
 *        [--idle-s 120] [--check-bands] [--no-backup] [--force] [--exclude-km R]
 *   vbsd --cache base_cache.json --decide LAT LON        # print the decision and exit
 *
 * --force and --exclude-km are for field tests only (docs/vbs/field-test.md), never production.
 */
#include <arpa/inet.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <poll.h>
#include <signal.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "vbs/vbs_basecache.h"
#include "vbs/vbs_basemeta.h"
#include "vbs/vbs_c2wcal.h"
#include "vbs/vbs_cluster.h"
#include "vbs/vbs_feed.h"
#include "vbs/vbs_gga.h"
#include "vbs/vbs_policy.h"
#include "vbs/vbs_session.h"
#include "vbs/vbs_tile.h"

#define MAXCLIENTS	256
#define MAXSESSIONS	64
#define MAXFEEDS	256
#define REQBUF		2048
#define MAXALLOW	16

static volatile sig_atomic_t stop, reload;

static void on_signal(int sig) {
	if (sig == SIGHUP)
		reload = 1;
	else
		stop = 1;
}

static double now_unix(void) {
	struct timeval tv;

	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1e-6;
}

static void logmsg(const char *fmt, ...) {
	char stamp[32];
	time_t t = time(NULL);
	va_list ap;

	strftime(stamp, sizeof stamp, "%Y-%m-%d %H:%M:%S", gmtime(&t));
	fprintf(stderr, "%s vbsd ", stamp);
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fputc('\n', stderr);
	fflush(stderr);
}

/* "off", "on" or "on (EC)": the constellations carrying the dispersive correction (D167) */
static const char *disp_str(int dispersive, int iono_sys) {
	static char buf[24];

	if (!dispersive)
		return "off";
	if (!iono_sys)
		return "on";
	snprintf(buf, sizeof buf, "on (%s%s%s%s%s)", iono_sys & SYS_GPS ? "G" : "",
	    iono_sys & SYS_GLO ? "R" : "", iono_sys & SYS_GAL ? "E" : "",
	    iono_sys & SYS_CMP ? "C" : "", iono_sys & SYS_QZS ? "J" : "");
	return buf;
}

/* ---- state -------------------------------------------------------------------------------- */

struct feedref {
	struct vbs_feed *feed;
	char mp[64];
	int cache_idx;
	int refs;
	double idle_since, opened;
	int pos_checked;		/* the streamed 1005/1006 has been compared to the cache */
	/*
	 * --c2w-cal (vbs_c2wcal.h): this base's own 2W-2X, its trailing estimate, and the verdict
	 * last logged (+2, 0 = never), so a change is said once. NULL unless the flag is on and the
	 * base is a UM980.
	 */
	struct vbs_c2wcal_base *c2w;
	int c2w_logged;
};

struct tilesess {
	int active;
	struct vbs_tile tile;
	char name[16];
	int staid;
	double xyz[3];
	struct vbs_session *sess;
	int nsta, master;
	char mp[VBS_SESSION_MAXSTA][64];
	int cache_idx[VBS_SESSION_MAXSTA];
	struct feedref *feed[VBS_SESSION_MAXSTA];
	struct feedref *pt;		/* pass-through: the nearest live base, served until converged */
	int family, dispersive, prefer_l2c, glonass, iono_sys;
	int nclients;
	int was_converged;
	double idle_since, created;
	long rebuilds;
};

struct client {
	int fd;
	int state;			/* 0 request, 1 streaming */
	int ntrip1;
	char in[REQBUF];
	int inlen;
	double lat, lon, h;
	int have_pos;
	struct vbs_tile tile;
	int have_tile;
	struct tilesess *sess;
	double connected, last_gga;
	char peer[64];
	long dropped;			/* packets not sent: socket full */
	long bytes;
	/*
	 * The unsent tail of a partial write. A non-blocking send() on a filling socket returns a
	 * short count well before it returns EAGAIN, so without this every slow rover was cut on
	 * "send error" and the drop policy below never ran. Flushed first on the next send; while
	 * it is pending, new data is dropped whole, so the stream stays frame-aligned.
	 */
	int pendlen;
	unsigned char pend[VBS_SESSION_OUTBUF];
};

struct daemon {
	/* config */
	const char *cache_path, *caster_host, *user, *password;
	/*
	 * --coords: the operator's reference-coordinate document, applied over the cache at startup
	 * (D146). --coords-max-m bounds how far a base may be moved by it; the file this was built
	 * for carries two entries 160 and 215 km out, and a mountpoint mismatch would look exactly
	 * the same, so an unbounded override is not an option. 5 m is far beyond any real survey
	 * difference and far short of a teleport.
	 */
	const char *coords_path;
	double coords_max_m;
	int caster_port, listen_port, check_bands;
	double hysteresis_m, stale_s, idle_s, feed_idle_s, connect_grace_s;
	int max_aux;
	int min_aux;		/* cluster admission floor; 0 = the engine's default */
	/*
	 * --no-backup: the engine as measured before "2 + 1 backup" -- no reserve in the policy or
	 * the cluster, and a fit floor of 3 auxiliaries per satellite. For running the rest of this
	 * build while the reserve's leave-one-out results are still out.
	 */
	int no_backup;
	/*
	 * Field-test mode, never production (docs/vbs/field-test.md). --exclude-km hides every base
	 * within that distance of the tile centre, as if the area had a coverage gap that size;
	 * --force serves wherever a homogeneous Mosaic or UM980 cluster can form (VBS_POLICY_FORCE),
	 * and each session's log line keeps what production would have answered.
	 */
	double exclude_km;
	int force;
	/*
	 * --c2w-cal run:SEC | l2c:SEC (D167): each UM980 base's GPS L2 P(Y) code brought to its own
	 * L2C level, from a trailing SEC-second estimate of its own 2W-2X, once per base epoch before
	 * the epoch reaches any session -- so every session sees the same corrected observations, as
	 * the offline pipeline does. Gated per firmware build (c2w_groups). Off by default; the
	 * acausal `day` mode of vbs_pipeline cannot run live and is not offered.
	 */
	int c2w_mode;			/* 0 off, 1 constant, 2 l2c */
	double c2w_win;
	struct vbs_c2wcal_groups c2w_groups;
	/*
	 * Who may connect. Rovers are authenticated by the millipede caster, not here (spec §6b:
	 * "NTRIP client credentials like any other client of the millipede caster"), so this port
	 * should see only that caster. `--allow` lists the addresses or prefixes it may come from;
	 * empty means anyone, which is right for a test and wrong for a deployment.
	 */
	struct {
		unsigned char addr[16];
		int family, bits;
	} allow[MAXALLOW];
	int nallow;
	long denied;
	/* state */
	struct vbs_basecache bc;
	struct vbs_policy_base *pb;
	nav_t nav;
	struct feedref feeds[MAXFEEDS];
	int nfeeds;
	struct tilesess sess[MAXSESSIONS];
	struct client cl[MAXCLIENTS];
	int ncl;
	int lfd;
	double last_housekeeping, last_status, started;
	long served, refused, connections;
	long refused_by_rule[8];
};

/* ---- access ------------------------------------------------------------------------------- */

/* "1.2.3.4", "1.2.3.0/24", "2001:db8::/32" */
static int allow_add(struct daemon *d, const char *spec) {
	char buf[80], *slash;
	int bits = -1;

	if (d->nallow >= MAXALLOW)
		return -1;
	snprintf(buf, sizeof buf, "%s", spec);
	if ((slash = strchr(buf, '/')) != NULL) {
		*slash = '\0';
		bits = atoi(slash + 1);
	}
	if (inet_pton(AF_INET, buf, d->allow[d->nallow].addr) == 1) {
		d->allow[d->nallow].family = AF_INET;
		d->allow[d->nallow].bits = bits < 0 || bits > 32 ? 32 : bits;
	} else if (inet_pton(AF_INET6, buf, d->allow[d->nallow].addr) == 1) {
		d->allow[d->nallow].family = AF_INET6;
		d->allow[d->nallow].bits = bits < 0 || bits > 128 ? 128 : bits;
	} else
		return -1;
	d->nallow++;
	return 0;
}

static int prefix_match(const unsigned char *a, const unsigned char *b, int bits) {
	int full = bits / 8, rest = bits % 8;

	if (memcmp(a, b, full) != 0)
		return 0;
	if (rest == 0)
		return 1;
	return ((a[full] ^ b[full]) & (0xff << (8 - rest))) == 0;
}

static int allowed(const struct daemon *d, const struct sockaddr_in6 *sa) {
	const unsigned char *v4 = NULL, *v6 = NULL;
	int i;

	if (d->nallow == 0)
		return 1;
	if (sa->sin6_family == AF_INET6) {
		v6 = sa->sin6_addr.s6_addr;
		/* a v4 peer on the dual-stack socket arrives as ::ffff:a.b.c.d */
		if (IN6_IS_ADDR_V4MAPPED(&sa->sin6_addr))
			v4 = v6 + 12;
	} else
		v4 = (const unsigned char *)&((const struct sockaddr_in *)sa)->sin_addr;
	for (i = 0; i < d->nallow; i++) {
		if (d->allow[i].family == AF_INET && v4 && prefix_match(d->allow[i].addr, v4, d->allow[i].bits))
			return 1;
		if (d->allow[i].family == AF_INET6 && v6 && prefix_match(d->allow[i].addr, v6, d->allow[i].bits))
			return 1;
	}
	return 0;
}

/* ---- feeds -------------------------------------------------------------------------------- */

static struct feedref *feed_get(struct daemon *d, int cache_idx, double now) {
	int i;
	struct vbs_feed_opt fo;

	for (i = 0; i < d->nfeeds; i++)
		if (d->feeds[i].feed && d->feeds[i].cache_idx == cache_idx) {
			d->feeds[i].refs++;
			return &d->feeds[i];
		}
	for (i = 0; i < MAXFEEDS; i++)
		if (d->feeds[i].feed == NULL)
			break;
	if (i == MAXFEEDS)
		return NULL;
	memset(&fo, 0, sizeof fo);
	fo.host = d->caster_host;
	fo.port = d->caster_port;
	fo.user = d->user;
	fo.password = d->password;
	fo.stale_s = d->stale_s;
	d->feeds[i].feed = vbs_feed_open(&fo, d->bc.e[cache_idx].mp, &d->nav);
	if (d->feeds[i].feed == NULL)
		return NULL;
	snprintf(d->feeds[i].mp, sizeof d->feeds[i].mp, "%s", d->bc.e[cache_idx].mp);
	d->feeds[i].cache_idx = cache_idx;
	d->feeds[i].refs = 1;
	d->feeds[i].opened = now;
	d->feeds[i].c2w_logged = 0;
	if (d->c2w_mode && d->bc.e[cache_idx].family == VBS_FAM_UM980) {
		if (d->feeds[i].c2w == NULL)
			d->feeds[i].c2w = malloc(sizeof *d->feeds[i].c2w);
		if (d->feeds[i].c2w != NULL)
			vbs_c2wcal_base_init(d->feeds[i].c2w, d->c2w_mode, d->c2w_win);
	}
	if (i >= d->nfeeds)
		d->nfeeds = i + 1;
	logmsg("feed %s opened (%s)", d->feeds[i].mp, d->bc.e[cache_idx].rec);
	return &d->feeds[i];
}

static void feed_put(struct feedref *f, double now) {
	if (f == NULL)
		return;
	if (--f->refs <= 0) {
		f->refs = 0;
		f->idle_since = now;
	}
}

/* ---- sessions ----------------------------------------------------------------------------- */

static double dist_m(const double *a, const double *b) {
	return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
	    (a[2] - b[2]) * (a[2] - b[2]));
}

/*
 * The base to pass through until the tile converges: the nearest one that is actually
 * delivering. Open feeds that are healthy come first (the cluster's own members are opening
 * anyway, so within a second or two one of them is); failing that, the nearest base the live
 * view still admits, opened now. The registry's `live` is an opinion -- a base can answer 404 on
 * the caster with `live: true` -- so a choice made from it alone is re-made every housekeeping
 * tick until the tile converges.
 */
static struct feedref *choose_passthrough(struct daemon *d, const double *xyz, double now,
    struct feedref *avoid) {
	int i, best = -1;
	double bestd = 0.0;

	for (i = 0; i < d->nfeeds; i++) {
		struct feedref *f = &d->feeds[i];
		double dd;

		if (f->feed == NULL || f == avoid || !vbs_feed_healthy(f->feed, now))
			continue;
		dd = dist_m(xyz, d->bc.e[f->cache_idx].xyz);
		if (d->exclude_km > 0.0 && dd <= d->exclude_km * 1000.0)
			continue;		/* hidden by the field test, as in a real gap */
		if (best < 0 || dd < bestd) {
			best = i;
			bestd = dd;
		}
	}
	if (best >= 0) {
		d->feeds[best].refs++;
		return &d->feeds[best];
	}
	for (i = 0; i < d->bc.n; i++) {
		double dd;

		if (!d->pb[i].live || d->bc.e[i].stream_ok < 0)
			continue;
		if (avoid && avoid->cache_idx == i)
			continue;
		dd = dist_m(xyz, d->bc.e[i].xyz);
		if (d->exclude_km > 0.0 && dd <= d->exclude_km * 1000.0)
			continue;
		if (best < 0 || dd < bestd) {
			best = i;
			bestd = dd;
		}
	}
	return best >= 0 ? feed_get(d, best, now) : NULL;
}

static void session_detach_all(struct daemon *d, struct tilesess *ts, int close_clients);

static void session_free(struct daemon *d, struct tilesess *ts, double now) {
	int k;

	if (!ts->active)
		return;
	vbs_session_free(ts->sess);
	ts->sess = NULL;
	for (k = 0; k < ts->nsta; k++)
		feed_put(ts->feed[k], now);
	feed_put(ts->pt, now);
	logmsg("tile %s: session closed", ts->name);
	memset(ts, 0, sizeof *ts);
}

/*
 * The whole serve decision for one tile: the rules, then the engine's own gate (frame
 * homogeneity) and geometry (hull, min_aux, baselines). Fills `xyz`, `dec`, and on success
 * `gated`/`ngated` and `cl`. Returns 1 to serve, 0 refused with the rule and reason in `dec`.
 *
 * Both the server and --decide go through here. They used to diverge: --decide stopped after
 * vbs_policy_decide(), so it answered "do the rules admit this position?" while the server
 * answers "can a cluster actually be formed?". It therefore printed SERVE for positions the
 * server then refused -- e.g. a cross-border tile where the rules find four homogeneous
 * receivers but one of them sits in a different reference frame, so the gate drops it and the
 * cluster is one auxiliary short. The tool must not be a weaker test than the server.
 */
static int decide_once(struct daemon *d, const double *xyz, int no_backup, int force,
    struct vbs_decision *dec, struct vbs_base *gated, int *ngated, struct vbs_cluster *cl);

/*
 * --exclude-km: the bases the field test hides from a tile -- every one in the live view within
 * that distance of the tile centre -- nearest first. `idx` has room for d->bc.n entries.
 */
static int hidden_bases(const struct daemon *d, const double *xyz, int *idx) {
	int i, j, n = 0;

	if (d->exclude_km <= 0.0)
		return 0;
	for (i = 0; i < d->bc.n; i++) {
		double km = dist_m(xyz, d->bc.e[i].xyz) / 1000.0;

		if (!d->pb[i].live || km > d->exclude_km)
			continue;
		for (j = n++; j > 0 && dist_m(xyz, d->bc.e[idx[j - 1]].xyz) / 1000.0 > km; j--)
			idx[j] = idx[j - 1];
		idx[j] = i;
	}
	return n;
}

/* The rules on the network a field test leaves visible: hidden bases drop out for this call only. */
static void policy_visible(struct daemon *d, const double *xyz, int flags, struct vbs_decision *dec) {
	int *hid = NULL, nh = 0, k;

	if (d->exclude_km > 0.0) {
		if ((hid = (int *)malloc(sizeof *hid * (d->bc.n > 0 ? d->bc.n : 1))) == NULL) {
			memset(dec, 0, sizeof *dec);
			snprintf(dec->rule, sizeof dec->rule, "R0");
			snprintf(dec->reason, sizeof dec->reason, "out of memory hiding bases");
			return;
		}
		nh = hidden_bases(d, xyz, hid);
		for (k = 0; k < nh; k++)
			d->pb[hid[k]].live = 0;
	}
	vbs_policy_decide(xyz, d->pb, d->bc.n, flags, dec);
	for (k = 0; k < nh; k++)
		d->pb[hid[k]].live = 1;		/* hidden_bases() only takes live ones */
	free(hid);
}

static int decide_tile(struct daemon *d, const struct vbs_tile *tile, double h, double *xyz,
    int force, struct vbs_decision *dec, struct vbs_base *gated, int *ngated,
    struct vbs_cluster *cl) {
	struct vbs_decision dec2;
	struct vbs_base gated2[VBS_POLICY_MAXCLUSTER];
	struct vbs_cluster cl2;
	int ngated2, k;

	vbs_tile_xyz(tile, h, xyz);
	if (decide_once(d, xyz, d->no_backup, force, dec, gated, ngated, cl))
		return 1;
	/*
	 * "2 + 1 backup" must never serve less than the engine without it. The reserve can make a
	 * short family eligible, and when that family then fails the frame gate the tile is refused
	 * although the other family, chosen without the reserve, would have served (T1003C0977:
	 * three UM980s, two of them in CHTRS95, displaced a Mosaic cluster of five). So a refusal
	 * with the reserve is asked again without it, and that answer is kept if it serves.
	 */
	if (d->no_backup || !decide_once(d, xyz, 1, force, &dec2, gated2, &ngated2, &cl2))
		return 0;
	*dec = dec2;
	*cl = cl2;
	*ngated = ngated2;
	for (k = 0; k < ngated2; k++)
		gated[k] = gated2[k];
	return 1;
}

static int decide_once(struct daemon *d, const double *xyz, int no_backup, int force,
    struct vbs_decision *dec, struct vbs_base *gated, int *ngated, struct vbs_cluster *cl) {
	struct vbs_cluster_opt copt;
	gtime_t t0 = {0, 0.0};
	int i, n, r;

	*ngated = 0;
	policy_visible(d, xyz, (d->check_bands ? VBS_POLICY_CHECK_BANDS : 0) |
	    (no_backup ? VBS_POLICY_NO_BACKUP : 0) | (force ? VBS_POLICY_FORCE : 0), dec);
	if (!dec->serve)
		return 0;

	n = dec->ncluster;
	for (i = 0; i < n; i++) {
		const struct vbs_basecache_entry *e = &d->bc.e[dec->cluster[i]];

		vbs_base_init(&gated[i], e->mp, e->frame, e->xyz);
		vbs_base_check(&gated[i], d->bc.e[dec->cluster[0]].frame, t0, NULL);
	}
	*ngated = n;
	vbs_cluster_default_opt(&copt);
	if (d->max_aux > 0)
		copt.max_aux = d->max_aux;
	if (d->min_aux > 0)
		copt.min_aux = d->min_aux;
	if (no_backup)
		copt.backup_min_aux = 0;
	copt.max_extrapolation_km = 0.0;	/* the rules already bound the distance */
	/* the policy may offer a backup from beyond its radius: never as an auxiliary */
	copt.aux_max_vbs_km = vbs_policy_cluster_radius_km();
	if ((r = vbs_cluster_select(xyz, gated, n, NULL, &copt, cl)) != VBS_CLUSTER_OK) {
		dec->serve = 0;
		snprintf(dec->rule, sizeof dec->rule, "R0");
		snprintf(dec->reason, sizeof dec->reason, "cluster: %s", vbs_cluster_reason_str(r));
		return 0;
	}
	return 1;
}

/*
 * What a field-test decision differs by, in one line: the answer production would give for the
 * same tile and visible network (cluster stage included), and the bases hidden from it. Empty
 * outside test mode.
 */
static void test_describe(struct daemon *d, const struct vbs_tile *tile, double h, char *buf,
    int len) {
	struct vbs_decision prod;
	struct vbs_base gated[VBS_POLICY_MAXCLUSTER];
	struct vbs_cluster cl;
	double xyz[3];
	int ngated, served, *hid, nh, k, o = 0;

	buf[0] = '\0';
	if (d->force) {
		served = decide_tile(d, tile, h, xyz, 0, &prod, gated, &ngated, &cl);
		o = snprintf(buf, len, "production %s %s: %s", prod.rule, served ? "SERVE" : "404",
		    prod.reason);
	}
	if (d->exclude_km <= 0.0 || o >= len - 1)
		return;
	vbs_tile_xyz(tile, h, xyz);
	if ((hid = (int *)malloc(sizeof *hid * (d->bc.n > 0 ? d->bc.n : 1))) == NULL)
		return;
	nh = hidden_bases(d, xyz, hid);
	o += snprintf(buf + o, len - o, "%shidden %d within %.0f km", o ? "; " : "", nh, d->exclude_km);
	for (k = 0; k < nh && o < len - 1; k++)
		o += snprintf(buf + o, len - o, "%s%s %.1f", k ? ", " : ": ", d->bc.e[hid[k]].mp,
		    dist_m(xyz, d->bc.e[hid[k]].xyz) / 1000.0);
	free(hid);
}

/*
 * Decide for a tile and, if the rules serve it, form the cluster. Returns 1 with `ts` built,
 * 0 with the rule in `dec` when refused, -1 on a resource failure.
 */
static int session_build(struct daemon *d, struct tilesess *ts, const struct vbs_tile *tile,
    struct vbs_decision *dec, double now) {
	struct vbs_base gated[VBS_POLICY_MAXCLUSTER];
	struct vbs_cluster cl;
	struct vbs_session_station sta[VBS_SESSION_MAXSTA];
	struct vbs_session_opt so;
	double h, xyz[3];
	int i, k, ngated, nprimary;

	/* the live view: a base whose stream has been tried is judged by its stream */
	for (i = 0; i < d->bc.n; i++) {
		int tried = d->bc.e[i].last_epoch_unix > 0.0 || d->bc.e[i].stream_ok < 0;

		d->pb[i].live = !d->bc.e[i].untrusted &&
		    (tried ? (d->bc.e[i].stream_ok > 0) : d->bc.e[i].live);
	}
	if (vbs_basecache_height_at(&d->bc, tile->clat, tile->clon,
	    vbs_policy_cluster_radius_km(), &h) < 0) {
		memset(dec, 0, sizeof *dec);
		snprintf(dec->rule, sizeof dec->rule, "R0");
		snprintf(dec->reason, sizeof dec->reason, "no base within reach");
		return 0;
	}
	if (!decide_tile(d, tile, h, xyz, d->force, dec, gated, &ngated, &cl))
		return 0;

	memset(ts, 0, sizeof *ts);
	ts->tile = *tile;
	vbs_tile_name(tile, ts->name, sizeof ts->name);
	ts->staid = vbs_tile_staid(tile);
	memcpy(ts->xyz, xyz, sizeof xyz);
	ts->family = dec->family;
	ts->dispersive = dec->dispersive;
	ts->prefer_l2c = dec->prefer_l2c;
	ts->iono_sys = dec->iono_sys;
	ts->glonass = dec->glonass;
	/* master first, then the auxiliaries the cluster module kept */
	ts->nsta = 0;
	ts->cache_idx[ts->nsta] = dec->cluster[cl.master];
	ts->nsta++;
	for (k = 0; k < cl.naux && ts->nsta < VBS_SESSION_MAXSTA; k++)
		ts->cache_idx[ts->nsta++] = dec->cluster[cl.aux[k]];
	nprimary = ts->nsta;
	/* the "2 + 1 backup" reserve, last, flagged so the interpolator holds it back */
	if (cl.backup >= 0 && ts->nsta < VBS_SESSION_MAXSTA)
		ts->cache_idx[ts->nsta++] = dec->cluster[cl.backup];
	ts->master = 0;
	for (k = 0; k < ts->nsta; k++) {
		const struct vbs_basecache_entry *e = &d->bc.e[ts->cache_idx[k]];

		sta[k].backup = k >= nprimary;	/* sta[] is not zeroed: set it for every entry */
		snprintf(ts->mp[k], sizeof ts->mp[k], "%s", e->mp);
		snprintf(sta[k].name, sizeof sta[k].name, "%.31s", e->mp);
		memcpy(sta[k].xyz, e->xyz, sizeof sta[k].xyz);
		snprintf(sta[k].rectype, sizeof sta[k].rectype, "%.31s", e->rec);
		if ((ts->feed[k] = feed_get(d, ts->cache_idx[k], now)) == NULL) {
			for (i = 0; i < k; i++)
				feed_put(ts->feed[i], now);
			return -1;
		}
	}
	vbs_session_default_opt(&so);
	so.sys_mask = SYS_GPS | SYS_GAL | SYS_CMP | (ts->glonass ? SYS_GLO : 0);
	so.dispersive = ts->dispersive;
	so.prefer_l2c = ts->prefer_l2c;
	so.iono_sys_mask = ts->iono_sys;
	if (d->no_backup)
		so.interp_min_aux = 3;
	so.staid = ts->staid;
	ts->sess = vbs_session_new(xyz, sta, ts->nsta, ts->master, &so, &d->nav);
	if (ts->sess == NULL) {
		for (k = 0; k < ts->nsta; k++)
			feed_put(ts->feed[k], now);
		return -1;
	}
	if (vbs_feed_has_position(ts->feed[0]->feed))
		vbs_session_set_master_station(ts->sess, vbs_feed_station(ts->feed[0]->feed));
	ts->pt = choose_passthrough(d, xyz, now, NULL);
	if (ts->pt == NULL)
		ts->pt = feed_get(d, ts->cache_idx[0], now);
	ts->active = 1;
	ts->created = now;
	ts->idle_since = now;
	/*
	 * No cluster mixes reference and registry coordinates (D148): every base the coordinate
	 * document does not position is untrusted, so it never reaches a cluster.
	 */
	logmsg("tile %s (staid %d, %.4f %.4f h %.0f m): %s, master %s at %.1f km, %d aux, "
	    "dispersive %s, L2C %s, GLONASS %s, %s hull, leverage %.2f; pass-through from %s until "
	    "converged",
	    ts->name, ts->staid, tile->clat, tile->clon, h, vbs_policy_family_str(ts->family),
	    ts->mp[0], cl.master_dist_km, ts->nsta - 1, disp_str(ts->dispersive, ts->iono_sys),
	    ts->prefer_l2c ? "preferred" : "not preferred", ts->glonass ? "on" : "off", cl.inside_hull ? "inside" : "OUTSIDE", cl.leverage,
	    ts->pt ? ts->pt->mp : "(none)");
	if (d->force || d->exclude_km > 0.0) {
		char why[640];

		test_describe(d, tile, h, why, sizeof why);
		logmsg("tile %s TEST: %s; %s", ts->name, dec->reason, why);
	}
	return 1;
}

static struct tilesess *session_find(struct daemon *d, const struct vbs_tile *t) {
	int i;

	for (i = 0; i < MAXSESSIONS; i++)
		if (d->sess[i].active && d->sess[i].tile.row == t->row && d->sess[i].tile.col == t->col)
			return &d->sess[i];
	return NULL;
}

static struct tilesess *session_slot(struct daemon *d) {
	int i;

	for (i = 0; i < MAXSESSIONS; i++)
		if (!d->sess[i].active)
			return &d->sess[i];
	return NULL;
}

/* ---- clients ------------------------------------------------------------------------------ */

static int send_all(struct client *c, const void *buf, int n) {
	ssize_t w;

	if (c->pendlen > 0) {
		w = send(c->fd, c->pend, c->pendlen, MSG_DONTWAIT | MSG_NOSIGNAL);
		if (w < 0) {
			if (errno != EAGAIN && errno != EWOULDBLOCK)
				return -1;
			w = 0;
		}
		c->bytes += w;
		if (w < c->pendlen) {
			memmove(c->pend, c->pend + w, c->pendlen - w);
			c->pendlen -= (int)w;
			c->dropped++;
			return 0;
		}
		c->pendlen = 0;
	}
	w = send(c->fd, buf, n, MSG_DONTWAIT | MSG_NOSIGNAL);
	if (w == n) {
		c->bytes += n;
		return 0;
	}
	if (w < 0) {
		if (errno != EAGAIN && errno != EWOULDBLOCK)
			return -1;
		/* real-time data: a rover that cannot keep up gets the next epoch, not a backlog */
		c->dropped++;
		return 0;
	}
	if (n - w > (ssize_t)sizeof c->pend)
		return -1;
	c->bytes += w;
	memcpy(c->pend, (const unsigned char *)buf + w, n - w);
	c->pendlen = (int)(n - w);
	return 0;
}

static void client_close(struct daemon *d, struct client *c, const char *why) {
	if (c->fd < 0)
		return;
	if (c->sess) {
		c->sess->nclients--;
		if (c->sess->nclients <= 0) {
			c->sess->nclients = 0;
			c->sess->idle_since = now_unix();
		}
	}
	logmsg("client %s: closed (%s), %ld bytes, %ld dropped, %.0f s", c->peer, why, c->bytes,
	    c->dropped, now_unix() - c->connected);
	close(c->fd);
	memset(c, 0, sizeof *c);
	c->fd = -1;
	(void)d;
}

static void session_detach_all(struct daemon *d, struct tilesess *ts, int close_clients) {
	int i;

	for (i = 0; i < MAXCLIENTS; i++)
		if (d->cl[i].fd >= 0 && d->cl[i].sess == ts) {
			if (close_clients)
				client_close(d, &d->cl[i], "tile no longer servable");
			else
				d->cl[i].sess = NULL;
		}
	ts->nclients = 0;
}

static void client_attach(struct daemon *d, struct client *c, struct tilesess *ts) {
	unsigned char buf[1024];
	int n;

	if (c->sess == ts)
		return;
	if (c->sess) {
		c->sess->nclients--;
		if (c->sess->nclients <= 0) {
			c->sess->nclients = 0;
			c->sess->idle_since = now_unix();
		}
	}
	c->sess = ts;
	ts->nclients++;
	/* a late joiner needs the station messages before any observation is useful */
	if (vbs_session_converged(ts->sess))
		n = vbs_session_station_msgs(ts->sess, buf, sizeof buf);
	else
		n = ts->pt ? vbs_feed_station_msgs(ts->pt->feed, buf, sizeof buf) : 0;
	if (n > 0)
		send_all(c, buf, n);
	(void)d;
}

static void respond(struct client *c, const char *status, const char *extra) {
	char hdr[1024];

	if (c->ntrip1 && !strncmp(status, "200", 3))
		snprintf(hdr, sizeof hdr, "ICY 200 OK\r\n\r\n");
	else
		snprintf(hdr, sizeof hdr,
		    "HTTP/1.1 %s\r\nNtrip-Version: Ntrip/2.0\r\nServer: NTRIP vbsd/0.1\r\n"
		    "%s%sConnection: close\r\n\r\n",
		    status, !strncmp(status, "200", 3) ? "Content-Type: gnss/data\r\n" : "Content-Length: 0\r\n",
		    extra ? extra : "");
	send_all(c, hdr, (int)strlen(hdr));
}

static const char *SOURCETABLE =
	"SOURCETABLE 200 OK\r\nServer: NTRIP vbsd/0.1\r\nContent-Type: text/plain\r\n\r\n"
	"STR;<LAT>/<LON>/<ALT>;VBS by position;RTCM 3.3;1005(10),1033(10),1230(10),1077(1),1087(1),"
	"1097(1),1127(1);2;GPS+GLO+GAL+BDS;VBS;XX;0.00;0.00;1;1;vbsd;none;B;N;3000;\r\n"
	"ENDSOURCETABLE\r\n";

/* Place (or move) a client by position. Returns 1 served, 0 refused (client closed), -1 error. */
static int client_place(struct daemon *d, struct client *c, double now, int first) {
	struct vbs_tile t;
	struct tilesess *ts;
	struct vbs_decision dec;
	char extra[320];
	int moved;

	moved = vbs_tile_update(c->have_tile ? &c->tile : NULL, c->lat, c->lon, d->hysteresis_m, &t);
	if (!moved && !first)
		return 1;
	if (!moved)
		t = c->tile;
	if ((ts = session_find(d, &t)) == NULL) {
		int r;

		if ((ts = session_slot(d)) == NULL) {
			if (first)
				respond(c, "503 Service Unavailable", "X-VBS-Reason: tile limit reached\r\n");
			client_close(d, c, "no session slot");
			return -1;
		}
		r = session_build(d, ts, &t, &dec, now);
		if (r <= 0) {
			int ri = dec.rule[1] >= '0' && dec.rule[1] <= '5' ? dec.rule[1] - '0' : 6;

			d->refused++;
			d->refused_by_rule[ri]++;
			snprintf(extra, sizeof extra, "X-VBS-Rule: %s\r\nX-VBS-Reason: %s\r\n", dec.rule, dec.reason);
			logmsg("client %s: %.4f %.4f -> 404 %s (%s)", c->peer, c->lat, c->lon, dec.rule, dec.reason);
			if (first)
				respond(c, "404 Not Found", extra);
			client_close(d, c, "refused");
			return 0;
		}
	}
	c->tile = t;
	c->have_tile = 1;
	if (first) {
		respond(c, "200 OK", NULL);
		c->state = 1;
		d->served++;
	}
	client_attach(d, c, ts);
	logmsg("client %s: %.4f %.4f -> tile %s (staid %d, %s)%s", c->peer, c->lat, c->lon, ts->name,
	    ts->staid, vbs_policy_family_str(ts->family), first ? "" : " [moved]");
	return 1;
}

static void client_request(struct daemon *d, struct client *c, double now) {
	char *eol, *line, *p, path[256] = "";
	int ntrip2 = 0;
	struct vbs_gga g;
	int have_gga = 0;

	c->in[c->inlen] = '\0';
	if ((p = strstr(c->in, "\r\n\r\n")) == NULL && (p = strstr(c->in, "\n\n")) == NULL) {
		if (c->inlen >= REQBUF - 1)
			client_close(d, c, "request too long");
		return;
	}
	d->connections++;
	/* request line */
	line = c->in;
	eol = strpbrk(line, "\r\n");
	if (eol)
		*eol = '\0';
	if (sscanf(line, "GET %255s", path) != 1) {
		respond(c, "400 Bad Request", NULL);
		client_close(d, c, "not a GET");
		return;
	}
	/* headers */
	if (eol) {
		char *h = eol + 1;

		while (h < p) {
			char *he;

			while (*h == '\r' || *h == '\n')
				h++;
			he = strpbrk(h, "\r\n");
			if (he)
				*he = '\0';
			if (!strncasecmp(h, "Ntrip-Version:", 14))
				ntrip2 = 1;
			else if (!strncasecmp(h, "Ntrip-GGA:", 10) && vbs_gga_parse(h + 10, &g) == 0)
				have_gga = 1;
			if (he == NULL)
				break;
			h = he + 1;
		}
	}
	c->ntrip1 = !ntrip2;
	if (!strcmp(path, "/") || path[0] == '\0') {
		send_all(c, SOURCETABLE, (int)strlen(SOURCETABLE));
		client_close(d, c, "sourcetable");
		return;
	}
	/* /LAT/LON/ALT; the slot is reused, so an omitted altitude must not inherit the last client's */
	c->h = 0.0;
	if (sscanf(path, "/%lf/%lf/%lf", &c->lat, &c->lon, &c->h) < 2) {
		if (have_gga) {
			c->lat = g.lat_deg;
			c->lon = g.lon_deg;
			c->h = g.h_ellipsoidal;
		} else {
			respond(c, "404 Not Found", "X-VBS-Reason: mountpoint must be /LAT/LON/ALT\r\n");
			client_close(d, c, "bad mountpoint");
			return;
		}
	}
	if (!isfinite(c->h))
		c->h = 0.0;
	/* sscanf reads "nan": every comparison with it is false, so test finiteness explicitly */
	if (!isfinite(c->lat) || !isfinite(c->lon) || fabs(c->lat) > 90.0 || fabs(c->lon) > 360.0) {
		respond(c, "404 Not Found", "X-VBS-Reason: position out of range\r\n");
		client_close(d, c, "bad position");
		return;
	}
	c->have_pos = 1;
	c->inlen = 0;
	client_place(d, c, now, 1);
}

static void client_gga(struct daemon *d, struct client *c, double now) {
	char *nl;

	c->in[c->inlen] = '\0';
	while ((nl = strchr(c->in, '\n')) != NULL) {
		struct vbs_gga g;
		int len = (int)(nl - c->in) + 1;

		*nl = '\0';
		if (vbs_gga_parse(c->in, &g) == 0) {
			c->lat = g.lat_deg;
			c->lon = g.lon_deg;
			if (g.has_height)
				c->h = g.h_ellipsoidal;
			c->last_gga = now;
			if (client_place(d, c, now, 0) <= 0)
				return;
		}
		memmove(c->in, c->in + len, c->inlen - len);
		c->inlen -= len;
		c->in[c->inlen] = '\0';
	}
	if (c->inlen >= REQBUF - 1)
		c->inlen = 0;		/* garbage without newlines: discard */
}

static void client_read(struct daemon *d, struct client *c, double now) {
	ssize_t n = recv(c->fd, c->in + c->inlen, REQBUF - 1 - c->inlen, 0);

	if (n == 0) {
		client_close(d, c, "eof");
		return;
	}
	if (n < 0) {
		if (errno == EAGAIN || errno == EWOULDBLOCK || errno == EINTR)
			return;
		client_close(d, c, "read error");
		return;
	}
	c->inlen += (int)n;
	if (c->state == 0)
		client_request(d, c, now);
	else
		client_gga(d, c, now);
}

static void accept_client(struct daemon *d, double now) {
	struct sockaddr_in6 sa;
	socklen_t sl = sizeof sa;
	int fd = accept(d->lfd, (struct sockaddr *)&sa, &sl), i, one = 1;
	char host[64] = "?";

	if (fd < 0)
		return;
	if (!allowed(d, &sa)) {
		if (sa.sin6_family == AF_INET6)
			inet_ntop(AF_INET6, &sa.sin6_addr, host, sizeof host);
		if (d->denied++ % 100 == 0)
			logmsg("connection from %s refused: not in --allow (%ld so far)", host, d->denied);
		close(fd);
		return;
	}
	for (i = 0; i < MAXCLIENTS; i++)
		if (d->cl[i].fd < 0)
			break;
	if (i == MAXCLIENTS) {
		close(fd);
		return;
	}
	fcntl(fd, F_SETFL, fcntl(fd, F_GETFL) | O_NONBLOCK);
	setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &one, sizeof one);
	if (sa.sin6_family == AF_INET6)
		inet_ntop(AF_INET6, &sa.sin6_addr, host, sizeof host);
	else
		inet_ntop(AF_INET, &((struct sockaddr_in *)&sa)->sin_addr, host, sizeof host);
	memset(&d->cl[i], 0, sizeof d->cl[i]);
	d->cl[i].fd = fd;
	d->cl[i].connected = now;
	snprintf(d->cl[i].peer, sizeof d->cl[i].peer, "%s#%d", host, fd);
	if (i >= d->ncl)
		d->ncl = i + 1;
}

/* ---- feed dispatch and housekeeping ------------------------------------------------------- */

struct dispatch {
	struct daemon *d;
	struct feedref *fr;
	double now;
};

static void on_message(struct vbs_feed *feed, int type, const unsigned char *frame, int len, void *arg) {
	struct dispatch *dp = (struct dispatch *)arg;
	struct daemon *d = dp->d;
	int i;

	(void)feed;
	(void)type;
	for (i = 0; i < MAXSESSIONS; i++) {
		struct tilesess *ts = &d->sess[i];

		if (ts->active && ts->pt == dp->fr && !vbs_session_converged(ts->sess))
			vbs_session_passthrough(ts->sess, frame, len);
	}
}

/*
 * --c2w-cal: one base epoch, corrected once for every session (vbs_c2wcal_base_epoch). Returns the
 * corrected record count in buf, or 0 to use the epoch as it came (flag off, not a UM980, verdict
 * not "distinct", estimate not settled). The build is the cache's RTCM 1033 string, or the stream's
 * own 1033 when the cache has none; a base of unknown build is judged on its own samples.
 */
static int c2w_correct(struct daemon *d, struct feedref *fr, const obsd_t *obs, int n,
    obsd_t *buf) {
	const char *fw = d->bc.e[fr->cache_idx].fw;
	double off;
	int v, nout;

	if (fr->c2w == NULL)
		return 0;
	if (fw[0] == '\0' && fr->feed != NULL)
		fw = vbs_feed_station(fr->feed)->recver;
	nout = vbs_c2wcal_base_epoch(fr->c2w, &d->c2w_groups, vbs_policy_fw_build(fw), obs, n, buf,
	    MAXOBS, &v, &off);
	if (v >= 0 && v + 2 != fr->c2w_logged) {
		double mean;

		vbs_c2wcal_probe_verdict(&fr->c2w->probe, &mean);
		logmsg("c2w: %s fw \"%s\" (build %ld): verdict %s, own 2W-2X %+.3f m over %ld",
		    fr->mp, fw, vbs_policy_fw_build(fw), v ? "distinct, correcting" :
		    "RELABELLED, left as is", mean, fr->c2w->probe.n_both);
		fr->c2w_logged = v + 2;
	}
	return nout;
}

static void on_epoch(struct vbs_feed *feed, const obsd_t *obs, int n, void *arg) {
	static obsd_t c2wbuf[MAXOBS];
	struct dispatch *dp = (struct dispatch *)arg;
	struct daemon *d = dp->d;
	int i, k;

	(void)feed;
	if ((k = c2w_correct(d, dp->fr, obs, n, c2wbuf)) > 0) {
		obs = c2wbuf;
		n = k;
	}
	d->bc.e[dp->fr->cache_idx].last_epoch_unix = dp->now;
	d->bc.e[dp->fr->cache_idx].stream_ok = 1;
	for (i = 0; i < MAXSESSIONS; i++) {
		struct tilesess *ts = &d->sess[i];

		if (!ts->active)
			continue;
		for (k = 0; k < ts->nsta; k++)
			if (ts->feed[k] == dp->fr)
				vbs_session_push(ts->sess, k, obs, n, dp->now);
	}
}

static void housekeeping(struct daemon *d, double now) {
	int i, k;

	/* feed health -> cache live view */
	for (i = 0; i < d->nfeeds; i++) {
		struct feedref *fr = &d->feeds[i];
		struct vbs_basecache_entry *e;

		if (fr->feed == NULL)
			continue;
		e = &d->bc.e[fr->cache_idx];
		/*
		 * The base's own 1005 against the registry. The sessions use the registry (the
		 * coordinate the network publishes); a disagreement is either a stale registry or a
		 * mis-surveyed base, and D56 found stations 7-13 cm off. Either way it lands in the
		 * corrections, so it is said once, loudly, with the numbers.
		 */
		if (!fr->pos_checked && vbs_feed_has_position(fr->feed)) {
			const sta_t *st = vbs_feed_station(fr->feed);
			double dm = dist_m(st->pos, e->xyz);

			fr->pos_checked = 1;
			if (dm > 0.05)
				logmsg("COORDINATE: %s streams 1005 %.3f %.3f %.3f, registry %.3f %.3f %.3f: "
				    "%.3f m apart (sessions use the registry)", fr->mp, st->pos[0], st->pos[1],
				    st->pos[2], e->xyz[0], e->xyz[1], e->xyz[2], dm);
			else if (dm > 0.01)
				logmsg("coordinate: %s streams 1005 %.3f m from the registry", fr->mp, dm);
		}
		if (vbs_feed_healthy(fr->feed, now)) {
			e->stream_ok = 1;
		} else if (e->last_epoch_unix > 0.0 || now - fr->opened > d->connect_grace_s) {
			if (e->stream_ok >= 0 && (e->stream_ok == 1 || e->last_epoch_unix == 0.0))
				logmsg("feed %s: no epochs for %.0f s, marked down", fr->mp,
				    e->last_epoch_unix > 0.0 ? now - e->last_epoch_unix : now - fr->opened);
			e->stream_ok = -1;	/* tried and failing */
		}
		if (fr->refs <= 0 && now - fr->idle_since > d->feed_idle_s) {
			logmsg("feed %s: idle, closed", fr->mp);
			vbs_feed_close(fr->feed);
			free(fr->c2w);
			memset(fr, 0, sizeof *fr);
		}
	}
	/* sessions: rebuild on a lost base, retire when idle */
	for (i = 0; i < MAXSESSIONS; i++) {
		struct tilesess *ts = &d->sess[i];
		int lost = 0;

		if (!ts->active)
			continue;
		for (k = 0; k < ts->nsta; k++)
			if (d->bc.e[ts->cache_idx[k]].stream_ok < 0)
				lost++;
		if (!vbs_session_converged(ts->sess) && ts->pt != NULL &&
		    !vbs_feed_healthy(ts->pt->feed, now) && now - ts->pt->opened > 3.0) {
			struct feedref *alt = choose_passthrough(d, ts->xyz, now, ts->pt);

			if (alt != NULL && alt != ts->pt) {
				unsigned char sm[1024];
				int nn, j;

				logmsg("tile %s: pass-through %s is silent, switching to %s", ts->name,
				    ts->pt->mp, alt->mp);
				feed_put(ts->pt, now);
				ts->pt = alt;
				nn = vbs_feed_station_msgs(alt->feed, sm, sizeof sm);
				for (j = 0; j < d->ncl && nn > 0; j++)
					if (d->cl[j].fd >= 0 && d->cl[j].sess == ts && d->cl[j].state == 1)
						send_all(&d->cl[j], sm, nn);
			} else if (alt != NULL)
				feed_put(alt, now);
		}
		if (lost) {
			struct tilesess nts;
			struct vbs_decision dec;
			struct vbs_tile tile = ts->tile;
			int r, j;

			logmsg("tile %s: %d base(s) down, re-forming the cluster", ts->name, lost);
			r = session_build(d, &nts, &tile, &dec, now);
			if (r == 1) {
				/* swap in place: clients stay attached, the station id is the tile's */
				int ncl = ts->nclients;

				vbs_session_free(ts->sess);
				for (k = 0; k < ts->nsta; k++)
					feed_put(ts->feed[k], now);
				feed_put(ts->pt, now);
				nts.nclients = ncl;
				nts.rebuilds = ts->rebuilds + 1;
				nts.created = ts->created;
				*ts = nts;
				for (j = 0; j < MAXCLIENTS; j++)
					if (d->cl[j].fd >= 0 && d->cl[j].sess == ts) {
						unsigned char buf[1024];
						int n = ts->pt ? vbs_feed_station_msgs(ts->pt->feed, buf, sizeof buf) : 0;

						if (n > 0)
							send_all(&d->cl[j], buf, n);
					}
			} else {
				logmsg("tile %s: no longer servable (%s: %s), closing %d client(s)",
				    ts->name, dec.rule, dec.reason, ts->nclients);
				session_detach_all(d, ts, 1);
				session_free(d, ts, now);
			}
			continue;
		}
		if (ts->nclients <= 0 && now - ts->idle_since > d->idle_s)
			session_free(d, ts, now);
	}
	/* clients that never sent a request */
	for (i = 0; i < d->ncl; i++)
		if (d->cl[i].fd >= 0 && d->cl[i].state == 0 && now - d->cl[i].connected > 30.0)
			client_close(d, &d->cl[i], "no request");
}

/*
 * SIGHUP: re-read the cache. Sessions and feeds hold indices into it, so entries are updated
 * in place by mountpoint and new ones appended; a mountpoint that disappeared is marked not
 * live, and any session using it is re-formed by the next housekeeping. Indices never move.
 */
static void reload_cache(struct daemon *d) {
	struct vbs_basecache nb;
	struct vbs_policy_base *pb;
	int updated = 0, added = 0, gone = 0;

	if (vbs_basecache_load(d->cache_path, &nb) <= 0) {
		logmsg("reload: %s unreadable, keeping the current cache", d->cache_path);
		return;
	}
	/* the startup repositioning (D146), applied again: the merge copies coordinates as given */
	if (d->coords_path != NULL &&
	    vbs_basecache_apply_coords(&nb, d->coords_path, d->coords_max_m, NULL, NULL) < 0) {
		logmsg("reload: --coords %s could not be applied, keeping the current cache", d->coords_path);
		vbs_basecache_free(&nb);
		return;
	}
	if (d->coords_path != NULL)
		vbs_basecache_require_ref(&nb);
	if (vbs_basecache_merge(&d->bc, &nb, &updated, &added, &gone) < 0 ||
	    (pb = (struct vbs_policy_base *)realloc(d->pb, d->bc.n * sizeof *d->pb)) == NULL) {
		logmsg("reload: out of memory, keeping the current cache");
		vbs_basecache_free(&nb);
		return;
	}
	d->pb = pb;
	vbs_basecache_policy_view(&d->bc, 0, d->pb);
	vbs_basecache_free(&nb);
	logmsg("reload: cache %s (generated %s): %d bases, %d updated, %d added, %d gone",
	    d->cache_path, d->bc.generated, d->bc.n, updated, added, gone);
}

static void status(struct daemon *d, double now) {
	int i, ns = 0, nc = 0, nf = 0;

	for (i = 0; i < MAXSESSIONS; i++)
		if (d->sess[i].active)
			ns++;
	for (i = 0; i < d->ncl; i++)
		if (d->cl[i].fd >= 0)
			nc++;
	for (i = 0; i < d->nfeeds; i++)
		if (d->feeds[i].feed)
			nf++;
	logmsg("status: up %.0f s, %d tiles, %d clients, %d feeds; served %ld, refused %ld "
	    "(R0 %ld R1 %ld R2s %ld R3 %ld R4 %ld R5 %ld)", now - d->started, ns, nc, nf, d->served,
	    d->refused, d->refused_by_rule[0], d->refused_by_rule[1], d->refused_by_rule[2],
	    d->refused_by_rule[3], d->refused_by_rule[4], d->refused_by_rule[5]);
	{
		int ne[4] = {0, 0, 0, 0}, ng = 0;

		for (i = 0; i < d->nav.n; i++) {
			int sys = d->nav.eph[i].sat > 0 ? satsys(d->nav.eph[i].sat, NULL) : 0;

			if (sys == SYS_GPS) ne[0]++;
			else if (sys == SYS_GAL) ne[2]++;
			else if (sys == SYS_CMP) ne[3]++;
		}
		for (i = 0; i < d->nav.ng; i++)
			if (d->nav.geph[i].sat > 0)
				ng++;
		logmsg("  shared nav: ephemerides GPS %d, GLO %d, GAL %d, BDS %d", ne[0], ng, ne[2], ne[3]);
	}
	for (i = 0; i < MAXSESSIONS; i++) {
		struct tilesess *ts = &d->sess[i];
		struct vbs_session_stats st;

		if (!ts->active)
			continue;
		vbs_session_get_stats(ts->sess, &st);
		logmsg("  tile %s: %d client(s), %s, master %s, in %ld out %ld incomplete %ld no-corr %ld, "
		    "%.1f sat/epoch, %ld bytes (%ld passed through), rebuilds %ld", ts->name, ts->nclients,
		    st.converged ? "VBS" : (ts->pt ? ts->pt->mp : "no pass-through"), ts->mp[0],
		    st.epochs_in, st.epochs_out, st.epochs_incomplete, st.epochs_no_corr,
		    st.epochs_out ? (double)st.sats_out / st.epochs_out : 0.0, st.bytes_out,
		    st.passthrough_bytes, ts->rebuilds);
		logmsg("    per system G/R/E/C: in %ld/%ld/%ld/%ld, aux corr %ld/%ld/%ld/%ld, "
		    "tile corr %ld/%ld/%ld/%ld, out %ld/%ld/%ld/%ld",
		    st.in_sys[0], st.in_sys[1], st.in_sys[2], st.in_sys[3],
		    st.aux_corr_sys[0], st.aux_corr_sys[1], st.aux_corr_sys[2], st.aux_corr_sys[3],
		    st.corr_sys[0], st.corr_sys[1], st.corr_sys[2], st.corr_sys[3],
		    st.out_sys[0], st.out_sys[1], st.out_sys[2], st.out_sys[3]);
		logmsg("    common bands G/R/E/C: %d/%d/%d/%d; no position G/R/E/C %ld/%ld/%ld/%ld; filter "
		    "rejections elev %ld noslot %ld code %ld phase %ld geom %ld snr %ld",
		    st.bands_sys[0], st.bands_sys[1], st.bands_sys[2], st.bands_sys[3],
		    st.nopos_sys[0], st.nopos_sys[1], st.nopos_sys[2], st.nopos_sys[3],
		    st.rej[0], st.rej[1], st.rej[2], st.rej[3], st.rej[4], st.rej[5]);
		logmsg("    first filter G/R/E/C: sats %d/%d/%d/%d, mean ion sigma %.3f/%.3f/%.3f/%.3f m, "
		    "mean worst amb sigma %.2f/%.2f/%.2f/%.2f cyc, longest run %ld/%ld/%ld/%ld epochs",
		    st.fsat_sys[0], st.fsat_sys[1], st.fsat_sys[2], st.fsat_sys[3],
		    st.fion_sys[0], st.fion_sys[1], st.fion_sys[2], st.fion_sys[3],
		    st.famb_sys[0], st.famb_sys[1], st.famb_sys[2], st.famb_sys[3],
		    st.ftrack_sys[0], st.ftrack_sys[1], st.ftrack_sys[2], st.ftrack_sys[3]);
		{
			char w[3][8] = {"-", "-", "-"};
			int q;

			int sy[3] = {0, 0, 0};

			for (q = 0; q < 3; q++)
				if (st.worst_sat[q] > 0) {
					satno2id(st.worst_sat[q], w[q]);
					sy[q] = satsys(st.worst_sat[q], NULL);
				}
			logmsg("    most-slipping master satellites: %s[#%d sys 0x%02x] x%ld (el %.0f), "
			    "%s[#%d sys 0x%02x] x%ld (el %.0f), %s[#%d sys 0x%02x] x%ld (el %.0f)",
			    w[0], st.worst_sat[0], sy[0], st.worst_n[0], st.worst_el[0],
			    w[1], st.worst_sat[1], sy[1], st.worst_n[1], st.worst_el[1],
			    w[2], st.worst_sat[2], sy[2], st.worst_n[2], st.worst_el[2]);
		}
		logmsg("    master phases G/R/E/C: %ld/%ld/%ld/%ld, LLI lock-loss %ld/%ld/%ld/%ld, half-cycle %ld/%ld/%ld/%ld; "
		    "slips LLI/GF/MW/GAP/NEW G %ld/%ld/%ld/%ld/%ld R %ld/%ld/%ld/%ld/%ld E %ld/%ld/%ld/%ld/%ld C %ld/%ld/%ld/%ld/%ld",
		    st.phase_sys[0], st.phase_sys[1], st.phase_sys[2], st.phase_sys[3],
		    st.lli1_sys[0], st.lli1_sys[1], st.lli1_sys[2], st.lli1_sys[3],
		    st.lli2_sys[0], st.lli2_sys[1], st.lli2_sys[2], st.lli2_sys[3],
		    st.slip_sys[0][0], st.slip_sys[0][1], st.slip_sys[0][2], st.slip_sys[0][3], st.slip_sys[0][4],
		    st.slip_sys[1][0], st.slip_sys[1][1], st.slip_sys[1][2], st.slip_sys[1][3], st.slip_sys[1][4],
		    st.slip_sys[2][0], st.slip_sys[2][1], st.slip_sys[2][2], st.slip_sys[2][3], st.slip_sys[2][4],
		    st.slip_sys[3][0], st.slip_sys[3][1], st.slip_sys[3][2], st.slip_sys[3][3], st.slip_sys[3][4]);
	}
}

/* ---- main loop ---------------------------------------------------------------------------- */

static int listen_on(int port) {
	int fd = socket(AF_INET6, SOCK_STREAM, 0), one = 1, zero = 0;
	struct sockaddr_in6 sa;

	if (fd < 0)
		return -1;
	setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, &one, sizeof one);
	setsockopt(fd, IPPROTO_IPV6, IPV6_V6ONLY, &zero, sizeof zero);
	memset(&sa, 0, sizeof sa);
	sa.sin6_family = AF_INET6;
	sa.sin6_addr = in6addr_any;
	sa.sin6_port = htons((unsigned short)port);
	if (bind(fd, (struct sockaddr *)&sa, sizeof sa) < 0 || listen(fd, 64) < 0) {
		close(fd);
		return -1;
	}
	fcntl(fd, F_SETFL, fcntl(fd, F_GETFL) | O_NONBLOCK);
	return fd;
}

static int run(struct daemon *d) {
	static struct pollfd pfd[MAXCLIENTS + MAXFEEDS + 1];
	static int kind[MAXCLIENTS + MAXFEEDS + 1], idx[MAXCLIENTS + MAXFEEDS + 1];
	unsigned char obuf[VBS_SESSION_OUTBUF];
	int i;

	d->started = now_unix();
	d->last_housekeeping = d->last_status = d->started;
	if ((d->lfd = listen_on(d->listen_port)) < 0) {
		logmsg("cannot listen on port %d: %s", d->listen_port, strerror(errno));
		return 1;
	}
	logmsg("listening on %d (%s); bases from %s:%d; cache %s (%d bases, generated %s); R2s %s; "
	    "2 + 1 backup %s",
	    d->listen_port, d->nallow ? "restricted by --allow" : "open to anyone",
	    d->caster_host, d->caster_port, d->cache_path, d->bc.n, d->bc.generated,
	    d->check_bands ? "on" : "off", d->no_backup ? "off" : "on");
	if (d->force || d->exclude_km > 0.0)
		logmsg("TEST MODE: force %s, exclude-km %.1f -- for field tests only, never run this "
		    "in production", d->force ? "on" : "off", d->exclude_km);
	while (!stop) {
		int nfd = 0, k, r;
		double now = now_unix();

		pfd[nfd].fd = d->lfd;
		pfd[nfd].events = POLLIN;
		kind[nfd] = 0;
		nfd++;
		for (i = 0; i < d->ncl; i++)
			if (d->cl[i].fd >= 0) {
				pfd[nfd].fd = d->cl[i].fd;
				pfd[nfd].events = POLLIN;
				kind[nfd] = 1;
				idx[nfd] = i;
				nfd++;
			}
		for (i = 0; i < d->nfeeds; i++) {
			int fd;

			if (d->feeds[i].feed == NULL)
				continue;
			fd = vbs_feed_fd(d->feeds[i].feed);
			if (fd < 0) {
				struct dispatch dp = {d, &d->feeds[i], now};
				struct vbs_feed_sink fs = {on_epoch, on_message, &dp};

				vbs_feed_pump(d->feeds[i].feed, now, &fs);	/* reconnect timer */
				continue;
			}
			pfd[nfd].fd = fd;
			pfd[nfd].events = vbs_feed_events(d->feeds[i].feed);
			kind[nfd] = 2;
			idx[nfd] = i;
			nfd++;
		}
		r = poll(pfd, nfd, 200);
		if (r < 0 && errno != EINTR) {
			logmsg("poll: %s", strerror(errno));
			break;
		}
		now = now_unix();
		for (i = 0; i < nfd && r > 0; i++) {
			if (!pfd[i].revents)
				continue;
			if (kind[i] == 0)
				accept_client(d, now);
			else if (kind[i] == 1)
				client_read(d, &d->cl[idx[i]], now);
			else {
				struct dispatch dp = {d, &d->feeds[idx[i]], now};
				struct vbs_feed_sink fs = {on_epoch, on_message, &dp};

				if (d->feeds[idx[i]].feed)
					vbs_feed_pump(d->feeds[idx[i]].feed, now, &fs);
			}
		}
		/* run every session and hand its output to its clients */
		for (i = 0; i < MAXSESSIONS; i++) {
			struct tilesess *ts = &d->sess[i];
			int n;

			if (!ts->active)
				continue;
			vbs_session_run(ts->sess, now);
			if (vbs_session_converged(ts->sess) != ts->was_converged) {
				ts->was_converged = vbs_session_converged(ts->sess);
				logmsg("tile %s: %s", ts->name, ts->was_converged ?
				    "converged, serving the synthesis" : "lost convergence, back to pass-through");
				if (!ts->was_converged && ts->pt) {
					unsigned char sm[1024];
					int nn = vbs_feed_station_msgs(ts->pt->feed, sm, sizeof sm);

					for (k = 0; k < d->ncl && nn > 0; k++)
						if (d->cl[k].fd >= 0 && d->cl[k].sess == ts && d->cl[k].state == 1)
							send_all(&d->cl[k], sm, nn);
				}
			}
			while ((n = vbs_session_output(ts->sess, obuf, sizeof obuf)) > 0)
				for (k = 0; k < d->ncl; k++)
					if (d->cl[k].fd >= 0 && d->cl[k].sess == ts && d->cl[k].state == 1)
						if (send_all(&d->cl[k], obuf, n) < 0)
							client_close(d, &d->cl[k], "send error");
		}
		if (reload) {
			reload = 0;
			reload_cache(d);
		}
		if (now - d->last_housekeeping >= 1.0) {
			housekeeping(d, now);
			d->last_housekeeping = now;
		}
		if (now - d->last_status >= 60.0) {
			status(d, now);
			d->last_status = now;
		}
	}
	logmsg("shutting down");
	for (i = 0; i < d->ncl; i++)
		if (d->cl[i].fd >= 0)
			client_close(d, &d->cl[i], "shutdown");
	for (i = 0; i < MAXSESSIONS; i++)
		if (d->sess[i].active)
			session_free(d, &d->sess[i], now_unix());
	for (i = 0; i < d->nfeeds; i++)
		if (d->feeds[i].feed)
			vbs_feed_close(d->feeds[i].feed);
	close(d->lfd);
	return 0;
}

static void usage(void) {
	fprintf(stderr,
	    "usage: vbsd --cache FILE --caster HOST[:PORT] [--user U --password P] [--listen PORT]\n"
	    "            [--allow IP[/bits],...] [--hysteresis-m M] [--stale-s S] [--idle-s S]\n"
	    "            [--no-check-bands] [--max-aux N] [--min-aux N] [--no-backup]\n"
	    "            [--force] [--exclude-km R]       (field tests only, never production)\n"
	    "            [--coords FILE [--coords-max-m M]]\n"
	    "            [--c2w-cal run:SEC|l2c:SEC]          (UM980 L2 P(Y) self-calibration, D167)\n"
    "                                                  (SIGHUP re-reads the cache)\n"
	    "       vbsd --cache FILE --decide LAT LON\n");
}

int main(int argc, char **argv) {
	static struct daemon d;
	int i;
	double dlat = 0.0, dlon = 0.0;
	int decide_only = 0;

	memset(&d, 0, sizeof d);
	d.coords_max_m = 5.0;
	d.caster_port = 2101;
	d.listen_port = 2102;
	d.hysteresis_m = 300.0;
	d.stale_s = 30.0;
	d.idle_s = 120.0;
	d.feed_idle_s = 60.0;
	d.connect_grace_s = 20.0;
	d.max_aux = 0;
	d.min_aux = 0;
	d.check_bands = 1;		/* R2s: a station measured below three GPS bands is refused */
	for (i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--cache") && i + 1 < argc)
			d.cache_path = argv[++i];
		else if (!strcmp(argv[i], "--coords") && i + 1 < argc)
			d.coords_path = argv[++i];
		else if (!strcmp(argv[i], "--coords-max-m") && i + 1 < argc)
			d.coords_max_m = atof(argv[++i]);
		else if (!strcmp(argv[i], "--caster") && i + 1 < argc) {
			char *colon;

			d.caster_host = argv[++i];
			if ((colon = strrchr(d.caster_host, ':')) != NULL) {
				*colon = '\0';
				d.caster_port = atoi(colon + 1);
			}
		} else if (!strcmp(argv[i], "--user") && i + 1 < argc)
			d.user = argv[++i];
		else if (!strcmp(argv[i], "--password") && i + 1 < argc)
			d.password = argv[++i];
		else if (!strcmp(argv[i], "--listen") && i + 1 < argc)
			d.listen_port = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--hysteresis-m") && i + 1 < argc)
			d.hysteresis_m = atof(argv[++i]);
		else if (!strcmp(argv[i], "--stale-s") && i + 1 < argc)
			d.stale_s = atof(argv[++i]);
		else if (!strcmp(argv[i], "--idle-s") && i + 1 < argc)
			d.idle_s = atof(argv[++i]);
		else if (!strcmp(argv[i], "--max-aux") && i + 1 < argc)
			d.max_aux = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--min-aux") && i + 1 < argc)
			d.min_aux = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--check-bands"))
			d.check_bands = 1;
		else if (!strcmp(argv[i], "--no-check-bands"))
			d.check_bands = 0;
		else if (!strcmp(argv[i], "--no-backup"))
			d.no_backup = 1;
		else if (!strcmp(argv[i], "--c2w-cal") && i + 1 < argc) {
			const char *mv = argv[++i];

			d.c2w_win = atof(mv + 4);
			if (!strncmp(mv, "run:", 4) && d.c2w_win > 0.0)
				d.c2w_mode = 1;
			else if (!strncmp(mv, "l2c:", 4) && d.c2w_win > 0.0)
				d.c2w_mode = 2;
			else {
				fprintf(stderr, "bad --c2w-cal %s (run:SEC or l2c:SEC)\n", mv);
				return 2;
			}
		}
		else if (!strcmp(argv[i], "--force"))
			d.force = 1;
		else if (!strcmp(argv[i], "--exclude-km") && i + 1 < argc)
			d.exclude_km = atof(argv[++i]);
		else if (!strcmp(argv[i], "--allow") && i + 1 < argc) {
			char *tok, *sv, *list = argv[++i];

			for (tok = strtok_r(list, ",", &sv); tok; tok = strtok_r(NULL, ",", &sv))
				if (allow_add(&d, tok) < 0) {
					fprintf(stderr, "bad --allow entry: %s\n", tok);
					return 1;
				}
		}
		else if (!strcmp(argv[i], "--decide") && i + 2 < argc) {
			decide_only = 1;
			dlat = atof(argv[++i]);
			dlon = atof(argv[++i]);
		} else {
			usage();
			return 1;
		}
	}
	/*
	 * The coordinate rule (2026-10-02): a base absent from the coordinate JSON is untrustworthy
	 * and the VBS ignores it. Without the JSON every base is absent, so a server started
	 * without one would serve from bases nobody vouches for. --decide is a diagnostic and still
	 * runs, but says so.
	 */
	if (!decide_only && d.coords_path == NULL) {
		fprintf(stderr, "vbsd: --coords is required: bases absent from the coordinate document "
		    "are untrusted and ignored, so without one there is nothing to serve from\n");
		return 2;
	}
	if (d.force || d.exclude_km > 0.0)
		fprintf(stderr, "vbsd: TEST MODE%s%s -- for field tests only, never run this in production\n",
		    d.force ? ", forcing a VBS wherever a homogeneous cluster forms" : "",
		    d.exclude_km > 0.0 ? ", hiding the bases near each tile" : "");
	if (decide_only && d.coords_path == NULL)
		fprintf(stderr, "vbsd: no --coords: every base is UNTRUSTED under the coordinate rule; "
		    "this decision uses registry coordinates and guessed frames\n");
	if (d.cache_path == NULL || (!decide_only && d.caster_host == NULL)) {
		usage();
		return 1;
	}
	if (vbs_basecache_load(d.cache_path, &d.bc) <= 0)
		return 1;
	/*
	 * Reposition the cache from the reference solutions before anything reads it
	 * (D146). This has to happen here, between the load and the policy view, because every
	 * later consumer -- cluster formation, the sessions, the 1005 the VBS emits -- takes the
	 * cache's coordinate and none of them look again.
	 *
	 * Measured worth: on a Mosaic rung whose stations were already correct to under 9 mm,
	 * applying these coordinates halved the rover's harmful distance and cut its harmful
	 * episodes from 11 to 5 on the worst day. It is not a repair of bad stations; it is a
	 * better number in front of an integer decision taken at a quarter of a cycle.
	 */
	if (d.coords_path != NULL) {
		int refused = 0;
		double worst = 0.0;

		if (vbs_basecache_apply_coords(&d.bc, d.coords_path, d.coords_max_m,
		    &refused, &worst) < 0) {
			fprintf(stderr, "vbsd: --coords %s could not be applied\n", d.coords_path);
			return 1;
		}
		if (refused > 0)
			logmsg("coords: %d entries refused; see stderr. The bound is --coords-max-m "
			    "%.1f and it exists because the reference file can carry a mountpoint "
			    "mismatch, which would teleport a base rather than nudge it.",
			    refused, d.coords_max_m);
		logmsg("coords: cache repositioned from %s, worst shift %.3f m", d.coords_path, worst);
		logmsg("coords: %d of %d bases not positioned by the document are ignored "
		    "(absent from the JSON means untrusted)",
		    vbs_basecache_require_ref(&d.bc), d.bc.n);
	}
	d.pb = (struct vbs_policy_base *)calloc(d.bc.n, sizeof *d.pb);
	vbs_basecache_policy_view(&d.bc, 0, d.pb);
	if (vbs_feed_nav_init(&d.nav) < 0)
		return 1;
	for (i = 0; i < MAXCLIENTS; i++)
		d.cl[i].fd = -1;

	if (decide_only) {
		struct vbs_base gated[VBS_POLICY_MAXCLUSTER];
		struct vbs_cluster cl;
		struct vbs_tile t;
		struct vbs_decision dec;
		double h, xyz[3];
		char name[16];
		int ngated, served;

		vbs_tile_from_llh(dlat, dlon, &t);
		if (vbs_basecache_height_at(&d.bc, t.clat, t.clon, vbs_policy_cluster_radius_km(), &h) < 0)
			h = 0.0;
		served = decide_tile(&d, &t, h, xyz, d.force, &dec, gated, &ngated, &cl);
		printf("tile %s staid %d centre %.5f %.5f h %.0f m\n", vbs_tile_name(&t, name, sizeof name),
		    vbs_tile_staid(&t), t.clat, t.clon, h);
		printf("%s %s: %s (nearest base %.1f km)\n", dec.rule, served ? "SERVE" : "404",
		    dec.reason, dec.nearest_km);
		/*
		 * The frame and the per-base verdict are printed because the common refusal is a
		 * cross-border tile: the rules find enough homogeneous receivers, one of them is in
		 * another reference frame, the gate drops it, and the cluster is an auxiliary short.
		 * Without the frame column that reads as an inexplicable "too few auxiliaries".
		 */
		for (i = 0; i < dec.ncluster; i++) {
			const struct vbs_basecache_entry *e = &d.bc.e[dec.cluster[i]];
			char role[48];
			int k, isaux = 0;

			if (i < ngated && gated[i].reason != VBS_BASE_OK)
				snprintf(role, sizeof role, "GATED: %s",
				    vbs_base_reason_str(gated[i].reason));
			else if (served) {
				for (k = 0; k < cl.naux; k++)
					if (cl.aux[k] == i)
						isaux = 1;
				snprintf(role, sizeof role, "%s",
				    i == cl.master ? "master" : isaux ? "aux" :
				    i == cl.backup ? "backup (2 + 1)" : "unused");
			} else
				snprintf(role, sizeof role, "%s", i < ngated ? "usable" : "candidate");
			printf("  %-8s %-10s %-24s %s\n", e->mp, e->frame, e->rec, role);
		}
		if (served)
			printf("cluster: master %s, %d auxiliaries, longest baseline %.1f km, %s\n",
			    d.bc.e[dec.cluster[cl.master]].mp, cl.naux, cl.max_baseline_km,
			    cl.inside_hull ? "inside hull" : "outside the hull");
		if (served)
			printf("recipe: dispersive %s, L2C %s, GLONASS %s\n",
			    disp_str(dec.dispersive, dec.iono_sys),
			    dec.prefer_l2c ? "preferred" : "not preferred", dec.glonass ? "on" : "off");
		/*
		 * The two numbers vbs.max_extrapolation_km is calibrated against (D17 item 3), on a
		 * line of their own so that nothing parsing the cluster line above has to change.
		 */
		if (served)
			printf("geometry: extrapolation %.1f km, leverage %.2f, master %.1f km\n",
			    cl.extrapolation_km, cl.leverage, cl.master_dist_km);
		if (served && cl.backup >= 0)
			printf("backup: %s, %.1f km from the VBS, %.1f km from the master -- held in "
			    "reserve, fitted only when a satellite's primaries fall short\n",
			    d.bc.e[dec.cluster[cl.backup]].mp, cl.backup_dist_km,
			    cl.backup_baseline_km);
		if (d.force || d.exclude_km > 0.0) {
			char why[640];

			test_describe(&d, &t, h, why, sizeof why);
			printf("test: %s\n", why);
		}
		return served ? 0 : 2;
	}

	signal(SIGINT, on_signal);
	signal(SIGTERM, on_signal);
	signal(SIGHUP, on_signal);
	signal(SIGPIPE, SIG_IGN);
	return run(&d);
}

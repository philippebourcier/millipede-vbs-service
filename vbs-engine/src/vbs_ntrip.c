#include <errno.h>
#include <fcntl.h>
#include <netdb.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <poll.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "vbs/vbs_ntrip.h"

/*
 * NTRIP client with verbatim recording -- see the header for the contract.
 */

#define BUFSZ			8192
#define DEFAULT_RECONNECT_S	10
#define HDRMAX			4096

enum state {
	ST_IDLE,		/* waiting out a reconnect delay */
	ST_CONNECTING,
	ST_REQUEST,		/* connected, request not yet sent */
	ST_HEADER,		/* request sent, reading the response header */
	ST_STREAM
};

struct vbs_ntrip {
	char *host, *mountpoint, *user, *password, *record_dir;
	int port, reconnect_delay_s;

	int fd;
	enum state state;
	time_t retry_at;

	char hdr[HDRMAX];
	size_t hdrlen;

	FILE *data_fp, *idx_fp;
	long long data_off;

	struct vbs_ntrip_stats stats;
};

static char *dup_or_null(const char *s) {
	return s != NULL ? strdup(s) : NULL;
}

/* RFC 4648 base64, for the Basic credentials. Unused against an open caster, which needs none. */
static char *base64(const char *in) {
	static const char *tbl = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	size_t n = strlen(in), i, j = 0;
	char *out = (char *)malloc(4 * ((n + 2) / 3) + 1);

	if (out == NULL)
		return NULL;
	for (i = 0; i < n; i += 3) {
		unsigned v = (unsigned char)in[i] << 16;

		v |= (i + 1 < n ? (unsigned char)in[i + 1] : 0) << 8;
		v |= (i + 2 < n ? (unsigned char)in[i + 2] : 0);
		out[j++] = tbl[(v >> 18) & 63];
		out[j++] = tbl[(v >> 12) & 63];
		out[j++] = i + 1 < n ? tbl[(v >> 6) & 63] : '=';
		out[j++] = i + 2 < n ? tbl[v & 63] : '=';
	}
	out[j] = '\0';
	return out;
}

static void path_for(const struct vbs_ntrip *this, const char *ext, char *out, size_t n) {
	snprintf(out, n, "%s/%s%s", this->record_dir, this->mountpoint, ext);
}

static int open_recording(struct vbs_ntrip *this) {
	char path[1024];

	if (this->record_dir == NULL || this->data_fp != NULL)
		return 0;
	mkdir(this->record_dir, 0755);		/* fine if it already exists */
	path_for(this, ".rtcm3", path, sizeof path);
	if ((this->data_fp = fopen(path, "ab")) == NULL)
		return -1;
	this->data_off = ftello(this->data_fp);
	path_for(this, ".idx", path, sizeof path);
	if ((this->idx_fp = fopen(path, "ab")) == NULL) {
		fclose(this->data_fp);
		this->data_fp = NULL;
		return -1;
	}
	/*
	 * The index is text, one line per read: arrival time and the offset the bytes begin at.
	 * Text because a recording that cannot be inspected with head(1) six months later is a
	 * recording nobody will trust, and the volume is trivial -- one short line per read
	 * against a 20 kbit/s stream.
	 */
	if (this->data_off == 0)
		fprintf(this->idx_fp, "# vbs_ntrip record v1: <unix_sec>.<usec> <offset> <len>\n");
	return 0;
}

/*
 * Append one read's worth of payload to the recording and hand it to the sink.
 *
 * The index entry is written *before* the sink runs, so a recording stays complete even if
 * the consumer throws: a capture that has to be repeated costs hours of wall clock against a
 * live network that cannot be rewound.
 */
static int record_and_emit(struct vbs_ntrip *this, const void *buf, size_t n,
	void (*sink)(const char *, const void *, size_t, void *), void *arg) {
	struct timeval tv;

	gettimeofday(&tv, NULL);
	if (this->data_fp != NULL) {
		if (fwrite(buf, 1, n, this->data_fp) != n) {
			this->stats.errors++;
			return -1;
		}
		fprintf(this->idx_fp, "%lld.%06ld %lld %zu\n",
		    (long long)tv.tv_sec, (long)tv.tv_usec, this->data_off, n);
		this->data_off += (long long)n;
	}
	this->stats.bytes += (long long)n;
	if (this->stats.first_byte == 0)
		this->stats.first_byte = tv.tv_sec;
	this->stats.last_byte = tv.tv_sec;
	if (sink != NULL)
		sink(this->mountpoint, buf, n, arg);
	return (int)n;
}

static void close_stream(struct vbs_ntrip *this, time_t now) {
	if (this->fd >= 0) {
		close(this->fd);
		this->fd = -1;
		this->stats.disconnects++;
	}
	this->state = ST_IDLE;
	this->hdrlen = 0;
	this->retry_at = now + this->reconnect_delay_s;
	if (this->data_fp != NULL)
		fflush(this->data_fp);
	if (this->idx_fp != NULL)
		fflush(this->idx_fp);
}

struct vbs_ntrip *vbs_ntrip_open(const struct vbs_ntrip_opt *opt) {
	struct vbs_ntrip *this = (struct vbs_ntrip *)calloc(1, sizeof *this);

	if (this == NULL)
		return NULL;
	this->host = dup_or_null(opt->host);
	this->mountpoint = dup_or_null(opt->mountpoint);
	this->user = dup_or_null(opt->user);
	this->password = dup_or_null(opt->password);
	this->record_dir = dup_or_null(opt->record_dir);
	this->port = opt->port > 0 ? opt->port : 2101;
	this->reconnect_delay_s = opt->reconnect_delay_s > 0 ?
	    opt->reconnect_delay_s : DEFAULT_RECONNECT_S;
	this->fd = -1;
	this->state = ST_IDLE;
	this->retry_at = 0;
	if (this->host == NULL || this->mountpoint == NULL || open_recording(this) < 0) {
		vbs_ntrip_close(this);
		return NULL;
	}
	return this;
}

void vbs_ntrip_close(struct vbs_ntrip *this) {
	if (this == NULL)
		return;
	if (this->fd >= 0)
		close(this->fd);
	if (this->data_fp != NULL)
		fclose(this->data_fp);
	if (this->idx_fp != NULL)
		fclose(this->idx_fp);
	free(this->host);
	free(this->mountpoint);
	free(this->user);
	free(this->password);
	free(this->record_dir);
	free(this);
}

int vbs_ntrip_fd(const struct vbs_ntrip *this) {
	return this->fd;
}

short vbs_ntrip_events(const struct vbs_ntrip *this) {
	return this->state == ST_CONNECTING ? POLLOUT : POLLIN;
}

const char *vbs_ntrip_mountpoint(const struct vbs_ntrip *this) {
	return this->mountpoint;
}

void vbs_ntrip_get_stats(const struct vbs_ntrip *this, struct vbs_ntrip_stats *out) {
	*out = this->stats;
}

static int start_connect(struct vbs_ntrip *this, time_t now) {
	struct addrinfo hints, *res = NULL, *ai;
	char portstr[16];
	int fd = -1, on = 1;

	memset(&hints, 0, sizeof hints);
	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	snprintf(portstr, sizeof portstr, "%d", this->port);
	if (getaddrinfo(this->host, portstr, &hints, &res) != 0) {
		this->stats.errors++;
		this->retry_at = now + this->reconnect_delay_s;
		return -1;
	}
	for (ai = res; ai != NULL; ai = ai->ai_next) {
		if ((fd = socket(ai->ai_family, ai->ai_socktype | SOCK_NONBLOCK,
		    ai->ai_protocol)) < 0)
			continue;
		setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &on, sizeof on);
		if (connect(fd, ai->ai_addr, ai->ai_addrlen) == 0) {
			this->state = ST_REQUEST;
			break;
		}
		if (errno == EINPROGRESS) {
			this->state = ST_CONNECTING;
			break;
		}
		close(fd);
		fd = -1;
	}
	freeaddrinfo(res);
	if (fd < 0) {
		this->stats.errors++;
		this->retry_at = now + this->reconnect_delay_s;
		return -1;
	}
	this->fd = fd;
	this->hdrlen = 0;
	return 0;
}

static int send_request(struct vbs_ntrip *this, time_t now) {
	char req[1024], auth[512] = "";
	int n;

	if (this->user != NULL && this->password != NULL) {
		char raw[256], *b64;

		snprintf(raw, sizeof raw, "%s:%s", this->user, this->password);
		if ((b64 = base64(raw)) != NULL) {
			snprintf(auth, sizeof auth, "Authorization: Basic %s\r\n", b64);
			free(b64);
		}
	}
	n = snprintf(req, sizeof req,
	    "GET /%s HTTP/1.1\r\n"
	    "Host: %s\r\n"
	    "Ntrip-Version: Ntrip/2.0\r\n"
	    "User-Agent: NTRIP vbs-engine/0.1\r\n"
	    "%s"
	    "Connection: close\r\n"
	    "\r\n",
	    this->mountpoint, this->host, auth);
	if (write(this->fd, req, n) != n) {
		this->stats.errors++;
		close_stream(this, now);
		return -1;
	}
	this->state = ST_HEADER;
	return 0;
}

/*
 * Accept the response. NTRIP v2 answers "HTTP/1.x 200"; NTRIP v1 answers the bare
 * "ICY 200 OK". Both are in the wild and both are fine -- anything else is not.
 */
static int header_ok(const char *hdr) {
	return strstr(hdr, "ICY 200") == hdr || strstr(hdr, " 200 ") != NULL;
}

static int read_header(struct vbs_ntrip *this, time_t now,
	void (*sink)(const char *, const void *, size_t, void *), void *arg) {
	char *end;
	ssize_t n;

	n = read(this->fd, this->hdr + this->hdrlen, sizeof this->hdr - this->hdrlen - 1);
	if (n == 0) {
		close_stream(this, now);
		return -1;
	}
	if (n < 0) {
		if (errno == EAGAIN || errno == EWOULDBLOCK || errno == EINTR)
			return 0;
		this->stats.errors++;
		close_stream(this, now);
		return -1;
	}
	this->hdrlen += (size_t)n;
	this->hdr[this->hdrlen] = '\0';

	if ((end = strstr(this->hdr, "\r\n\r\n")) == NULL) {
		if (this->hdrlen >= sizeof this->hdr - 1) {
			this->stats.errors++;
			close_stream(this, now);
			return -1;
		}
		return 0;
	}
	if (!header_ok(this->hdr)) {
		this->stats.errors++;
		close_stream(this, now);
		return -1;
	}
	this->state = ST_STREAM;
	this->stats.connects++;

	/* Payload can arrive in the same read as the header; it must not be dropped. */
	end += 4;
	{
		size_t rest = this->hdrlen - (size_t)(end - this->hdr);
		char tmp[HDRMAX];

		this->hdrlen = 0;
		if (rest > 0) {
			memcpy(tmp, end, rest);
			return record_and_emit(this, tmp, rest, sink, arg);
		}
	}
	return 0;
}

int vbs_ntrip_pump(struct vbs_ntrip *this, time_t now,
	void (*sink)(const char *, const void *, size_t, void *), void *arg) {
	char buf[BUFSZ];
	ssize_t n;

	switch (this->state) {
	case ST_IDLE:
		if (now < this->retry_at)
			return 0;
		return start_connect(this, now) < 0 ? -1 : 0;

	case ST_CONNECTING: {
		int err = 0;
		socklen_t len = sizeof err;

		if (getsockopt(this->fd, SOL_SOCKET, SO_ERROR, &err, &len) < 0 || err != 0) {
			this->stats.errors++;
			close_stream(this, now);
			return -1;
		}
		this->state = ST_REQUEST;
	}
	/* fall through */
	case ST_REQUEST:
		return send_request(this, now) < 0 ? -1 : 0;

	case ST_HEADER:
		return read_header(this, now, sink, arg);

	case ST_STREAM:
		n = read(this->fd, buf, sizeof buf);
		if (n == 0) {
			close_stream(this, now);
			return -1;
		}
		if (n < 0) {
			if (errno == EAGAIN || errno == EWOULDBLOCK || errno == EINTR)
				return 0;
			this->stats.errors++;
			close_stream(this, now);
			return -1;
		}
		return record_and_emit(this, buf, (size_t)n, sink, arg);
	}
	return 0;
}

long vbs_ntrip_replay(const char *record_dir, const char *mountpoint, double speed,
	void (*sink)(const char *, const void *, size_t, void *), void *arg) {
	char path[1024], line[256];
	FILE *data_fp, *idx_fp;
	double prev = 0.0;
	long chunks = 0;

	snprintf(path, sizeof path, "%s/%s.rtcm3", record_dir, mountpoint);
	if ((data_fp = fopen(path, "rb")) == NULL)
		return -1;
	snprintf(path, sizeof path, "%s/%s.idx", record_dir, mountpoint);
	if ((idx_fp = fopen(path, "r")) == NULL) {
		fclose(data_fp);
		return -1;
	}
	while (fgets(line, sizeof line, idx_fp) != NULL) {
		double t;
		long long off;
		long len;
		char *buf;

		if (line[0] == '#')
			continue;
		if (sscanf(line, "%lf %lld %ld", &t, &off, &len) != 3 || len <= 0)
			continue;
		if ((buf = (char *)malloc((size_t)len)) == NULL)
			break;
		if (fseeko(data_fp, (off_t)off, SEEK_SET) != 0 ||
		    fread(buf, 1, (size_t)len, data_fp) != (size_t)len) {
			free(buf);
			break;
		}
		if (speed > 0.0 && prev > 0.0 && t > prev) {
			double wait = (t - prev) / speed;
			struct timespec ts;

			ts.tv_sec = (time_t)wait;
			ts.tv_nsec = (long)((wait - (double)ts.tv_sec) * 1e9);
			nanosleep(&ts, NULL);
		}
		prev = t;
		if (sink != NULL)
			sink(mountpoint, buf, (size_t)len, arg);
		free(buf);
		chunks++;
	}
	fclose(idx_fp);
	fclose(data_fp);
	return chunks;
}

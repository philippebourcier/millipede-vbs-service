/*
 * vbs_record -- record live NTRIP streams verbatim, for spec §8's replay requirement.
 *
 * Records a whole cluster simultaneously, which is the point: network processing needs the
 * same epochs at every station, so streams captured one after another are not a cluster.
 *
 * One poll loop over every stream rather than a thread each. A few dozen sockets moving
 * 20 kbit/s apiece is not a threading problem, and the caster this feeds is careful about
 * exactly this (`docs/vbs/DECISIONS.md` D18: "no reason not to keep good habits").
 *
 * Usage:
 *   vbs_record --dir <out> [--host localhost] [--port 2101] [--minutes N] <mountpoint>...
 *
 * Writes <out>/<mountpoint>.rtcm3 (raw bytes) and <out>/<mountpoint>.idx (arrival index).
 * Both append, so an interrupted capture can be resumed.
 */
#include <errno.h>
#include <fcntl.h>
#include <poll.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "vbs/vbs_ntrip.h"

#define MAXSTREAM	32

static volatile sig_atomic_t stop;

/*
 * Refuse to share an output directory with another recorder.
 *
 * Both files are opened append-only, so two processes writing the same mountpoint do not
 * fail -- they interleave, and the result looks like a recording until you check it. That
 * happened: a launch that appeared to have failed was in fact running, a second was started
 * into the same directory, and four hours of capture were silently doubled up. The tell is
 * that the last index entry's offset+length is half the file size, since the index of one
 * writer cannot account for the other's bytes.
 *
 * flock() rather than a lockfile with O_EXCL: the kernel drops it when the process dies, so
 * a killed recorder never leaves a stale lock behind for the next one to puzzle over.
 */
static int lock_dir(const char *dir) {
	char path[1024];
	int fd;

	mkdir(dir, 0755);
	snprintf(path, sizeof path, "%s/.recording.lock", dir);
	if ((fd = open(path, O_CREAT | O_RDWR, 0644)) < 0)
		return -1;
	if (flock(fd, LOCK_EX | LOCK_NB) < 0) {
		close(fd);
		return -1;
	}
	return fd;			/* held for the lifetime of the process */
}

static void on_signal(int sig) {
	(void)sig;
	stop = 1;
}

static void report(struct vbs_ntrip **s, int n, time_t started) {
	long long total = 0;
	int i, live = 0;
	time_t now = time(NULL);

	for (i = 0; i < n; i++) {
		struct vbs_ntrip_stats st;

		vbs_ntrip_get_stats(s[i], &st);
		total += st.bytes;
		if (st.last_byte != 0 && now - st.last_byte <= 30)
			live++;
	}
	fprintf(stderr, "[%4ld min] %d/%d streaming, %.1f MB total\n",
	    (long)((now - started) / 60), live, n, total / 1048576.0);
	for (i = 0; i < n; i++) {
		struct vbs_ntrip_stats st;

		vbs_ntrip_get_stats(s[i], &st);
		fprintf(stderr, "    %-10s %8.2f MB  %ld connect(s) %ld drop(s) %ld error(s)%s\n",
		    vbs_ntrip_mountpoint(s[i]), st.bytes / 1048576.0,
		    st.connects, st.disconnects, st.errors,
		    (st.last_byte == 0 || now - st.last_byte > 30) ? "  [SILENT]" : "");
	}
}

int main(int argc, char **argv) {
	struct vbs_ntrip *streams[MAXSTREAM];
	struct pollfd pfd[MAXSTREAM];
	const char *dir = NULL, *host = "localhost", *mp[MAXSTREAM];
	int port = 2101, minutes = 0, n = 0, i;
	time_t started, last_report;

	for (i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--dir") && i + 1 < argc)
			dir = argv[++i];
		else if (!strcmp(argv[i], "--host") && i + 1 < argc)
			host = argv[++i];
		else if (!strcmp(argv[i], "--port") && i + 1 < argc)
			port = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--minutes") && i + 1 < argc)
			minutes = atoi(argv[++i]);
		else if (argv[i][0] == '-') {
			fprintf(stderr, "unknown option %s\n", argv[i]);
			return 1;
		} else if (n < MAXSTREAM)
			mp[n++] = argv[i];
	}
	if (dir == NULL || n == 0) {
		fprintf(stderr, "usage: %s --dir <out> [--host H] [--port P] [--minutes N] "
		    "<mountpoint>...\n", argv[0]);
		return 1;
	}
	/* Before opening anything: no two recorders may share an output directory. */
	if (lock_dir(dir) < 0) {
		fprintf(stderr, "%s is already being recorded into by another process "
		    "(or is not writable) -- refusing to interleave\n", dir);
		return 1;
	}
	for (i = 0; i < n; i++) {
		struct vbs_ntrip_opt opt;

		memset(&opt, 0, sizeof opt);
		opt.host = host;
		opt.port = port;
		opt.mountpoint = mp[i];
		opt.record_dir = dir;
		if ((streams[i] = vbs_ntrip_open(&opt)) == NULL) {
			fprintf(stderr, "cannot open %s\n", mp[i]);
			return 1;
		}
	}

	signal(SIGINT, on_signal);
	signal(SIGTERM, on_signal);
	signal(SIGPIPE, SIG_IGN);

	started = last_report = time(NULL);
	fprintf(stderr, "recording %d stream(s) from %s:%d into %s%s\n",
	    n, host, port, dir, minutes ? "" : " (until interrupted)");

	while (!stop) {
		time_t now = time(NULL);
		int nfd = 0, idx[MAXSTREAM];

		if (minutes > 0 && now - started >= minutes * 60)
			break;

		for (i = 0; i < n; i++) {
			int fd = vbs_ntrip_fd(streams[i]);

			if (fd < 0) {
				/* Idle: give it a pump so its reconnect timer can fire. */
				vbs_ntrip_pump(streams[i], now, NULL, NULL);
				continue;
			}
			pfd[nfd].fd = fd;
			pfd[nfd].events = vbs_ntrip_events(streams[i]);
			pfd[nfd].revents = 0;
			idx[nfd] = i;
			nfd++;
		}
		if (poll(pfd, nfd, 1000) < 0 && errno != EINTR) {
			perror("poll");
			break;
		}
		now = time(NULL);
		for (i = 0; i < nfd; i++)
			if (pfd[i].revents)
				vbs_ntrip_pump(streams[idx[i]], now, NULL, NULL);

		if (now - last_report >= 300) {
			report(streams, n, started);
			last_report = now;
		}
	}

	report(streams, n, started);
	for (i = 0; i < n; i++)
		vbs_ntrip_close(streams[i]);
	return 0;
}

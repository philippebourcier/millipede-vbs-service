/*
 * vbs_xcheck -- the same receiver, two decoders: does MRTKLIB's live RTCM decode agree with
 * RENAG's RINEX conversion of the same stream?
 *
 * D19/D20 left this open in Phase 0 and Phase 1 never needed it, because every measurement ran
 * on RINEX. Phase 2 runs on the live decode, and the first live cluster produced no BeiDou
 * corrections while the Phase 1 emits of the same receiver family were full of them -- with
 * the Melbourne-Wubbena slip detector firing on BDS every few epochs live and never offline.
 * Either the decoders disagree about a BDS observable, or they disagree about which signal is
 * which, and this tool says which by comparing them satellite by satellite, band by band, at
 * the same epochs.
 *
 * Usage: vbs_xcheck <record_dir> <mountpoint> <rinex.rnx> [max_epochs]
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_feed.h"
#include "vbs/vbs_signal.h"

#define MAXEP	9000

struct epoch {
	gtime_t t;
	obsd_t obs[MAXOBS];
	int n;
};

static struct epoch *live;
static int nlive;

static void on_epoch(struct vbs_feed *feed, const obsd_t *obs, int n, void *arg) {
	(void)feed;
	(void)arg;
	if (nlive >= MAXEP)
		return;
	live[nlive].t = obs[0].time;
	live[nlive].n = n > MAXOBS ? MAXOBS : n;
	memcpy(live[nlive].obs, obs, live[nlive].n * sizeof *obs);
	nlive++;
}

struct ctx {
	struct vbs_feed *feed;
	double clock;
};

static void sink(const char *mp, const void *buf, size_t n, void *arg) {
	struct ctx *c = (struct ctx *)arg;
	struct vbs_feed_sink fs = {on_epoch, NULL, NULL};

	(void)mp;
	c->clock += 0.001;
	vbs_feed_input(c->feed, (const unsigned char *)buf, (int)n, c->clock, &fs);
}

static int sysidx(int sat) {
	switch (satsys(sat, NULL)) {
	case SYS_GPS: return 0;
	case SYS_GLO: return 1;
	case SYS_GAL: return 2;
	case SYS_CMP: return 3;
	default: return -1;
	}
}

static const char *SYSN[4] = {"GPS", "GLO", "GAL", "BDS"};

int main(int argc, char **argv) {
	static obs_t rnx;
	static nav_t nav;
	struct ctx c;
	gtime_t t0 = {0, 0.0};
	int max_ep = argc > 4 ? atoi(argv[4]) : MAXEP, i, j, r = 0;
	/* per system, per band: matched, |dP| sum/max, |dL| sum/max (cycles, after removing the
	 * per-satellite constant), code-attribute disagreements, live-only, rinex-only */
	long matched[4][10], live_only[4][10], rnx_only[4][10], attr[4][10];
	double dP[4][10], dPmax[4][10], dL[4][10], dLmax[4][10];
	long total_live[4], total_rnx[4];
	long nonint_l[4][10];

	if (argc < 4) {
		fprintf(stderr, "usage: %s <record_dir> <mountpoint> <rinex.rnx> [max_epochs]\n", argv[0]);
		return 1;
	}
	memset(matched, 0, sizeof matched); memset(live_only, 0, sizeof live_only);
	memset(rnx_only, 0, sizeof rnx_only); memset(attr, 0, sizeof attr);
	memset(dP, 0, sizeof dP); memset(dPmax, 0, sizeof dPmax);
	memset(dL, 0, sizeof dL); memset(dLmax, 0, sizeof dLmax);
	memset(total_live, 0, sizeof total_live); memset(total_rnx, 0, sizeof total_rnx);
	memset(nonint_l, 0, sizeof nonint_l);

	live = (struct epoch *)calloc(MAXEP, sizeof *live);
	vbs_feed_nav_init(&nav);
	c.feed = vbs_feed_open_offline(argv[2], &nav, 30.0);
	c.clock = 0.0;
	if (vbs_ntrip_replay(argv[1], argv[2], 0.0, sink, &c) <= 0) {
		fprintf(stderr, "no recording for %s in %s\n", argv[2], argv[1]);
		return 1;
	}
	printf("live decode: %d epochs from %s/%s\n", nlive, argv[1], argv[2]);

	memset(&rnx, 0, sizeof rnx);
	if (readrnxt(argv[3], 1, t0, t0, 0.0, "", &rnx, NULL, NULL) < 0 || rnx.n <= 0) {
		fprintf(stderr, "cannot read %s\n", argv[3]);
		return 1;
	}
	printf("rinex: %d records from %s\n", rnx.n, argv[3]);

	/* walk the live epochs; find the same epoch in the RINEX */
	for (i = 0, r = 0; i < nlive && i < max_ep; i++) {
		int rs, re, k;

		while (r < rnx.n && timediff(rnx.data[r].time, live[i].t) < -1e-3)
			r++;
		if (r >= rnx.n)
			break;
		if (timediff(rnx.data[r].time, live[i].t) > 1e-3)
			continue;			/* epoch missing from the RINEX */
		rs = r;
		re = r;
		while (re < rnx.n && fabs(timediff(rnx.data[re].time, live[i].t)) <= 1e-3)
			re++;
		/* per live record, find the same satellite in the RINEX epoch */
		for (k = 0; k < live[i].n; k++) {
			const obsd_t *lo = &live[i].obs[k];
			const obsd_t *ro = NULL;
			int sx = sysidx(lo->sat), f;

			if (sx < 0)
				continue;
			for (j = rs; j < re; j++)
				if (rnx.data[j].sat == lo->sat) {
					ro = &rnx.data[j];
					break;
				}
			for (f = 0; f < NFREQ + NEXOBS; f++) {
				int lb, rb, g, found = 0;

				if (lo->code[f] == 0 || lo->P[f] == 0.0)
					continue;
				lb = vbs_signal_band(lo->code[f]);
				if (lb <= 0 || lb >= 10)
					continue;
				total_live[sx]++;
				if (ro == NULL) {
					live_only[sx][lb]++;
					continue;
				}
				for (g = 0; g < NFREQ + NEXOBS; g++) {
					if (ro->code[g] == 0 || ro->P[g] == 0.0)
						continue;
					rb = vbs_signal_band(ro->code[g]);
					if (rb != lb)
						continue;
					found = 1;
					matched[sx][lb]++;
					if (ro->code[g] != lo->code[f])
						attr[sx][lb]++;
					{
						double p = fabs(lo->P[f] - ro->P[g]);
						double l = lo->L[f] - ro->L[g];
						double lf = fabs(l - floor(l + 0.5));

						dP[sx][lb] += p;
						if (p > dPmax[sx][lb])
							dPmax[sx][lb] = p;
						/* a converter may carry a different integer; the fractional
						 * part is what must agree */
						dL[sx][lb] += lf;
						if (lf > dLmax[sx][lb])
							dLmax[sx][lb] = lf;
						if (lf > 0.01)
							nonint_l[sx][lb]++;
					}
					break;
				}
				if (!found)
					live_only[sx][lb]++;
			}
			/* bands the RINEX has that the live decode did not */
			if (ro != NULL)
				for (f = 0; f < NFREQ + NEXOBS; f++) {
					int rb, g, found = 0;

					if (ro->code[f] == 0 || ro->P[f] == 0.0)
						continue;
					rb = vbs_signal_band(ro->code[f]);
					if (rb <= 0 || rb >= 10)
						continue;
					total_rnx[sx]++;
					for (g = 0; g < NFREQ + NEXOBS; g++)
						if (lo->code[g] != 0 && lo->P[g] != 0.0 && vbs_signal_band(lo->code[g]) == rb)
							found = 1;
					if (!found)
						rnx_only[sx][rb]++;
				}
		}
	}

	printf("\n%-4s %-5s %8s %9s %9s %10s %9s %9s %9s %9s\n", "sys", "band", "matched",
	    "attr!=", "live-only", "rinex-only", "|dP|mean", "|dP|max", "frac|dL|", "dL>.01");
	for (i = 0; i < 4; i++) {
		int b;

		for (b = 1; b < 10; b++) {
			long m = matched[i][b];

			if (m == 0 && live_only[i][b] == 0 && rnx_only[i][b] == 0)
				continue;
			printf("%-4s %-5d %8ld %9ld %9ld %10ld %9.4f %9.4f %9.4f %9ld\n", SYSN[i], b, m,
			    attr[i][b], live_only[i][b], rnx_only[i][b],
			    m ? dP[i][b] / m : 0.0, dPmax[i][b], m ? dL[i][b] / m : 0.0, nonint_l[i][b]);
		}
	}
	printf("\nobservables: live %ld/%ld/%ld/%ld, rinex %ld/%ld/%ld/%ld (G/R/E/C)\n",
	    total_live[0], total_live[1], total_live[2], total_live[3],
	    total_rnx[0], total_rnx[1], total_rnx[2], total_rnx[3]);
	return 0;
}

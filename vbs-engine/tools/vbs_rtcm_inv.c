/*
 * vbs_rtcm_inv -- what is actually inside a recorded RTCM3 stream.
 *
 * Answers three questions that shape the engine and that the RINEX archive cannot answer:
 *
 *   1. Does MRTKLIB decode this base's live stream at all, and which message types arrive?
 *      D19/D20 left open whether MRTKLIB reads Galileo E5b correctly from a real-time stream
 *      as opposed to from RINEX -- the RINEX path came from RENAG's own converter, so it
 *      never exercised MRTKLIB's RTCM decoder.
 *   2. Which signals does each constellation carry? `01-design-study.md` §3.1 concluded
 *      "BeiDou is single-frequency on this network", but that was measured on F9P-dominated
 *      data, and an F9P looks for B2I, which BDS-3 does not transmit. A Unicore UM98x tracks
 *      B1I/B3I/B2a, so the conclusion may be an artefact of the receiver population.
 *   3. How many bands per satellite -- which is what selects the §4.3 three-frequency cascade
 *      over the dual-frequency wide-lane path, per baseline.
 *
 * Usage: vbs_rtcm_inv <file.rtcm3> [max_epochs]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rtklib.h"

#include "vbs/vbs_signal.h"

#define VBS_NSYS	7

static const int SYSLIST[VBS_NSYS] = {SYS_GPS, SYS_GLO, SYS_GAL, SYS_CMP, SYS_QZS, SYS_SBS, SYS_IRN};
static const char *SYSNAME[VBS_NSYS] = {"GPS", "GLO", "GAL", "BDS", "QZS", "SBS", "IRN"};

static int sysindex(int sys) {
	int i;

	for (i = 0; i < VBS_NSYS; i++)
		if (SYSLIST[i] == sys)
			return i;
	return -1;
}

int main(int argc, char **argv) {
	static rtcm_t rtcm;
	FILE *fp;
	long nobs = 0, nepoch = 0;
	long sigcount[VBS_NSYS][VBS_NBAND_MAX];
	/* Histogram of "how many distinct bands did this satellite carry this epoch". */
	long nbandhist[VBS_NSYS][VBS_NBAND_MAX];
	long codeseen[VBS_NSYS][MAXCODE + 1];
	int max_epochs = argc > 2 ? atoi(argv[2]) : 0;
	int c, i, j;

	if (argc < 2) {
		fprintf(stderr, "usage: %s <file.rtcm3> [max_epochs]\n", argv[0]);
		return 1;
	}
	memset(sigcount, 0, sizeof sigcount);
	memset(nbandhist, 0, sizeof nbandhist);
	memset(codeseen, 0, sizeof codeseen);

	if ((fp = fopen(argv[1], "rb")) == NULL) {
		fprintf(stderr, "cannot open %s\n", argv[1]);
		return 1;
	}
	init_rtcm(&rtcm);

	while ((c = fgetc(fp)) != EOF) {
		int ret = input_rtcm3(&rtcm, (uint8_t)c);

		if (ret == 0 || ret == -1)
			continue;
		if (ret != 1)			/* 1 == an observation set was completed */
			continue;
		nepoch++;
		for (i = 0; i < rtcm.obs.n; i++) {
			const obsd_t *o = &rtcm.obs.data[i];
			int si = sysindex(satsys(o->sat, NULL));
			int bands[VBS_NBAND_MAX] = {0}, nb = 0;

			if (si < 0)
				continue;
			nobs++;
			for (j = 0; j < NFREQ + NEXOBS; j++) {
				int band;

				if (o->code[j] == CODE_NONE)
					continue;
				if (o->L[j] == 0.0 && o->P[j] == 0.0)
					continue;
				if (o->code[j] <= MAXCODE)
					codeseen[si][o->code[j]]++;
				band = vbs_signal_band(o->code[j]);
				if (band <= 0 || band >= VBS_NBAND_MAX)
					continue;
				sigcount[si][band]++;
				if (!bands[band]) {
					bands[band] = 1;
					nb++;
				}
			}
			if (nb < VBS_NBAND_MAX)
				nbandhist[si][nb]++;
		}
		if (max_epochs > 0 && nepoch >= max_epochs)
			break;
	}
	fclose(fp);

	printf("%s: %ld observation epochs, %ld satellite-epochs\n", argv[1], nepoch, nobs);
	/* nmsg3[] indexes 1..299 as types 1001..1299, and 300..329 as 4070..4099. */
	printf("\nRTCM message types:\n");
	for (i = 1; i < 400; i++)
		if (rtcm.nmsg3[i])
			printf("  %4d  %8u\n",
			    i < 300 ? 1000 + i : 4070 + (i - 300), rtcm.nmsg3[i]);
	if (rtcm.nmsg3[0])
		printf("  other %8u\n", rtcm.nmsg3[0]);

	printf("\nSignals per constellation (satellite-epochs carrying each band):\n");
	for (i = 0; i < VBS_NSYS; i++) {
		long tot = 0;

		for (j = 1; j < VBS_NBAND_MAX; j++)
			tot += sigcount[i][j];
		if (tot == 0)
			continue;
		printf("  %-4s", SYSNAME[i]);
		for (j = 1; j < VBS_NBAND_MAX; j++)
			if (sigcount[i][j])
				printf("  L%d:%-8ld", j, sigcount[i][j]);
		printf("\n        codes:");
		for (j = 0; j <= MAXCODE; j++)
			if (codeseen[i][j])
				printf(" %s", code2obs((uint8_t)j));
		printf("\n");
	}

	printf("\nBands per satellite-epoch (this is what picks dual- vs triple-frequency, §4.3):\n");
	for (i = 0; i < VBS_NSYS; i++) {
		long tot = 0;

		for (j = 0; j < VBS_NBAND_MAX; j++)
			tot += nbandhist[i][j];
		if (tot == 0)
			continue;
		printf("  %-4s", SYSNAME[i]);
		for (j = 1; j < VBS_NBAND_MAX; j++)
			if (nbandhist[i][j])
				printf("  %d band:%5.1f%%", j, 100.0 * nbandhist[i][j] / tot);
		printf("\n");
	}
	free_rtcm(&rtcm);
	return 0;
}

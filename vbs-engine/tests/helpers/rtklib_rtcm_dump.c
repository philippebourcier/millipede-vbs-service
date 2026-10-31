/*
 * Independent RTCM3 decoder for the round-trip test.
 *
 * Spec S9 Phase 1 asks for an "encoding round-trip with an independent decoder". The engine
 * encodes with MRTKLIB, so decoding with MRTKLIB would let a shared bug cancel itself out.
 * This helper links RTKLIB-EX (rtklibexplorer/demo5) instead -- a separate codebase, built
 * with its own ABI (NFREQ=3 vs MRTKLIB's 5), which is exactly why it has to be a separate
 * binary rather than another object in the test: the two libraries export the same symbol
 * names with incompatible struct layouts.
 *
 * RTKLIB-EX's convbin could do the decoding instead, but it can only write RINEX, whose
 * printed resolution is ~1 mm -- coarser than the MSM7 carrier-phase LSB (0.14 mm) the test
 * needs to resolve. This dumps full double precision so the tolerance being checked is the
 * one D17 actually specifies.
 *
 * Usage: rtklib_rtcm_dump <file.rtcm3> <YYYY/MM/DD> <HH:MM:SS>
 * Output: one line per observable -- "<satid> <obscode> <P metres> <L cycles>"
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rtklib.h"

int main(int argc, char **argv) {
	static rtcm_t rtcm;
	double ep[6];
	FILE *fp;
	char id[8];
	int ret, i, j;

	if (argc < 4) {
		fprintf(stderr, "usage: %s <file.rtcm3> <YYYY/MM/DD> <HH:MM:SS>\n", argv[0]);
		return 2;
	}
	if (sscanf(argv[2], "%lf/%lf/%lf", ep, ep + 1, ep + 2) != 3 ||
	    sscanf(argv[3], "%lf:%lf:%lf", ep + 3, ep + 4, ep + 5) != 3) {
		fprintf(stderr, "cannot parse time\n");
		return 2;
	}
	if (!init_rtcm(&rtcm))
		return 1;

	/*
	 * MSM messages carry only a time of week, so the decoder needs an approximate epoch to
	 * resolve the week. Without it RTKLIB falls back to the wall clock and the decoded
	 * epoch lands in whatever week today happens to be.
	 */
	rtcm.time = epoch2time(ep);

	if ((fp = fopen(argv[1], "rb")) == NULL) {
		fprintf(stderr, "cannot open %s\n", argv[1]);
		free_rtcm(&rtcm);
		return 1;
	}
	while ((ret = input_rtcm3f(&rtcm, fp)) != -2) {
		if (ret != 1)		/* 1 = observation data complete */
			continue;
		for (i = 0; i < rtcm.obs.n; i++) {
			satno2id(rtcm.obs.data[i].sat, id);
			for (j = 0; j < NFREQ + NEXOBS; j++) {
				if (rtcm.obs.data[i].code[j] == CODE_NONE)
					continue;
				/*
				 * The WHOLE LLI byte is printed, not just bit 0. It was masked to
				 * bit 0 here, which made the half-cycle-ambiguity flag (bit 1)
				 * invisible to every test that used this dumper -- so a VBS that
				 * dropped DF420 would have looked correct, and a rover fixing
				 * integers on a half-cycle-ambiguous phase can land half a cycle
				 * out with nothing downstream able to tell.
				 *
				 * LLI is printed because it is what a rover acts on: RTKLIB's
				 * decoder sets it from the MSM Lock Time Indicator via
				 * lossoflock(), so a reset that survives encode and decode shows
				 * up here as LLI bit 0. Without it this dumper could confirm the
				 * observables round-trip while saying nothing about whether a
				 * discontinuity was ever signalled (spec §5.3).
				 */
				char ts[64];

				time2str(rtcm.obs.data[i].time, ts, 0);
				/*
				 * LLI and the epoch are appended, not prepended: test_zero_baseline
				 * parses the first four fields and must keep working. LLI is here
				 * because it is what a rover acts on -- RTKLIB's decoder sets it
				 * from the MSM Lock Time Indicator via lossoflock(), so a reset that
				 * survives encode and decode appears here as bit 0. Without it this
				 * dumper could confirm the observables round-trip while saying
				 * nothing about whether a discontinuity was signalled at all (§5.3).
				 */
				printf("%s %s %.6f %.6f %d %s\n", id,
				    code2obs(rtcm.obs.data[i].code[j]),
				    rtcm.obs.data[i].P[j], rtcm.obs.data[i].L[j],
				    rtcm.obs.data[i].LLI[j], ts);
			}
		}
	}
	fclose(fp);
	free_rtcm(&rtcm);
	return 0;
}

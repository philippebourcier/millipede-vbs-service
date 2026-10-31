#include <stdio.h>
#include <string.h>

#include "vbs/vbs_signal.h"

/*
 * Signal canonicalisation: (constellation, frequency band), never the RINEX attribute.
 * See the header and docs/vbs/DECISIONS.md D4 for why.
 */

/*
 * Band digit of a CODE_* value, or 0 if unknown.
 *
 * code2obs() returns the RINEX observation code without the leading type letter, e.g. "2X"
 * for L2C(M+L) or "7Q" for E5b pilot. Its first character is the band digit, which is what
 * identifies the physical carrier. An unknown or empty code yields 0.
 */
int vbs_signal_band(uint8_t code) {
	const char *obs = code2obs(code);

	if (obs == NULL || obs[0] < '1' || obs[0] > '9')
		return 0;
	return obs[0] - '0';
}

const char *vbs_signal_str(int sys, int band) {
	static char buf[8];
	char c;

	switch (sys) {
	case SYS_GPS: c = 'G'; break;
	case SYS_GLO: c = 'R'; break;
	case SYS_GAL: c = 'E'; break;
	case SYS_CMP: c = 'C'; break;
	case SYS_QZS: c = 'J'; break;
	case SYS_SBS: c = 'S'; break;
	case SYS_IRN: c = 'I'; break;
	default: c = '?'; break;
	}
	snprintf(buf, sizeof buf, "%c%d", c, band);
	return buf;
}

/*
 * Record which CODE_* a station reports for each (sys, band), from one epoch.
 *
 * Indexed [system index][band digit]. A station may in principle report two attributes of the
 * same band in different slots (e.g. both L2W and L2X, which 6 dataset stations do); we keep
 * the first one seen, since for band-level matching either is equivalent and the engine only
 * needs one observable per band. Which one it was is preserved in the signal set so the
 * attribute-mismatch check stays honest.
 */
#define VBS_NSYS	7
#define VBS_NBAND	10

static int sys_index(int sys) {
	switch (sys) {
	case SYS_GPS: return 0;
	case SYS_GLO: return 1;
	case SYS_GAL: return 2;
	case SYS_CMP: return 3;
	case SYS_QZS: return 4;
	case SYS_SBS: return 5;
	case SYS_IRN: return 6;
	}
	return -1;
}

static const int sys_of_index[VBS_NSYS] = {
	SYS_GPS, SYS_GLO, SYS_GAL, SYS_CMP, SYS_QZS, SYS_SBS, SYS_IRN
};

/*
 * Fill table[sysidx][band] with the first CODE_* observed for that (sys, band), scanning
 * every frequency slot of every satellite in the epoch. Scanning slots rather than trusting
 * a fixed slot convention is the point: the band comes from the code, not from where the
 * reader happened to put it.
 */
static void collect_codes(const obsd_t *obs, int n, int sys_mask,
	uint8_t table[VBS_NSYS][VBS_NBAND]) {
	int i, j, sys, si, band;

	memset(table, 0, VBS_NSYS * VBS_NBAND * sizeof(uint8_t));

	for (i = 0; i < n; i++) {
		sys = satsys(obs[i].sat, NULL);
		if (!(sys & sys_mask))
			continue;
		if ((si = sys_index(sys)) < 0)
			continue;

		for (j = 0; j < NFREQ + NEXOBS; j++) {
			if (obs[i].code[j] == CODE_NONE)
				continue;
			/* An observable with neither phase nor code is not usable. */
			if (obs[i].L[j] == 0.0 && obs[i].P[j] == 0.0)
				continue;
			band = vbs_signal_band(obs[i].code[j]);
			if (band <= 0 || band >= VBS_NBAND)
				continue;
			if (table[si][band] == 0)
				table[si][band] = obs[i].code[j];
		}
	}
}

int vbs_sigset_common(const obsd_t *a, int na, const obsd_t *b, int nb,
	int sys_mask, struct vbs_sigset *out) {
	uint8_t ta[VBS_NSYS][VBS_NBAND], tb[VBS_NSYS][VBS_NBAND];
	int si, band;

	collect_codes(a, na, sys_mask, ta);
	collect_codes(b, nb, sys_mask, tb);

	out->n = 0;
	for (si = 0; si < VBS_NSYS; si++) {
		for (band = 1; band < VBS_NBAND; band++) {
			if (ta[si][band] == 0 || tb[si][band] == 0)
				continue;
			if (out->n >= VBS_MAX_SIGNALS)
				return -1;
			struct vbs_signal *s = &out->sig[out->n++];
			s->sys = sys_of_index[si];
			s->band = band;
			s->code_a = ta[si][band];
			s->code_b = tb[si][band];
			s->attr_mismatch = (s->code_a != s->code_b);
		}
	}
	return out->n;
}

int vbs_sigset_nbands(const struct vbs_sigset *this, int sys) {
	int i, n = 0;

	for (i = 0; i < this->n; i++)
		if (this->sig[i].sys == sys)
			n++;
	return n;
}

double vbs_signal_band_freq(int sys, int band) {
	switch (sys) {
	case SYS_GPS:
		switch (band) {
		case VBS_BAND_L1: return FREQ1;		/* L1 C/A, L1C */
		case VBS_BAND_L2: return FREQ2;
		case VBS_BAND_L5: return FREQ5;
		}
		/*
		 * GPS has no band 6, and must not borrow QZSS's. Sharing one case with QZSS made
		 * this answer FREQ6 for GPS, which cost nothing while no combination in the
		 * cascade used band 6 -- and the moment one did (D76) it manufactured four GPS
		 * steps on a band no GPS satellite transmits. They found no satellites and so
		 * changed no result, but they were attempted every epoch and counted in the
		 * statistics, and a frequency table that answers for a signal that does not exist
		 * is a trap left armed for the next caller.
		 */
		return 0.0;
	case SYS_QZS:
		switch (band) {
		case VBS_BAND_L1: return FREQ1;
		case VBS_BAND_L2: return FREQ2;
		case VBS_BAND_L5: return FREQ5;
		case VBS_BAND_L6: return FREQ6;		/* LEX/CLAS, QZSS only */
		}
		return 0.0;
	case SYS_GAL:
		switch (band) {
		case VBS_BAND_L1: return FREQ1;		/* E1 */
		case VBS_BAND_L5: return FREQ5;		/* E5a */
		case VBS_BAND_L6: return FREQ6;		/* E6 */
		case VBS_BAND_L7: return FREQ7;		/* E5b */
		case VBS_BAND_L8: return FREQ8;		/* E5 a+b */
		}
		return 0.0;
	case SYS_GLO:
		/*
		 * The k = 0 nominal, NOT any real satellite's frequency: GLONASS satellites sit at
		 * 1602 + k*0.5625 MHz and 1246 + k*0.4375 MHz. This is here so the cascade can be
		 * *ordered* -- the G1-G2 wide lane is nominally 0.842 m, which is where it belongs
		 * among the other constellations' combinations -- and for nothing else.
		 *
		 * Anything computing a wavelength that a measurement depends on must call
		 * vbs_signal_sat_freq(), which takes the satellite and the ephemerides. Using this
		 * value for a GLONASS ambiguity would be wrong by up to 8 MHz, about half a per
		 * cent of the wide lane, which is a quarter of a cycle over a few hundred cycles.
		 */
		switch (band) {
		case VBS_BAND_L1: return FREQ1_GLO;
		case VBS_BAND_L2: return FREQ2_GLO;
		}
		return 0.0;
	case SYS_CMP:
		switch (band) {
		case VBS_BAND_L1: return FREQ1;		/* B1C, 1575.42 -- not B1I */
		case VBS_BAND_L2: return FREQ1_CMP;	/* B1I, 1561.098 */
		case VBS_BAND_L5: return FREQ5;		/* B2a */
		case VBS_BAND_L6: return FREQ3_CMP;	/* B3I */
		case VBS_BAND_L7: return FREQ2_CMP;	/* B2b */
		case VBS_BAND_L8: return FREQ8;		/* B2 a+b */
		}
		return 0.0;
	}
	return 0.0;
}

/*
 * Carrier frequency of one *satellite* on one band, in Hz, or 0 if it does not carry it.
 *
 * `vbs_signal_band_freq` answers per constellation, which is right for every CDMA system and
 * wrong for GLONASS: its satellites are separated by frequency, not by code, sitting at
 * 1602 + k*0.5625 MHz on G1 and 1246 + k*0.4375 MHz on G2 with the channel number k carried in
 * the broadcast ephemeris. Anything that needs a wavelength for a specific satellite -- a
 * double-differenced ambiguity, a phase translation in the synthesis -- has to use this one.
 * `vbs_signal_band_freq` remains correct for the genuinely per-constellation uses, such as
 * ordering the cascade, where a nominal wavelength is what is wanted.
 *
 * The channel lookup is delegated to the library's `sat2freq`, which reads `geph[].frq` and is
 * already exercised by its positioning path; reimplementing the arithmetic here would be a
 * second place to get it wrong. `nav` must carry GLONASS ephemerides, so a nav file without
 * them yields 0 and the caller must treat that as "no frequency" rather than as a value --
 * the same in-band-sentinel trap as D21's phase_m() and the tropmodel() case in vbs_interp.
 */
double vbs_signal_sat_freq(int sat, int band, const nav_t *nav) {
	int sys = satsys(sat, NULL);

	if (sys != SYS_GLO)
		return vbs_signal_band_freq(sys, band);
	if (nav == NULL)
		return 0.0;
	/* G1 and G2 only; GLONASS band 3 (L3OC, CDMA) is not FDMA and is not handled here. */
	switch (band) {
	case VBS_BAND_L1: return sat2freq(sat, CODE_L1C, nav);
	case VBS_BAND_L2: return sat2freq(sat, CODE_L2C, nav);
	}
	return 0.0;
}

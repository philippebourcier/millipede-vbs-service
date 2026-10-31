#ifndef __VBS_SIGNAL_H__
#define __VBS_SIGNAL_H__

#include <stdint.h>

#include "rtklib.h"

/*
 * Canonical signal identity for the VBS engine.
 *
 * docs/vbs/DECISIONS.md D4: signals are matched between the two ends of a baseline by
 * (constellation, frequency band) -- never by the RINEX tracking-attribute string. On this
 * network 54 of 76 sampled stations label GPS L2C as "L2X" and 17 label it "L2L"/"L2S", from
 * identical u-blox F9P hardware; matching on the attribute letter reports no common GPS L2 on
 * a third of baselines and would discard them for a labelling difference with no physical
 * meaning.
 *
 * The band here is the RINEX 3 band digit (the first character of an observation code:
 * "2X" -> 2), which identifies the physical carrier unambiguously within a constellation.
 * It is deliberately NOT MRTKLIB's frequency slot index: that index is a positioning-engine
 * convention which, for Galileo, assigns slot 1 to E5a, while every base on this network
 * streams E5b (D19). Deriving the band from the code itself is independent of whatever slot
 * a library chose to park the observation in.
 */

/* Band digits seen on this network. Values are the RINEX band digit itself. */
#define VBS_BAND_L1	1	/* GPS L1 / GAL E1 / BDS B1 / GLO G1 */
#define VBS_BAND_L2	2	/* GPS L2 / BDS B1I(2I) / GLO G2 */
#define VBS_BAND_L5	5	/* GPS L5 / GAL E5a / BDS B2a */
#define VBS_BAND_L6	6	/* GAL E6 / BDS B3 */
#define VBS_BAND_L7	7	/* GAL E5b / BDS B2b */
#define VBS_BAND_L8	8	/* GAL E5(a+b) */

#define VBS_MAX_SIGNALS	16	/* per baseline, across all constellations */
#define VBS_NBAND_MAX	10	/* exclusive upper bound on band digits (1..9) */

/*
 * One canonical signal common to both ends of a baseline.
 *
 * code_a/code_b are the CODE_* values each end actually reports for this band. They differ
 * whenever the two ends label the same physical signal differently (the L2X vs L2L case);
 * attr_mismatch records that, because mixing tracking-attribute variants across a baseline
 * could in principle leave a constant receiver-dependent phase offset, which would make the
 * double-difference ambiguity non-integer. All receivers involved here are the same model, so
 * the expectation is zero -- but D4 requires that to be measured in the leave-one-out
 * residuals rather than assumed, and it cannot be measured if it is not first recorded.
 */
struct vbs_signal {
	int sys;		/* SYS_GPS, SYS_GAL, SYS_CMP, ... */
	int band;		/* VBS_BAND_*, the RINEX band digit */
	uint8_t code_a;		/* CODE_* reported at end A */
	uint8_t code_b;		/* CODE_* reported at end B */
	int attr_mismatch;	/* code_a != code_b (same band, different attribute) */
};

/* The common signal set of one baseline. */
struct vbs_sigset {
	int n;
	struct vbs_signal sig[VBS_MAX_SIGNALS];
};

/*
 * Band digit of a CODE_* value, or 0 if the code is unknown/empty.
 */
int vbs_signal_band(uint8_t code);

/*
 * Short human-readable name of a (sys, band) pair, e.g. "G2", "E7", "C2".
 * Returns a pointer to a static per-call buffer: for logging, not for keeping.
 */
const char *vbs_signal_str(int sys, int band);

/*
 * Build the common signal set of a baseline from one epoch observed at each end.
 *
 * sys_mask is an OR of SYS_* to consider (the engine runs GPS+Galileo, D5). Signals present
 * at only one end are not included -- a double difference needs both. Returns the number of
 * common signals, or -1 if out would overflow.
 */
int vbs_sigset_common(const obsd_t *a, int na, const obsd_t *b, int nb,
	int sys_mask, struct vbs_sigset *out);

/*
 * Nominal carrier frequency of a (constellation, band digit) pair, in Hz, or 0 if that
 * constellation does not use that band.
 *
 * Deliberately keyed on the band digit rather than on a library frequency-slot index, for the
 * same reason as the rest of this module (D4): MRTKLIB's slot 1 is E5a for Galileo while every
 * F9P on this network streams E5b, so a slot index says nothing reliable about which carrier
 * is meant. BeiDou is where this matters most -- band 2 is B1I at 1561.098 MHz while band 1 is
 * B1C at 1575.42 MHz, and the two are 14 MHz apart on the same constellation.
 *
 * GLONASS is absent on purpose: its FDMA frequencies are per-satellite, so a (system, band)
 * pair does not determine one. D12 excludes GLONASS permanently, which is what makes it safe
 * for this function to have no answer for it.
 */
double vbs_signal_band_freq(int sys, int band);

/*
 * Number of distinct bands the set has for one constellation -- i.e. whether this baseline
 * can be processed dual- or triple-frequency for that system (spec S4.3).
 */
int vbs_sigset_nbands(const struct vbs_sigset *this, int sys);

/*
 * Carrier frequency of one satellite on one band, in Hz, or 0.
 *
 * Correct for GLONASS, where `vbs_signal_band_freq` cannot be: FDMA puts each satellite on its
 * own frequency and the channel number comes from the broadcast ephemeris, so `nav` is
 * required. For every other constellation this is just `vbs_signal_band_freq`.
 *
 * Use this wherever the wavelength of a *particular satellite* matters -- a double-differenced
 * ambiguity, a phase translation -- and `vbs_signal_band_freq` only where a per-constellation
 * nominal is what is meant.
 */
double vbs_signal_sat_freq(int sat, int band, const nav_t *nav);

#endif /* __VBS_SIGNAL_H__ */

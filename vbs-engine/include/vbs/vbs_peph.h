#ifndef __VBS_PEPH_H__
#define __VBS_PEPH_H__

#include <stdio.h>

#include "rtklib.h"

/*
 * Precise orbits in place of the broadcast ones -- for measurement, not for serving (D163).
 *
 * Spec §4.7 says broadcast ephemerides are enough: the orbit error reaches a rover only through
 * the change of line of sight between two receivers, |dr|·b/ρ, about 2.5 mm per metre of orbit
 * error at 50 km. That argument is what an SSR stream such as Galileo HAS would be bought
 * against, so this module lets the tools measure it instead. A run with SP3 differs from the
 * broadcast run by the satellite positions and nothing else:
 *
 *  - a precise position replaces the broadcast one only where both exist, so the satellite set,
 *    the health flags, the clocks and the transmission times stay the broadcast run's (satposs
 *    derives the transmission time from the broadcast clock in both modes);
 *  - a satellite the SP3 lacks keeps its broadcast position and is counted. Dropping it would
 *    change the sample set the two runs are compared on, which is the D32 failure;
 *  - SP3 positions are the centre of mass and broadcast ones the antenna phase centre. The
 *    satellite ANTEX closes that gap, and the gap is about a metre -- the size of the effect
 *    being measured -- so loading SP3 without it is refused.
 *
 * The statistics say what the substitution changed, per constellation: the SP3 - broadcast
 * difference, whose radial mean exposes an antenna-offset mismatch (it would be a constant), and
 * its projection on e_a - e_b for the two receivers the caller names, which is the error the
 * synthesis carried.
 *
 * Not reentrant (one static work area): tools only, never vbsd.
 */

#define VBS_PEPH_MAXFILE	8
#define VBS_PEPH_NSYS		4	/* GPS, GLONASS, Galileo, BeiDou */

struct vbs_peph_stats {
	long replaced[VBS_PEPH_NSYS];
	long kept[VBS_PEPH_NSYS];	/* no SP3 position: the broadcast one stayed */
	double sum_d2[VBS_PEPH_NSYS];	/* |SP3 - broadcast|², m² */
	double sum_radial[VBS_PEPH_NSYS];
	double sum_proj2[VBS_PEPH_NSYS];	/* (SP3 - broadcast)·(e_a - e_b), m² */
	double max_proj[VBS_PEPH_NSYS];
};

/*
 * Read the SP3 files and the satellite antennas into a nav that already holds the broadcast
 * ephemerides (GLONASS antenna offsets need its frequency channels). Returns the number of SP3
 * epochs, or -1 with the reason on stderr.
 */
int vbs_peph_load(nav_t *nav, const char *const *sp3, int nsp3, const char *atx);

/*
 * Call right after satposs(..., EPHOPT_BRDC, rs, ...) on the same records: overwrites rs with
 * the precise positions and velocities. a and b are the two receivers whose differenced geometry
 * the statistics are projected on; st may be NULL.
 */
void vbs_peph_substitute(const nav_t *nav, const obsd_t *obs, int n, double *rs,
	const double *a, const double *b, struct vbs_peph_stats *st);

/* One "# sp3 ..." line per constellation that had any position. */
void vbs_peph_report(FILE *fp, const struct vbs_peph_stats *st, const char *geometry);

#endif

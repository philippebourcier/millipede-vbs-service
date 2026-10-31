/*
 * The GGA parser, on sentences as rovers actually send them: with and without a checksum, with
 * a bad checksum (refused), southern and western hemispheres, no height, a different talker,
 * junk before the '$' (the "Ntrip-GGA: " header prefix), and a sentence that is not GGA.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "vbs/vbs_gga.h"

static int nfail;

static void check(int ok, const char *what) {
	printf("  %-66s %s\n", what, ok ? "ok" : "FAIL");
	if (!ok)
		nfail++;
}

int main(void) {
	struct vbs_gga g;

	printf("=== GGA ===\n");
	/* a real-looking Paris sentence, checksum computed for it */
	check(vbs_gga_parse("$GPGGA,123519,4851.3960,N,00221.1320,E,4,12,0.9,35.0,M,44.6,M,,*4C", &g) == 0 ||
	    vbs_gga_parse("$GPGGA,123519,4851.3960,N,00221.1320,E,4,12,0.9,35.0,M,44.6,M,,", &g) == 0,
	    "a well-formed GGA parses");
	vbs_gga_parse("$GPGGA,123519,4851.3960,N,00221.1320,E,4,12,0.9,35.0,M,44.6,M,,", &g);
	check(fabs(g.lat_deg - 48.8566) < 1e-4 && fabs(g.lon_deg - 2.3522) < 1e-4, "lat/lon in degrees");
	check(g.has_height && fabs(g.h_ellipsoidal - 79.6) < 1e-6, "height is orthometric + geoid separation");
	check(g.fix_quality == 4 && g.nsat == 12, "fix quality and satellite count");

	check(vbs_gga_parse("$GNGGA,000000.00,3352.1234,S,15112.5678,W,1,08,1.0,10.0,M,20.0,M,,*00", &g) != 0,
	    "a wrong checksum is refused");
	check(vbs_gga_parse("$GNGGA,000000.00,3352.1234,S,15112.5678,W,1,08,1.0,10.0,M,20.0,M,,", &g) == 0 &&
	    g.lat_deg < 0 && g.lon_deg < 0 && fabs(g.lat_deg + 33.86872) < 1e-4,
	    "south and west are negative, GN talker accepted");
	check(vbs_gga_parse("Ntrip-GGA: $GPGGA,123519,4851.3960,N,00221.1320,E,1,05,2.0,,M,,M,,", &g) == 0 &&
	    !g.has_height, "header prefix skipped; empty height reported as absent");
	check(vbs_gga_parse("$GPRMC,123519,A,4807.038,N,01131.000,E,022.4,084.4,230394,003.1,W*6A", &g) != 0,
	    "RMC is not GGA");
	check(vbs_gga_parse("$GPGGA,123519,9951.3960,N,00221.1320,E,1,05,2.0,,M,,M,,", &g) != 0,
	    "latitude 99 degrees is refused");
	check(vbs_gga_parse("", &g) != 0 && vbs_gga_parse(NULL, &g) != 0, "empty and NULL refused");
	/* strtod() reads "nan" and "inf"; a position built from them passed every range check */
	check(vbs_gga_parse("$GPGGA,123519,nan,N,00221.1320,E,1,05,2.0,,M,,M,,", &g) != 0,
	    "a NaN latitude is refused");
	check(vbs_gga_parse("$GPGGA,123519,4851.3960,N,inf,E,1,05,2.0,,M,,M,,", &g) != 0,
	    "an infinite longitude is refused");
	check(vbs_gga_parse("$GPGGA,123519,4851.3960,N,00221.1320,E,1,05,2.0,nan,M,44.6,M,,", &g) == 0 &&
	    !g.has_height, "a NaN height is reported as absent");
	check(vbs_gga_parse("$GPGGA,123519,4851.3960,N,00221.1320,E,1,05,2.0,1e308,M,1e308,M,,", &g) == 0 &&
	    !g.has_height, "a height that overflows to infinity is reported as absent");
	/* the range check looked only at whole degrees: 90 deg 52 min and 180 deg 30 min both passed */
	check(vbs_gga_parse("$GPGGA,123519,9052.0000,S,00221.1320,E,1,05,2.0,,M,,M,,", &g) != 0,
	    "latitude 90 deg 52 min is refused");
	check(vbs_gga_parse("$GPGGA,123519,4851.3960,N,18030.0000,W,1,05,2.0,,M,,M,,", &g) != 0,
	    "longitude 180 deg 30 min is refused");
	check(vbs_gga_parse("$GPGGA,123519,9000.0000,N,18000.0000,E,1,05,2.0,,M,,M,,", &g) == 0 &&
	    g.lat_deg == 90.0 && g.lon_deg == 180.0, "exactly 90 N and 180 E still parse");

	printf("%s: %d failure(s)\n", nfail ? "FAILED" : "PASSED", nfail);
	return nfail != 0;
}

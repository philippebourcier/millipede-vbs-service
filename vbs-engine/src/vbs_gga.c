#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_gga.h"

/* ddmm.mmmm -> degrees; -1 on a malformed field */
static int dm_to_deg(const char *s, int lon, double *out) {
	double v, deg, min;
	char *end;

	if (s == NULL || *s == '\0')
		return -1;
	v = strtod(s, &end);
	if (end == s || !isfinite(v) || v < 0.0)
		return -1;
	deg = floor(v / 100.0);
	min = v - deg * 100.0;
	/* the assembled value, not just the degrees: "9052" is 90 deg 52 min */
	if (min >= 60.0 || deg + min / 60.0 > (lon ? 180.0 : 90.0))
		return -1;
	*out = deg + min / 60.0;
	return 0;
}

int vbs_gga_parse(const char *line, struct vbs_gga *out) {
	const char *dollar, *star;
	char buf[256], *f[20], *p;
	int nf = 0;
	double lat, lon;

	memset(out, 0, sizeof *out);
	if (line == NULL || (dollar = strchr(line, '$')) == NULL)
		return -1;
	/* $GPGGA, $GNGGA, $GLGGA ... : talker varies, sentence is what matters */
	if (strlen(dollar) < 10 || dollar[1] != 'G' || strncmp(dollar + 3, "GGA,", 4) != 0)
		return -1;
	snprintf(buf, sizeof buf, "%s", dollar + 1);
	/* checksum, when present: XOR of everything between $ and * */
	if ((star = strchr(buf, '*')) != NULL) {
		unsigned x = 0, want;

		for (p = buf; p < star; p++)
			x ^= (unsigned char)*p;
		if (sscanf(star + 1, "%2x", &want) == 1 && want != x)
			return -1;
		buf[star - buf] = '\0';
	} else {
		/* strip CR/LF */
		buf[strcspn(buf, "\r\n")] = '\0';
	}
	for (p = buf; nf < 20 && p != NULL; nf++) {
		f[nf] = p;
		p = strchr(p, ',');
		if (p)
			*p++ = '\0';
	}
	if (nf < 10)
		return -1;
	/* 0 talker+GGA, 1 time, 2 lat, 3 N/S, 4 lon, 5 E/W, 6 quality, 7 nsat, 8 hdop, 9 alt, 10 M,
	 * 11 geoid sep, 12 M, 13 age, 14 station */
	if (dm_to_deg(f[2], 0, &lat) < 0 || dm_to_deg(f[4], 1, &lon) < 0)
		return -1;
	if (f[3][0] == 'S')
		lat = -lat;
	else if (f[3][0] != 'N')
		return -1;
	if (f[5][0] == 'W')
		lon = -lon;
	else if (f[5][0] != 'E')
		return -1;
	out->lat_deg = lat;
	out->lon_deg = lon;
	out->fix_quality = atoi(f[6]);
	out->nsat = atoi(f[7]);
	if (f[9][0] != '\0') {
		double alt = strtod(f[9], NULL), sep = 0.0;

		if (nf > 11 && f[11][0] != '\0')
			sep = strtod(f[11], NULL);
		out->h_ellipsoidal = alt + sep;
		out->has_height = isfinite(out->h_ellipsoidal);
	}
	return 0;
}

#ifndef __VBS_GGA_H__
#define __VBS_GGA_H__

/*
 * NMEA GGA, the one sentence the VBS caster reads from a rover (spec §6b: the millipede caster
 * forwards the client's GGA on the upstream connection, and the VBS caster moves the rover
 * between tiles from it).
 *
 * Only the position is used. Height is returned ELLIPSOIDAL: GGA carries orthometric height
 * (field 9) and geoid separation (field 11) separately, and the tile grid, the geometry and the
 * troposphere model all want the sum. A parser that returned field 9 alone would put every
 * rover ~50 m low in France and nobody would notice until the vertical drifted.
 */

struct vbs_gga {
	double lat_deg, lon_deg;
	double h_ellipsoidal;		/* metres; 0 if the sentence carried no height */
	int fix_quality;		/* field 6: 0 none, 1 GPS, 2 DGPS, 4 RTK fixed, 5 float, ... */
	int nsat;
	int has_height;
};

/*
 * Parse one line (any leading junk before the '$' is skipped; the NMEA checksum is verified
 * when present). Returns 0 on success, -1 if the line is not a usable GGA.
 */
int vbs_gga_parse(const char *line, struct vbs_gga *out);

#endif /* __VBS_GGA_H__ */

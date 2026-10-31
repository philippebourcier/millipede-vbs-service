#ifndef __VBS_BASECACHE_H__
#define __VBS_BASECACHE_H__

/*
 * The base cache: what the VBS caster knows about every physical base before it has seen a
 * byte of its stream (spec §6b, §6c).
 *
 * Generated from the network registry and a RINEX-derived site audit (see
 * examples/base_cache.json), and read here. Format, one object per base:
 *
 *   {"version": 1, "generated": "...", "bases": [
 *     {"mp": "BASE1", "marker": "X00100EXA", "rec": "UBL-F9P", "ant": "BEIBT800S",
 *      "frame": "RGF93", "xyz": [4582081.2, 55183.6, 4421736.6],
 *      "gps_bands": "12", "live": true}, ...]}
 *
 * `gps_bands` is the set of GPS bands the station was OBSERVED to deliver, as digits from
 * {1,2,5}; absent or empty means never audited, which the policy treats as unknown rather than
 * as bad. `frame` is the datum its coordinates are in -- clusters never mix frames (Phase 2
 * decision 0). `live` is the registry's opinion; the daemon overrides it from stream health.
 *
 * The parser is deliberately a small purpose-built one rather than a JSON library: the engine
 * has no other JSON to read, json-c is not available on the compute host, and the format is
 * flat enough that a full library would be the larger dependency.
 */

#include "vbs/vbs_policy.h"

#define VBS_BC_MAXSTR	48

/* Applied coordinate shifts past this are named individually in the log: a centimetre is the
 * normal case, a metre means the registry was wrong rather than stale. */
#define VBS_BC_COORD_LOUD_M	0.50

struct vbs_basecache_entry {
	char mp[VBS_BC_MAXSTR];
	char marker[VBS_BC_MAXSTR];
	char rec[VBS_BC_MAXSTR];
	char ant[VBS_BC_MAXSTR];
	char frame[VBS_BC_MAXSTR];
	/*
	 * Receiver firmware, from RTCM 1033's Receiver Firmware Version field as the caster records
	 * it. Empty means unknown, and unknown is NOT a finding -- the same convention gps_bands
	 * follows. D159 found a firmware regression that no other field exposes: builds 11833 and
	 * 13504 declare a GPS L2C column and never fill it while tracking P(Y) on a fifth of
	 * satellite-records, on the same configuration that build 17548 runs flawlessly.
	 *
	 * This is deliberately the LIVE string rather than a cached verdict, so a station that gets
	 * upgraded stops being refused as soon as the caster's JSON reflects the new firmware and
	 * the cache is re-read (SIGHUP). Never memoise a decision derived from it.
	 */
	char fw[VBS_BC_MAXSTR];
	double xyz[3];
	double llh[3];			/* geodetic, rad/rad/m, derived */
	int gps_bands;			/* count of distinct bands among 1/2/5; 0 = unknown */
	/*
	 * Where this entry's coordinate came from: 0 the registry, 1 the reference document
	 * (D148). It matters because the reference document can be PARTIAL -- the first one
	 * delivered covered 1 079 of 1 520 bases -- and a cluster half repositioned and half not is
	 * internally inconsistent by the few millimetres that D146 measured to change integers.
	 * Between two *consistent* sets, the better one wins; between a mixture, neither does.
	 */
	int coord_from_ref;
	/*
	 * Not in the reference document, or refused by it. The coordinate rule (2026-10-02): a base
	 * absent from the coordinate JSON is untrustworthy and is ignored by the VBS entirely --
	 * no cluster, no density count, no "covered" verdict, no pass-through, no tile height.
	 * Set by vbs_basecache_require_ref(); never cleared by stream state.
	 */
	int untrusted;
	int family;			/* VBS_FAM_* */
	int live;
	/* Runtime state, owned by the daemon: */
	int stream_ok;			/* delivering epochs recently */
	double last_epoch_unix;		/* wall time of the last epoch seen */
};

struct vbs_basecache {
	struct vbs_basecache_entry *e;
	int n;
	char generated[64];
};

/* Load from a JSON file. Returns the number of bases, or -1 with a message on stderr. */
int vbs_basecache_load(const char *path, struct vbs_basecache *out);
void vbs_basecache_free(struct vbs_basecache *bc);

/* Index by mountpoint, or -1. */
int vbs_basecache_find(const struct vbs_basecache *bc, const char *mp);

/*
 * Merge a freshly loaded cache `nb` into the live `bc` without moving any index: sessions and
 * feeds hold indices into `bc`. Entries are updated in place by mountpoint -- coordinates,
 * frame, receiver, firmware, bands, family, liveness and where the coordinate came from -- and
 * new mountpoints are appended. An entry absent from `nb` stays, marked not live.
 *
 * Apply any reference coordinates to `nb` BEFORE merging. The merge copies coordinates as
 * they are, so merging an un-repositioned `nb` would silently put every repositioned base back
 * on its registry coordinate (and leave coord_from_ref claiming otherwise).
 *
 * Returns 0, or -1 if growing `bc` failed (then `bc` is unchanged). Counts are optional.
 */
int vbs_basecache_merge(struct vbs_basecache *bc, const struct vbs_basecache *nb,
	int *n_updated, int *n_added, int *n_gone);

/*
 * The cache as the policy wants it. `pb` must have room for bc->n entries. Uses the runtime
 * `stream_ok` flag when `use_stream_state` is set, else the registry's `live`.
 */
int vbs_basecache_policy_view(const struct vbs_basecache *bc, int use_stream_state,
	struct vbs_policy_base *pb);

/*
 * Ellipsoidal height for a point, from the bases within `radius_km` by inverse-distance
 * weighting. A DEM would be the proper source (D15, D37); until one is wired in, the network's
 * own station heights are the best information available and, unlike a rover's GGA, they make
 * the tile height a pure function of position and cache. Returns 0 on success, -1 if no base
 * is within reach.
 */
int vbs_basecache_height_at(const struct vbs_basecache *bc, double lat_deg, double lon_deg,
	double radius_km, double *h_out);

/*
 * Replace cached coordinates from an external reference solution (D146).
 *
 * The registry's positions are what the network publishes; they are not what a precise
 * point-positioning solution says the antennas are. Measured against precise point-positioning solutions,
 * the two agree to a median of 5 mm horizontal and 3 mm vertical over 648 stations -- and
 * applying the difference anyway **halves the harmful distance a rover travels**, on clusters
 * where no station was more than 9 mm out. Seven millimetres is 4 % of an L1 cycle, five
 * stations' worth of it feeds one integer decision taken at 0.25 cycles, and the effect shows
 * up as fewer wrong integers rather than smaller errors: harmful episodes 11 -> 5 on the day it
 * mattered most.
 *
 * `path` is the operator's reference document, read as it ships -- the same file production
 * will hold, not a converted copy. Its `sources` object is keyed by mountpoint, and two fields
 * beyond the coordinate are used rather than discarded:
 *
 *   calculation_quality   EXPLORATORY_GLOBAL_BASELINE entries are refused. They are the ones
 *                         computed without a national transformation, and they sit a median of
 *                         0.66 m from the registry -- far worse than the coordinate they would
 *                         replace. The file says so; believing it is free.
 *   reference_frame       copied into the entry, so the existing "clusters never mix frames"
 *                         check finally has real data to work from. It found nothing before
 *                         because the registry does not record a frame at all.
 *
 * `max_shift_m` is not optional and there is no sensible "unlimited". The reference file this
 * was built for contains two entries 160 and 215 km from the registry, which are wrong
 * somewhere -- a bad transformation, a mismatched mountpoint, a typo -- and silently teleporting
 * a base is far worse than using a slightly stale coordinate. An override beyond the bound is
 * refused, counted and logged, never applied.
 *
 * Returns the number applied, or -1 with a message on stderr. `n_refused` and `worst_applied_m`
 * are optional out-parameters for the caller's log line.
 */
/*
 * After vbs_basecache_apply_coords(): mark every entry the document did not position --
 * absent from it, refused for distance, or exploratory -- as untrusted. Returns how many.
 */
int vbs_basecache_require_ref(struct vbs_basecache *bc);

int vbs_basecache_apply_coords(struct vbs_basecache *bc, const char *path, double max_shift_m,
	int *n_refused, double *worst_applied_m);

#endif /* __VBS_BASECACHE_H__ */

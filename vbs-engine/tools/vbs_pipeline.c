/*
 * vbs_pipeline -- the whole chain on replayed data: cluster -> netekf x N -> ambres ->
 * interp -> synth, spec §4.1 through §4.6 joined up.
 *
 * Until now every module has been exercised on its own. This is where they have to agree, and
 * it is the workhorse for §8's validation:
 *
 *   --mode loo    synthesize a VBS at the target's own coordinates and double-difference
 *                 against the target's real observations. With corrections applied, the same
 *                 measurement `03-error-budget.md` made with corrections forced to zero --
 *                 so the two are directly comparable, and the difference is what the network
 *                 actually removed.
 *
 *   --mode emit   synthesize at an arbitrary point and write a RINEX 3 observation file
 *                 (MRTKLIB's outrnxobs*), for the rover run of §8 step 4. D30 requires
 *                 sweeping the rover-to-VBS offset up to the grid cell's half-diagonal, which
 *                 is what --offset is for. The header's APPROX POSITION is the VBS's
 *                 *catalogue* (tide-free) coordinate, as the 1005 would carry (D23), so a
 *                 rover engine reading the base position from the header gets the same
 *                 convention a physical base gives it.
 *
 *   --mode identity  the standing-rule check for the rover harness (D40): a VBS placed at the
 *                 master with corrections zeroed must give rover results identical to
 *                 pass-through against the master, epoch by epoch. So this mode emits *every*
 *                 master satellite in *every* epoch with zero corrections and the master's own
 *                 LLI, through the same RINEX writer -- no filters, no validity gating, no
 *                 convergence delay, no §5.3 re-stamping. `--no-corrections` is deliberately
 *                 *not* that: it keeps the gated satellite set so the paired residual
 *                 decomposition compares identical samples, which for a rover means ~12
 *                 satellites instead of ~30 and a re-entry flag on every gate crossing (D46).
 *
 * The target is excluded from the network in both modes: that is the whole point of
 * leave-one-out, and passing its observations to the filter would make the result meaningless
 * in a way the numbers would not reveal.
 *
 * Usage:
 *   vbs_pipeline --nav <nav.rnx> --vbs <x,y,z> [--mode loo|emit] [--offset <e,n,u>]
 *                [--target <file,x,y,z>] [--out <file.obs>]
 *                [--sp3 <orbits.sp3> ... --sat-atx <igs20.atx>]
 *                --base <file,x,y,z> --base <file,x,y,z> ...
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_ambres.h"
#include "vbs/vbs_basemeta.h"
#include "vbs/vbs_c2wcal.h"
#include "vbs/vbs_cluster.h"
#include "vbs/vbs_ddres.h"
#include "vbs/vbs_interp.h"
#include "vbs/vbs_ionopost.h"
#include "vbs/vbs_lock.h"
#include "vbs/vbs_netekf.h"
#include "vbs/vbs_peph.h"
#include "vbs/vbs_policy.h"
#include "vbs/vbs_preproc.h"
#include "vbs/vbs_synth.h"

#define MAXBASE		12
#define MAXPAIR		(MAXBASE * (MAXBASE - 1) / 2)
#define ELEV_MASK_DEG	15.0

struct station {
	char path[640];
	double xyz[3];
	/*
	 * The station's own RINEX header, kept for the receiver type. GLONASS admission depends
	 * on it: FDMA inter-channel biases cancel in a double difference only when both ends
	 * have matching frequency-dependent hardware delays, which in practice means the same
	 * receiver model (D90).
	 */
	sta_t sta;
	obs_t obs;
	int cursor;
	struct vbs_netekf *ekf;			/* NULL for the master and the target */
	struct vbs_amb_fix fix[256];
	int nfix;
	/*
	 * This epoch's slice, captured once. slice_at() advances the station's cursor, so a
	 * second call for the same epoch returns nothing -- which silently starved the
	 * auxiliary-auxiliary filters and made loop closure check zero triangles while
	 * reporting 100 % coverage. Same failure as the target cursor in the paired
	 * decomposition; the fix is the same: read each station's epoch once and share it.
	 */
	const obsd_t *ep;
	int nep;
	struct vbs_corr corr[MAXOBS];
	double sd_I[MAXOBS], sd_T[MAXOBS];
	int ncorr;
	/*
	 * One preprocessor per station, over that station's whole stream (D45). Every filter
	 * containing this station reads its arc ids from here, so a slip resets the ambiguities
	 * in all of them in the same epoch, and the §5.3 output signalling reads the master's.
	 */
	struct vbs_preproc *pp;
	long arc[MAXSAT];
};

static struct station bases[MAXBASE];
/* Which bases the cluster actually selected: master plus the auxiliaries that passed the gate. */
static int used[MAXBASE];
static struct station target;
static nav_t nav;
static int nbase;

static int parse_site(const char *arg, struct station *s) {
	char buf[600], *p;
	int i;

	snprintf(buf, sizeof buf, "%s", arg);
	if ((p = strchr(buf, ',')) == NULL)
		return -1;
	*p++ = '\0';
	snprintf(s->path, sizeof s->path, "%s", buf);
	for (i = 0; i < 3; i++) {
		char *q = strchr(p, ',');

		if (i < 2 && q == NULL)
			return -1;
		if (q != NULL)
			*q++ = '\0';
		s->xyz[i] = atof(p);
		p = q;
	}
	return 0;
}

static int parse_xyz(const char *arg, double *v) {
	return sscanf(arg, "%lf,%lf,%lf", &v[0], &v[1], &v[2]) == 3 ? 0 : -1;
}

/* Records belonging to the epoch at *i; advances *i past them. */
static int slice(const obs_t *obs, int *i, const obsd_t **out) {
	int start = *i, n = 0;

	if (start >= obs->n)
		return 0;
	*out = &obs->data[start];
	while (start + n < obs->n &&
	    fabs(timediff(obs->data[start + n].time, obs->data[start].time)) <= 1e-3)
		n++;
	*i = start + n;
	return n;
}

/*
 * --c2w-cal (see main()). Every UM980 base supplied feeds its build's pooled gate -- the more
 * units of a build, the firmer its verdict -- and every UM980 member of the selected cluster
 * (used[]) then takes its own correction, or none of them does.
 */
/*
 * --c2w-raw-out: the master's observations as they came, for the synthesis alone. The correction
 * then reaches the network estimate (filters, ambiguities, ionosphere) but not the code the VBS
 * emits, so a rover tracking the same `2W` as the master keeps differencing like against like.
 */
static obs_t c2w_raw_master;
#define C2W_SYNTH_EP(ep)	(c2w_raw_master.data != NULL ? \
	c2w_raw_master.data + ((ep) - bases[master].obs.data) : (ep))

static int c2w_apply_cluster(int mode, double win_s, int partial, int master, int raw_out) {
	static struct vbs_c2wcal_groups grp;
	static struct vbs_c2wcal_probe own[MAXBASE];
	static struct vbs_c2wcal_est est;
	int ok[MAXBASE], i, nwant = 0, nok = 0;

	memset(ok, 0, sizeof ok);
	vbs_c2wcal_groups_reset(&grp);
	for (i = 0; i < nbase; i++) {
		struct vbs_c2wcal_probe *gp;

		vbs_c2wcal_probe_reset(&own[i]);
		if (vbs_policy_family(bases[i].sta.rectype) != VBS_FAM_UM980)
			continue;
		vbs_c2wcal_probe_add_epoch(&own[i], bases[i].obs.data, bases[i].obs.n);
		if ((gp = vbs_c2wcal_group(&grp, vbs_policy_fw_build(bases[i].sta.recver))) != NULL)
			vbs_c2wcal_probe_add_epoch(gp, bases[i].obs.data, bases[i].obs.n);
	}
	for (i = 0; i < nbase; i++) {
		const struct vbs_c2wcal_probe *gp;
		long build;
		int v;
		double mean;

		ok[i] = 0;
		if (!used[i] || vbs_policy_family(bases[i].sta.rectype) != VBS_FAM_UM980)
			continue;
		nwant++;
		build = vbs_policy_fw_build(bases[i].sta.recver);
		gp = vbs_c2wcal_group(&grp, build);
		v = vbs_c2wcal_probe_verdict(gp != NULL ? gp : &own[i], NULL);
		vbs_c2wcal_probe_verdict(&own[i], &mean);
		ok[i] = v == 1 && own[i].n_both >= VBS_C2WCAL_MIN_SAMPLES;
		nok += ok[i];
		printf("# c2w B%d %s fw \"%s\" build %ld: own 2W-2X %+.3f m over %ld (%.1f %% identical), "
		    "%ld L2C-less; %s verdict %s -> %s\n", i, bases[i].sta.name[0] ?
		    bases[i].sta.name : "?", bases[i].sta.recver, build, mean, own[i].n_both,
		    own[i].n_both ? 100.0 * own[i].n_ident / own[i].n_both : 0.0, own[i].n_py_only,
		    gp != NULL ? "build" : "own", v == 1 ? "distinct" : v == 0 ? "RELABELLED" :
		    "undecided", ok[i] ? "correct" : "cannot correct");
	}
	if (raw_out && master >= 0 && ok[master]) {
		c2w_raw_master.n = c2w_raw_master.nmax = bases[master].obs.n;
		c2w_raw_master.data = malloc(sizeof(obsd_t) * (size_t)bases[master].obs.n);
		if (c2w_raw_master.data == NULL)
			return -1;
		memcpy(c2w_raw_master.data, bases[master].obs.data,
		    sizeof(obsd_t) * (size_t)bases[master].obs.n);
	}
	if (nok < nwant && !partial) {
		printf("# c2w: withheld from the whole cluster -- %d of %d UM980 members can take it "
		    "(a cluster corrected in part is inconsistent, D148)\n", nok, nwant);
		return 0;
	}
	for (i = 0; i < nbase; i++) {
		const obsd_t *ep;
		int cur = 0, n;
		long nep = 0, nskip = 0, n_l2c = 0, n_py = 0;
		double b = 0.0, bmin = 1e9, bmax = -1e9, bsum = 0.0;

		if (!ok[i])
			continue;
		vbs_c2wcal_est_init(&est, mode == 1 && win_s <= 0.0 ? 0.0 : win_s);
		if (win_s <= 0.0) {		/* acausal: one constant from everything loaded */
			while ((n = slice(&bases[i].obs, &cur, &ep)) > 0)
				vbs_c2wcal_est_add_epoch(&est, ep[0].time, ep, n);
			if (!vbs_c2wcal_est_value(&est, bases[i].obs.data[0].time, &b, NULL)) {
				printf("# c2w B%d: no estimate\n", i);
				continue;
			}
			cur = 0;
		}
		while ((n = slice(&bases[i].obs, &cur, &ep)) > 0) {
			obsd_t *mp = (obsd_t *)ep;	/* our own array: the slice is const only by type */
			int a1 = 0, a2 = 0;

			if (win_s > 0.0) {
				vbs_c2wcal_est_add_epoch(&est, ep[0].time, ep, n);
				if (!vbs_c2wcal_est_value(&est, ep[0].time, &b, NULL)) {
					nskip++;
					continue;
				}
			}
			if (mode == 1)
				vbs_c2wcal_apply(mp, n, b);
			else {
				vbs_c2wcal_apply_l2c(mp, n, b, &a1, &a2);
				n_l2c += a1;
				n_py += a2;
			}
			nep++;
			bsum += b;
			if (b < bmin)
				bmin = b;
			if (b > bmax)
				bmax = b;
		}
		printf("# c2w B%d: %s, %ld epochs corrected (%ld before the estimate settled), "
		    "offset mean %+.3f range %+.3f..%+.3f m", i, mode == 1 ? "constant" : "l2c",
		    nep, nskip, nep ? bsum / nep : 0.0, nep ? bmin : 0.0, nep ? bmax : 0.0);
		if (mode == 2)
			printf(", %ld codes from L2C, %ld L2C-less offset", n_l2c, n_py);
		printf("\n");
	}
	return 0;
}

/*
 * How many distinct (constellation, band) signals a station carries.
 *
 * This is the `nsignal` argument `vbs_cluster_select` has always accepted and that this tool
 * has always passed as NULL -- so the master score's "prefer a richer master" term
 * (score = -distance + 0.5*nsignal) was identically zero and master selection was pure
 * proximity. `vbs_cluster.h` states why it matters: "the master's signal set is the VBS's
 * signal set, since nothing can be synthesized that the master did not observe", and every
 * network baseline runs through the master, so a cluster whose nearest station is its least
 * capable receiver is capped at that station's bands throughout.
 *
 * MIX_F9P_UM9 is exactly that: an F9P 0.9 km nearer the VBS than a UM980 became master of a
 * mixed cluster, capping every baseline to GPS L1/L2 plus Galileo E1/E5b -- no extra-wide lane
 * anywhere, on a cluster containing two UM980s.
 *
 * Counted by asking `vbs_sigset_common` for the station against itself, which reuses the tested
 * band mapping and usable-observation rules rather than restating them. Unioned over several
 * epochs, because a band carried by only some satellites -- GPS L5 on about two thirds of them
 * -- can be missing from any single epoch. The station's own cursor is untouched: the main loop
 * has not started, and advancing it here would silently skip its first epochs.
 */
#define NSIG_EPOCHS	30

static int station_nsignal(const struct station *s, int sys_mask) {
	struct { int sys, band; } seen[VBS_MAX_SIGNALS * 2];
	struct vbs_sigset set;
	const obsd_t *ep;
	int cursor = 0, k, n, nseen = 0, i, j;

	for (k = 0; k < NSIG_EPOCHS; k++) {
		if ((n = slice(&s->obs, &cursor, &ep)) <= 0)
			break;
		if (vbs_sigset_common(ep, n, ep, n, sys_mask, &set) < 0)
			continue;
		for (i = 0; i < set.n; i++) {
			for (j = 0; j < nseen; j++)
				if (seen[j].sys == set.sig[i].sys &&
				    seen[j].band == set.sig[i].band)
					break;
			if (j == nseen && nseen < (int)(sizeof seen / sizeof seen[0])) {
				seen[nseen].sys = set.sig[i].sys;
				seen[nseen].band = set.sig[i].band;
				nseen++;
			}
		}
	}
	return nseen;
}

/*
 * Which bands a station actually delivers, as a bit per RINEX band digit.
 *
 * Separate from station_nsignal() above because the question is the opposite one: that counts
 * signals to score a master, this asks whether a band is *absent*. A negative claim needs a
 * longer window, since a single epoch can miss GPS L5 -- about two thirds of the constellation
 * carries it -- without the station lacking the band. Galileo settles it either way: E5a is on
 * every Galileo satellite, so a station tracking Galileo at all and still reporting no band 5
 * over several minutes is not a sampling artefact.
 *
 * The header is not evidence here. Of 68 u-blox F9P stations in this network, 43 declare a
 * band 5 observation code in SYS / # / OBS TYPES and none of them ever populates it, so a check
 * written against the declaration would have accused two thirds of the F9P fleet of a fault
 * their chip makes impossible.
 */
#define NBAND_EPOCHS	300

static unsigned station_band_mask(const struct station *s, int sys) {
	struct vbs_sigset set;
	const obsd_t *ep;
	int cursor = 0, k, n, i;
	unsigned mask = 0;

	for (k = 0; k < NBAND_EPOCHS; k++) {
		if ((n = slice(&s->obs, &cursor, &ep)) <= 0)
			break;
		if (vbs_sigset_common(ep, n, ep, n, sys, &set) < 0)
			continue;
		for (i = 0; i < set.n; i++)
			if (set.sig[i].band > 0 && set.sig[i].band < VBS_NBAND_MAX)
				mask |= 1u << set.sig[i].band;
	}
	return mask;
}

/*
 * Receivers that publish L2C under the L2W label.
 *
 * A Unicore UM98x emits GPS L2C twice: honestly as `2X`, and again in a `2W` column where a
 * Septentrio or a Trimble puts genuine semi-codeless L2P(Y). Measured on twelve stations of this
 * network across both firmware lines, the `2W` C/N0 sits 0 to 1.4 dB below L1 where real
 * semi-codeless tracking costs about 6 dB.
 *
 * This matters because the engine picks one slot per band by code, so a UM980 and a Mosaic
 * paired on their first L2 slot both present `2W` and the pairing looks clean. It is not: those
 * are different observables with different hardware delays, and D136 measured such baselines
 * carrying two to five times the per-satellite dispersive error of matched ones. Naming the
 * models is narrow on purpose -- a wrong entry here produces a false accusation against a base
 * owner -- and it lists only what this network runs and has been checked against.
 */
static int rectype_relabels_l2w(const char *rectype) {
	return strncmp(rectype, "UNICORE UM98", 12) == 0;
}

/*
 * Which observation CODE each end of a baseline actually feeds the filter, per band.
 *
 * The filter differences one slot per band per station, chosen by a rule that lives inside it.
 * Nothing outside could ask what it chose, and that blind spot cost real time: a UM980 publishes
 * L2C twice, honestly as `2X` and again as `2W`, a Mosaic publishes genuine semi-codeless `2W`
 * and L2C as `2L`, and some stations of both families publish only one of the pair. Whether a
 * given baseline ends up differencing 2W against 2W, 2X against 2L, or 2W against 2L was
 * invisible, and it is not a small thing: mixed-code baselines carry two to five times the
 * per-satellite dispersive error of matched ones (D136).
 *
 * The code is a property of the station's column set, not of the epoch, so no time alignment
 * between the two stations is needed or attempted -- each is sampled on its own cursor and the
 * resulting per-band codes compared. The station's own cursor is untouched: the main loop has
 * not started.
 *
 * What this reports is the LABEL the two ends use, which is the only thing the engine can see.
 * A UM980's `2W` carries L2C on the satellites that transmit it and semi-codeless on the rest,
 * so matching labels do not by themselves prove matching physics. A mismatch, though, is
 * certain: those two ends are differencing observables the RINEX itself says are different.
 */
#define PAIR_EPOCHS	60
#define PAIR_NSYS	7
#define PAIR_NBAND	10

static int pair_sysidx(int sys) {
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

static const int pair_sysof[PAIR_NSYS] = {
	SYS_GPS, SYS_GLO, SYS_GAL, SYS_CMP, SYS_QZS, SYS_SBS, SYS_IRN
};

/* table[sysidx][band] = the code the filter would take, CODE_NONE if the band is absent. */
static void station_band_codes(const struct station *s, int prefer_l2c,
	uint8_t table[PAIR_NSYS][PAIR_NBAND]) {
	const obsd_t *ep;
	int cursor = 0, k, n, q, band, si;

	for (si = 0; si < PAIR_NSYS; si++)
		for (band = 0; band < PAIR_NBAND; band++)
			table[si][band] = CODE_NONE;
	for (k = 0; k < PAIR_EPOCHS; k++) {
		if ((n = slice(&s->obs, &cursor, &ep)) <= 0)
			break;
		for (q = 0; q < n; q++) {
			if ((si = pair_sysidx(satsys(ep[q].sat, NULL))) < 0)
				continue;
			for (band = 1; band < PAIR_NBAND; band++) {
				uint8_t c;

				if (table[si][band] != CODE_NONE)
					continue;
				c = vbs_netekf_band_code(&ep[q], band, prefer_l2c);
				if (c != CODE_NONE)
					table[si][band] = c;
			}
		}
	}
}

static void report_signal_pairing(int master_idx, const struct vbs_cluster *cl, int prefer_l2c) {
	uint8_t tm[PAIR_NSYS][PAIR_NBAND], ta[PAIR_NSYS][PAIR_NBAND];
	int ai, si, band, nmis = 0;

	station_band_codes(&bases[master_idx], prefer_l2c, tm);
	for (ai = 0; ai < cl->naux; ai++) {
		const char *rm = bases[master_idx].sta.rectype;
		const char *ra = bases[cl->aux[ai]].sta.rectype;

		station_band_codes(&bases[cl->aux[ai]], prefer_l2c, ta);
		/*
		 * The case the code comparison below cannot see: both ends say 2W and only one of
		 * them means it. Reported separately because it is the opposite of a mismatch --
		 * the labels agree, which is exactly why it stayed invisible.
		 */
		if (tm[0][2] == CODE_L2W && ta[0][2] == CODE_L2W &&
		    rectype_relabels_l2w(rm) != rectype_relabels_l2w(ra)) {
			printf("# signal pairing WARNING: baseline B%d differences G2 as 2W at both "
			    "ends, but %s publishes L2C under that label while %s publishes genuine "
			    "semi-codeless L2P(Y). The labels match and the observables do not. "
			    "D136 measured such pairs carrying two to five times the per-satellite "
			    "dispersive error of a matched baseline.\n",
			    cl->aux[ai],
			    rectype_relabels_l2w(rm) ? "the master" : "this auxiliary",
			    rectype_relabels_l2w(rm) ? "this auxiliary" : "the master");
			nmis++;
		}
		for (si = 0; si < PAIR_NSYS; si++)
			for (band = 1; band < PAIR_NBAND; band++) {
				if (tm[si][band] == CODE_NONE || ta[si][band] == CODE_NONE)
					continue;
				if (tm[si][band] == ta[si][band])
					continue;
				printf("# signal pairing: baseline B%d differences %s as %s at the "
				    "master and %s at the auxiliary -- the two ends report different "
				    "codes for this band\n",
				    cl->aux[ai], vbs_signal_str(pair_sysof[si], band),
				    code2obs(tm[si][band]), code2obs(ta[si][band]));
				nmis++;
			}
	}
	if (nmis == 0)
		printf("# signal pairing: every baseline differences matching codes on every band\n");
}

/*
 * Receivers whose chip tracks a third frequency.
 *
 * Deliberately short. The table's only job is to decide whether a missing band is the chip's
 * doing or the site's, and a wrong entry turns that into a false accusation pointed at a base
 * owner, so it lists only the models this network runs and only those checked against its own
 * observations: every SEPT MOSAIC-X5 (58 of 58), every Trimble NetR9 and Alloy (5 of 5) and 77
 * of 80 UNICORE UM980 deliver band 5. The F9P is absent on purpose and is the control -- it is
 * L1/L2 by construction, 67 of 68 deliver no band 5, and nothing here should ever fire on
 * hardware that never could.
 */
static int rectype_triple_freq(const char *rectype) {
	static const char *const triple[] = {
		"UNICORE UM98",		/* UM980, UM982 */
		"SEPT MOSAIC",
		"TRIMBLE NETR9",
		"TRIMBLE ALLOY",
	};
	size_t t;

	for (t = 0; t < sizeof triple / sizeof triple[0]; t++)
		if (strncmp(rectype, triple[t], strlen(triple[t])) == 0)
			return 1;
	return 0;
}

/*
 * "G1,G2" or "E1,E7", or "none", into buf.
 *
 * Named per constellation rather than by the bare band digit, because the digit is ambiguous to
 * a reader in exactly the place this line is read: Galileo band 7 is E5b, not anything called
 * L7, and a warning about a missing L5 that lists "L1,L7" invites the conclusion that L7 is
 * some higher frequency that made it through. vbs_signal_str() already holds the mapping.
 */
static const char *band_list(int sys, unsigned mask, char *buf, size_t len) {
	size_t off = 0;
	int b;

	buf[0] = '\0';
	for (b = 1; b < VBS_NBAND_MAX; b++)
		if ((mask & (1u << b)) && off + 4 < len)
			off += (size_t)snprintf(buf + off, len - off, "%s%s", off ? "," : "",
			    vbs_signal_str(sys, b));
	return buf[0] ? buf : "none";
}

/*
 * Warn where a station's delivered bands disagree with its declared receiver.
 *
 * Both directions are the site owner's to fix, and neither is visible to a rover -- which sees
 * only a VBS quietly missing a frequency, with no way to tell that it was ever meant to be
 * there. That invisibility is the whole reason to say it here.
 *
 *   - A triple-frequency receiver delivering no band 5. The chip tracks it, so the cause lies
 *     downstream of the chip: a dual-band antenna, or a receiver output set that omits L5/E5a.
 *     The observations cannot separate those two, and one case argues against the reflex
 *     answer -- one such station carries a BEIBT800S, the same antenna as a station that does
 *     deliver L5 -- so the line names both causes and leaves the diagnosis to the owner.
 *
 *   - A receiver delivering bands its declared model cannot track. One station is registered as
 *     UBL-F9P and delivers GPS L5, Galileo E6 and BeiDou B3; an F9P tracks none of the three.
 *     The metadata is wrong, and that is not cosmetic here: the serving rules key on chip type,
 *     so a triple-frequency station registered as an F9P is refused by the all-F9P rule for a
 *     limitation it does not have.
 *
 * A station carrying neither GPS nor Galileo is skipped rather than reported: with no band-5
 * capable constellation tracked there is nothing to conclude.
 */
static void report_band_capability(void) {	/* over the file-scope `bases`/`nbase` */
	char gb[32], eb[32];
	int i;

	for (i = 0; i < nbase; i++) {
		const char *rec = bases[i].sta.rectype[0] ? bases[i].sta.rectype : "";
		const char *nm = bases[i].sta.name[0] ? bases[i].sta.name : "?";
		unsigned g = station_band_mask(&bases[i], SYS_GPS);
		unsigned e = station_band_mask(&bases[i], SYS_GAL);
		int has5 = ((g | e) & (1u << VBS_BAND_L5)) != 0;

		if (g == 0 && e == 0)
			continue;
		if (rec[0] == '\0') {
			printf("# site warning: B%d (%s) declares no receiver type; chip capability "
			    "cannot be checked\n", i, nm);
			continue;
		}
		if (rectype_triple_freq(rec) && !has5)
			printf("# site warning: B%d (%s, %s) delivers no L5/E5a -- GPS %s, GAL %s. "
			    "The chip tracks three frequencies, so this is a site limit, not a chip "
			    "limit: a dual-band antenna, or a receiver output set that omits L5/E5a. "
			    "The base owner can fix it; until then this station caps the VBS at two "
			    "frequencies whenever it is master.\n",
			    i, nm, rec, band_list(SYS_GPS, g, gb, sizeof gb), band_list(SYS_GAL, e, eb, sizeof eb));
		else if (!rectype_triple_freq(rec) && has5)
			printf("# site warning: B%d (%s) is registered as %s but delivers L5/E5a "
			    "(GPS %s, GAL %s), which that model cannot track. The receiver type in "
			    "the station metadata is wrong, and the serving rules key on chip type.\n",
			    i, nm, rec, band_list(SYS_GPS, g, gb, sizeof gb), band_list(SYS_GAL, e, eb, sizeof eb));
	}
}

/*
 * A band the station declares, and almost never sends.
 *
 * Presence is not coverage. `station_band_mask` sets a band's bit if it appears once in 300
 * epochs, which is right for its question ("can this site do L5 at all?") and useless for this
 * one. Four stations in this network publish a GPS L2 column and populate it on 12.6 to 22.2 %
 * of epochs while their neighbours run 89 to 100 %. Such a
 * station is single-frequency on GPS most of the time, and it is silent about it: the header
 * promises L2, the cluster admits it, and the network filter simply gets very little to work
 * with. Three measurements in this project were taken across such a baseline before anyone
 * noticed, and the diagnosis they produced was wrong (D138).
 *
 * GPS band 2 only, deliberately. Every GPS satellite transmits L2, so a shortfall there is the
 * receiver's doing and the claim is exact. GPS L5 is carried by about two thirds of the
 * constellation and Galileo by all of it, so the same threshold applied across bands would
 * accuse healthy stations. The 60 % line sits in open space: the four bad stations are all
 * below 23 % and every healthy one above 88 %.
 */
#define L2COV_EPOCHS	300
#define L2COV_MIN	0.60

static void report_band_coverage(void) {	/* over the file-scope `bases`/`nbase` */
	int i;

	for (i = 0; i < nbase; i++) {
		const obsd_t *ep;
		const char *nm = bases[i].sta.name[0] ? bases[i].sta.name : "?";
		const char *rec = bases[i].sta.rectype[0] ? bases[i].sta.rectype : "?";
		long n1 = 0, n2 = 0;
		int cursor = 0, k, n, q;

		for (k = 0; k < L2COV_EPOCHS; k++) {
			if ((n = slice(&bases[i].obs, &cursor, &ep)) <= 0)
				break;
			for (q = 0; q < n; q++) {
				if (satsys(ep[q].sat, NULL) != SYS_GPS)
					continue;
				if (vbs_netekf_band_code(&ep[q], 1, 0) != CODE_NONE)
					n1++;
				if (vbs_netekf_band_code(&ep[q], 2, 0) != CODE_NONE)
					n2++;
			}
		}
		if (n1 < 100 || n2 == 0)
			continue;
		if ((double)n2 / (double)n1 < L2COV_MIN)
			printf("# site warning: B%d (%s, %s) carries GPS L2 on only %.1f %% of the "
			    "satellite-epochs it carries GPS L1 (%ld of %ld). Every GPS satellite "
			    "transmits L2, so this is a receiver output or tracking setting, not the "
			    "sky. This station is single-frequency on GPS most of the time and every "
			    "baseline through it is weakened accordingly.\n",
			    i, nm, rec, 100.0 * (double)n2 / (double)n1, n2, n1);
	}
}

/* Advance a station to time t; returns its epoch slice, or 0 if it has none at t. */
static int slice_at(struct station *s, gtime_t t, const obsd_t **out) {
	const obsd_t *ep;
	int save, n;

	for (;;) {
		save = s->cursor;
		if ((n = slice(&s->obs, &s->cursor, &ep)) <= 0)
			return 0;
		if (timediff(ep[0].time, t) < -1e-3)
			continue;		/* behind: keep walking */
		if (fabs(timediff(ep[0].time, t)) <= 1e-3) {
			*out = ep;
			return n;
		}
		s->cursor = save;		/* ahead: put it back for a later epoch */
		return 0;
	}
}

/* rnxopt_t's system index: {GPS,GLO,GAL,QZS,SBS,CMP,IRN}. */
static int rnx_sysidx(int sys) {
	switch (sys) {
	case SYS_GPS: return 0;
	case SYS_GLO: return 1;
	case SYS_GAL: return 2;
	case SYS_QZS: return 3;
	case SYS_SBS: return 4;
	case SYS_CMP: return 5;
	case SYS_IRN: return 6;
	}
	return -1;
}

/* Add one observation type (e.g. "C1C") to a system's list if not already present. */
static void rnx_add_tobs(rnxopt_t *o, int si, const char *type) {
	int k;

	if (si < 0)
		return;
	for (k = 0; k < o->nobs[si]; k++)
		if (strcmp(o->tobs[si][k], type) == 0)
			return;
	if (o->nobs[si] >= MAXOBSTYPE)
		return;
	snprintf(o->tobs[si][o->nobs[si]], sizeof o->tobs[si][0], "%s", type);
	o->nobs[si]++;
}

/*
 * Build the RINEX header options from the *master's* observations, restricted to the systems
 * the VBS carries.
 *
 * Not from the first synthesized epoch, which is what the first version did: that epoch held
 * four Galileo satellites -- convergence and the 0.05 m validity gate had admitted nothing
 * else yet -- so the header declared a Galileo-only file, and every GPS observation in the
 * 1 678 epochs that followed was written under a system with no observation types. The
 * observation types the VBS can carry are exactly the master's (§4.6: nothing can be
 * synthesized that the master did not observe) -- but NOT every one is present in every
 * epoch: a Block IIR satellite carries L2W and no L2C, a BeiDou-2 satellite B2I and no B2b,
 * and a first epoch with none of them in view produced a header without those types, after
 * which the writer silently dropped every such observation (D52: T60/T100/T150 identities
 * failed by exactly those records). The types are therefore collected from the master's
 * entire record set, which is in memory before the first VBS epoch is written.
 */
static void rnx_init(rnxopt_t *o, const obsd_t *synth, int n, const double *apppos,
	gtime_t t, int sys_mask) {
	int i, j;

	memset(o, 0, sizeof *o);
	o->rnxver = 304;
	o->ts = t;
	o->tstart = t;		/* whichever field this MRTKLIB writes as TIME OF FIRST OBS */
	o->navsys = 0;
	snprintf(o->prog, sizeof o->prog, "vbs_pipeline");
	snprintf(o->runby, sizeof o->runby, "vbs-engine");
	snprintf(o->marker, sizeof o->marker, "VBS");
	snprintf(o->rec[1], sizeof o->rec[1], "VIRTUAL");
	snprintf(o->ant[1], sizeof o->ant[1], "VIRTUAL         NONE");
	for (i = 0; i < 3; i++)
		o->apppos[i] = apppos[i];
	for (i = 0; i < n; i++) {
		int sys = satsys(synth[i].sat, NULL), si = rnx_sysidx(sys);
		char type[8];

		if (si < 0 || !(sys & sys_mask))
			continue;
		o->navsys |= sys;
		for (j = 0; j < NFREQ + NEXOBS; j++) {
			if (synth[i].code[j] == CODE_NONE)
				continue;
			snprintf(type, sizeof type, "C%s", code2obs(synth[i].code[j]));
			rnx_add_tobs(o, si, type);
			snprintf(type, sizeof type, "L%s", code2obs(synth[i].code[j]));
			rnx_add_tobs(o, si, type);
			snprintf(type, sizeof type, "S%s", code2obs(synth[i].code[j]));
			rnx_add_tobs(o, si, type);
		}
	}
	/* Code mask: every code the VBS carries, per system. */
	for (i = 0; i < n; i++) {
		int sys = satsys(synth[i].sat, NULL), si = rnx_sysidx(sys);

		if (si < 0 || !(sys & sys_mask))
			continue;
		for (j = 0; j < NFREQ + NEXOBS; j++)
			if (synth[i].code[j] != CODE_NONE && synth[i].code[j] < MAXCODE)
				o->mask[si][synth[i].code[j] - 1] = 1;
	}
	o->obstype = OBSTYPE_PR | OBSTYPE_CP | OBSTYPE_SNR;
	o->freqtype = FREQTYPE_ALL;
}

/*
 * The float and integer of (sat, ref, comb) on one baseline, re-referenced exactly as
 * vbs_ambres_loop_check does it (N(s,r) = N(s,k) - N(r,k)), so that a misclosure can be shown
 * next to the floats the three fixes were taken from.
 */
static int fix_lookup(const struct vbs_amb_fix *list, int n, int sat, int ref, int comb_id,
	double *value, double *integer, double *sigma) {
	const struct vbs_amb_fix *fs = NULL, *fr = NULL;
	int i;

	for (i = 0; i < n; i++)
		if (list[i].sat == sat && list[i].comb_id == comb_id)
			fs = &list[i];
	if (fs == NULL)
		return 0;
	if (fs->refsat == ref) {
		*value = fs->value; *integer = fs->integer; *sigma = fs->sigma;
		return 1;
	}
	for (i = 0; i < n; i++)
		if (list[i].sat == ref && list[i].comb_id == comb_id)
			fr = &list[i];
	if (fr == NULL || fr->refsat != fs->refsat)
		return 0;
	*value = fs->value - fr->value;
	*integer = fs->integer - fr->integer;
	*sigma = sqrt(fs->sigma * fs->sigma + fr->sigma * fr->sigma);
	return 1;
}

/* Copy records, leaving out satellites whose MW scatter at this station exceeds the gate. */
static int mw_gate_records(const struct vbs_preproc *pp, const obsd_t *in, int n, double gate,
	obsd_t *out, long *ngated) {
	int i, m = 0;

	for (i = 0; i < n; i++) {
		double sg = vbs_preproc_mw_sigma(pp, in[i].sat);

		if (!isnan(sg) && sg > gate) {
			(*ngated)++;
			continue;
		}
		out[m++] = in[i];
	}
	return m;
}

int main(int argc, char **argv) {
	static double rs[MAXOBS * 6], dts[MAXOBS * 2], var[MAXOBS];
	static double rs_sat[MAXSAT * 6];	/* the same geometry for every filter, by satellite */
	static int svh[MAXOBS];
	static struct vbs_ambres_stats ar;
	static struct vbs_interp_aux iaux[MAXBASE];
	static struct vbs_corr vcorr[MAXOBS];
	static struct vbs_interp_quality vqual[MAXOBS];
	static obsd_t synth[MAXOBS];
	static double azel[MAXSAT * 2];
	static struct vbs_slip mslip[MAXOBS];
	int nmslip = 0;
	/*
	 * §5.3 signalling on the output. Case 1 comes from the master's station-level
	 * preprocessor; case 2 is derived here by comparing each master baseline's fixed
	 * integers against the previous epoch's -- a changed integer, or an ambiguity that was
	 * reset (arc change) and is no longer fixed, is a discontinuity the rover must be told
	 * about. Cases 3 and 4 are vbs_lock's own.
	 */
	struct vbs_lock *lock = NULL;
	struct vbs_lock_opt lopt;
	struct vbs_lock_stats lock_st;
	static rnxopt_t ropt;
	int rnx_header_written = 0;
#define NCOMB	16
	static double last_fix[MAXBASE][MAXSAT][NCOMB];
	/*
	 * Misclosure *events*. A wrong fix is held for the life of the constraint, so per-epoch
	 * counts multiply one event by however many epochs it persists -- T30's 21 % extra-wide
	 * rate was, in its first 28 minutes, two events. An event begins when a (triangle,
	 * satellite, combination) miscloses after not misclosing at the previous epoch, or when
	 * the misclosure's integer value changes.
	 */
	static long mis_last_epoch[MAXPAIR][MAXSAT][NCOMB];
	static double mis_last_int[MAXPAIR][MAXSAT][NCOMB];
	long mis_onsets[VBS_AMB_NSTEP] = {0, 0, 0};
	long mis_onsets_sys[VBS_AMB_NSTEP][4] = {{0}};	/* by constellation: GPS, GAL, BDS, other */
	static long mis_onsets_sat[VBS_AMB_NSTEP][MAXSAT];	/* by satellite, for the PRN listing */
	static long last_arc_m[MAXSAT];
	static int refixed[MAXSAT], nrefixed, mslip_sat[MAXOBS], nmslip_sat;
	const char *navpath = NULL, *mode = "loo", *outpath = NULL;
	/*
	 * --sp3 (repeatable) / --sat-atx: precise orbits in place of the broadcast ones for every
	 * filter and the synthesis, the rover keeping its broadcast set -- what a server fed by an
	 * SSR stream such as Galileo HAS would serve (D163). See vbs_peph.h for why the two runs
	 * differ by the positions alone.
	 */
	const char *sp3path[VBS_PEPH_MAXFILE], *sat_atx = NULL;
	int nsp3 = 0;
	struct vbs_peph_stats sp3st;
	double vbs_xyz[3] = {0}, offset[3] = {0}, vbs_eff[3];
	int have_vbs = 0, have_target = 0, zero_corr = 0, no_ar = 0, i, k, master = -1;
	double min_ratio = 0.0;
	int paired = 0, no_tropo = 0, no_iono = 0, loop_closure = 0;
	/*
	 * --iono-sys GREC: apply the dispersive correction only to these constellations, the others
	 * served as --no-iono serves them (D167). The UM980's code-levelled ionosphere carries
	 * per-satellite constants of 0.15-0.35 m on GPS L2 and 0.5-2.3 m on GLONASS where a Mosaic's
	 * carry a tenth of that (D128); Galileo and BeiDou have neither a P(Y) column nor FDMA. If
	 * the damage is confined to those two, the clean constellations can still carry the network.
	 */
	int iono_sys_mask = 0;		/* 0: every constellation */
	/*
	 * --nondisp / --no-nondisp: apply the *estimated* non-dispersive residual (the
	 * interpolated m_w*dZTD plus the troposphere gradients) on top of the deterministic
	 * model term. -1 leaves the library default.
	 *
	 * Until now only --no-nondisp existed, and the library default is 0, so the flag was a
	 * no-op and the residual could not be switched on from here at all: the filter estimated
	 * it, vbs_netekf put the gradient into the emitted dT (D69), vbs_interp interpolated it,
	 * and then it was multiplied by zero. D33 item 1 asks for it to be applied only where it
	 * measurably reduces the error -- which needs it to be reachable to measure.
	 */
	int nondisp = -1;
	/*
	 * --corr-max-sigma: the interpolation's gate on whether a correction is offered at all
	 * (vbs_interp_opt.max_sigma, default 0.05 m). It had no flag, so the one parameter that
	 * decides between withholding a correction and offering a doubtful one was unreachable
	 * from here -- the same shape of gap as `apply_nondispersive` before D78.
	 *
	 * It is the lever D86 identifies for long range: M80 at 62 km fixes 42.8 % of epochs with
	 * 51.81 % of them wrong, because corrections pass this gate while still carrying error,
	 * whereas a cluster whose corrections mostly fail it fixes rarely and is right. Beyond
	 * ~45 km integrity is worth more than availability, so the gate should tighten with range.
	 */
	double corr_max_sigma = -1.0;
	int ewl_dump = 0, ewl_dumped = 0;
	/*
	 * --corr-dump: the per-satellite, per-auxiliary-baseline correction table a
	 * bisect needs. A single T30 baseline converges on its own (vbs_netrun) while the T30
	 * *network* produces no usable dispersive correction, so the question is what the three
	 * baselines are saying about the same satellite at the same epoch, and which of them the
	 * 0.05 m sigma gate is throwing away. One line per (epoch, baseline, satellite), with the
	 * pre-gate correction and its sigma, so the gate's effect is visible rather than implied.
	 */
	int corr_dump = 0;
	long corr_dump_max = 0, corr_dump_n = 0;
	double ar_corr_time_s = 0.0;
	double snr_mask = 0.0;		/* --snr-mask dB-Hz: D51 C/N0 mask on every filter; 0 = off */
	int prefer_l2c = 0;		/* --prefer-l2c: D53, GPS L2 from L2C rather than semi-codeless L2W */
	/*
	 * --amb-dump: every fix with the float it came from. The float network on the Unicore
	 * cases is healthy while its integers are not, which is the signature of a fractional
	 * phase bias rather than of noise: a DD float sitting near a half integer is rounded one
	 * way on one baseline and the other way on the next, and both pass the ratio test. The
	 * fractional part of (float - integer) says so directly -- a peak at 0 is noise, a peak
	 * at +-0.5 is a bias, and nothing else in the pipeline can tell them apart.
	 */
	int amb_dump = 0;
	long amb_dump_max = 0, amb_dump_n = 0;
	/*
	 * --base-resid i,E,N,U: how far base i's declared position was measured to sit from
	 * where its neighbours put it (D56). The engine cannot measure this itself -- in
	 * production it comes from the leave-one-out integrity monitor (D59); here it comes from
	 * an external coordinate check. A base outside the network's published acceptance gate is
	 * then refused by vbs_base_check and never enters a cluster.
	 */
	double base_resid[MAXBASE][3] = {{0}};
	int base_resid_known[MAXBASE] = {0};
	/*
	 * --l2-bias i,METRES: base i's own inter-signal code bias on GPS L2, measured as the
	 * median of C2W - C2X over a day (D159, l2pair.py), and subtracted from its L2 pseudorange
	 * before anything reads it.
	 *
	 * Why this and not --dd-datum. --dd-datum removes the constant of each BASELINE's estimated
	 * dispersive term, which also destroys any genuine common-mode ionosphere on that baseline.
	 * This removes a constant belonging to one STATION, which is instrumental by construction
	 * and cannot contain ionosphere. D155 refuted the first; the second has never been tried.
	 *
	 * Measured across 79 two-column stations the bias spans 3.09 m (-1.85 to +1.24), the same
	 * order as D136's four metres of per-baseline datum, and drifts 0.15-0.36 m over 44 days --
	 * a 15:1 ratio, so it is dominated by a per-station constant. The expected residual after
	 * removal is therefore ~0.2 m, which is still 28x the 7 mm D146 showed changes integers.
	 */
	double l2_bias[MAXBASE] = {0};
	int l2_bias_known[MAXBASE] = {0};
	/*
	 * --c2w-cal MODE: each UM980 base's own GPS L2 P(Y) offset (2W-2X), measured from its own
	 * data rather than read from a table (vbs_c2wcal.h). D160's table was one snapshot applied to
	 * days 64 days apart, and D161 measured units drifting metres in days, so here nothing is
	 * configured:
	 *   day      one constant per base, the mean over everything loaded (acausal: the best a
	 *            fresh constant can do, not something a server can run);
	 *   run:SEC  the mean over a trailing SEC window, re-read every epoch (causal, what vbsd runs);
	 *   l2c:SEC  the L2C code itself in the P(Y) slot wherever both are tracked, the trailing
	 *            estimate only on L2C-less satellites (SEC 0: the acausal day mean there).
	 * Gated by vbs_c2wcal's bit-identity test pooled per firmware build, and withheld from the
	 * whole cluster unless every UM980 member can take it (--c2w-partial lifts that): a cluster
	 * corrected in part is internally inconsistent, D148's lesson, and D160's U29 was exactly that.
	 */
	int c2w_mode = 0, c2w_partial = 0, c2w_raw_out = 0;	/* mode 0 off, 1 constant, 2 l2c */
	double c2w_win = 0.0;
	/*
	 * D66: the slant-ionosphere random walk. The default, 1e-6 m2/s, is 5.5 mm per 30 s --
	 * defensible at 30 km and far too tight at 100, where the differential ionosphere moves
	 * much faster than that. A model that constrains the ionosphere harder than physics
	 * allows pushes the real variation into the ambiguities, which is the ionosphere acting
	 * as an input to resolution instead of an output of it.
	 */
	double q_ion = 0.0, var_ion0 = 0.0;
	/*
	 * --ion-from-fixes (D66): take the dispersive correction from the resolved integers via
	 * the geometry-free combination, and make "the ambiguities are fixed" the validity test
	 * instead of "the ionosphere state's variance came below 5 cm". The old test cannot be
	 * met at network range and is not what the literature asks of a correction.
	 */
	int ion_from_fixes = 0, if_nl = 0, if_min_epochs = 0;
	/*
	 * -1 means "leave the module's own default alone". Writing the CLI variable in
	 * unconditionally silently overrode a default the library had just been given, so
	 * switching gradients on by default changed nothing and the measurement said so.
	 */
	int est_gradients = -1;
	double grad_sigma = 0.0, max_mapg = 0.0, iono_scale = 0.0;
	double iono_gate = 0.0;		/* --iono-gate <m>: vbs_interp_opt.iono_gate_m, D165 */
	/* --iono-l5, --iono-gps-scale, --iono-model, --iono-model-check, --out-drop,
	 * --iono-shrink: vbs_ionopost */
	struct vbs_ionopost_opt ipo;
	struct vbs_ionopost_stats ips;
	int anchor_master = -1, max_aux = 0, method = -1;
	double corr_len = 0.0, sig_field = 0.0;
	double if_max_sigma = 0.0, if_sys = -1.0, if_max_frac = 0.0;
	long nrep_total = 0, ncorr_total = 0;
	int ar_max_step = -1;		/* --ar-max-step ewl|wl|nl: stop the cascade short (diagnostic) */
	/*
	 * --cascade-order table restores the pre-D74 row order, so the reordering can be
	 * measured against what it replaced from this same binary. -1 leaves the library
	 * default (wavelength order), which is the correct behaviour.
	 */
	int cascade_by_lambda = -1;
	int use_band6 = -1;		/* --no-band6: drop Galileo E6 / BeiDou B3I (measurement) */
	/*
	 * --glonass / --no-glonass. -1 decides from the cluster: GLONASS is admitted only when
	 * every station in it reports the same receiver type.
	 *
	 * GLONASS is FDMA, so each satellite sits on its own frequency and the receiver's
	 * hardware delay is frequency-dependent. In a double difference between two receivers
	 * those delays cancel exactly only if they match; between different models they leave an
	 * inter-channel bias that makes the ambiguity non-integer, and this network mixes u-blox
	 * F9P, Unicore UM980 and Septentrio Mosaic-X5. Calibrating the bias per receiver pair is
	 * the standard remedy and is not attempted here, so the gate is homogeneity instead.
	 *
	 * Float GLONASS corrections are not a fallback: D24 measured float slant ionosphere at
	 * 30-370 mm against an interpolation budget of 14-125 mm, so a mixed cluster should drop
	 * GLONASS rather than serve it badly.
	 */
	int use_glonass = -1;
	/*
	 * Feed the per-station signal counts to master selection, so its documented
	 * "prefer a richer master" term is non-zero (D80). **Default off: measured harmful.**
	 *
	 * The defect it addresses is real -- the term was always multiplied by a NULL nsignal, so
	 * selection was pure proximity -- but preferring the richer master outright is the wrong
	 * criterion, and MIX_F9P_UM9 shows why: promoting a UM980 (11 signals) over an F9P
	 * (6) among three F9P auxiliaries took the rover from 78.8 % fixed / 2.49 % false to
	 * 32.2 % / 21.94 %, with V95 8.9 -> 84.8 cm.
	 *
	 * The master's signal set is the VBS's ceiling, but a correction for a signal needs the
	 * *auxiliaries* to carry it as well. Raising the master above its auxiliaries therefore
	 * lets the VBS synthesize bands no baseline can correct, and those reach the rover with
	 * absent or unconverged corrections. The right score is over the signals common to the
	 * master and the auxiliaries it will actually use, which is a different computation and
	 * needs its own measurement.
	 *
	 * 1 restores the richest-master behaviour for that measurement; --master-nearest is a
	 * no-op while this is 0, kept so the flag still reads as an explicit choice.
	 */
	int master_by_signals = 0;
	int dd_datum = 0;		/* --dd-datum: D55, remove each auxiliary's ionosphere datum before the plane fit */
	int no_bds = 0, no_band5 = 0, no_gal = 0;	/* D53 diagnostics: network without BeiDou / band 5 / Galileo */
	int no_band2 = 0;			/* --no-band2: GPS L2 out of the network filter, ionosphere from L1/L5 */
	double sig_code = 0.0;			/* --sig-code: zenith code sigma of the network filter, m (0 = default) */
	double sig_code_gps = 0.0, q_ion_gps = 0.0;	/* --sig-code-gps, --q-ion-gps: D168 */
	int no_wl12_sys = 0;			/* --no-wl12[-gps|-glo]: D129, skip that system's L1-L2 wide lane.
					 * The split flags are the control: dropping it for GLONASS alone also
					 * removes GLONASS's ONLY wide lane, so "no-wl12 helps" and "GLONASS
					 * ambiguity resolution was the problem" are different claims. */
	int wl12_derive_sys = 0;		/* --wl12-derive: D131, N(1-2) = N(1-5) - N(2-5) */
	int one_ref_sys = 0;			/* --one-ref: D132, one reference satellite per constellation */
	double mw_gate_cyc = 0.0;	/* --mw-gate: D53, drop a station-satellite whose MW scatter over the arc exceeds this */
	long mw_gated = 0;
	static obsd_t gated_m[MAXOBS], gated_a[MAXBASE][MAXOBS];
	double out_elmask_deg = 0.0;	/* --out-elmask: D52, do not emit satellites below this at the VBS */
	double aux_max_vbs_km = 0.0;	/* --aux-max-vbs-km: the server sets the policy radius */
	long out_elmasked = 0;	/* --ar-corr-time: D48 search-side covariance inflation; 0 = off */
	int identity = 0;
	/*
	 * Nominal observation interval, for the preprocessor. D22 made the geometry-free slip
	 * threshold scale with the interval -- at 30 s the ionosphere has thirty times longer to
	 * move than at 1 s, and a fixed threshold either flags the ionospheric tail as slips or
	 * misses real ones -- so a 1 Hz run must say so. Defaults to the dataset's 30 s.
	 */
	double interval_s = 30.0;
	/*
	 * Auxiliary<->auxiliary filters, one per consecutive pair, each estimating that baseline
	 * *independently*. That independence is the whole point: if the third baseline's
	 * ambiguities were derived from the other two, the closure identity would hold by
	 * construction and the check would confirm nothing.
	 */
	struct vbs_netekf *aa_ekf[MAXPAIR];
	static struct vbs_amb_fix aa_fix[MAXPAIR][256];
	int aa_a[MAXPAIR], aa_b[MAXPAIR], aa_n[MAXPAIR], naa = 0;
	/*
	 * Per-satellite blame accounting. With every auxiliary pair filtered, each baseline sits
	 * in several triangles, so a baseline that is wrong fails *all* of its triangles while an
	 * innocent one fails only those sharing the guilty baseline. The baseline common to every
	 * failing triangle is therefore the culprit, and unfixing only it costs one satellite on
	 * one baseline instead of three -- which matters because unfixing drives availability
	 * (D30 item 5).
	 */
	static unsigned char fail_tri[MAXSAT][MAXPAIR + MAXBASE];
	static int fail_sats[MAXSAT], nfail_sats;
	struct vbs_amb_loop_stats loop_stats;
	long loop_unfixed = 0, loop_isolated = 0, fixes_total = 0;
	long mis_arc_mismatch = 0, mis_same_arc = 0, blame_master = 0, blame_auxaux = 0;
	/*
	 * Report-only by default: a check that is itself under suspicion must not unfix
	 * good ambiguities, or it biases the very measurement meant to judge it.
	 */
	int report_only = 1;
	int q2;
	struct vbs_interp_opt iopt;
	struct vbs_ddres *ddref = NULL;
	struct vbs_ambres_opt aopt;
	struct vbs_netekf_opt nopt;
	struct vbs_basemeta_opt bopt;
	struct vbs_ddres *dd = NULL;
	gtime_t t0 = {0, 0.0};
	long nepoch = 0, nsat_total = 0, nepoch_out = 0;
	long nskip_aux = 0, nskip_corr = 0, nsat_min = 9999, nsat_max = 0;
	long nsat_premask_total = 0;
	FILE *out = NULL;

	memset(&ipo, 0, sizeof ipo);
	memset(&ips, 0, sizeof ips);
	for (i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--nav") && i + 1 < argc)
			navpath = argv[++i];
		else if (!strcmp(argv[i], "--sp3") && i + 1 < argc && nsp3 < VBS_PEPH_MAXFILE)
			sp3path[nsp3++] = argv[++i];
		else if (!strcmp(argv[i], "--sat-atx") && i + 1 < argc)
			sat_atx = argv[++i];
		else if (!strcmp(argv[i], "--mode") && i + 1 < argc)
			mode = argv[++i];
		else if (!strcmp(argv[i], "--out") && i + 1 < argc)
			outpath = argv[++i];
		else if (!strcmp(argv[i], "--vbs") && i + 1 < argc)
			have_vbs = parse_xyz(argv[++i], vbs_xyz) == 0;
		else if (!strcmp(argv[i], "--no-corrections"))
			zero_corr = 1;
		else if (!strcmp(argv[i], "--interval") && i + 1 < argc)
			interval_s = atof(argv[++i]);
		else if (!strcmp(argv[i], "--loop-closure"))
			loop_closure = 1;
		else if (!strcmp(argv[i], "--aux-max-vbs-km") && i + 1 < argc)
			aux_max_vbs_km = atof(argv[++i]);
		else if (!strcmp(argv[i], "--out-elmask") && i + 1 < argc)
			out_elmask_deg = atof(argv[++i]);
		else if (!strcmp(argv[i], "--mw-gate") && i + 1 < argc)
			mw_gate_cyc = atof(argv[++i]);
		else if (!strcmp(argv[i], "--no-bds"))
			no_bds = 1;
		else if (!strcmp(argv[i], "--gps-only"))
			no_bds = no_gal = 1;
		else if (!strcmp(argv[i], "--one-ref"))
			one_ref_sys |= SYS_GPS | SYS_GAL | SYS_CMP | SYS_GLO;
		else if (!strcmp(argv[i], "--wl12-derive"))
			wl12_derive_sys |= SYS_GPS;
		else if (!strcmp(argv[i], "--no-wl12"))
			no_wl12_sys |= SYS_GPS | SYS_GLO;
		else if (!strcmp(argv[i], "--no-wl12-gps"))
			no_wl12_sys |= SYS_GPS;
		else if (!strcmp(argv[i], "--no-wl12-glo"))
			no_wl12_sys |= SYS_GLO;
		else if (!strcmp(argv[i], "--sig-code") && i + 1 < argc)
			sig_code = atof(argv[++i]);
		else if (!strcmp(argv[i], "--sig-code-gps") && i + 1 < argc)
			sig_code_gps = atof(argv[++i]);
		else if (!strcmp(argv[i], "--q-ion-gps") && i + 1 < argc)
			q_ion_gps = atof(argv[++i]);
		else if (!strcmp(argv[i], "--no-band2"))
			no_band2 = 1;
		else if (!strcmp(argv[i], "--no-band5"))
			no_band5 = 1;
		else if (!strcmp(argv[i], "--l2-bias") && i + 1 < argc) {
			int bi2 = -1;
			double v2 = 0.0;

			i++;
			if (sscanf(argv[i], "%d,%lf", &bi2, &v2) != 2 ||
			    bi2 < 0 || bi2 >= MAXBASE) {
				fprintf(stderr, "bad --l2-bias %s\n", argv[i]);
				return 2;
			}
			l2_bias[bi2] = v2;
			l2_bias_known[bi2] = 1;
		}
		else if (!strcmp(argv[i], "--c2w-cal") && i + 1 < argc) {
			const char *mv = argv[++i];

			if (!strcmp(mv, "day")) {
				c2w_mode = 1;
				c2w_win = 0.0;
			} else if (!strncmp(mv, "run:", 4) && atof(mv + 4) > 0.0) {
				c2w_mode = 1;
				c2w_win = atof(mv + 4);
			} else if (!strncmp(mv, "l2c:", 4) && atof(mv + 4) >= 0.0) {
				c2w_mode = 2;
				c2w_win = atof(mv + 4);
			} else {
				fprintf(stderr, "bad --c2w-cal %s (day, run:SEC, l2c:SEC)\n", mv);
				return 2;
			}
		}
		else if (!strcmp(argv[i], "--c2w-partial"))
			c2w_partial = 1;
		else if (!strcmp(argv[i], "--c2w-raw-out"))
			c2w_raw_out = 1;
		else if (!strcmp(argv[i], "--base-resid") && i + 1 < argc) {
			int bi;
			double e2, n2, u2;

			if (sscanf(argv[++i], "%d,%lf,%lf,%lf", &bi, &e2, &n2, &u2) != 4 ||
			    bi < 0 || bi >= MAXBASE) {
				fprintf(stderr, "bad --base-resid %s\n", argv[i]);
				return 1;
			}
			base_resid[bi][0] = e2;
			base_resid[bi][1] = n2;
			base_resid[bi][2] = u2;
			base_resid_known[bi] = 1;
		} else if (!strcmp(argv[i], "--amb-dump") && i + 1 < argc) {
			amb_dump = 1;
			amb_dump_max = atol(argv[++i]);
		} else if (!strcmp(argv[i], "--ar-max-step") && i + 1 < argc) {
			const char *v = argv[++i];

			ar_max_step = !strcmp(v, "ewl") ? VBS_AMB_STEP_EWL :
			    !strcmp(v, "wl") ? VBS_AMB_STEP_WL : VBS_AMB_STEP_NL;
		} else if (!strcmp(argv[i], "--cascade-order") && i + 1 < argc) {
			const char *v = argv[++i];

			cascade_by_lambda = strcmp(v, "table") != 0;
		} else if (!strcmp(argv[i], "--master-nearest")) {
			master_by_signals = 0;
		} else if (!strcmp(argv[i], "--master-richest")) {
			master_by_signals = 1;
		} else if (!strcmp(argv[i], "--glonass")) {
			use_glonass = 1;
		} else if (!strcmp(argv[i], "--no-glonass")) {
			use_glonass = 0;
		} else if (!strcmp(argv[i], "--no-band6")) {
			use_band6 = 0;
		} else if (!strcmp(argv[i], "--dd-datum"))
			dd_datum = 1;
		else if (!strcmp(argv[i], "--ion-from-fixes"))
			ion_from_fixes = 1;
		else if (!strcmp(argv[i], "--anchor-master"))
			anchor_master = 1;
		else if (!strcmp(argv[i], "--iono-scale") && i + 1 < argc)
			iono_scale = atof(argv[++i]);
		else if (!strcmp(argv[i], "--iono-gate") && i + 1 < argc)
			iono_gate = atof(argv[++i]);
		else if (!strcmp(argv[i], "--collocation"))
			method = VBS_INTERP_METHOD_COLLOCATION;
		else if (!strcmp(argv[i], "--plane-fit"))
			method = VBS_INTERP_METHOD_PLANE;
		else if (!strcmp(argv[i], "--corr-len") && i + 1 < argc)
			corr_len = atof(argv[++i]);
		else if (!strcmp(argv[i], "--sig-field") && i + 1 < argc)
			sig_field = atof(argv[++i]);
		else if (!strcmp(argv[i], "--max-aux") && i + 1 < argc)
			max_aux = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--gradients"))
			est_gradients = 1;
		else if (!strcmp(argv[i], "--no-gradients"))
			est_gradients = 0;
		else if (!strcmp(argv[i], "--no-anchor-master"))
			anchor_master = 0;
		else if (!strcmp(argv[i], "--grad-sigma") && i + 1 < argc)
			grad_sigma = atof(argv[++i]);
		else if (!strcmp(argv[i], "--max-mapg") && i + 1 < argc)
			max_mapg = atof(argv[++i]);
		else if (!strcmp(argv[i], "--if-narrow-lane"))
			if_nl = 1;
		else if (!strcmp(argv[i], "--if-min-epochs") && i + 1 < argc)
			if_min_epochs = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--corr-max-sigma") && i + 1 < argc)
			corr_max_sigma = atof(argv[++i]);
		else if (!strcmp(argv[i], "--if-max-sigma") && i + 1 < argc)
			if_max_sigma = atof(argv[++i]);
		else if (!strcmp(argv[i], "--if-sys") && i + 1 < argc)
			if_sys = atof(argv[++i]);
		else if (!strcmp(argv[i], "--if-max-frac") && i + 1 < argc)
			if_max_frac = atof(argv[++i]);
		else if (!strcmp(argv[i], "--q-ion") && i + 1 < argc)
			q_ion = atof(argv[++i]);
		else if (!strcmp(argv[i], "--var-ion0") && i + 1 < argc)
			var_ion0 = atof(argv[++i]);
		else if (!strcmp(argv[i], "--prefer-l2c"))
			prefer_l2c = 1;
		else if (!strcmp(argv[i], "--snr-mask") && i + 1 < argc)
			snr_mask = atof(argv[++i]);
		else if (!strcmp(argv[i], "--ar-corr-time") && i + 1 < argc)
			ar_corr_time_s = atof(argv[++i]);
		else if (!strcmp(argv[i], "--corr-dump") && i + 1 < argc) {
			corr_dump = 1;
			corr_dump_max = atol(argv[++i]);
		} else if (!strcmp(argv[i], "--ewl-dump"))
			ewl_dump = 1;
		else if (!strcmp(argv[i], "--loop-unfix"))
			report_only = 0;
		else if (!strcmp(argv[i], "--no-iono"))
			no_iono = 1;
		else if (!strcmp(argv[i], "--iono-sys") && i + 1 < argc) {
			const char *c;

			for (c = argv[++i]; *c; c++)
				iono_sys_mask |= *c == 'G' ? SYS_GPS : *c == 'R' ? SYS_GLO :
				    *c == 'E' ? SYS_GAL : *c == 'C' ? SYS_CMP : *c == 'J' ? SYS_QZS : 0;
			if (iono_sys_mask == 0) {
				fprintf(stderr, "bad --iono-sys %s (letters of GREC)\n", argv[i]);
				return 2;
			}
		}
		else if (!strcmp(argv[i], "--iono-l5") && i + 1 < argc) {
			if ((ipo.need_band5_sys = vbs_ionopost_parse_sys(argv[++i])) <= 0) {
				fprintf(stderr, "bad --iono-l5 %s (letters of GREC)\n", argv[i]);
				return 2;
			}
		} else if (!strcmp(argv[i], "--iono-gps-scale") && i + 1 < argc)
			ipo.scale_gps = atof(argv[++i]);
		else if (!strcmp(argv[i], "--iono-model") && i + 1 < argc) {
			/* LETTERS[:plane|idw|quad] */
			char buf[32], *colon;

			snprintf(buf, sizeof buf, "%s", argv[++i]);
			if ((colon = strchr(buf, ':')) != NULL) {
				*colon++ = '\0';
				ipo.model_mode = !strcmp(colon, "idw") ? VBS_IONOPOST_IDW :
				    !strcmp(colon, "quad") ? VBS_IONOPOST_QUAD :
				    !strcmp(colon, "plane") ? VBS_IONOPOST_PLANE : -1;
			}
			if ((ipo.model_sys = vbs_ionopost_parse_sys(buf)) <= 0 || ipo.model_mode < 0) {
				fprintf(stderr, "bad --iono-model %s (GREC[:plane|idw|quad])\n", argv[i]);
				return 2;
			}
		} else if (!strcmp(argv[i], "--iono-model-check") && i + 1 < argc)
			ipo.model_check_m = atof(argv[++i]);
		else if (!strcmp(argv[i], "--iono-shrink") && i + 1 < argc)
			ipo.shrink_m = atof(argv[++i]);
		else if (!strcmp(argv[i], "--out-drop") && i + 1 < argc) {
			if ((ipo.drop_sys = vbs_ionopost_parse_sys(argv[++i])) <= 0) {
				fprintf(stderr, "bad --out-drop %s (letters of GREC)\n", argv[i]);
				return 2;
			}
		}
		else if (!strcmp(argv[i], "--no-tropo"))
			no_tropo = 1;
		else if (!strcmp(argv[i], "--nondisp"))
			nondisp = 1;
		else if (!strcmp(argv[i], "--no-nondisp"))
			nondisp = 0;
		else if (!strcmp(argv[i], "--paired"))
			paired = 1;
		else if (!strcmp(argv[i], "--no-ar"))
			no_ar = 1;
		else if (!strcmp(argv[i], "--ratio") && i + 1 < argc)
			min_ratio = atof(argv[++i]);
		else if (!strcmp(argv[i], "--offset") && i + 1 < argc)
			parse_xyz(argv[++i], offset);
		else if (!strcmp(argv[i], "--target") && i + 1 < argc)
			have_target = parse_site(argv[++i], &target) == 0;
		else if (!strcmp(argv[i], "--base") && i + 1 < argc && nbase < MAXBASE) {
			if (parse_site(argv[++i], &bases[nbase]) == 0)
				nbase++;
		} else {
			fprintf(stderr, "unknown or malformed argument: %s\n", argv[i]);
			return 1;
		}
	}
	if (navpath == NULL || !have_vbs || nbase < 2) {
		fprintf(stderr, "usage: %s --nav <nav.rnx> --vbs <x,y,z> [--mode loo|emit] "
		    "[--offset <e,n,u>] [--target <file,x,y,z>] [--out <f>] [--ar-corr-time <s>] "
		    "--base <file,x,y,z> ...\n", argv[0]);
		return 1;
	}

	memset(&nav, 0, sizeof nav);
	/*
	 * §5.3 output gaps are judged against the *nominal* interval. The default 10 s gap
	 * threshold is a 1 Hz number: at 30 s it declared every epoch a gap and stamped LLI on
	 * 100 % of the emitted phases, which turned the rover's solution code-only -- and the
	 * 30 s rover figures were the harness identity's only, so nothing caught it (D48).
	 */
	vbs_lock_default_opt(&lopt);
	if (3.0 * interval_s > lopt.max_gap_s)
		lopt.max_gap_s = 3.0 * interval_s;
	if ((lock = vbs_lock_new(&lopt)) == NULL) {
		fprintf(stderr, "vbs_lock_new failed\n");
		return 1;
	}

	for (i = 0; i < nbase; i++) {
		memset(&bases[i].obs, 0, sizeof bases[i].obs);
		if (readrnxt(bases[i].path, i + 1, t0, t0, 0.0, "", &bases[i].obs, &nav,
		    &bases[i].sta) < 0 || bases[i].obs.n <= 0) {
			fprintf(stderr, "cannot read %s\n", bases[i].path);
			return 1;
		}
		if (l2_bias_known[i] && l2_bias[i] != 0.0) {
			int kb, jb, nfix2 = 0;

			for (kb = 0; kb < bases[i].obs.n; kb++) {
				obsd_t *o = &bases[i].obs.data[kb];

				if (satsys(o->sat, NULL) != SYS_GPS)
					continue;
				for (jb = 0; jb < NFREQ + NEXOBS; jb++) {
					if (o->code[jb] == CODE_NONE || o->P[jb] == 0.0)
						continue;
					if (vbs_signal_band(o->code[jb]) != 2)
						continue;
					o->P[jb] -= l2_bias[i];
					nfix2++;
				}
			}
			printf("  base %d (%s): GPS L2 code bias %+.3f m removed from %d observables\n",
			    i, bases[i].sta.name[0] ? bases[i].sta.name : "?", l2_bias[i], nfix2);
		}
	}
	if (have_target) {
		memset(&target.obs, 0, sizeof target.obs);
		if (readrnxt(target.path, MAXBASE + 1, t0, t0, 0.0, "", &target.obs, &nav,
		    NULL) < 0 || target.obs.n <= 0) {
			fprintf(stderr, "cannot read target %s\n", target.path);
			return 1;
		}
	}
	if (readrnxt(navpath, 1, t0, t0, 0.0, "", NULL, &nav, NULL) < 0) {
		fprintf(stderr, "cannot read nav %s\n", navpath);
		return 1;
	}
	/*
	 * readrnxt() returns 0, not -1, for a path that does not exist. With no ephemerides every
	 * satellite position fails and the run "succeeds" with zero VBS epochs: a day-251 sweep
	 * whose tree lacked its navigation file reported 162 clean exits and no output.
	 */
	if (nav.n + nav.ng <= 0) {
		fprintf(stderr, "no broadcast ephemerides in %s\n", navpath);
		return 1;
	}
	memset(&sp3st, 0, sizeof sp3st);
	if (nsp3 > 0) {
		int ne = vbs_peph_load(&nav, sp3path, nsp3, sat_atx);

		if (ne < 0)
			return 1;
		printf("# precise orbits: %d SP3 epochs from %d file(s), satellite antennas from %s\n",
		    ne, nsp3, sat_atx);
	}

	/*
	 * The constellations in play, needed by the signal count below. `nopt` is not built
	 * until after cluster selection, so the mask is derived here from the same flags.
	 */
	int nopt_sys_mask = SYS_GPS | SYS_GAL | SYS_CMP;

	if (no_bds)
		nopt_sys_mask &= ~SYS_CMP;
	if (no_gal)
		nopt_sys_mask &= ~SYS_GAL;

	/*
	 * GLONASS admission (D90). Admitted only when every station in the cluster reports the
	 * same receiver type, because its FDMA inter-channel biases cancel in a double
	 * difference only between matching hardware. `--glonass` forces it on for measurement,
	 * `--no-glonass` off.
	 *
	 * A VBS that carries no GLONASS makes the rover's GLONASS unusable -- RTK needs the base
	 * observation for every satellite the rover wants -- so a user moved from a real base to
	 * a virtual one loses a constellation. That is the reason to serve it, and it is why the
	 * gate is on the cluster rather than a global default.
	 */
	{
		int homogeneous = nbase > 0, i2;

		for (i2 = 1; i2 < nbase; i2++)
			if (strcmp(bases[i2].sta.rectype, bases[0].sta.rectype) != 0)
				homogeneous = 0;
		if (use_glonass > 0 || (use_glonass < 0 && homogeneous))
			nopt_sys_mask |= SYS_GLO;
		printf("# GLONASS: %s (cluster receivers %s", 
		    (nopt_sys_mask & SYS_GLO) ? "ADMITTED" : "excluded",
		    homogeneous ? "match" : "differ");
		for (i2 = 0; i2 < nbase && i2 < MAXBASE; i2++)
			printf("%s B%d=%s", i2 ? "," : ": ", i2,
			    bases[i2].sta.rectype[0] ? bases[i2].sta.rectype : "?");
		printf(")%s\n", use_glonass >= 0 ? " [forced by flag]" : "");
	}

	report_band_capability();
	report_band_coverage();

	/*
	 * Cluster selection over the supplied bases. The target is never among them: in
	 * leave-one-out it is what the network is being asked to predict, and including it
	 * would produce an excellent result that means nothing.
	 */
	{
		struct vbs_base gated[MAXBASE];
		struct vbs_cluster cl;
		int r;

		int nsig[MAXBASE];

		for (i = 0; i < nbase; i++) {
			char nm[16];

			snprintf(nm, sizeof nm, "B%d", i);
			vbs_base_init(&gated[i], nm, "RGF93", bases[i].xyz);
			if (base_resid_known[i])
				vbs_base_set_coord_resid(&gated[i], base_resid[i]);
			vbs_base_check(&gated[i], "RGF93", t0, NULL);
			nsig[i] = station_nsignal(&bases[i], nopt_sys_mask);
			if (gated[i].reason != VBS_BASE_OK)
				printf("# base B%d refused: %s\n", i,
				    vbs_base_reason_str(gated[i].reason));
		}
		struct vbs_cluster_opt copt;

		vbs_cluster_default_opt(&copt);
		if (max_aux > 0)
			copt.max_aux = max_aux;
		copt.aux_max_vbs_km = aux_max_vbs_km;	/* as the server: a far reserve is never an auxiliary */
		if ((r = vbs_cluster_select(vbs_xyz, gated, nbase,
		    master_by_signals ? nsig : NULL, &copt, &cl)) !=
		    VBS_CLUSTER_OK) {
			fprintf(stderr, "no cluster: %s\n", vbs_cluster_reason_str(r));
			return 1;
		}
		master = cl.master;
		/*
		 * Honour the selection, not just its master. Until this was added the pipeline
		 * read `cl.master` and then used *every* supplied base as an auxiliary, so no
		 * gate verdict ever excluded one: a base refused for a frame mismatch (D13), for
		 * staleness, or for a bad coordinate (D56) still contributed a filter, a
		 * correction and a place in the interpolation. The symptom was a gated run whose
		 * rover output was identical to the ungated one to the digit, on a fresh emit.
		 */
		for (i = 0; i < nbase; i++)
			used[i] = 0;
		used[cl.master] = 1;
		for (i = 0; i < cl.naux; i++)
			used[cl.aux[i]] = 1;
		/*
		 * The "2 + 1 backup" reserve: 2, not 1, so every loop that runs a station's filter
		 * still runs it (its ambiguities must be fixed before it is ever needed), while the
		 * interpolator is told to hold it back.
		 */
		if (cl.backup >= 0)
			used[cl.backup] = 2;
		if (master_by_signals) {
			printf("# master by signal count:");
			for (i = 0; i < nbase; i++)
				printf(" B%d=%d%s", i, nsig[i],
				    i == cl.master ? "*" : "");
			printf("   (* chosen; score = -km + 0.5*nsignal)\n");
		}
		printf("# cluster: master B%d at %.1f km, %d auxiliaries, %s hull, "
		    "extrapolation %.1f km, leverage %.2f\n",
		    master, cl.master_dist_km, cl.naux,
		    cl.inside_hull ? "inside" : "OUTSIDE", cl.extrapolation_km, cl.leverage);
		if (cl.backup >= 0)
			printf("# backup: B%d at %.1f km from the VBS, %.1f km from the master, "
			    "held in reserve\n", cl.backup, cl.backup_dist_km, cl.backup_baseline_km);
		/*
		 * The master-to-auxiliary baseline lengths, which are what the network filter's
		 * ambiguity resolution actually faces -- and which nothing reported until now.
		 * The distance above is master-to-VBS, a different quantity serving a different
		 * purpose: it sets how far the corrections are carried, while these set how hard
		 * they are to produce in the first place. Two clusters can share a master-to-VBS
		 * distance and differ by 20 km in their shortest network baseline, so a sweep
		 * over the former alone cannot attribute what it measures.
		 */
		{
			double bl, lo = 1e30, hi = 0.0, sum = 0.0;
			int a, c;

			for (i = 0; i < cl.naux; i++) {
				a = cl.aux[i];
				bl = 0.0;
				for (c = 0; c < 3; c++) {
					double dx = bases[a].xyz[c] - bases[master].xyz[c];

					bl += dx * dx;
				}
				bl = sqrt(bl) / 1000.0;
				if (bl < lo)
					lo = bl;
				if (bl > hi)
					hi = bl;
				sum += bl;
			}
			if (cl.naux > 0)
				printf("# network baselines (master to auxiliary): "
				    "%.1f km shortest, %.1f km longest, %.1f km mean\n",
				    lo, hi, sum / cl.naux);
		}
		report_signal_pairing(master, &cl, prefer_l2c);
	}

	vbs_interp_default_opt(&iopt);
	if (dd_datum)
		iopt.datum = VBS_INTERP_DATUM_COMMON;
	if (iono_scale > 0.0)
		iopt.iono_scale = iono_scale;
	iopt.iono_gate_m = iono_gate;
	if (anchor_master >= 0)
		iopt.anchor_master = anchor_master;
	if (method >= 0)
		iopt.method = method;
	if (corr_len > 0.0)
		iopt.corr_len_km = corr_len;
	if (sig_field > 0.0)
		iopt.sig_field_m = sig_field;
	if (nondisp >= 0)
		iopt.apply_nondispersive = nondisp;
	if (corr_max_sigma > 0.0)
		iopt.max_sigma = corr_max_sigma;
	vbs_netekf_default_opt(&nopt);
	nopt.min_snr_dbhz = snr_mask;
	nopt.prefer_l2c = prefer_l2c;
	if (est_gradients >= 0)
		nopt.est_gradients = est_gradients;
	if (grad_sigma > 0.0)
		nopt.var_grad0 = grad_sigma * grad_sigma;
	if (max_mapg > 0.0)
		nopt.max_mapg = max_mapg;
	if (q_ion > 0.0)
		nopt.q_ion = q_ion;
	if (var_ion0 > 0.0)
		nopt.var_ion0 = var_ion0 * var_ion0;
	nopt.drop_band5 = no_band5;
	nopt.drop_band2 = no_band2;
	if (sig_code > 0.0)
		nopt.sig_code = sig_code;
	nopt.sig_code_gps = sig_code_gps;
	nopt.q_ion_gps = q_ion_gps;
	if (no_bds)
		nopt.sys_mask &= ~SYS_CMP;
	if (no_gal)
		nopt.sys_mask &= ~SYS_GAL;
	if (nopt_sys_mask & SYS_GLO)
		nopt.sys_mask |= SYS_GLO;
	else
		nopt.sys_mask &= ~SYS_GLO;
	vbs_ambres_default_opt(&aopt);
	/*
	 * D48: at 1 Hz the formal ambiguity covariance is optimistic by roughly (residual
	 * correlation time / interval); inflating it at the search alone kept availability and
	 * cut narrow-lane misclosure events 179 -> 29 on E1-0. Provisional: the correlation
	 * time is asserted, not measured. Exactly 1 (no effect) at or above that interval.
	 */
	if (ar_corr_time_s > 0.0 && interval_s > 0.0 && interval_s < ar_corr_time_s)
		aopt.q_inflate = ar_corr_time_s / interval_s;
	if (min_ratio > 0.0)
		aopt.min_ratio = min_ratio;
	if (ar_max_step >= 0)
		aopt.max_step = ar_max_step;
	if (cascade_by_lambda >= 0)
		aopt.cascade_by_lambda = cascade_by_lambda;
	if (use_band6 >= 0)
		aopt.use_band6 = use_band6;
	aopt.skip_wl12_sys = no_wl12_sys;
	aopt.wl12_derive_sys = wl12_derive_sys;
	aopt.one_ref_sys = one_ref_sys;
	aopt.if_narrow_lane = if_nl;
	if (if_min_epochs > 0)
		aopt.if_min_epochs = if_min_epochs;
	if (if_max_sigma > 0.0)
		aopt.if_max_sigma = if_max_sigma;
	if (if_sys >= 0.0)
		aopt.if_sys_m = if_sys;
	if (if_max_frac > 0.0)
		aopt.if_max_frac = if_max_frac;
	vbs_basemeta_default_opt(&bopt);
	for (i = 0; i < nbase; i++) {
		if ((bases[i].pp = vbs_preproc_new(interval_s, 3.0)) == NULL) {
			fprintf(stderr, "cannot create preprocessor for base %d\n", i);
			return 1;
		}
		memset(bases[i].arc, 0, sizeof bases[i].arc);
	}
	for (i = 0; i < nbase; i++) {
		if (i == master || !used[i])
			continue;
		bases[i].ekf = vbs_netekf_new(bases[master].xyz, bases[i].xyz, &nopt);
		if (bases[i].ekf == NULL) {
			fprintf(stderr, "cannot create filter for base %d\n", i);
			return 1;
		}
	}

	memset(&loop_stats, 0, sizeof loop_stats);
	for (i = 0; i < MAXBASE; i++) {
		int a3, c3;

		for (a3 = 0; a3 < MAXSAT; a3++)
			for (c3 = 0; c3 < NCOMB; c3++)
				last_fix[i][a3][c3] = NAN;
	}
	if (lock == NULL) {
		fprintf(stderr, "cannot create lock tracker\n");
		return 1;
	}
	if (loop_closure) {
		int a2, b2;

		/*
		 * Every auxiliary pair, not a chain. With four stations each baseline then sits in
		 * two or more triangles, which is what makes the culprit identifiable; a chain
		 * leaves most baselines in exactly one triangle, where a failure can only say
		 * "one of these three".
		 *
		 * Identical configuration to the master baselines -- same vbs_netekf_new(NULL)
		 * defaults, same vbs_ambres options. If these filters converged less well that
		 * would be a setup difference to remove, not a property to accept, so there is
		 * deliberately nothing to configure differently.
		 */
		for (a2 = 0; a2 < nbase; a2++) {
			if (a2 == master || !used[a2])
				continue;
			for (b2 = a2 + 1; b2 < nbase; b2++) {
				if (b2 == master || !used[b2] || naa >= MAXPAIR)
					continue;
				aa_ekf[naa] = vbs_netekf_new(bases[a2].xyz, bases[b2].xyz, &nopt);
				if (aa_ekf[naa] == NULL) {
					fprintf(stderr, "cannot create aux-aux filter\n");
					return 1;
				}
				aa_a[naa] = a2;
				aa_b[naa] = b2;
				naa++;
			}
		}
		printf("# loop closure: %d auxiliary-auxiliary baseline(s), all pairs\n", naa);
	}

	/* The VBS point, plus D30 item 2's rover-to-VBS offset in local ENU. */
	{
		double pos[3], denu[3], dr[3];

		ecef2pos(vbs_xyz, pos);
		denu[0] = offset[0];
		denu[1] = offset[1];
		denu[2] = offset[2];
		enu2ecef(pos, denu, dr);
		for (i = 0; i < 3; i++)
			vbs_eff[i] = vbs_xyz[i] + dr[i];
	}

	identity = !strcmp(mode, "identity");
	if (!strcmp(mode, "loo")) {
		if (!have_target) {
			fprintf(stderr, "--mode loo needs --target\n");
			return 1;
		}
		dd = vbs_ddres_new(ELEV_MASK_DEG);
		/*
		 * The reference accumulator is fed from the *same* epoch, the same satellites
		 * and the same synthesis call, with only the correction values zeroed. Running
		 * it as a second process could not guarantee that, and D32 records what happens
		 * when the two sample sets are allowed to differ.
		 */
		if (paired)
			ddref = vbs_ddres_new(ELEV_MASK_DEG);
	}
	if (outpath != NULL && (out = fopen(outpath, "w")) == NULL) {
		fprintf(stderr, "cannot write %s\n", outpath);
		return 1;
	}

	if (c2w_raw_out && mw_gate_cyc > 0.0) {
		fprintf(stderr, "--c2w-raw-out cannot be combined with --mw-gate\n");
		return 2;
	}
	if (c2w_mode && c2w_apply_cluster(c2w_mode, c2w_win, c2w_partial, master, c2w_raw_out) < 0)
		return 1;

	/* --- the epoch loop ------------------------------------------------------------- */
	for (;;) {
		const obsd_t *mep, *aep, *tep, *mep_full;
		int nm_full;
		double mxyz[3], vxyz[3], vpos[3];
		int nm, na, nt = 0, nvalid = 0, ncl = 0, nsynth;
		gtime_t t;

		if ((nm = slice(&bases[master].obs, &bases[master].cursor, &mep)) <= 0)
			break;
		t = mep[0].time;
		if (nm > MAXOBS)
			nm = MAXOBS;
		nepoch++;
		/*
		 * D53 MW-scatter gate: the preprocessor keeps seeing every record (so a gated
		 * satellite's statistics keep updating and it can come back); the filters, the
		 * geometry and the output see the gated set. Applied here, before satposs, so
		 * every record-ordered array below refers to the same records. The gate uses
		 * the statistics up to the previous epoch.
		 */
		mep_full = mep;
		nm_full = nm;
		if (mw_gate_cyc > 0.0) {
			nm = mw_gate_records(bases[master].pp, mep, nm, mw_gate_cyc, gated_m, &mw_gated);
			mep = gated_m;
			if (nm <= 0)
				continue;
		}

		satposs(t, mep, nm, &nav, EPHOPT_BRDC, rs, dts, var, svh);
		if (nsp3 > 0)
			vbs_peph_substitute(&nav, mep, nm, rs, bases[master].xyz, vbs_eff, &sp3st);
		vbs_satpos_by_sat(mep, nm, rs, rs_sat);

		if (identity) {
			/*
			 * Exactly what test_zero_baseline verifies at the observation level --
			 * synthesis at the master with ncorr = 0 reproduces the master -- but taken
			 * all the way through the RINEX writer to a rover. The master's LLI is kept
			 * verbatim: the point is to isolate the output path, not to judge the
			 * engine's slip flagging against the receiver's.
			 */
			vbs_basemeta_tide_xyz(bases[master].xyz, t, &bopt, mxyz);
			nsynth = vbs_synth_epoch(mep, nm, rs, mxyz, mxyz, NULL, 0, &nav, synth);
			if (nsynth <= 0)
				continue;
			nepoch_out++;
			nsat_total += nsynth;
			if (out != NULL) {
				if (!rnx_header_written) {
					/*
					 * nopt_sys_mask: the constellation set belongs to one
					 * variable, not to each call site. Spelled out inline
					 * here and at the emit header below, the RINEX writer
					 * omitted GLONASS from both the header's SYS/#/OBS TYPES
					 * and the records, so satellites that were synthesized
					 * *with valid corrections* never reached the file --
					 * the fifth place in this tool where a duplicated
					 * constellation list silently discarded work already done
					 * (D90).
					 */
					rnx_init(&ropt, bases[master].obs.data, bases[master].obs.n,
					    bases[master].xyz, t, nopt_sys_mask);
					if (!outrnxobsh(out, &ropt, &nav)) {
						fprintf(stderr, "cannot write RINEX header\n");
						return 1;
					}
					rnx_header_written = 1;
				}
				outrnxobsb(out, &ropt, synth, nsynth, 0);
			}
			continue;
		}

		/* The master's station-level arcs, read once and shared by every filter. */
		nmslip = vbs_preproc_epoch(bases[master].pp, mep_full, nm_full, &nav, mslip, MAXOBS);
		vbs_preproc_arc_ids(bases[master].pp, bases[master].arc);
		nmslip_sat = 0;
		for (k = 0; k < nmslip && k < MAXOBS; k++)
			mslip_sat[nmslip_sat++] = mslip[k].sat;
		nrefixed = 0;

		/*
		 * Tide-corrected instantaneous positions (D23). The master and the VBS get the
		 * same treatment from the same model at the same epoch; applying it to one and
		 * not the other would turn a term that cancels into a real error.
		 */
		vbs_basemeta_tide_xyz(bases[master].xyz, t, &bopt, mxyz);
		vbs_basemeta_tide_xyz(vbs_eff, t, &bopt, vxyz);
		ecef2pos(vxyz, vpos);

		for (i = 0; i < nbase; i++) {
			if (i == master || !used[i])
				continue;
			bases[i].ep = NULL;
			bases[i].nep = 0;
			if ((na = slice_at(&bases[i], t, &aep)) <= 0) {
				bases[i].ncorr = 0;
				continue;
			}
			vbs_preproc_epoch(bases[i].pp, aep, na, &nav, NULL, 0);
			vbs_preproc_arc_ids(bases[i].pp, bases[i].arc);
			if (mw_gate_cyc > 0.0) {
				na = mw_gate_records(bases[i].pp, aep, na, mw_gate_cyc, gated_a[i],
				    &mw_gated);
				aep = gated_a[i];
			}
			bases[i].ep = aep;
			bases[i].nep = na;
			vbs_netekf_epoch(bases[i].ekf, t, mep, nm, rs_sat, dts, aep, na, &nav,
			    bases[master].arc, bases[i].arc);
			if (!no_ar)
				/*
				 * nopt_sys_mask, not a second copy of the constellation list.
				 * This call used to spell the mask out inline, which meant the
				 * resolver silently ran on a different set of constellations
				 * from the filter: GLONASS entered the filter's sigset, was
				 * never offered to the cascade, and its float corrections then
				 * failed the interpolation's sigma gate -- so GLONASS was
				 * synthesized, dropped for want of a valid correction, and no
				 * GLONASS reached the rover at all (D90).
				 */
				vbs_ambres_resolve(bases[i].ekf, &aopt, nopt_sys_mask, &ar,
				    bases[i].fix, 256, &bases[i].nfix);
			/*
			 * 0.05 m: below the 114-341 mm that newly-risen, low-elevation
			 * satellites sit at before their ambiguities are ever attempted, and
			 * well above the 16-20 mm a converged one reaches. D30 item 5 requires
			 * that those satellites reach neither interp nor the output.
			 */
			if (amb_dump) {
				int q4;

				for (q4 = 0; q4 < bases[i].nfix && amb_dump_n < amb_dump_max; q4++) {
					const struct vbs_amb_fix *f = &bases[i].fix[q4];
					char sid[8], rid[8];

					satno2id(f->sat, sid);
					satno2id(f->refsat, rid);
					printf("AMB %s B%d %s %s %d %d %.4f %.4f %.4f\n",
					    time_str(t, 0), i, sid, rid, f->comb_id, f->step,
					    f->value, f->integer, f->sigma);
					amb_dump_n++;
				}
			}
			if (ion_from_fixes) {
				int nc, nrep, q5, keep = 0;

				/*
				 * Widest gate the function accepts, so every tracked satellite
				 * appears and the fixed-integer readout decides, not the state's
				 * variance. Then keep only what the integers actually produced.
				 */
				nc = vbs_netekf_corrections(bases[i].ekf, bases[i].corr, MAXOBS,
				    sqrt(nopt.var_ion0) * 0.999);
				vbs_netekf_correction_sigmas(bases[i].ekf, bases[i].corr, nc,
				    bases[i].sd_I, bases[i].sd_T);
				/*
				 * Union, not replacement. A satellite whose integers are fixed gets
				 * the exact geometry-free ionosphere; one whose integers are not yet
				 * fixed but whose state has genuinely converged keeps the old path.
				 * Requiring the integers alone cut coverage from 15 satellites an
				 * epoch to 3, because the carrier step this cascade still uses is the
				 * one the ionosphere-free method exists to replace.
				 */
				for (q5 = 0; q5 < nc; q5++)
					bases[i].corr[q5].valid =
					    bases[i].sd_I[q5] <= 0.05 ? 2 : 0;
				nrep = vbs_ambres_ion_from_fixes(bases[i].ekf, bases[i].fix,
				    bases[i].nfix, bases[i].corr, nc, bases[i].sd_I);
				for (q5 = 0; q5 < nc; q5++) {
					if (!bases[i].corr[q5].valid ||
					    bases[i].sd_T[q5] > 0.05)
						continue;
					bases[i].corr[q5].valid = 1;
					bases[i].corr[keep] = bases[i].corr[q5];
					bases[i].sd_I[keep] = bases[i].sd_I[q5];
					bases[i].sd_T[keep] = bases[i].sd_T[q5];
					keep++;
				}
				bases[i].ncorr = keep;
				nrep_total += nrep;
				ncorr_total += nc;
			} else {
				bases[i].ncorr = vbs_netekf_corrections(bases[i].ekf, bases[i].corr,
				    MAXOBS, 0.05);
				vbs_netekf_correction_sigmas(bases[i].ekf, bases[i].corr,
				    bases[i].ncorr, bases[i].sd_I, bases[i].sd_T);
			}
			if (corr_dump && corr_dump_n < corr_dump_max) {
				/*
				 * The same call with the gate open, so the line says both what the
				 * filter believes and whether 0.05 m let it through: a satellite
				 * missing from the table is one the filter has no state for at all,
				 * which is a different failure from one it gates away.
				 */
				static struct vbs_corr all[MAXOBS];
				static double asd_I[MAXOBS], asd_T[MAXOBS];
				/*
				 * Just under the prior: vbs_netekf_corrections() refuses a gate at
				 * or above sqrt(var_ion0), on the ground that such a threshold
				 * cannot tell an informed state from an untouched one. For the dump
				 * that refusal would return nothing at all, so ask for the widest
				 * gate the function accepts and let the columns speak.
				 */
				int na2 = vbs_netekf_corrections(bases[i].ekf, all, MAXOBS,
				    sqrt(nopt.var_ion0) * 0.999);
				int q3;

				vbs_netekf_correction_sigmas(bases[i].ekf, all, na2, asd_I, asd_T);
				for (q3 = 0; q3 < na2 && corr_dump_n < corr_dump_max; q3++) {
					char sid[8];

					satno2id(all[q3].sat, sid);
					printf("CORR %s %s B%d %.4f %.4f %.4f %.4f %d %d\n",
					    time_str(t, 0), sid, i, all[q3].dI, asd_I[q3],
					    all[q3].dT, asd_T[q3], all[q3].valid,
					    asd_I[q3] <= 0.05 && asd_T[q3] <= 0.05);
					corr_dump_n++;
				}
			}

			memcpy(iaux[ncl].xyz, bases[i].xyz, sizeof iaux[ncl].xyz);
			iaux[ncl].backup = used[i] == 2;	/* iaux is static: set it every time */
			iaux[ncl].corr = bases[i].corr;
			iaux[ncl].sd_I = bases[i].sd_I;
			iaux[ncl].sd_T = bases[i].sd_T;
			iaux[ncl].ncorr = bases[i].ncorr;
			ncl++;
		}
		/*
		 * §5.3 case 2: which satellites had a network ambiguity re-fixed to a different
		 * integer, or lost, on any master baseline this epoch. Compared per
		 * (baseline, satellite, combination) against the previous epoch.
		 */
		for (i = 0; i < nbase; i++) {
			int q;

			if (i == master || bases[i].ekf == NULL)
				continue;
			for (q = 0; q < bases[i].nfix; q++) {
				const struct vbs_amb_fix *f = &bases[i].fix[q];
				double *slot;

				if (f->sat <= 0 || f->sat > MAXSAT || f->comb_id < 0 ||
				    f->comb_id >= NCOMB)
					continue;
				slot = &last_fix[i][f->sat - 1][f->comb_id];
				if (!isnan(*slot) && *slot != f->integer) {
					int r, dup = 0;

					for (r = 0; r < nrefixed; r++)
						if (refixed[r] == f->sat)
							dup = 1;
					if (!dup && nrefixed < MAXSAT)
						refixed[nrefixed++] = f->sat;
				}
				*slot = f->integer;
			}
		}
		/* An arc change at the master resets every baseline: the ambiguity is lost. */
		for (k = 0; k < MAXSAT; k++) {
			if (bases[master].arc[k] != last_arc_m[k]) {
				int r, dup = 0;

				for (r = 0; r < nrefixed; r++)
					if (refixed[r] == k + 1)
						dup = 1;
				if (!dup && nrefixed < MAXSAT && last_arc_m[k] != 0)
					refixed[nrefixed++] = k + 1;
				/* Forget the integers of an arc that no longer exists. */
				for (i = 0; i < nbase; i++) {
					int c;

					for (c = 0; c < NCOMB; c++)
						last_fix[i][k][c] = NAN;
				}
			}
			last_arc_m[k] = bases[master].arc[k];
		}

		/*
		 * Loop closure (D33 item 2). Each auxiliary-auxiliary baseline is filtered
		 * independently, then each triangle M-A-B is checked: N_AB must equal
		 * N_MB - N_MA. A misclosure proves a wrong fix somewhere in the triangle, and
		 * since one triangle cannot say which, all three are unfixed.
		 */
		nfail_sats = 0;
		if (loop_closure)
			for (k = 0; k < naa; k++) {
				struct vbs_amb_fix bad[64];
				int nbad, q;

				if (bases[aa_a[k]].ep == NULL || bases[aa_b[k]].ep == NULL)
					continue;
				/*
				 * One geometry for the whole network: rs_sat is indexed by satellite, so
				 * the aux-aux filter reads exactly the positions the master baselines
				 * used (D42, D49).
				 */
				vbs_netekf_epoch(aa_ekf[k], t, bases[aa_a[k]].ep,
				    bases[aa_a[k]].nep, rs_sat, dts,
				    bases[aa_b[k]].ep, bases[aa_b[k]].nep, &nav,
				    bases[aa_a[k]].arc, bases[aa_b[k]].arc);
				vbs_ambres_resolve(aa_ekf[k], &aopt, nopt_sys_mask, &ar,
				    aa_fix[k], 256, &aa_n[k]);

				nbad = vbs_ambres_loop_check(
				    bases[aa_a[k]].fix, bases[aa_a[k]].nfix,
				    bases[aa_b[k]].fix, bases[aa_b[k]].nfix,
				    aa_fix[k], aa_n[k], 0.25, &loop_stats, bad, 64);
				/*
				 * Record which baselines this failing triangle implicates, rather
				 * than unfixing now: the culprit only becomes visible once every
				 * triangle has voted.
				 */
				for (q = 0; q < nbad; q++) {
					int sat = bad[q].sat;

					if (sat <= 0 || sat > MAXSAT)
						continue;
					if (bad[q].comb_id >= 0 && bad[q].comb_id < NCOMB &&
					    bad[q].step >= 0 && bad[q].step < VBS_AMB_NSTEP) {
						long *le = &mis_last_epoch[k][sat - 1][bad[q].comb_id];
						double *li = &mis_last_int[k][sat - 1][bad[q].comb_id];
						double mi = bad[q].integer;

						if (*le != nepoch - 1 || *li != mi) {
							int sy = satsys(sat, NULL);

							mis_onsets[bad[q].step]++;
							mis_onsets_sat[bad[q].step][sat - 1]++;
							mis_onsets_sys[bad[q].step][sy == SYS_GPS ? 0 :
							    sy == SYS_GAL ? 1 : sy == SYS_CMP ? 2 : 3]++;
						}
						*le = nepoch;
						*li = mi;
					}
					/*
					 * D47's bounded check: for the first extra-wide-lane misclosures,
					 * show the three floats beside the three integers. Each float
					 * near its own integer means the *floats* disagree across the
					 * baselines (a state defect); a float far from the integer it
					 * was rounded to means the search accepted what it should not.
					 */
					if (ewl_dump && bad[q].step == VBS_AMB_STEP_EWL && ewl_dumped < 40) {
						double va, ia, sa, vb, ib, sb, vc, ic, sc;
						char tbuf[40], ida[8];

						if (fix_lookup(bases[aa_a[k]].fix, bases[aa_a[k]].nfix, sat,
						    bad[q].refsat, bad[q].comb_id, &va, &ia, &sa) &&
						    fix_lookup(bases[aa_b[k]].fix, bases[aa_b[k]].nfix, sat,
						    bad[q].refsat, bad[q].comb_id, &vb, &ib, &sb) &&
						    fix_lookup(aa_fix[k], aa_n[k], sat,
						    bad[q].refsat, bad[q].comb_id, &vc, &ic, &sc)) {
							time2str(t, tbuf, 0);
							satno2id(sat, ida);
							printf("EWL %s %s ref%3d tri%d  MA %+9.3f->%+6.0f (%.2f)"
							    "  MB %+9.3f->%+6.0f (%.2f)  AB %+9.3f->%+6.0f (%.2f)"
							    "  mis float %+7.3f int %+5.0f\n",
							    tbuf, ida, bad[q].refsat, k, va, ia, sa, vb, ib, sb,
							    vc, ic, sc, vc - (vb - va), ic - (ib - ia));
							ewl_dumped++;
						}
					}
					/*
					 * Do the three baselines agree about which *arc* this
					 * satellite is on? Each runs its own slip detection over
					 * its own epoch set, so one can have reset after a slip
					 * while another has not. Differencing across that boundary
					 * gives an integer equal to the slip size -- which looks
					 * exactly like a wrong fix and is not one.
					 */
					/*
					 * With arcs shared per station this is an assertion: the
					 * three filters read the same station ids, so they cannot
					 * disagree unless the sharing is broken. Any count here is
					 * a bug in the plumbing, not a property of the data.
					 */
					{
						long ma_m, ma_a, mb_m, mb_b, ab_a, ab_b;
						int ok = vbs_netekf_station_arcs(
						    bases[aa_a[k]].ekf, sat, &ma_m, &ma_a)
						    && vbs_netekf_station_arcs(
						    bases[aa_b[k]].ekf, sat, &mb_m, &mb_b)
						    && vbs_netekf_station_arcs(
						    aa_ekf[k], sat, &ab_a, &ab_b);

						if (!ok || ma_m != mb_m || ma_a != ab_a ||
						    mb_b != ab_b)
							mis_arc_mismatch++;
						else
							mis_same_arc++;
					}
					if (!fail_tri[sat - 1][MAXPAIR + aa_a[k]] &&
					    !fail_tri[sat - 1][MAXPAIR + aa_b[k]] &&
					    !fail_tri[sat - 1][k] && nfail_sats < MAXSAT)
						fail_sats[nfail_sats++] = sat;
					fail_tri[sat - 1][MAXPAIR + aa_a[k]]++;
					fail_tri[sat - 1][MAXPAIR + aa_b[k]]++;
					fail_tri[sat - 1][k]++;
				}
			}

		/*
		 * Blame. For each satellite that failed somewhere, the baseline implicated by
		 * *every* failing triangle is the culprit: an innocent baseline is only dragged in
		 * by the triangles it shares with the guilty one, so it cannot appear in all of
		 * them unless it is the guilty one. If no single baseline is common to all -- which
		 * includes the case of a lone failing triangle, where nothing can be isolated --
		 * every implicated baseline is unfixed, as before.
		 */
		for (q2 = 0; q2 < nfail_sats; q2++) {
			int sat = fail_sats[q2], nmax = 0, culprit = -1, nculprit = 0, b;

			for (b = 0; b < MAXPAIR + MAXBASE; b++)
				if (fail_tri[sat - 1][b] > nmax)
					nmax = fail_tri[sat - 1][b];
			for (b = 0; b < MAXPAIR + MAXBASE; b++)
				if (fail_tri[sat - 1][b] == nmax) {
					culprit = b;
					nculprit++;
				}
			if (nmax >= 2 && nculprit == 1) {
				/*
				 * Which kind of baseline gets the blame? If it falls mostly on the
				 * auxiliary-auxiliary filters, the fault is in the check's own
				 * machinery rather than in the corrections it is auditing.
				 */
				if (culprit >= MAXPAIR)
					blame_master++;
				else
					blame_auxaux++;
				if (!report_only) {
					if (culprit >= MAXPAIR)
						vbs_netekf_unfix_sat(
						    bases[culprit - MAXPAIR].ekf, sat);
					else
						vbs_netekf_unfix_sat(aa_ekf[culprit], sat);
					loop_unfixed++;
				}
				loop_isolated++;
			} else {
				for (b = 0; b < MAXPAIR + MAXBASE; b++) {
					if (!fail_tri[sat - 1][b])
						continue;
					if (b >= MAXPAIR)
						blame_master++;
					else
						blame_auxaux++;
					if (!report_only) {
						if (b >= MAXPAIR)
							vbs_netekf_unfix_sat(
							    bases[b - MAXPAIR].ekf, sat);
						else
							vbs_netekf_unfix_sat(aa_ekf[b], sat);
						loop_unfixed++;
					}
				}
			}
			memset(fail_tri[sat - 1], 0, sizeof fail_tri[sat - 1]);
		}

		/*
		 * Coverage: what fraction of the fixes made on master baselines this epoch belong
		 * to a triangle that was actually checked. A check that rarely runs gives the same
		 * false assurance as a guard that never fires, so the denominator is reported
		 * alongside the misclosure rate rather than left implicit.
		 */
		if (loop_closure)
			for (i = 0; i < nbase; i++)
				if (i != master && bases[i].ekf != NULL)
					fixes_total += bases[i].nfix;

		if (ncl < 3) {
			nskip_aux++;
			continue;
		}

		/*
		 * Satellite az/el at the VBS, keyed by satellite number, for §4.5's computed
		 * troposphere model term. Built from the same satellite positions the filters
		 * used, so the geometry cannot disagree between the two.
		 */
		memset(azel, 0, sizeof azel);
		for (k = 0; k < nm; k++) {
			double ev[3], ae[2];

			if (mep[k].sat <= 0 || mep[k].sat > MAXSAT)
				continue;
			if (geodist(rs + k * 6, vxyz, ev) <= 0.0)
				continue;
			satazel(vpos, ev, ae);
			azel[(mep[k].sat - 1) * 2] = ae[0];
			azel[(mep[k].sat - 1) * 2 + 1] = ae[1];
		}

		/*
		 * --no-tropo withholds the troposphere model term as well, by giving interp no
		 * elevations. That leaves the dispersive correction alone, which makes an exact
		 * identity available: the ionosphere-free combination of the applied correction
		 * must be *zero*, because ΔI is applied as μ_i = (f_ref/f_i)² and the iono-free
		 * combination is built to annihilate exactly that. Any non-zero result means the
		 * μ_i used in synthesis disagree with the band frequencies used in the
		 * combination -- and would have leaked into every alpha figure in D34.
		 */
		nvalid = vbs_interp_epoch(t, mxyz, vxyz, iaux, ncl, &iopt,
		    no_tropo ? NULL : azel, vcorr, vqual, MAXOBS);
		if (nvalid <= 0) {
			nskip_corr++;
			continue;
		}
		/*
		 * D52: the satellites that stopped a rover fixing for whole ten-minute windows
		 * on E1-0 were the low ones (E36 at 19°, E32 at 21°, E09 at 17°); excluding one
		 * of them at the rover made the window fix. A correction's formal sigma does not
		 * see the bias a low satellite carries, so the emit itself masks by elevation at
		 * the VBS, above the rover's own 15°. Policy knob, measured, not a default yet.
		 */
		{
			/*
			 * Count what is actually offered to the rover, not what was interpolated.
			 *
			 * The mask clears `valid` without touching `nvalid`, and the satellites-per-
			 * epoch statistic below accumulated `nvalid`, so it reported the *pre-mask*
			 * count and was identical at 20 and 25 degrees. That made it useless for the
			 * one question it was being used to answer -- whether a mask that helps is
			 * trading satellite count for correction quality, and therefore whether
			 * adding a constellation would restore the count it costs.
			 */
			int nmasked = 0, nserved;

			if (out_elmask_deg > 0.0)
				for (k = 0; k < nvalid; k++) {
					int sat = vcorr[k].sat;

					if (sat >= 1 && sat <= MAXSAT && vcorr[k].valid &&
					    azel[(sat - 1) * 2 + 1] < out_elmask_deg * D2R) {
						vcorr[k].valid = 0;
						out_elmasked++;
						nmasked++;
					}
				}
			nserved = nvalid - nmasked;
			nsat_total += nserved;
			if (nserved < nsat_min)
				nsat_min = nserved;
			if (nserved > nsat_max)
				nsat_max = nserved;
			nsat_premask_total += nvalid;
		}
		nepoch_out++;

		/*
		 * The uncorrected reference, for a comparison that means something.
		 *
		 * The corrections are zeroed *here* -- after interpolation, not before it -- so
		 * that the epochs kept, the satellites kept and the validity decisions are
		 * byte-identical to the corrected run. Comparing a corrected run against a
		 * separately-generated uncorrected one does not measure the corrections: it also
		 * measures the difference between the two epoch sets, and since the epochs the
		 * corrected run drops are the hard ones, that difference flatters the
		 * corrections. The first comparison made here had exactly that defect (1679
		 * epochs against 2637).
		 */
		if (zero_corr)
			for (k = 0; k < nvalid; k++)
				vcorr[k].dI = vcorr[k].dT = 0.0;

		/*
		 * The target's epoch is taken once and shared. slice_at() advances the station's
		 * cursor, so calling it in both the reference and the corrected branch consumed
		 * the epoch in the first and left the second with nothing -- which surfaced as
		 * the paired comparison refusing itself on mismatched sample counts. The refusal
		 * was the design working: mispaired series would have produced a plausible
		 * decomposition of nonsense.
		 */
		nt = have_target ? slice_at(&target, t, &tep) : 0;

		if (ddref != NULL && nt > 0) {
			static struct vbs_corr zcorr[MAXOBS];
			static obsd_t zsynth[MAXOBS];
			int nz;

			memcpy(zcorr, vcorr, nvalid * sizeof *vcorr);
			for (k = 0; k < nvalid; k++)
				zcorr[k].dI = zcorr[k].dT = 0.0;
			nz = vbs_synth_epoch(C2W_SYNTH_EP(mep), nm, rs, mxyz, vxyz, zcorr, nvalid, &nav,
			    zsynth);
			if (nz > 0)
				vbs_ddres_add_epoch(ddref, zsynth, nz, tep, nt, rs, vxyz, &nav);
		}

		/*
		 * The mirror identity: with the dispersive correction withheld, the
		 * *geometry-free* combination of what remains must be zero, because ΔT is applied
		 * identically on every band and L_a - L_b annihilates anything band-independent.
		 * Together with --no-tropo's iono-free identity this pins both components to the
		 * combinations they are supposed to live in.
		 */
		if (no_iono)
			for (k = 0; k < nvalid; k++)
				vcorr[k].dI = 0.0;
		else if (vbs_ionopost_active(&ipo)) {
			ipo.sys_mask = iono_sys_mask;
			vbs_ionopost_apply(vcorr, nvalid, azel, vpos, mep, nm, &nav, &ipo, &ips);
		} else if (iono_sys_mask)
			for (k = 0; k < nvalid; k++)
				if (!(satsys(vcorr[k].sat, NULL) & iono_sys_mask))
					vcorr[k].dI = 0.0;
		if (no_iono && ipo.drop_sys)
			for (k = 0; k < nvalid; k++)
				if (satsys(vcorr[k].sat, NULL) & ipo.drop_sys)
					vcorr[k].valid = 0;

		nsynth = vbs_synth_epoch(C2W_SYNTH_EP(mep), nm, rs, mxyz, vxyz, vcorr, nvalid, &nav,
		    synth);
		if (nsynth <= 0)
			continue;

		/* §5.3: stamp the output's lock-time signalling. In place, on `synth`. */
		vbs_lock_update(lock, t, mep, nm, synth, nsynth, vcorr, nvalid,
		    mslip_sat, nmslip_sat, refixed, nrefixed);

		if (dd != NULL && nt > 0)
			vbs_ddres_add_epoch(dd, synth, nsynth, tep, nt, rs, vxyz, &nav);
		if (out != NULL) {
			if (!rnx_header_written) {
				/*
				 * The header carries the tide-free catalogue coordinate, not the
				 * tide-displaced synthesis point: the same convention a 1005 gives a
				 * rover (D23), and what makes the rover's own tidal displacement
				 * cancel against the VBS's.
				 */
				rnx_init(&ropt, bases[master].obs.data, bases[master].obs.n,
				    vbs_eff, t, nopt_sys_mask);
				if (!outrnxobsh(out, &ropt, &nav)) {
					fprintf(stderr, "cannot write RINEX header\n");
					return 1;
				}
				rnx_header_written = 1;
			}
			outrnxobsb(out, &ropt, synth, nsynth, 0);
		}
	}

	/*
	 * D30 item 5: partial ambiguity resolution is the design assumption, so the satellite
	 * count that survives it is a headline number, not a footnote -- a VBS carrying too few
	 * satellites is useless to a rover however good each correction is.
	 */
	printf("# %ld master epochs: %ld produced a VBS, %ld skipped for <3 auxiliaries, "
	    "%ld for no converged correction\n", nepoch, nepoch_out, nskip_aux, nskip_corr);
	if (loop_closure)
		printf("# loop closure: %ld checked, %ld misclosed (%.2f%%), %ld unfixed, "
		    "%ld unusable (not fixed on all three or not re-referenceable), "
		    "worst %.2f cycles; coverage %.1f%% of fixes; %ld culprit(s) isolated\n",
		    loop_stats.checked, loop_stats.misclosed,
		    loop_stats.checked ? 100.0 * loop_stats.misclosed / loop_stats.checked : 0.0,
		    loop_unfixed, loop_stats.skipped_ref, loop_stats.worst_misclosure,
		    fixes_total ? 100.0 * loop_stats.checked / fixes_total : 0.0, loop_isolated);
	if (loop_closure) {
		static const char *SN[VBS_AMB_NSTEP] = {"extra-wide", "wide", "narrow"};
		static const char *BN[VBS_LOOP_NBIN] = {"closed", "1 cyc", "2 cyc",
		    ">2 cyc", "NON-INT"};
		int b;

		printf("# misclosure attribution: %ld on consistent station arcs, %ld arc "
		    "mismatches (must be 0 with per-station arcs -- any count is a plumbing bug)\n",
		    mis_same_arc, mis_arc_mismatch);
		printf("# blame: %ld master-auxiliary, %ld auxiliary-auxiliary%s\n",
		    blame_master, blame_auxaux, report_only ? "  (report-only: nothing unfixed)" : "");
		printf("# narrow-lane tail (>2.5 cyc): %ld where the wide lane ALSO misclosed, "
		    "%ld where it did not\n", loop_stats.nl_tail_with_wl,
		    loop_stats.nl_tail_without_wl);
		printf("# narrow-lane 1-cycle:        %ld with wide-lane misclosure, %ld without\n",
		    loop_stats.nl_one_with_wl, loop_stats.nl_one_without_wl);
		printf("# wide-lane tail (>2.5 cyc):   %ld where the EXTRA-wide lane also misclosed, "
		    "%ld where it did not\n", loop_stats.wl_tail_with_ewl,
		    loop_stats.wl_tail_without_ewl);
		printf("# narrow-lane tail (>2.5 cyc): %ld with extra-wide misclosure, %ld without\n",
		    loop_stats.nl_tail_with_ewl, loop_stats.nl_tail_without_ewl);
		printf("# misclosure events (onsets, not epochs): EWL %ld  WL %ld  NL %ld\n",
		    mis_onsets[VBS_AMB_STEP_EWL], mis_onsets[VBS_AMB_STEP_WL],
		    mis_onsets[VBS_AMB_STEP_NL]);
		printf("# misclosure events by constellation (GPS/GAL/BDS/other): EWL %ld/%ld/%ld/%ld  "
		    "WL %ld/%ld/%ld/%ld  NL %ld/%ld/%ld/%ld\n",
		    mis_onsets_sys[0][0], mis_onsets_sys[0][1], mis_onsets_sys[0][2], mis_onsets_sys[0][3],
		    mis_onsets_sys[1][0], mis_onsets_sys[1][1], mis_onsets_sys[1][2], mis_onsets_sys[1][3],
		    mis_onsets_sys[2][0], mis_onsets_sys[2][1], mis_onsets_sys[2][2], mis_onsets_sys[2][3]);
		{
			int st, rank;

			for (st = 0; st < VBS_AMB_NSTEP; st++) {
				printf("# misclosure events by PRN, %s, top 12:",
				    st == VBS_AMB_STEP_EWL ? "EWL" : st == VBS_AMB_STEP_WL ? "WL" : "NL");
				for (rank = 0; rank < 12; rank++) {
					int best = -1, s2;
					char id[8];

					for (s2 = 0; s2 < MAXSAT; s2++)
						if (mis_onsets_sat[st][s2] > 0 &&
						    (best < 0 || mis_onsets_sat[st][s2] > mis_onsets_sat[st][best]))
							best = s2;
					if (best < 0)
						break;
					satno2id(best + 1, id);
					printf(" %s:%ld", id, mis_onsets_sat[st][best]);
					mis_onsets_sat[st][best] = -mis_onsets_sat[st][best];	/* consumed */
				}
				printf("\n");
			}
		}
		printf("# output elevation mask %.0f deg: %ld satellite-epochs withheld; MW gate %.2f cyc: %ld records dropped\n",
		    out_elmask_deg, out_elmasked, mw_gate_cyc, mw_gated);
		printf("# misclosure histogram by cascade level\n");
		printf("#   %-11s", "step");
		for (b = 0; b < VBS_LOOP_NBIN; b++)
			printf(" %9s", BN[b]);
		printf("\n");
		for (k = 0; k < VBS_AMB_NSTEP; k++) {
			printf("    %-11s", SN[k]);
			for (b = 0; b < VBS_LOOP_NBIN; b++)
				printf(" %9ld", loop_stats.hist[k][b]);
			printf("\n");
		}
	}
	vbs_lock_get_stats(lock, &lock_st);
	printf("# §5.3 lock resets on the output: %ld master slip, %ld ambiguity re-fix/loss, "
	    "%ld correction jump, %ld re-entry, %ld arc start\n",
	    lock_st.resets[0], lock_st.resets[1], lock_st.resets[2], lock_st.resets[3],
	    lock_st.resets[4]);
	if (out_elmask_deg > 0.0 && nepoch_out > 0)
		printf("# output elevation mask %.0f deg: %.1f of %.1f satellites per epoch "
		    "survive it (%ld withheld in total)\n", out_elmask_deg,
		    (double)nsat_total / nepoch_out, (double)nsat_premask_total / nepoch_out,
		    out_elmasked);
	printf("# satellites per VBS epoch: %.1f mean, %ld min, %ld max\n",
	    nepoch_out ? (double)nsat_total / nepoch_out : 0.0,
	    nepoch_out ? nsat_min : 0, nsat_max);
	if (nsp3 > 0)
		vbs_peph_report(stdout, &sp3st, "master-VBS");
	/*
	 * What the cascade actually fixed, per step. Added because three separate experiments
	 * turned out to be measuring a step that never fires: the dump of a six-hour run holds
	 * 2.9 million accepted fixes and **not one** of them is a narrow lane, so the shipped
	 * configuration resolves wide lanes and stops. Nothing in the report said so.
	 */
	printf("# cascade fixes: EWL %ld/%ld attempted, WL %ld/%ld, NL %ld/%ld\n",
	    ar.fixed[VBS_AMB_STEP_EWL], ar.attempts[VBS_AMB_STEP_EWL],
	    ar.fixed[VBS_AMB_STEP_WL], ar.attempts[VBS_AMB_STEP_WL],
	    ar.fixed[VBS_AMB_STEP_NL], ar.attempts[VBS_AMB_STEP_NL]);
	/*
	 * D66's ionosphere-from-fixes, reported. The counter existed and was never printed, so
	 * whether the feature did anything was invisible from the outside -- and it turns out to
	 * do nothing at all on some clusters, which nobody could have seen. A capability with no
	 * observable effect is the defect class D90 named and counted five times in this engine;
	 * the cheapest guard against a sixth is to make every such path report its own count.
	 */
	if (vbs_ionopost_active(&ipo))
		printf("# ionosphere post: %ld model fits (mean residual %.3f m), %ld refused; %ld "
		    "corrections filled, %ld replaced, %ld kept; %ld withheld for no band 5; %ld shrunk "
		    "(mean weight %.3f)\n",
		    ips.fits, ips.fits ? ips.fit_rms_sum / ips.fits : 0.0, ips.fit_refused, ips.filled,
		    ips.replaced, ips.kept, ips.band5_withheld, ips.shrunk,
		    ips.shrunk ? ips.shrink_w_sum / ips.shrunk : 0.0);
	if (iono_gate > 0.0) {
		long chk = 0, drop = 0, wh = 0;

		vbs_interp_gate_stats(&chk, &drop, &wh);
		printf("# ionosphere gate at %.3f m: %ld satellite fits checked, %ld auxiliaries "
		    "dropped, %ld satellites withheld (%.2f %%)\n", iono_gate, chk, drop, wh,
		    chk ? 100.0 * wh / chk : 0.0);
	}
	if (ion_from_fixes) {
		long no_ref = 0, no_pair = 0, okc = 0;

		vbs_ambres_ion_diag(&no_ref, &no_pair, &okc);
		printf("# ionosphere from fixed integers: %ld of %ld corrections replaced (%.1f %%)\n",
		    nrep_total, ncorr_total,
		    ncorr_total ? 100.0 * nrep_total / ncorr_total : 0.0);
		printf("#   skipped: %ld for no usable reference satellite, %ld for no "
		    "wide-lane/narrow-lane pair against it\n", no_ref, no_pair);
	}

	if (dd != NULL) {
		struct vbs_dd_result res;
		int n = vbs_ddres_finish(dd, &res);

			printf("# leave-one-out residuals, corrections %s\n",
		    zero_corr ? "ZEROED (reference)" : "applied");
		for (i = 0; i < n && i < res.nstat; i++) {
			const struct vbs_dd_stat *s = &res.stat[i];
			const char *cn = s->comb == VBS_DD_GEOMFREE ? "geom-free" :
			    s->comb == VBS_DD_IONOFREE ? "iono-free" : "band";

			printf("  %-4s %-10s band %d: n=%-7d rms=%8.4f m  p95=%8.4f  max=%8.4f\n",
			    s->sys == SYS_GPS ? "GPS" : s->sys == SYS_GAL ? "GAL" : "BDS",
			    cn, s->band, s->nsample, s->rms, s->p95, s->max);
		}
		if (ddref != NULL) {
			struct vbs_dd_result refres;
			struct vbs_dd_pair pr[VBS_DD_MAX_STATS];
			int np;

			vbs_ddres_finish(ddref, &refres);
			/*
			 * Pairing by index assumes the two runs segment arcs identically. They
			 * need not: vbs_ddres breaks an arc on a data-driven DD jump
			 * (DD_SLIP_THRESHOLD_M), and the corrected and zeroed series carry
			 * different values, so a jump can cross the threshold in one and not the
			 * other. Different segmentation means different per-arc means, which puts
			 * step functions into c and would corrupt every alpha in the decomposition.
			 */
			for (i = 0; i < refres.nstat && i < res.nstat; i++)
				if (refres.stat[i].narc != res.stat[i].narc)
					printf("# ARC MISMATCH %s %d: ref %d arcs, test %d\n",
					    refres.stat[i].sys == SYS_GPS ? "GPS" :
					    refres.stat[i].sys == SYS_GAL ? "GAL" : "BDS",
					    refres.stat[i].comb, refres.stat[i].narc,
					    res.stat[i].narc);
			np = vbs_ddres_pair(ddref, dd, pr, VBS_DD_MAX_STATS);
			if (np < 0)
				printf("# paired comparison refused: sample counts differ\n");
			else {
				printf("# paired decomposition: what each correction removes "
				    "against what it adds\n");
				printf("#   %-4s %-10s %7s %9s %9s %9s %7s %9s %9s %9s\n",
				    "sys", "comb", "n", "ref_rms", "corr_rms", "applied",
				    "alpha", "best_rms", "c_r2_sin", "c_r2_map");
				for (i = 0; i < np; i++) {
					const char *cn = pr[i].comb == VBS_DD_GEOMFREE ?
					    "geom-free" : pr[i].comb == VBS_DD_IONOFREE ?
					    "iono-free" : "band";

					if (pr[i].comb == VBS_DD_BAND)
						continue;
					printf("    %-4s %-10s %7d %9.4f %11.3e %9.4f %7.2f %9.4f %9.3f %9.3f\n",
					    pr[i].sys == SYS_GPS ? "GPS" :
					    pr[i].sys == SYS_GAL ? "GAL" : "BDS", cn,
					    pr[i].nsample, pr[i].rms_ref, pr[i].rms_corr,
					    pr[i].rms_test, pr[i].alpha, pr[i].rms_best,
					    pr[i].r2_sin_corr, pr[i].r2_map_corr);
				}
			}
			{
				struct vbs_dd_shape sh[VBS_DD_MAX_STATS];
				int ns = vbs_ddres_shape(ddref, sh, VBS_DD_MAX_STATS);

				printf("# elevation shape of the UNCORRECTED residual "
				    "(D35 item 1): sin(el) is antenna-like, m(el) troposphere-like\n");
				printf("#   %-4s %-10s %7s %8s %8s %10s %10s %8s\n", "sys",
				    "comb", "n", "r2_sin", "r2_map", "coef_sin", "coef_map",
				    "reg_corr");
				for (i = 0; i < ns; i++) {
					const char *cn = sh[i].comb == VBS_DD_GEOMFREE ?
					    "geom-free" : sh[i].comb == VBS_DD_IONOFREE ?
					    "iono-free" : "band";

					if (sh[i].comb == VBS_DD_BAND)
						continue;
					printf("    %-4s %-10s %7d %8.3f %8.3f %10.4f %10.4f %8.2f\n",
					    sh[i].sys == SYS_GPS ? "GPS" :
					    sh[i].sys == SYS_GAL ? "GAL" : "BDS", cn,
					    sh[i].nsample, sh[i].r2_sin, sh[i].r2_map,
					    sh[i].coef_sin, sh[i].coef_map,
					    sh[i].regressor_corr);
				}
			}
			vbs_ddres_free(ddref);
		}
		vbs_ddres_free(dd);
	}
	if (out != NULL)
		fclose(out);
	for (i = 0; i < nbase; i++) {
		vbs_netekf_free(bases[i].ekf);
		vbs_preproc_free(bases[i].pp);
		free(bases[i].obs.data);
	}
	for (i = 0; i < naa; i++)
		vbs_netekf_free(aa_ekf[i]);
	vbs_lock_free(lock);
	if (have_target)
		free(target.obs.data);
	freenav(&nav, 0xFF);
	return 0;
}

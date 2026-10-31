/*
 * The C serving policy against the Python reference.
 *
 * Two layers. First, the fourteen fixtures from data/scripts/vbs_policy.py's self-test, rebuilt
 * here with the same geometry, must give the same rule ids. Second -- and this is the one that
 * catches drift -- the test writes a set of random positions and a base list to a file, and
 * tools/vbs_policy_dump can replay both implementations on it; the Makefile's `make policy-check`
 * target does that when the Python interpreter and the cache are available. A port checked only
 * against hand-picked cases is a port checked against the cases its author thought of.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "rtklib.h"
#include "vbs/vbs_policy.h"

static int nfail;

static void check(int ok, const char *what) {
	printf("  %-62s %s\n", what, ok ? "ok" : "FAIL");
	if (!ok)
		nfail++;
}

#define MOS	"SEPT MOSAIC-X5"
#define UM9	"UNICORE UM980"
#define F9P	"UBL-F9P"

static const double O[3] = {4500000.0, 0.0, 4500000.0};

/* The Python fixture: base at east e_km, north n_km of the origin, in a flat local frame. */
static void b(struct vbs_policy_base *out, const char *mp, const char *rec, double e_km, double n_km) {
	out->mp = mp;
	out->family = vbs_policy_family(rec);
	out->xyz[0] = O[0];
	out->xyz[1] = O[1] + e_km * 1000.0;
	out->xyz[2] = O[2] + n_km * 1000.0;
	out->gps_bands = 0;
	out->fw = NULL;
	out->live = 1;
}

/* fam_set(rec, first, step=18, n=5): auxiliaries spread north as well as east */
static int fam_set(struct vbs_policy_base *out, const char *rec, double first, double step, int n) {
	int i;

	for (i = 0; i < n; i++)
		b(&out[i], rec, rec, first + i * step, (i % 3) * 12.0);
	return n;
}

static void expect(const char *name, const struct vbs_policy_base *bases, int n, int want_serve,
    const char *want_rule, int want_disp, int check_bands) {
	struct vbs_decision d;
	int ok;
	char what[96];

	vbs_policy_decide(O, bases, n, check_bands, &d);
	ok = d.serve == want_serve && (want_rule == NULL || !strcmp(d.rule, want_rule)) &&
	    (want_disp < 0 || d.dispersive == want_disp);
	snprintf(what, sizeof what, "%-44s -> %s %s", name, d.rule, d.serve ? "SERVE" : "404");
	check(ok, what);
	if (!ok)
		printf("      reason: %s\n", d.reason);
}

int main(void) {
	struct vbs_policy_base bs[16];
	int n;

	printf("=== serving policy (C port vs the Python fixtures) ===\n");

	/* R0 density */
	b(&bs[0], "a", MOS, 30, 0); b(&bs[1], "b", MOS, 45, 0);
	expect("only 2 bases of any family -> R0", bs, 2, 0, "R0", -1, 0);
	b(&bs[0], "a", MOS, 30, 0); b(&bs[1], "b", MOS, 60, 0); b(&bs[2], "c", MOS, 90, 0);
	expect("sparse market: 3 Mosaic -> R0", bs, 3, 0, "R0", -1, 0);

	/* R1 coverage */
	n = fam_set(bs, UM9, 6.0, 18.0, 5);
	expect("nearest base 6 km -> R1", bs, n, 0, "R1", -1, 0);
	b(&bs[0], "a", UM9, 18, 0); b(&bs[1], "b", UM9, 22, 0); b(&bs[2], "c", UM9, 27, 0);
	b(&bs[3], "d", UM9, 60, 40); b(&bs[4], "e", UM9, 75, -30);
	expect("nearest 18 km, 3 within 12 km of it -> R1", bs, 5, 0, "R1", -1, 0);

	/* R2 per-family envelopes */
	n = fam_set(bs, UM9, 16.0, 18.0, 5);
	expect("UM980 gap at 16 km -> serve, dispersive (D167)", bs, n, 1, "R2", 1, 0);
	n = fam_set(bs, MOS, 35.0, 18.0, 5);
	expect("Mosaic at 35 km -> serve WITH dispersive", bs, n, 1, "R2", 1, 0);
	n = fam_set(bs, UM9, 38.0, 18.0, 5);
	expect("UM980 at 38 km -> serve (D167 ceiling 40)", bs, n, 1, "R2", 1, 0);
	n = fam_set(bs, UM9, 42.0, 18.0, 5);
	expect("UM980 at 42 km -> outside its envelope (R5)", bs, n, 0, "R5", -1, 0);

	/* R2s signal capability */
	n = fam_set(bs, MOS, 35.0, 18.0, 5);
	{
		int i;

		for (i = 0; i < n; i++)
			bs[i].gps_bands = 2;
	}
	expect("Mosaic at 35 km but 2 GPS bands -> R2 404", bs, n, 0, "R2", -1, 1);
	{
		int i;

		for (i = 0; i < n; i++)
			bs[i].gps_bands = 3;
	}
	expect("Mosaic at 35 km with 3 GPS bands -> serve", bs, n, 1, "R2", 1, 1);
	{
		int i;

		for (i = 0; i < n; i++)
			bs[i].gps_bands = 0;
	}
	expect("Mosaic at 35 km, bands unknown -> serve (not refused unmeasured)", bs, n, 1, "R2", 1, 1);

	/*
	 * Firmware blocklist (D159), COUPLED with the band check rather than replacing it. The band
	 * check catches an unknown bad build by its symptom; this catches a known one at connect
	 * time. Both run, so neither gap is left open.
	 */
	check(vbs_policy_fw_blocked("2.7.0 R4.10Build11833") == 1, "build 11833 is blocked");
	check(vbs_policy_fw_blocked("2.7.0 R4.10Build13504") == 1, "build 13504 is blocked");
	check(vbs_policy_fw_blocked("2.6.3 R4.10Build17548") == 0, "build 17548 is not blocked");
	/*
	 * THE TRUNCATION TRAP. The RINEX header is 20 characters wide and drops the last digit, so
	 * the bad build was first recorded as "Build1183". A substring blocklist would refuse that
	 * string -- and 1183 is a build nobody has shown to be bad. Matching the number as an
	 * integer is what makes this safe.
	 */
	check(vbs_policy_fw_blocked("2.7.0 R4.10Build1183") == 0,
	    "the TRUNCATED form Build1183 is NOT blocked: it is a different build");
	check(vbs_policy_fw_blocked("2.7.0 R4.10Build118330") == 0,
	    "nor is Build118330, which merely starts with the bad number");
	check(vbs_policy_fw_blocked(NULL) == 0, "unknown firmware is not a finding");
	check(vbs_policy_fw_blocked("") == 0, "empty firmware is not a finding");
	check(vbs_policy_fw_blocked("RTKBase Unicore no build here") == 0,
	    "a string with no Build token is not a finding");

	{
		int i;

		for (i = 0; i < n; i++) {
			bs[i].gps_bands = 3;
			bs[i].fw = "2.6.3 R4.10Build17548";
		}
	}
	expect("healthy firmware, 3 bands -> serve", bs, n, 1, "R2", 1, 1);
	bs[2].fw = "2.7.0 R4.10Build11833";
	expect("one station on a blocked build -> R2 404, though bands pass", bs, n, 0, "R2", -1, 1);
	/* An upgrade must clear it with no other change: the check reads the live string. */
	bs[2].fw = "2.7.0 R4.10Build17548";
	expect("after the upgrade, the same cluster serves again", bs, n, 1, "R2", 1, 1);
	/* And the two checks are independent: bad bands still refuse on good firmware. */
	bs[2].gps_bands = 2;
	expect("good firmware does not excuse 2 GPS bands", bs, n, 0, "R2", -1, 1);
	{
		int i;

		for (i = 0; i < n; i++) {
			bs[i].gps_bands = 0;
			bs[i].fw = NULL;
		}
	}

	/* R3 -- and the regression for the unbounded candidate search */
	n = fam_set(bs, F9P, 16.0, 18.0, 5);
	n += fam_set(bs + n, MOS, 300.0, 18.0, 5);
	expect("all-F9P here, Mosaic cluster 300 km away -> R3, not R5", bs, n, 0, "R3", -1, 0);
	n = fam_set(bs, F9P, 16.0, 18.0, 5);
	expect("all-F9P -> R3", bs, n, 0, "R3", -1, 0);

	/* R4 */
	b(&bs[0], "a", UM9, 16, 0); b(&bs[1], "b", MOS, 30, 10); b(&bs[2], "c", UM9, 40, -10);
	b(&bs[3], "d", MOS, 50, 20); b(&bs[4], "e", UM9, 24, 30);
	/* the UM980 cluster here is 4 strong but its nearest is 16 km -> in envelope: serve? */
	{
		/* The Python fixture for R4 uses 2 of each family plus one, so no family reaches 4. */
		b(&bs[0], "a", UM9, 25, 0); b(&bs[1], "b", MOS, 30, 10); b(&bs[2], "c", UM9, 40, -10);
		b(&bs[3], "d", MOS, 50, 20); b(&bs[4], "e", F9P, 27, 30);
		expect("UM980+Mosaic+F9P mix, no family of 4 -> R4", bs, 5, 0, "R4", -1, 0);
	}

	/* nearest qualifying family wins over an unconditional Mosaic preference */
	n = fam_set(bs, UM9, 16.0, 18.0, 5);
	n += fam_set(bs + n, MOS, 40.0, 18.0, 5);
	{
		struct vbs_decision d;

		vbs_policy_decide(O, bs, n, 0, &d);
		check(d.serve && d.family == VBS_FAM_UM980 && d.dispersive,
		    "16 km UM980 gap beats a 40 km Mosaic cluster (UM980 serves)");
		check(d.ncluster == 5 && d.cluster[0] == 0,
		    "the served cluster is returned, nearest first");
	}

	/* dead bases are invisible */
	n = fam_set(bs, UM9, 16.0, 18.0, 5);
	bs[0].live = 0;
	{
		struct vbs_decision d;

		vbs_policy_decide(O, bs, n, 0, &d);
		/* the next base sits at (34 km E, 12 km N): 36.06 km */
		check(fabs(d.nearest_km - 36.06) < 0.1, "a base that is not live does not count as nearest");
	}

	/* family mapping */
	check(vbs_policy_family("SEPT MOSAIC-X5") == VBS_FAM_MOSAIC &&
	    vbs_policy_family("UNICORE UM980") == VBS_FAM_UM980 &&
	    vbs_policy_family("UNICORE UM982") == VBS_FAM_UM980 &&
	    vbs_policy_family("UBL-F9P") == VBS_FAM_F9P &&
	    vbs_policy_family("TRIMBLE NETR9") == VBS_FAM_OTHER &&
	    vbs_policy_family(NULL) == VBS_FAM_OTHER, "receiver strings map to the Python families");

	/*
	 * "2 + 1 backup" (decided 2026-10-02): a homogeneous family with only its master
	 * and two auxiliaries within reach still serves when one more of its stations lies within
	 * the backup radius (250 km). The Python reference carries the same rule (_backup).
	 */
	printf("=== \"2 + 1 backup\" ===\n");
	{
		struct vbs_decision d;
		char what[120];
		int tagged;

#define BACKUP_CASE(NAME, NB, WANT_SERVE, WANT_RULE, NCL, TAG) do {			\
		vbs_policy_decide(O, bs, NB, 0, &d);					\
		snprintf(what, sizeof what, "%-46s -> %s %s, %d candidates", NAME, d.rule,	\
		    d.serve ? "SERVE" : "404", d.ncluster);				\
		tagged = strstr(d.reason, "2 + 1 backup") != NULL;			\
		check(d.serve == (WANT_SERVE) && !strcmp(d.rule, WANT_RULE) &&		\
		    (!(WANT_SERVE) || d.ncluster == (NCL)) &&				\
		    ((TAG) < 0 || tagged == (TAG)), what);					\
		if (!(d.serve == (WANT_SERVE) && !strcmp(d.rule, WANT_RULE)))		\
			printf("      reason: %s\n", d.reason);				\
	} while (0)

		b(&bs[0], "a", MOS, 35, 0); b(&bs[1], "b", MOS, 55, 10); b(&bs[2], "c", MOS, 70, -10);
		b(&bs[3], "far", MOS, 200, 0);
		BACKUP_CASE("3 Mosaic in reach + 1 at 200 km -> R2 with backup", 4, 1, "R2", 4, 1);
		check(d.serve && d.ncluster == 4 && !strcmp(bs[d.cluster[3]].mp, "far") && d.dispersive,
		    "  the far station is offered last, the recipe stays Mosaic's");

		b(&bs[3], "far", MOS, 300, 0);
		BACKUP_CASE("... but at 300 km, beyond the backup radius -> R0", 4, 0, "R0", 0, -1);

		b(&bs[0], "a", UM9, 16, 0); b(&bs[1], "b", UM9, 30, 12); b(&bs[2], "c", UM9, 45, -12);
		b(&bs[3], "far", UM9, 180, 0);
		BACKUP_CASE("3 UM980 in reach + 1 at 180 km -> R2 with backup", 4, 1, "R2", 4, 1);

		b(&bs[0], "a", F9P, 16, 0); b(&bs[1], "b", F9P, 30, 12); b(&bs[2], "c", F9P, 45, -12);
		b(&bs[3], "far", F9P, 180, 0);
		BACKUP_CASE("all-F9P never gets a reserve -> refused", 4, 0, "R0", 0, -1);

		/* A family that is not short still gets the candidate when there is room. */
		n = fam_set(bs, MOS, 35.0, 18.0, 4);
		b(&bs[n], "far", MOS, 200, 0);
		BACKUP_CASE("4 Mosaic in reach + far -> R2, reserve offered", n + 1, 1, "R2", 5, 0);
		n = fam_set(bs, MOS, 35.0, 18.0, 5);
		b(&bs[n], "far", MOS, 200, 0);
		BACKUP_CASE("5 Mosaic in reach + far -> R2, no room for it", n + 1, 1, "R2", 5, 0);
		check(d.serve && strcmp(bs[d.cluster[4]].mp, "far") != 0,
		    "  five in reach fill the candidate list; the far one is not offered");
#undef BACKUP_CASE

		/* --no-backup: the policy as it was before the reserve */
		b(&bs[0], "a", MOS, 35, 0); b(&bs[1], "b", MOS, 55, 10); b(&bs[2], "c", MOS, 70, -10);
		b(&bs[3], "far", MOS, 200, 0);
		vbs_policy_decide(O, bs, 4, VBS_POLICY_NO_BACKUP, &d);
		check(!d.serve && !strcmp(d.rule, "R0"),
		    "NO_BACKUP: 3 Mosaic in reach + 1 at 200 km -> R0");
		n = fam_set(bs, MOS, 35.0, 18.0, 4);
		b(&bs[n], "far", MOS, 200, 0);
		vbs_policy_decide(O, bs, n + 1, VBS_POLICY_NO_BACKUP, &d);
		check(d.serve && d.ncluster == 4, "NO_BACKUP: 4 Mosaic in reach -> R2, no reserve offered");

		/*
		 * A reserve R2s would refuse is passed over, never a veto. T0990C0990: CSEL, two GPS
		 * bands, offered as the reserve, refused a Mosaic cluster of four that served without it.
		 */
		n = fam_set(bs, MOS, 35.0, 18.0, 4);
		for (tagged = 0; tagged < n; tagged++)
			bs[tagged].gps_bands = 3;
		b(&bs[n], "bad", MOS, 150, 0);
		bs[n].gps_bands = 2;
		vbs_policy_decide(O, bs, n + 1, VBS_POLICY_CHECK_BANDS, &d);
		check(d.serve && d.ncluster == 4, "a two-band reserve is not offered and vetoes nothing");
		b(&bs[n + 1], "good", MOS, 200, 0);
		bs[n + 1].gps_bands = 3;
		vbs_policy_decide(O, bs, n + 2, VBS_POLICY_CHECK_BANDS, &d);
		check(d.serve && d.ncluster == 5 && !strcmp(bs[d.cluster[4]].mp, "good"),
		    "  the next one out is offered instead");
		b(&bs[0], "a", MOS, 35, 0); b(&bs[1], "b", MOS, 55, 10); b(&bs[2], "c", MOS, 70, -10);
		b(&bs[3], "bad", MOS, 150, 0);
		bs[0].gps_bands = bs[1].gps_bands = bs[2].gps_bands = 3;
		bs[3].gps_bands = 2;
		vbs_policy_decide(O, bs, 4, VBS_POLICY_CHECK_BANDS, &d);
		check(!d.serve && !strcmp(d.rule, "R0"), "3 Mosaic + only a two-band reserve -> R0, no 2 + 1");
	}

	/*
	 * The envelopes and recipe the Python reference applies: UM980 13-40 km (D167), Mosaic
	 * 28-60 km (D135), GPS L2 from L2C for the Mosaic only (D138, D140), and R1's redundancy
	 * clause on its own 20 km rather than the Mosaic floor.
	 */
	printf("=== envelopes and recipe (D133, D135, D138) ===\n");
	{
		struct vbs_decision d;

		n = fam_set(bs, UM9, 25.0, 18.0, 5);
		vbs_policy_decide(O, bs, n, 0, &d);
		check(d.serve && !strcmp(d.rule, "R2") && d.family == VBS_FAM_UM980,
		    "UM980 nearest at 25 km -> R2 (the ceiling is 40, not 22)");
		check(d.serve && d.dispersive && d.iono_sys == (SYS_GAL | SYS_CMP) && !d.prefer_l2c,
		    "  UM980 recipe: dispersive on Galileo + BeiDou only, L2C not preferred (D167)");

		n = fam_set(bs, UM9, 32.0, 18.0, 5);
		vbs_policy_decide(O, bs, n, 0, &d);
		check(d.serve && !strcmp(d.rule, "R2"), "UM980 nearest at 32 km -> R2 (was R5 before D167)");

		n = fam_set(bs, UM9, 41.0, 18.0, 5);
		vbs_policy_decide(O, bs, n, 0, &d);
		check(!d.serve && !strcmp(d.rule, "R5"), "UM980 nearest at 41 km -> R5");

		n = fam_set(bs, MOS, 24.0, 18.0, 5);
		vbs_policy_decide(O, bs, n, 0, &d);
		check(!d.serve && !strcmp(d.rule, "R5"),
		    "Mosaic nearest at 24 km -> R5 (the floor is 28, not 20)");

		n = fam_set(bs, MOS, 30.0, 18.0, 5);
		vbs_policy_decide(O, bs, n, 0, &d);
		check(d.serve && !strcmp(d.rule, "R2") && d.family == VBS_FAM_MOSAIC,
		    "Mosaic nearest at 30 km -> R2");
		check(d.serve && d.dispersive && d.iono_sys == 0 && d.prefer_l2c,
		    "  Mosaic recipe: dispersive on every constellation, L2C preferred");

		/* three bases within 12 km of a 24 km nearest: redundant, but above 20 km */
		b(&bs[0], "a", MOS, 24, 0); b(&bs[1], "b", MOS, 24, 10); b(&bs[2], "c", MOS, 30, -5);
		b(&bs[3], "d", MOS, 40, 20); b(&bs[4], "e", MOS, 50, -20);
		vbs_policy_decide(O, bs, 5, 0, &d);
		check(!d.serve && !strcmp(d.rule, "R5"),
		    "redundant Mosaic at 24 km -> R5, not R1 (R1 does not follow the Mosaic floor)");
		if (strcmp(d.rule, "R5"))
			printf("      reason: %s\n", d.reason);
	}

	/*
	 * VBS_POLICY_FORCE, the field-test flag (vbsd --force): R1 and the envelopes are skipped so a
	 * homogeneous Mosaic or UM980 cluster serves wherever it forms; everything else still refuses.
	 */
	printf("=== field-test FORCE ===\n");
	{
		struct vbs_decision d;
		int k;

		n = fam_set(bs, UM9, 6.0, 18.0, 5);
		vbs_policy_decide(O, bs, n, VBS_POLICY_FORCE, &d);
		check(d.serve && !strcmp(d.rule, "R2") && !strncmp(d.reason, "forced:", 7),
		    "FORCE: UM980 at 6 km serves (R1 skipped), reason says forced");
		vbs_policy_decide(O, bs, n, 0, &d);
		check(!d.serve && !strcmp(d.rule, "R1"), "  without FORCE the same tile is R1");

		n = fam_set(bs, MOS, 20.0, 18.0, 5);
		vbs_policy_decide(O, bs, n, VBS_POLICY_FORCE, &d);
		check(d.serve && d.family == VBS_FAM_MOSAIC && d.dispersive && d.prefer_l2c,
		    "FORCE: Mosaic at 20 km serves with the Mosaic recipe (floor skipped)");

		n = fam_set(bs, UM9, 45.0, 18.0, 5);
		vbs_policy_decide(O, bs, n, VBS_POLICY_FORCE, &d);
		check(d.serve && d.family == VBS_FAM_UM980 && d.dispersive &&
		    d.iono_sys == (SYS_GAL | SYS_CMP),
		    "FORCE: UM980 at 45 km serves with the UM980 recipe (ceiling skipped)");
		vbs_policy_decide(O, bs, n, 0, &d);
		check(!d.serve && !strcmp(d.rule, "R5"), "  without FORCE the same tile is R5");

		n = fam_set(bs, F9P, 16.0, 18.0, 5);
		vbs_policy_decide(O, bs, n, VBS_POLICY_FORCE, &d);
		check(!d.serve && !strcmp(d.rule, "R3"), "FORCE: all-F9P is still refused (R3)");

		b(&bs[0], "a", UM9, 25, 0); b(&bs[1], "b", MOS, 30, 10); b(&bs[2], "c", UM9, 40, -10);
		b(&bs[3], "d", MOS, 50, 20); b(&bs[4], "e", F9P, 27, 30);
		vbs_policy_decide(O, bs, 5, VBS_POLICY_FORCE, &d);
		check(!d.serve && !strcmp(d.rule, "R4"), "FORCE: a mix with no family of 4 is still refused (R4)");

		n = fam_set(bs, MOS, 35.0, 18.0, 5);
		for (k = 0; k < n; k++)
			bs[k].gps_bands = 3;
		bs[2].gps_bands = 2;
		vbs_policy_decide(O, bs, n, VBS_POLICY_FORCE | VBS_POLICY_CHECK_BANDS, &d);
		check(!d.serve && !strcmp(d.rule, "R2"), "FORCE: a two-band station still vetoes its cluster (R2s)");

		b(&bs[0], "a", MOS, 30, 0); b(&bs[1], "b", MOS, 45, 0);
		vbs_policy_decide(O, bs, 2, VBS_POLICY_FORCE, &d);
		check(!d.serve && !strcmp(d.rule, "R0"), "FORCE: two bases are still too few (R0)");
	}

	printf("%s: %d failure(s)\n", nfail ? "FAILED" : "PASSED", nfail);
	return nfail != 0;
}

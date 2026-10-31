#!/bin/sh
#
# --decide must answer the same question the server answers.
#
# It used to stop after vbs_policy_decide(), so it reported "SERVE" for positions the server
# then refused: the rules admit four homogeneous receivers, the engine's gate drops one of them
# for sitting in another reference frame, and the cluster is left an auxiliary short. The caster's
# author hit exactly that on a Belgian tile -- the tool said SERVE, the endpoint said
# "cluster: too few auxiliaries". Both paths now go through decide_tile().
#
# Run from vbs-engine/ (as `make test` does). Needs tools/vbsd. Checks 1-3 also need a full network cache in VBS_TEST_CACHE.
set -u
VBSD=tools/vbsd
CACHE=${VBS_TEST_CACHE:-}
fail=0

if [ ! -x "$VBSD" ]; then
	echo "SKIPPED: need $VBSD"
	exit 0
fi
if [ -n "$CACHE" ] && [ -f "$CACHE" ]; then

# 1. The reported case: a cross-frame tile. Must refuse, name the cluster reason, and say which
#    base was gated and why -- a bare "too few auxiliaries" is what made it unexplainable.
out=$("$VBSD" --cache "$CACHE" --decide 50.29252 5.70688 2>&1); rc=$?
for want in "404" "cluster: too few auxiliaries" "GATED: frame mismatch"; do
	case "$out" in
	*"$want"*) ;;
	*) echo "FAIL: cross-frame tile output missing '$want'"; fail=1 ;;
	esac
done
[ "$rc" = 2 ] || { echo "FAIL: cross-frame tile exit $rc, want 2"; fail=1; }
case "$out" in
*SERVE*) echo "FAIL: cross-frame tile reported SERVE; the server refuses it"; fail=1 ;;
esac

# 2. A tile that genuinely serves: the roles and the cluster summary must be reported.
out=$("$VBSD" --cache "$CACHE" --decide 45.5 2.0 2>&1); rc=$?
for want in "SERVE" "master" "aux" "cluster: master"; do
	case "$out" in
	*"$want"*) ;;
	*) echo "FAIL: serving tile output missing '$want'"; fail=1 ;;
	esac
done
[ "$rc" = 0 ] || { echo "FAIL: serving tile exit $rc, want 0"; fail=1; }

# 3. The invariant, over a spread of positions: SERVE implies a cluster was actually formed
#    (exit 0 and a cluster summary). A future divergence would break this without being noticed
#    on the two cases above.
for p in "50.29252 5.70688" "45.5 2.0" "42.5 1.0" "48.82 2.34" "47.5 -2.0" "43.5 7.0" "46.5 2.5"; do
	# shellcheck disable=SC2086
	out=$("$VBSD" --cache "$CACHE" --decide $p 2>&1); rc=$?
	case "$out" in
	*SERVE*)
		[ "$rc" = 0 ] || { echo "FAIL: $p says SERVE but exits $rc"; fail=1; }
		case "$out" in
		*"cluster: master"*) ;;
		*) echo "FAIL: $p says SERVE with no cluster formed"; fail=1 ;;
		esac ;;
	*)
		[ "$rc" = 2 ] || { echo "FAIL: $p refused but exits $rc"; fail=1; } ;;
	esac
done

# The coordinate rule (2026-10-02): a base absent from the coordinate document is untrusted and
# ignored. The production document is not redistributable, so the checks use the invented
# example, which names none of the cache's bases -- with it, nothing may be served anywhere.
out=$("$VBSD" --cache "$CACHE" --coords tests/coords_cases.json --decide 45.5 2.0 2>&1)
case "$out" in
*"1520 of 1520 bases not positioned"*|*"bases not positioned by the document are ignored"*) ;;
*) echo "FAIL: --coords does not report the untrusted bases"; fail=1 ;;
esac
case "$out" in
*SERVE*) echo "FAIL: served from bases the document does not name"; fail=1 ;;
esac
# ... and the server itself refuses to start without a document.
"$VBSD" --cache "$CACHE" --caster 127.0.0.1:1 --listen 0 >/dev/null 2>&1
rc=$?
[ "$rc" = 2 ] || { echo "FAIL: the daemon started without --coords (exit $rc)"; fail=1; }

else
	echo "SKIPPED: checks 1-3 (set VBS_TEST_CACHE to a full network cache)"
fi

# 4. "2 + 1 backup" never serves less than the engine without it. The reserve makes a short UM980
#    family eligible; the nearest UM980 is in another frame, so that cluster fails the gate. The
#    tile must still be served by the Mosaic cluster the decision without the reserve picks.
#    Shape of T1003C0977 (Sion): synthetic stations, since the real coordinates are not ours.
FX=$(mktemp -d)
cat > "$FX/cache.json" <<'JSON'
{"version":1,"generated":"fixture","bases":[
{"mp":"UA","marker":"T00000XXX","rec":"UNICORE UM980","ant":"","frame":"F2","xyz":[4400849.09,561510.227,4567673.823],"live":true},
{"mp":"UB","marker":"T00001XXX","rec":"UNICORE UM980","ant":"","frame":"F1","xyz":[4377600.273,526425.789,4593939.562],"live":true},
{"mp":"UC","marker":"T00002XXX","rec":"UNICORE UM980","ant":"","frame":"F1","xyz":[4446185.252,544861.338,4525873.945],"live":true},
{"mp":"UR","marker":"T00003XXX","rec":"UNICORE UM980","ant":"","frame":"F1","xyz":[4378019.129,718132.963,4567673.823],"live":true},
{"mp":"MA","marker":"T00004XXX","rec":"SEPT MOSAIC-X5","ant":"","frame":"F1","xyz":[4374965.319,536133.634,4595317.488],"live":true},
{"mp":"MB","marker":"T00005XXX","rec":"SEPT MOSAIC-X5","ant":"","frame":"F1","xyz":[4394310.25,583885.473,4571139.067],"live":true},
{"mp":"MC","marker":"T00006XXX","rec":"SEPT MOSAIC-X5","ant":"","frame":"F1","xyz":[4439723.866,538982.277,4532868.326],"live":true},
{"mp":"MD","marker":"T00007XXX","rec":"SEPT MOSAIC-X5","ant":"","frame":"F1","xyz":[4409946.566,484910.025,4567673.823],"live":true},
{"mp":"ME","marker":"T00008XXX","rec":"SEPT MOSAIC-X5","ant":"","frame":"F1","xyz":[4369097.401,576519.829,4596006.283],"live":true}]}
JSON
cat > "$FX/coords.json" <<'JSON'
{"sources":{
"UA":{"x_m":4400849.09,"y_m":561510.227,"z_m":4567673.823,"reference_frame":"F2","calculation_quality":"FINAL"},
"UB":{"x_m":4377600.273,"y_m":526425.789,"z_m":4593939.562,"reference_frame":"F1","calculation_quality":"FINAL"},
"UC":{"x_m":4446185.252,"y_m":544861.338,"z_m":4525873.945,"reference_frame":"F1","calculation_quality":"FINAL"},
"UR":{"x_m":4378019.129,"y_m":718132.963,"z_m":4567673.823,"reference_frame":"F1","calculation_quality":"FINAL"},
"MA":{"x_m":4374965.319,"y_m":536133.634,"z_m":4595317.488,"reference_frame":"F1","calculation_quality":"FINAL"},
"MB":{"x_m":4394310.25,"y_m":583885.473,"z_m":4571139.067,"reference_frame":"F1","calculation_quality":"FINAL"},
"MC":{"x_m":4439723.866,"y_m":538982.277,"z_m":4532868.326,"reference_frame":"F1","calculation_quality":"FINAL"},
"MD":{"x_m":4409946.566,"y_m":484910.025,"z_m":4567673.823,"reference_frame":"F1","calculation_quality":"FINAL"},
"ME":{"x_m":4369097.401,"y_m":576519.829,"z_m":4596006.283,"reference_frame":"F1","calculation_quality":"FINAL"}}}
JSON
out=$("$VBSD" --cache "$FX/cache.json" --coords "$FX/coords.json" --decide 46.0 7.0 2>&1); rc=$?
case "$out" in
*"R2 SERVE: homogeneous MOSAIC"*) ;;
*) echo "FAIL: the reserve displaced a Mosaic cluster that serves without it (exit $rc)"; fail=1 ;;
esac

# 5. Field-test flags. --exclude-km hides the bases near the tile: UA (22 km) must vanish from
#    the decision and be reported as hidden. --force serves and says so, and keeps what
#    production would have answered on the same visible network.
out=$("$VBSD" --cache "$FX/cache.json" --coords "$FX/coords.json" --exclude-km 25 --decide 46.0 7.0 2>&1)
case "$out" in
*"hidden 1 within 25 km: UA "*) ;;
*) echo "FAIL: --exclude-km 25 does not report UA as hidden"; fail=1 ;;
esac
case "$out" in
*"  UA "*) echo "FAIL: a hidden base appears in the decision"; fail=1 ;;
esac
out=$("$VBSD" --cache "$FX/cache.json" --coords "$FX/coords.json" --force --decide 46.0 7.0 2>&1)
for want in "TEST MODE" "SERVE: forced:" "test: production R2 SERVE"; do
	case "$out" in
	*"$want"*) ;;
	*) echo "FAIL: --force output missing '$want'"; fail=1 ;;
	esac
done
rm -rf "$FX"

if [ "$fail" = 0 ]; then
	echo "PASSED: --decide agrees with the serve decision"
else
	echo "TESTS FAILED"
fi
exit "$fail"

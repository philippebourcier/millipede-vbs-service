# millipede-vbs-service: specification

This document describes what the service does, the measured rules that decide where it serves,
the files it reads, and how to build, test and deploy it. The code comments cite decision numbers
(`D146`, `D167`, ...) and research notes from the study that produced the engine. Those notes are
not published; the comment next to each citation carries the result that matters. The maths and
physics are in [THEORY.md](THEORY.md).

## 1. What a VBS is, and when it helps

A rover doing RTK against one physical base suffers from the atmosphere and orbit errors that do
not cancel over the baseline. Past 15-20 km, ambiguity resolution becomes slow and, worse, wrong
integers are fixed with a centimetre-level sigma. A VBS uses a cluster of bases around the
rover: it resolves the network's own ambiguities, estimates the differential ionosphere and
troposphere along each baseline, interpolates them to a point next to the rover, and synthesises
the observations a base at that point would have recorded. The rover sees one short-baseline base.

What the measurements say about where that is worth it, in short:

- Below about **13 km** of the nearest base, direct RTK is already at the acceptance ceiling. A
  VBS adds nothing and is refused.
- Between 13 and 40-60 km, a VBS beats the nearest base **only for some receiver families and
  only on homogeneous clusters**. The failure mode outside the envelope is not a worse fix but
  a wrong one.
- Past 60 km neither direct RTK nor the VBS is usable.

Acceptance, relative to the alternative at the same place: fix rate at least 90 %, H95/V95 at
most 10 cm, false fix at most 1 %.

## 2. Architecture

```
rover --NTRIP--> millipede-caster --(GET /<lat>/<lon>/<alt> + GGA)--> vbsd --NTRIP--> millipede-caster
                  (NEAR + VBS window)                                    (pulls the cluster's bases)
```

- **millipede-caster** keeps its NEAR mountpoint. With the patch, a configured VBS mountpoint
  serves NEAR as usual when the nearest live base is closer than `vbs_min_km` or further than
  `vbs_max_km`. In between, the caster opens one upstream connection to `vbsd`, addressed by
  position, and forwards the rover's GGA on it. A 404 from `vbsd` means "no VBS here": the
  caster falls back to NEAR. The caster makes only this coarse distance test; every finer rule
  lives in `vbsd`, so rules can change without redeploying the caster.
- **vbsd** is a separate process and may run on another host. It subscribes to the physical
  bases it needs from the caster over plain NTRIP, runs one synthesis per occupied tile, and
  streams RTCM to every rover on that tile. It reads tile changes from the GGA itself, so the
  caster never reconnects when a rover moves.

### 2.1 Tiles

The world is divided once into cells of about 15 km (`VBS_TILE_KM`): rows equal in latitude,
columns sized per row. A tile is a pure function of position, so two rovers at the same place
share one synthesis. The RTCM station ID is derived from the tile so that every spatial
neighbour has a different one. A rover leaving a tile sees its base change, which RTCM signals
by the station ID changing. Tile changes use a hysteresis (`--hysteresis-m`, default 300 m).

### 2.2 Session life

1. A client asks for `/<lat>/<lon>/<alt>`. `vbsd` applies the serving rules (§4). If they
   refuse, it answers `404` with `X-VBS-Rule:` and `X-VBS-Reason:` headers.
2. Otherwise it forms the tile's cluster (a master and its auxiliaries) and opens the feeds it
   needs. Feeds are shared between tiles.
3. **Pass-through.** The network filters need about three minutes before the first correction
   passes its gate, and a rover must not sit on an empty stream. Until then the session relays
   the nearest live base's own RTCM verbatim, which is exactly what NEAR would have sent. When
   the first VBS epoch is ready it switches to the tile's 1005/1033/1230 and the synthesised
   observations. The rover sees one base change.
4. A feed that stops delivering (`--stale-s`) is dropped from the live view, and every session
   using it is re-formed from the remaining stations. A tile with no client for `--idle-s` is
   closed.

## 3. Processing pipeline (per tile, per epoch)

| Step | Module | What it does |
|---|---|---|
| Cluster | `vbs_cluster` | Master = nearest eligible base; three or more auxiliaries around the point; frames never mixed; "2 + 1 backup" (below) |
| Base check | `vbs_basemeta` | Declared coordinate against the reference document; solid-earth tide; refuse a base outside the acceptance gate |
| Pre-processing | `vbs_preproc`, `vbs_lock` | Cycle slips (geometry-free and Melbourne-Wübbena, plus MSM lock-time indicators), station arcs |
| Signals | `vbs_signal`, `vbs_c2wcal` | Band matching across receivers (L2 P(Y) against L2C), GLONASS FDMA frequencies |
| Network AR | `vbs_netekf`, `vbs_ambres`, `vbs_ddres` | Per-baseline EKF with receiver clocks, troposphere and horizontal gradients, double-differenced ambiguities by LAMBDA; wide lane, then narrow lane on the ionosphere-free combination for long baselines |
| Corrections | `vbs_netekf`, `vbs_ionopost` | Dispersive (ionosphere) and non-dispersive (troposphere, orbit) corrections per satellite and baseline |
| Interpolation | `vbs_interp` | Least-squares collocation anchored at the master (its own correction is zero by definition) |
| Synthesis | `vbs_synth`, `vbs_session` | Master observations translated to the VBS point: geometry, interpolated corrections, troposphere model |
| Output | `vbs_rtcm_out` | RTCM 3: 1005 (tile point, constant within a session), 1033 and 1230 from the master, MSM7 1077/1087/1097/1127 |
| Orbits | `vbs_peph` | Broadcast ephemerides from the feeds; precise orbits optional (tools only) |

Four choices carry most of the performance:

- horizontal troposphere gradients per baseline, tightly constrained, with the mapping
  function capped;
- the interpolation anchored at the master;
- collocation instead of a plane fit. The plane cannot carry curvature, and the 10 km tile
  only works with collocation;
- for long clusters, the narrow lane resolved on the ionosphere-free combination
  (`--if-narrow-lane`).

**"2 + 1 backup".** A family with only its master and two auxiliaries within 120 km still forms a
cluster when one more of its stations lies within 250 km. That station is processed but held in
reserve, and is fitted only for a satellite the primaries lose. `vbsd --no-backup` disables it.

## 4. The serving decision

`vbs_policy` decides per position, before any stream is opened. Rules are evaluated in order;
the first that fires decides. The receiver family comes from the base cache (`rec`): Septentrio
Mosaic-X5, Unicore UM980/UM982, u-blox F9P, other.

| Rule | Name | Fires when | Answer |
|---|---|---|---|
| R0 | density | fewer than 4 live, trusted bases within 120 km | 404 |
| R1 | covered | nearest base under 13 km, or 2 or more bases within 12 km of the nearest (below 20 km) | 404 |
| R2 | homogeneous | a single-family cluster whose nearest station is inside that family's envelope; the nearest qualifying family wins | **serve** |
| R2s | capability | a station of that cluster delivers fewer than 3 GPS bands | 404 |
| R3 | all-F9P | only an F9P cluster within reach | 404 |
| R4 | mixed | the nearest five within reach span more than one family | 404 |
| R5 | out of range | nearest base outside every envelope | 404 |

Envelopes and per-family recipe (`ENVELOPE` in `src/vbs_policy.c`):

| Cluster | Serve between | Dispersive correction | GPS L2 used |
|---|---|---|---|
| Homogeneous Unicore UM980 / UM982 | 13-40 km | on, Galileo + BeiDou only (`--iono-sys EC`) | first slot |
| Homogeneous Septentrio Mosaic-X5 | 28-60 km | on, all constellations | L2C (`--prefer-l2c`) |
| All u-blox F9P | never | n/a | n/a |
| Mixed families | never | n/a | n/a |

Why these, briefly:

- **UM980 GPS ionosphere.** The network ionosphere is reliable on Galileo and BeiDou and erratic
  on GPS, whose L2 P(Y) code carries a per-unit bias of metres that drifts over days.
  Correcting that bias from the base's own L2C was built (`--c2w-cal`) and tested; it does not
  help, and it stays off.
- **All-F9P.** Refused because it is erratic, not because it always loses: nothing separates
  the good cases from the bad ones in advance.
- **Mixed clusters.** On leave-one-out tests over four days, a mixed cluster was never better
  than the nearest base.

Firmware builds known to produce unusable observations are refused at connect time
(`fw_bad_builds`: UM980 builds 11833 and 13504 declare an L2C column and never fill it).

**Coordinates are a precondition.** A base absent from the coordinate document (§5.2), refused
by it, or marked exploratory is *untrusted*: it is ignored entirely, with no cluster, no density
count and no pass-through. `vbsd` refuses to start without `--coords`. Shifts of 7-9 mm in base
coordinates were measured to halve a rover's harmful distance, so registry coordinates are not
good enough.

`vbsd --decide LAT LON` prints the decision for a position, with each base's role, and exits.
It goes through the same code path as the server.

## 5. Input files

### 5.1 Base cache (`--cache`)

What `vbsd` knows about every base before it has seen a byte of its stream. One object per
base; `examples/base_cache.json` is a complete small network.

```json
{"version": 1, "generated": "2026-10-01T00:00:00Z", "bases": [
  {"mp": "UA", "marker": "T00000XXX", "rec": "UNICORE UM980", "ant": "", "frame": "F1",
   "xyz": [4400849.090, 561510.227, 4567673.823], "gps_bands": "125", "fw": "", "live": true}
]}
```

| Field | Meaning |
|---|---|
| `mp` | mountpoint on the caster; also the key into the coordinate document |
| `marker` | RINEX marker name, informational |
| `rec` | receiver type; selects the family (`MOSAIC`, `UM98`/`UNICORE`, `F9P`, other) |
| `ant` | antenna type, informational |
| `frame` | reference frame of `xyz`; clusters never mix frames |
| `xyz` | ECEF, metres |
| `gps_bands` | GPS bands the station was *observed* to deliver, digits of `{1,2,5}`; absent means unknown, which is not treated as bad |
| `fw` | receiver firmware (RTCM 1033 string); empty means unknown |
| `live` | the registry's opinion; overridden at run time by stream health |

`SIGHUP` re-reads the cache (and re-applies the coordinate document) without dropping a rover.

### 5.2 Coordinate document (`--coords`)

The reference solution for each base, keyed by mountpoint. `examples/coordinates.json` matches
the example cache.

```json
{"sources": {
  "UA": {"x_m": 4400849.090, "y_m": 561510.227, "z_m": 4567673.823,
         "reference_frame": "F1", "calculation_quality": "FINAL"}
}}
```

- A shift larger than `--coords-max-m` (default 5 m) is refused, counted and logged, never
  applied. A wrong mountpoint looks exactly like a large shift. A vertical-only difference is
  accepted as the registry's height being wrong.
- `calculation_quality: EXPLORATORY_GLOBAL_BASELINE` is refused.
- `reference_frame` replaces the cache's `frame`.
- A base the document does not position is untrusted (§4).

`tests/coords_cases.json` exercises each outcome (applied, past the bound, exploratory, absent).

## 6. `vbsd` options

```
vbsd --cache FILE --coords FILE --caster HOST[:PORT] [--user U --password P]
     [--listen PORT] [--allow IP[/bits],...] [--hysteresis-m M] [--stale-s S] [--idle-s S]
     [--no-check-bands] [--max-aux N] [--min-aux N] [--no-backup] [--coords-max-m M]
     [--c2w-cal run:SEC|l2c:SEC]
     [--force] [--exclude-km R]          field tests only, never in production
vbsd --cache FILE [--coords FILE] --decide LAT LON
```

- `--caster`: where the physical bases are pulled from.
- `--listen` (2102 in the unit file) and `--allow`: only the millipede caster should reach
  `vbsd`; rovers authenticate at the caster.
- `--force` serves wherever a homogeneous cluster forms, and `--exclude-km` hides the bases near
  each tile so a VBS can be compared with a real base on the same spot. Both print `TEST MODE`.
- The log goes to stderr; the systemd unit sends it to journald.

## 7. Deployment

### 7.1 millipede-caster

`vbs-engine/deploy/millipede-caster-vbs.patch` is one diff against the upstream commit named at
its top. It adds the `vbs:` configuration block, the serving hook, and their tests
(`caster/tests.c` window cases, `tests/test-vbs-serving.py`).

```sh
cd millipede-caster && patch -p1 < ../millipede-vbs-service/vbs-engine/deploy/millipede-caster-vbs.patch
```

Configuration (also in the patched `sample-config/caster.yaml`):

```yaml
vbs:
  - host:            vbscaster.lan     # where vbsd listens
    port:            2102
    vbs_min_km:      13                # the caster's coarse window; vbsd refines it
    vbs_max_km:      60
    hysteresis_m:    500
    vbs_mountpoint:  VBS               # rovers ask for this mountpoint and send GGA
```

The `vbsd` host is exempted from the per-IP connection quota, since it pulls one stream per base
in use.

### 7.2 vbsd

`vbs-engine/deploy/vbsd.service` is a systemd unit (user `vbs`, files under `/opt/vbs`). Set
`--caster` to your caster, install the cache and the coordinate document, then
`systemctl enable --now vbsd`. `systemctl reload vbsd` re-reads both files.

## 8. Build and test

Dependencies: a C11 compiler (gcc or clang; the tree builds with `-Werror` under both), LAPACK
and BLAS, and a built [MRTKLIB](https://github.com/h-shiono/MRTKLIB) tree (`build/libmrtklib.a`).

```sh
cd vbs-engine
make MRTKLIB_DIR=/path/to/MRTKLIB -j
make MRTKLIB_DIR=/path/to/MRTKLIB test
```

MRTKLIB's compile definitions (`NFREQ=5`, `NEXOBS=5`, `ENA*`) change `obsd_t`. The Makefile
mirrors them, and a mismatch gives silently misaligned structures, so keep them in sync with
MRTKLIB's CMakeLists. `.buildflags` rebuilds everything when the flags change.

`make test` runs the unit tests and `tests/test_decide.sh`. Without extra data:

- tests on synthetic data run fully;
- tests on real observations print `SKIPPED`. To run them, put RINEX and RTCM files in
  `vbs-engine/testdata/` under the names the tests give, with the matching coordinates;
- `test_zero_baseline` also needs RTKLIB (demo5) as an independent RTCM decoder
  (`RTKLIB_DIR`);
- checks 1-3 of `test_decide.sh` need a full network cache in `VBS_TEST_CACHE`.

## 9. Tools

| Tool | Purpose |
|---|---|
| `tools/vbsd` | the VBS caster |
| `tools/vbs_pipeline` | the offline pipeline on RINEX files: every engine option, writes the synthesised observations |
| `tools/vbs_loo` | leave-one-out: what the corrections must remove at a real station, from its neighbours |
| `tools/vbs_loo_select` | the server's own cluster choice for each station left out, to build leave-one-out cases from a cache |
| `tools/vbs_netrun` | the network EKF and the ambiguity cascade over one baseline: fix rate per step, ionosphere sigma |
| `tools/vbs_record` | record live NTRIP streams to files (`--host`, `--port`) |
| `tools/vbs_rtcm_inv` | what is inside a recorded RTCM 3 stream, and whether MRTKLIB decodes it |
| `tools/vbs_xcheck` | the live RTCM decode against a RINEX conversion of the same stream |
| `tools/vbs_tidecheck` | how much of the solid-earth, pole and ocean-loading tide survives the interpolation |

## 10. Known limits

- The envelopes above were measured on one network's receivers and geometry. Different
  densities or receivers may move them. The leave-one-out tools exist to re-measure them on
  your own data.
- Nothing at run time predicts a wrong fix. The rules are a prior, applied before a stream is
  opened.
- The narrow lane on the raw carrier is not used; the shipped constants keep it gated off.

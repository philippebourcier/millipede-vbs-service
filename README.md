# millipede-vbs-service

A Virtual Base Station (VBS) service for [millipede-caster](https://github.com/pbeyssac/millipede-caster).

A rover connects to the caster as it would for NEAR. When its nearest base is too far for good
direct RTK but a suitable cluster of bases surrounds it, the caster hands the connection to
`vbsd`. `vbsd` synthesises a base close to the rover, from a network solution over that cluster,
and streams it as ordinary RTCM 3. Where the rules say a VBS would not beat the nearest base,
`vbsd` answers 404 and the caster falls back to NEAR.

- `vbs-engine/`: the engine library, `vbsd` (the VBS caster), research tools and tests.
- `vbs-engine/deploy/millipede-caster-vbs.patch`: the caster side, one patch against upstream.
- `vbs-engine/examples/`: the two JSON files `vbsd` reads, for a small invented network.
- [`SPEC.md`](SPEC.md): how it works, the serving rules, the file formats, build and deployment.

## Quick start

```sh
# MRTKLIB (BSD-2) supplies RTCM, RINEX, ephemerides and LAMBDA
git clone https://github.com/h-shiono/MRTKLIB.git && (cd MRTKLIB && cmake -B build && cmake --build build -j)

git clone https://github.com/philippebourcier/millipede-vbs-service.git
cd millipede-vbs-service/vbs-engine
make MRTKLIB_DIR=../../MRTKLIB -j
make MRTKLIB_DIR=../../MRTKLIB test

# Ask the serving rules about a position, without starting the server
tools/vbsd --cache examples/base_cache.json --coords examples/coordinates.json --decide 46.0 7.0
```

Then patch and configure millipede-caster, and run `vbsd` against it: see SPEC.md §7.

## License

BSD 3-Clause, see [LICENSE](LICENSE). MRTKLIB is BSD-2-Clause and is not included.

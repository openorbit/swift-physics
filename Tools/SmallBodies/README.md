# Small-body POI extraction

Python 3 standard library only. Run from the repository root:

```sh
python3 Tools/SmallBodies/extract.py select
python3 Tools/SmallBodies/extract.py enrich --output Data/JPL/pois-enriched.json
python3 Tools/SmallBodies/extract.py kernels --catalog Data/JPL/pois-enriched.json --start 2026-01-01 --stop 2036-01-01
```

Coverage in the last command is an example, not an adopted simulation interval.
Select actual application coverage before requesting the full kernel collection.
`enrich` and `kernels` accept `--limit N` for smoke checks. `select` accepts
`--asteroids N`, `--trojans N` (0 disables the supplement), `--comets N`,
`--manifest FILE`, and `--output FILE`. Default selection is 1000 known-diameter
asteroids, a 50-Trojan minimum inclusion list, mandatory dwarf planets, and 100 comets.
All are unioned/deduplicated by resolved SBDB SPK-ID; extra mandatory objects can
increase the asteroid count. Top-1000 size ranking is among **defined SBDB effective
diameters**, not a complete or uncertainty-aware physical size ranking.

On 2026-10-13 the generated selection contains 1004 asteroid-category objects and
100 comets. All 50 Trojan supplement entries were already in the top 1000; the
catalog contains 146 TJN objects overall. Ceres, Pluto, Haumea, Eris and Makemake
are included explicitly. SBDB query diameters for the latter four are absent in
this snapshot; the script preserves null values, not inferred zeros or invented
sizes. Their physical models need independently sourced data. Pluto's existing
library provider remains a separate application choice; selection does not replace it.
Greek/Trojan L4/L5 camps are not assigned by the TJN class field.

`selection.json` is editorial input. Its comet seeds cover familiar historical
apparitions, spacecraft targets, meteor-stream parents and an interstellar comet.
Remaining places are filled round-robin among SBDB orbit families, ranked within
families by known M1 (smaller first), then diameter and designation. Missing M1 ranks
after known values. This is a reproducible starter selection, not an objective
"most interesting" score or prediction of observed peak brightness. Parent comets
are selected with fragments excluded. Edit the manifest to curate the result.
Borisov is resolved as primary designation 2019 Q4; 2I is retained as an editorial alias.

The catalog stores numeric query fields as JSON numbers, names/classifications as
strings, and unknowns as null. `diameter_sigma` remains source text because it can
be asymmetric (e.g. -1/+4). Original API responses are retained in cache envelopes.
Field units are recorded: km, km^3/s^2, g/cm^3, AU, degrees, JD TDB and hours. Rotation
periods are the source's synodic values; spin-pole text does not imply an orientation
provider. Do not treat SBDB TDB element epochs as MPC TT epochs. SBDB SPK IDs are
used as supplied; do not construct them from asteroid numbers or assume old offsets.

`enrich` retrieves full per-object SBDB responses with physical parameters, their
references/uncertainties and alternate designations. It retains those responses in
`sbdb_details`, preserving source distinctions and missing properties. The base
bulk catalog already includes available physical/orbital fields. An enrichment
summary records how many objects have details; a limited run is visibly partial.

Requests are sequential, with caching, timeout, transient-error retries and API
version checks (Query 1.0, SBDB 1.3, Horizons 1.2). Cache entries retain the request
URL, retrieval timestamp and original response SHA-256. The commands can be rerun
without repeated downloads; `--offline` requires cache hits and `--refresh` fetches
new snapshots. Global options precede the subcommand:

```sh
python3 Tools/SmallBodies/extract.py --offline select
python3 Tools/SmallBodies/extract.py --cache /tmp/jpl-cache select --asteroids 10 --trojans 0 --comets 10
```

Horizons JSON SPK content is decoded with strict base64/header validation. Its
returned ID is preserved independently of the requested ID/alias. Each completed
file has an entry in the kernel manifest containing requested TDB coverage, source
cache key and SHA-256. Actual descriptor coverage/origin/axes must still be inspected
by the library. Reruns preserve manifest progress; use a different output directory
for a different coverage interval. Failures are explicit rather than silently
omitting curated targets or substituting trajectories.

Horizons small-body SPKs typically use type 21, now supported natively. The Ceres
smoke export uses Sun-relative J2000 axes. An SSB result therefore also needs a
compatible Sun ephemeris link; heliocentric output needs an explicit origin-centred
frame. Loading these files does not automatically register bodies or frame bindings.
The library continues to enforce coverage, origin, axes and time-scale contracts.

Validation:

```sh
python3 -m unittest discover -s Tools/SmallBodies -p 'test_*.py'
swift test --filter Spice
```

The scripts were exercised against live JPL selection, two SBDB enrichments and a
Ceres SPK export, then replayed from cache offline. Python tests cover selection
union/diversity, unknown/asymmetric values, response schema failures, offline cache
behavior, identity checks, kernel decoding and manifest resumption/conflicts.

Sources: [SBDB Query API](https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html),
[filters](https://ssd-api.jpl.nasa.gov/doc/sbdb_filter.html),
[SBDB API](https://ssd-api.jpl.nasa.gov/doc/sbdb.html),
[Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html).

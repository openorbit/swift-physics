# Reproducibility and floating-point limits

The library supports repeatable, explicit model evaluation. It does not promise
bit-for-bit identical results across operating systems, CPU architectures, Swift
versions, optimization settings or library releases. Fix the configuration and
compare physical outputs using meaningful tolerances. Mathematical determinism,
floating-point repeatability and agreement with the real solar system are separate
properties.

## What is reproducible

For a fixed build, configured model, immutable assets and deterministic injected
dependencies, queries evaluate the same mathematical model at the requested epoch.
Built-in analytical and two-body queries have no evolving integration cursor.
Numerical propagators start from their configured initial state on each query and
use local integration storage. Ordered force composition and specified maneuver
schedules are part of that configuration. Callbacks must not read wall time,
random values, mutable global state or changing external data if repeatability is
required.

Exact-epoch batch queries and trajectory samples preserve their documented order
and duplicates. Ephemeris caching stores native states without interpolation or
time rounding; eviction and hit order change cost, not the model being evaluated,
assuming faithful cache storage and a deterministic provider. Neither a cache hit
nor repeated agreement certifies physical accuracy. Model fallback choice is also
an input: different assets, coverage or policies can select a different model at
the same epoch. Record the selected source and provenance when diagnosing a result.

Within one executable on an unchanged host, identical pure query paths are expected
to repeat their numerical results. There is no broader bitwise compatibility contract.
Do not base gameplay identity, persistent object IDs or multiplayer authority on
floating-point output hashes. Runtime UUID cache scopes, frame object identities,
Swift hash values and dictionary/set iteration order are not durable identifiers.
Use explicit ordered arrays for model composition and stable application IDs for
saved configuration.

## Sources of numerical variation

`Double` is used for physical coordinates and calculations; finite precision still
applies. Platform transcendental functions, compiler contraction/reassociation and
optimization can alter the last bits. Subtraction of large nearby coordinates can
lose precision already absent from the inputs, and ill-conditioned geometries can
amplify small changes. Iterative roots, nearly parabolic motion and short/long
transfer limits deserve case-specific checks.

Adaptive integration accepts or rejects steps using calculated local errors. Small
rounding changes can alter the step sequence, derivative counts and accumulated
error. Local-error tolerances do not certify global or physical trajectory accuracy.
Root and event searches can terminate at slightly different values within their
criteria. Near an event, coverage edge, fuel boundary, feasibility threshold or
work limit, small differences can also change a discrete status. Treat those
boundaries explicitly rather than assuming every platform chooses the same branch.

Planner ranking is only as stable as its costs. For example, the Lambert window
search sorts unique departure epochs and flight times, retains branch input order,
and breaks exactly equal ranking keys using candidate order after its documented
cost/time/departure keys. Near-equal costs are still distinct floating-point values;
rounding may change their ranking. Retain candidate sets, residuals and termination
status when comparing plans. A tie-breaking rule does not establish a unique physical
optimum or platform-independent winner.

`Epoch` and `Duration` retain whole seconds separately from fractional seconds.
Keep those components when transferring precise time; flattening a distant epoch
into one `Double` second count or Julian date loses fractional resolution. UTC label
conversion additionally depends on the selected leap-second table and its coverage;
TDB conversion depends on its model and version. A label without its time scale and
conversion inputs is insufficient to reproduce an instant.

`SimulationClock` accepts supplied durations and makes no wall-clock reads. Replaying
identical advance, seek, rate and pause commands with the same starting configuration
is repeatable under the fixed-build assumptions above. It anchors elapsed time to
limit drift, but floating-point duration accumulation and rate scaling do not promise
exact equality between arbitrary partitions of an elapsed interval. Rate/pause/seek
changes also reanchor the clock. For deterministic game replay, record the commands
and durations or drive explicit epochs from a fixed tick schedule; do not substitute
fresh wall-clock deltas. Fixed ticks do not establish cross-platform bitwise physics.

Rendering introduces a separate rounding boundary: subtract the moving origin in
Double, then scale and convert to Float. Re-rendering from native SI states avoids
feeding presentation quantization back into propagation. A render origin or scale
change alters Float values without changing the physical model.

## Record enough context

A reproducibility report or saved input description should capture:

- Library revision, Swift compiler, build/optimization settings, target architecture
  and operating system. Include application version and custom-model versions.
- Whole/fractional epoch components, declared time scales, leap-second table version,
  TDB model and coverage, and any Earth-orientation or rotation data.
- Stable body and frame descriptions, physical origins, parent/provider bindings,
  units, handedness, GM convention and numerical constant values with provenance.
- Provider/model versions, requested and effective accuracy, coefficient selection
  coverage and explicit solver tolerances, work limits and derivative policies.
- Initial translational/attitude states, masses, ordered force/torque models, guidance,
  maneuvers, event rules, waypoint constraints and planner candidate order/objective.
- Asset content digests and load order, SPICE bindings and kernel versions, DE file
  identity, and fallback policy. Paths alone do not identify immutable content.
- Query epochs or replay commands, results and diagnostics, selected model sources,
  event/termination statuses and comparison tolerances.

This is a capture policy, not a serialization schema or a guarantee that arbitrary
custom providers can be serialized. The <doc:SERIALIZATION> defines versioned saved
inputs. Published benchmark reports are measurements of their stated build and
workloads, not reproducibility promises for other hardware.

## Comparing runs

Use an absolute and relative numerical criterion, for example
`abs(actual - reference) <= absoluteTolerance + relativeTolerance * abs(reference)`.
Select tolerances for the units and physical purpose of each quantity. Relative
error alone is inappropriate near zero. For Cartesian states, compare position and
velocity in the same frame/origin and at the same atomic epoch, preferably using
vector norms. Compare rotations geometrically, allowing equivalent quaternion
signs and periodic angle representations rather than comparing raw components.

Check discrete identities, model selection and schema versions exactly. Compare
event times, patch residuals and numerical states using justified tolerances; near
thresholds, report the ambiguity instead of silently expanding tolerances until a
test passes. Cross-platform work counters and iteration counts are diagnostics,
not generally exact compatibility checks.

Distinguish regression against a fixed implementation, truncation against a fuller
model, convergence under tighter solver settings, and physical validation against
independent reference data. Agreement in one category does not prove the others.
The existing Swift Testing suite contains analytical/reference comparisons and
query-order checks, but no multi-platform bitwise certification. A platform change
should rerun those checks and inspect meaningful residuals, especially at model and
event boundaries.

See <doc:ACCURACY>, <doc:TimeGuide>,
<doc:THREAD_SAFETY>,
<doc:GAME_INTEGRATION#Cache-invalidation> and
<doc:GAME_PERFORMANCE> for related limitations.

Versioned saved states and reconstruction policy are documented in
<doc:SERIALIZATION>.

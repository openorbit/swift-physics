# Simulation saves and schema versions

`SimulationSnapshot` is a Codable, Sendable data envelope with format
`celmek.simulation` and integer `schemaVersion: 1`. It contains an optional clock
checkpoint, ordered Cartesian-state snapshots, ordered two-body input records and
ordered application-owned dependency descriptions. Use JSONEncoder/JSONDecoder or
another Codable format; the schema is independent of live object identity.

Version 1 defines SI position, velocity, acceleration and mass fields explicitly.
A missing acceleration remains unavailable. Dry and propellant mass must either both
be present or both absent. Body identifiers use either a signed NAIF code or a custom
namespace/name pair. Two-body records retain their initial state, central body, GM
value/convention/provenance and all solver options; restore never substitutes current
solver defaults.

Epochs store an atomic offset from JD 2451545.0 TT. Whole seconds are canonical decimal
strings; fractional seconds are finite values in [0,1). This avoids integer loss in
JSON tools that represent all numbers as Double. Keep the fraction as a binary64 value
or an exact round-tripping decimal; avoid flattening the epoch into a Julian-date
number. No UTC or TDB conversion takes place while restoring an epoch.

Clock checkpoints include current epoch, rate, pause state, optional inclusive
coverage, anchor epoch and elapsed real duration. Restoring recomputes the current
epoch from that anchor using the clock's arithmetic and checks exact consistency.
This preserves fixed-build continuation rather than introducing a new rounding
anchor. A checkpoint may therefore fail consistency checks on a different arithmetic
implementation; it is not a cross-platform bitwise replay guarantee. For a deliberate
new playback session, construct a new clock from the saved epoch/rate/pause/coverage
and explicitly accept reanchoring. A saved clock is not a log of previous commands.

## Saving and restoring

```swift
let record = try SavedTwoBodyInput(id: "ship-coast", frameID: "earth-inertial-v1",
                                  propagator: coast)
let state = try SavedCartesianState(id: "ship", frameID: "earth-inertial-v1",
                                    state: coast.state(at: clock.epoch))
let snapshot = try SimulationSnapshot(clock: clock, states: [state],
                                      twoBodyInputs: [record])
let encoder = JSONEncoder()
encoder.outputFormatting = [.prettyPrinted, .sortedKeys]
let data = try encoder.encode(snapshot)

let loaded = try JSONDecoder().decode(SimulationSnapshot.self, from: data)
// The application creates/validates its physical frame graph before binding IDs.
let restoredCoast = try loaded.twoBodyInputs[0].restore { id in
    guard id == "earth-inertial-v1" else { throw SaveError.missingFrame }
    return earthInertialFrame
}
let restoredClock = try loaded.clock?.restore()
```

`SaveError`, `coast`, `clock` and the frame binding in this snippet belong to the
application. An ID must identify the full axes/origin/provider configuration, not just
a display name. Resolve each ID to the same frame object within the reconstructed
graph. Restore performs no implicit frame conversion. A two-body frame must satisfy
the propagator's inertial contract; the caller must also bind the intended central
origin, units and axes. Frame identities, trees, closures, open readers and caches
are not serialized.

For spacecraft snapshots, supply `mass:` when constructing `SavedCartesianState`.
After loading, combine `restore(using:)` and `restoreMass()` into `SpacecraftState`.
This stores an instantaneous mass inventory; it does not recreate the force model,
fuel schedule or a numerical integrator's internal steps.

## External and custom inputs

`SimulationDependency` provides stable application IDs, kind/version strings and
string-valued parameters. Applications define and validate these records when
reconstructing frames, ephemerides, time models, assets, forces, attitude models,
maneuvers and planners. The library retains them as data and does not execute a kind
string, open paths or select a fallback during decoding. Unknown application kinds
or versions must fail in the application's reconstruction dispatcher.

Record explicit parameter units and ordered dependency IDs. For binary assets,
include content digests and ordered load manifests; paths alone are insufficient.
For model selection, include requested/effective accuracy, coverage, time model,
fallback order, constants and provenance. Extend an application-owned envelope with
typed records for numerical forces, torque/guidance callbacks, maneuvers, waypoint
plans and mission state. Do not claim a Cartesian snapshot alone reproduces those
models. Custom executable callbacks require application-defined reconstruction code.
See [reproducibility capture policy](REPRODUCIBILITY.md#record-enough-context).

This first library schema directly reconstructs Cartesian/mass states, clocks and
two-body inputs. It establishes the envelope/versioning policy for broader saves,
rather than silently attempting to serialize every provider protocol. Model-specific
saved configurations can be added in future schema versions or application envelopes.

## Validation and compatibility

The envelope decoder rejects unknown format/schema versions before decoding payloads.
It requires the state/input/dependency arrays, validates unique IDs within each array,
checks time/clock consistency and uses runtime constructors to validate vectors,
mass, GM, body identifiers and solver options. State and input ID namespaces are
separate. Frame IDs are external references, so the envelope does not require a frame
description in its dependency array. Actual restoration validates the supplied frame
and propagates binding/provider errors.

Standalone Codable record decoding creates data records; call their restore methods
before use. Envelope decoding validates its contained records without trusting saved
physical frame claims. Additional unknown JSON fields are ignored by Codable in
version 1; they must not be used for new required semantics. An incompatible field,
unit, interpretation or required model input requires a new schema version. Consumers
must not treat a dependency description as validated application configuration merely
because envelope decoding succeeded.

There is no automatic migration or best-effort loading of future versions. Implement
explicit, tested migrations into a supported schema, retaining the original save for
recovery. Keep application schema versions separate from the library envelope version.
Sorted JSON keys aid review but do not define canonical signing bytes. Round-tripped
physics inputs are independent of whitespace/key order; do not hash encoder output
as a portable model identity without a separately specified canonical format.

Bound input file size and application object counts before decoding untrusted saves;
the Codable decoder is not a resource-limited asset loader. Stage validation and all
bindings before replacing the live scene. Publish a new scene revision and invalidate
old derived buffers under the [cache policy](GAME_INTEGRATION.md#cache-invalidation).
Write saves atomically using the application's storage layer. Do not persist runtime
cache entries, UUID scopes or frame object addresses.

Swift Testing covers round trips, precise distant epochs, missing acceleration/mass,
clock continuation, invalid payloads, duplicate IDs, unknown versions, binding failures
and non-inertial two-body frame rejection. A hand-written
[version-1 fixture](Tests/CelMekTests/Fixtures/simulation-save-v1.json) pins the wire
format independently of encoder round trips. This is not a promise of bitwise numerical
compatibility across library/compiler/platform versions.

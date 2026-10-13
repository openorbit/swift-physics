# Common ephemeris interface

`EphemerisProvider` extends `StateProvider` for prescribed motion of one body relative
to one origin. Its `EphemerisMetadata` declares body and origin identifiers, native
frame identity, usable coverage in atomic `Epoch`s, coverage basis, model time scale,
derivative methods, and source/version provenance. Native results are SI
`CartesianState`s; velocity and acceleration use elapsed atomic seconds regardless
of the underlying model's time argument. Physical body properties remain separate.

`coverageBasis` distinguishes caller-selected intervals from model-documented
coverage. Neither tag alone is an accuracy guarantee. A frame must actually use the
declared origin and axes: assigning metadata does not convert coordinates or prove
that a caller's declaration is correct. Body centres and system barycentres remain
distinct identifiers. Origin/frame conversion requires the frame system or a model
adapter that genuinely performs the conversion.

## Adapting an existing state provider

```swift
import CelMek
import Math

let start = try Epoch(julianDate: 2451545, scale: .tt)
let end = try start.advanced(by: CelMek.Duration(seconds: 86400))
let metadata = try EphemerisMetadata(
  body: .earth, origin: .solarSystemBarycentre, frame: earth.frame,
  coverage: start...end, coverageBasis: .callerSelected,
  modelTimeScale: .tdb, velocityMethod: .analytic,
  accelerationMethod: .unavailable,
  provenance: PhysicalDataProvenance(
    source: "Built-in VSOP87E, complete coefficients; default two-term TDB adapter",
    version: "application snapshot 1"))
let ephemeris = try StateProviderEphemeris(provider: earth, metadata: metadata)
let state = try ephemeris.state(at: start)

let collection = try EphemerisCollection(providers: [ephemeris])
let earthState = try collection.state(of: .earth, at: start)
let tree = try CelestialFrames.makeTree()
let icrf = try ephemeris.state(at: start, in: CelestialFrames.icrf, using: tree)
```

This example adds a caller-selected interval to the existing VSOP87E state provider;
it does not establish new scientific coverage. Choose the intersection of all
model, coefficient and time-data coverage. A configured `VSOP87StateProvider` can
be wrapped the same way; record its selected TDB model in provenance and restrict
the interval accordingly. Dedicated model adapters are described below.

The adapter checks frame identity on construction and coverage before evaluation.
Every result must carry the requested epoch, the declared native frame, and the
advertised acceleration availability. Provider errors propagate. Missing acceleration
is `nil`; an advertised available acceleration must be present, and an unavailable
one must remain absent. Velocity cannot be marked unavailable because CartesianState
always contains it. Methods distinguish analytic derivatives, finite differences
with a positive interval, and supplied derivatives with a required description.
Metadata describes the derivative method; it does not calculate missing derivatives.

## Collections and custom implementations

`EphemerisCollection` is an immutable selection of one provider per body. It preserves
input order in `bodies` and rejects duplicate target IDs. Each body may have a different
origin, frame and coverage; inspect `metadata(for:)` before combining states. Unknown
bodies throw `unsupportedBody`, with no hidden provider fallback. `provider(for:)`
returns the selected provider for direct integration with existing state/frame APIs.
An empty collection is valid.

Implement `EphemerisProvider` directly for a new model. `frame` defaults to
`metadata.frame`; `state(at:)` must enforce that model's contract. The protocol's
`checkedState(at:)` checks custom implementations, and all collection queries use
that path. A custom implementation's direct `state(at:)` is its own responsibility.
Frame conversion through an existing StateProvider API preserves provider errors.

Collections do not implicitly cache or depend on query order. Optional caching is
described below. Collections do not make
mutable or file-backed providers thread-safe. Providers retain their own data and
concurrency contracts, and must avoid strongly retaining a frame tree that owns them.
The metadata interval includes queryable endpoints; numerical derivative stencils
must be accounted for when selecting that interval.

Swift Testing covers coverage endpoints and rejection before evaluation, custom
provider validation, absent bodies, duplicate IDs, epoch/frame and acceleration
mismatches, derivative metadata validation, and VSOP87 integration through ICRF.

## Dedicated analytic adapters

`VSOP87Ephemeris` evaluates the complete built-in VSOP87E series with analytic
velocity, unavailable acceleration, caller-selected coverage, and a required
`TDBTimeModel`. Its coverage must lie inside the selected time model's coverage.
The native frame remains `CelestialFrames.vsop87E`, with Solar System barycentric
origin. Earth denotes the geocentre; the outer planet solutions are exposed as
system barycentres rather than claiming separately resolved planet-centre motion.
Mercury, Venus and the Sun retain their body-centre identifiers. The model does not
resolve planetary satellites or provide barycentre-to-planet offsets.

`SphericalAnalyticEphemeris` adapts the two existing spherical approximations:

| Model | Target / origin | Native axes | Position units before conversion |
| --- | --- | --- | --- |
| `meeusMoon` | Moon / Earth | Mean ecliptic/equinox of date | Radians and metres |
| `goffinPluto` | Pluto system barycentre / Sun | Ecliptic/equinox J2000 | Radians and AU |

The lunar implementation named `elp2000_82b` is the truncated Meeus chapter 47
series, not the full ELP coefficient model. Its longitude retains the existing
constant light-time term; the adapter does not claim a fully geometric or apparent
observational ephemeris. The Pluto series does not resolve Pluto and Charon separately.
The adapters use TT as legacy ephemeris time, without claiming that TT and TDB are
identical. Provenance records the chosen model and time convention.

Cartesian velocity uses a five-point derivative on atomic seconds. Supply a positive
`velocityStep` explicitly. `positionCoverage` is caller-selected; returned metadata
coverage removes two steps at each end. Too short an interval is rejected. No velocity
is obtained by differentiating wrapped longitude directly. Acceleration remains nil.
The legacy functions still return their original spherical units and conventions.
The internal evaluators now accept centuries from a split TT date to avoid first
collapsing epochs into a large scalar Julian date when estimating derivatives.

### Registering the native frames

```swift
let tdb = try FairheadTDBModel()
let earthModel = try VSOP87Ephemeris(body: earth, coverage: coverage, tdbModel: tdb)
let sunModel = try VSOP87Ephemeris(body: sun, coverage: coverage, tdbModel: tdb)
let tree = try CelestialFrames.makeTree()
let earthICRF = try FrameConvertedStateProvider(earthModel, in: CelestialFrames.icrf, using: tree)
let sunICRF = try FrameConvertedStateProvider(sunModel, in: CelestialFrames.icrf, using: tree)
let frames = try tree.makeAnalyticEphemerisFrames(
  earth: earthICRF, sun: sunICRF,
  accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 60)))
let moonModel = try SphericalAnalyticEphemeris(
  model: .meeusMoon, frames: frames, positionCoverage: coverage,
  velocityStep: CelMek.Duration(seconds: 30))
let collection = try EphemerisCollection(providers: [earthModel, sunModel, moonModel])
let moonNative = try collection.state(of: .moon, at: epoch)
let moonICRF = try moonModel.state(at: epoch, in: CelestialFrames.icrf, using: tree)
```

`coverage` and `epoch` are caller-selected; choose an interior epoch with all required
neighbouring samples available. Retain the tree for the converted centre providers.
Earth and Sun inputs to frame construction must use ICRF identities. The geocentric
frame uses Newtonian Earth translation and an EME2000 mean-of-date edge; the solar
frame uses Newtonian Sun translation and fixed ecliptic J2000 axes.

The lunar axes use IAU 1976 precession and IAU 1980 mean obliquity, consistent with
the legacy mean-date convention. The Swift precession rotation is checked against
[ERFA's published pmat76 reference](https://github.com/liberfa/erfa/blob/v2.0.1/src/t_erfa_c.c),
with no C runtime dependency. Its angular derivatives use a five-point stencil with
a caller-configurable 60-second default step. No nutation is applied. Using a newer
precession model without refitting the legacy series would change this convention.

Frame conversion adds origin motion and rotating-axis velocity terms; it does not
relabel native states as ICRF. Coverage for converted queries is the intersection of
model, TDB and centre-provider coverage, further restricted by acceleration estimates
and frame derivative stencils. The native ephemeris metadata does not expand that
intersection. Physical accuracy versus modern DE ephemerides remains to be measured;
complete evaluation of these truncated legacy models does not make them full ELP or
modern Pluto ephemerides.

Tests preserve the legacy published lunar/Pluto cases, verify spherical-to-Cartesian
units, compare derivative steps and differentiated ICRF trajectories, check precession
against ERFA, and validate native coverage endpoints and collection integration.

## Explicit origin conversion and barycentres

`makeOriginConvertedEphemeris` subtracts a new-origin state from a target state.
Both inputs must declare the same original origin and have the exact same native
frame identity. The method registers a translated child frame with unchanged axes,
so its new origin is represented by a physical frame edge rather than a metadata
relabel. Position, velocity and available acceleration are subtracted in SI units.
The usable coverage is the closed intersection of the input intervals; disjoint
coverage fails before registering a frame.

```swift
// moonBarycentric and earthBarycentric already use the same barycentric frame.
let lunarGeocentric = try tree.makeOriginConvertedEphemeris(
  target: moonBarycentric, newOrigin: earthBarycentric,
  frameName: "Earth-centred lunar ephemeris",
  accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 60)))
let relative = try lunarGeocentric.state(at: epoch)
```

This operation requires actual states of both objects. It cannot derive Earth from
an Earth–Moon barycentre alone or resolve a planet from its system barycentre without
the missing relative motion and mass data. Transform differently oriented inputs
into a common frame first, then wrap the converted StateProviders with metadata
that accurately describes that frame's origin, coverage and derivative methods.
When translation has already been performed by a frame transform, update the declared
origin too; do not subtract the same origin a second time.

The native subtraction query requires no frame-tree lookup. Its acceleration is
available only when both inputs supply acceleration; neither zero nor estimated
acceleration is selected implicitly. The registered edge's `OriginAccelerationPolicy`
governs later frame transformations. With `requireProvided`, transformations throw
if new-origin acceleration is missing. An explicitly estimated or zero acceleration
for that edge does not manufacture acceleration in the recentered ephemeris.
Neighbour coverage for estimation further restricts converted queries, beyond the
native metadata interval. The edge checks custom ephemeris results before use.

### Constructing a barycentre

`BarycentreEphemeris` computes the constant-mass weighted combination of at least
two distinct constituent bodies. All constituents must share the same declared
origin and native frame. Positive finite relative masses and their provenance are
required; they are normalized internally, so a common mass scale is irrelevant.
The output body identity must differ from every constituent. Missing constituents
cannot be inferred, and duplicate identities are rejected.

```swift
let ratios = try PhysicalDataProvenance(source: "Selected Earth/Moon mass ratio", version: "application data 1")
let emb = try BarycentreEphemeris(
  bodyIdentifier: .earthMoonBarycentre,
  components: [
    BarycentreComponent(ephemeris: earthBarycentric, relativeMass: earthToMoonMassRatio, provenance: ratios),
    BarycentreComponent(ephemeris: moonBarycentric, relativeMass: 1, provenance: ratios)
  ],
  provenance: PhysicalDataProvenance(source: "Computed Earth–Moon barycentre", version: "1"))
```

Use a ratio compatible with the actual ephemeris solution; no ratio is inferred from
the bundled body catalogue. GM ratios can represent mass ratios only when their
coordinate conventions and source solutions are compatible. Ratios are constant:
time-varying masses would require extra velocity/acceleration terms and a different
provider. The combined field records component ephemeris and mass-ratio provenance.

Coverage is the intersection of every component's interval. Acceleration remains
unknown if any component lacks it. Analytic derivative metadata is retained only
when every input derivative is analytic; other available combinations are labelled
as supplied compositions. Floating-point cancellation still limits small relative
states computed from large barycentric coordinates, and extreme mass ratios may
underflow normalized weights.

Mixed model time arguments are permitted because all input states have already
converted derivatives into atomic seconds. `modelTimeScale` retains the target's
(or first component's) primary argument for compatibility; `modelTimeScales` records
the complete set of arguments involved. No time conversion occurs during subtraction
or weighted composition. Constituent providers remain responsible for their own
TT/TDB conversions and coverage. Physical consistency across different solutions
is a caller decision, not guaranteed by matching identifiers.

Both providers are immutable and introduce no query-order cache or hidden fallback.
Registered origin edges retain providers without retaining their tree; any existing
FrameConvertedStateProviders still hold their tree weakly. Retain the tree for
later conversion. Tests check all derivatives, frame round trips, mixed time metadata,
coverage and identity rejection, missing acceleration, and scale-independent weights.

Analytic adapters accept an explicit shared accuracy level; supported levels and their meaning are documented in [ACCURACY.md](ACCURACY.md).

## Batch queries and explicit caches

Ephemeris providers expose `states(at: [Epoch])`. Collections expose
`states(of:at:)` for one body and `states(for: [EphemerisRequest])` for mixed bodies.
Results preserve the input order and duplicate requests, with exact epochs. Empty
provider/request batches return an empty array; a one-body collection query still
validates that body's selection even when its epoch list is empty. Batches are
sequential and fail at the first error, without returning a partial array. Earlier
successful calls may have populated an explicitly selected cache.

```swift
let cache = try MemoryEphemerisCache(capacity: 512)
let cached = try CachedEphemerisProvider(provider: earthModel, cache: cache)
let samples = try cached.states(at: epochs)
let collection = try EphemerisCollection(providers: [cached, moonModel])
let states = try collection.states(for: [
  EphemerisRequest(body: .moon, epoch: epoch),
  EphemerisRequest(body: .earth, epoch: epoch)
])
// After replacing data dependencies while preserving metadata:
try cached.invalidate()
// Release all shared storage when changing scenes:
cache.removeAll()
```

`CachedEphemerisProvider` memoizes native Cartesian states at exact atomic epochs.
It does not round time, interpolate, sample neighbours, cache transforms, or change
accuracy settings. Construct a new wrapper for new metadata, body, origin, frame,
coverage, time model or accuracy configuration. Each wrapper receives its own opaque
scope; wrappers sharing storage cannot collide even when their bodies and epochs
match. Metadata, including requested/effective accuracy, is retained unchanged.
Scopes are runtime identities, not persistent model/version keys.

Queries check coverage before accessing storage. Both cache hits and computed
results must match the requested epoch, declared frame and acceleration availability.
Provider and custom-cache errors propagate. Failed evaluations are never stored.
The cache trusts the actual numerical values supplied by its storage implementation;
it cannot independently certify an externally injected state's physical correctness.

`MemoryEphemerisCache` requires a positive capacity and uses least-recently-used
eviction with O(1) lookup/update. Both lookups and overwrites refresh recency. Capacity
bounds total entries across all scopes sharing that cache. `count` reports current
entries. Storage mutations use a lock, and list links avoid ownership cycles. Scope
invalidation removes only that wrapper's entries; clearing storage removes all
scopes. Entries retain their states/frame identities until invalidation, eviction or
cache release.

Invalidation rotates the wrapper's generation before clearing its old storage.
Queries already in flight may complete with their original results but cannot store
them in the old generation after invalidation. Later queries use the new generation.
External calls to the shared memory cache's `removeAll()` clear storage but do not
rotate wrapper generations; use `cached.invalidate()` when changing model data.
If custom-cache removal throws, the wrapper still uses its new generation and the
error propagates, while stale storage remains the custom cache's responsibility.

Cache storage locking does not serialize native provider evaluations or make custom
providers thread-safe. Concurrent misses may evaluate the same epoch more than once.
The caller must respect the underlying provider's concurrency contract. Cache
protocol callbacks must not reenter their invoking wrapper. No unchecked Sendable
conformance or concurrency guarantee for frame identities is introduced.

Custom storage implements `EphemerisStateCache` using the full `EphemerisCacheKey`
(scope plus exact epoch). It must preserve native states faithfully and handle its
own capacity, persistence and synchronization. The built-in cache is intentionally
in-memory; file coefficient/record caches remain independent model concerns.

Query-order independence assumes the configured ephemeris is deterministic at a
given epoch. Hits and eviction can change cost, but not numerical results. Mutable
model/data dependencies require explicit invalidation; no data changes are detected
automatically. Tests check ordered and reversed batches, duplicate and subnanosecond
keys, LRU eviction, scope isolation, generation changes during evaluation, failures,
invalid cache hits, and mixed-body fail-fast queries.

## Native JPL DE441

`DE441StateProvider` adapts a URL-opened `DE441Ephemeris` to the shared SI
ephemeris contract with explicit target/centre, coverage, frame and TDB model.
It evaluates complete file coefficients for all accuracy hints and supplies
analytic velocity converted to atomic seconds. The reader has a configurable
bounded record cache and internally serialized random-access IO. See
[DE441.md](DE441.md) for supported binary profiles, split-date interpolation,
time-model limits and optional local-data tests.

## Explicit fallback tiers

`FallbackEphemeris(de441:model:kepler:)` selects **DE441 → analytic/fitted model
→ Kepler** at each requested epoch. Each tier is optional. Omit DE441 when its
optional asset is absent; an invalid asset is an error, not permission to hide a
failure. Selection uses declared coverage, including endpoints. Errors from a
selected source (including narrower dependency coverage, corrupt coefficients,
and state-contract violations) propagate. Accuracy hints control model evaluation,
not tier ordering.

The SPK overload requires `fileBackedOrder: .de441ThenSPK` or `.spkThenDE441`:

```swift
let ephemeris = try FallbackEphemeris(
    de441: optionalDEProvider,
    spk: optionalSPKProvider,
    fileBackedOrder: .de441ThenSPK,
    spkFallbackPolicy: .missingLinksAndGaps,
    model: optionalModelProvider,
    kepler: optionalKeplerProvider)
let result = try ephemeris.selection(at: epoch)
// result.tier distinguishes .de441, .spk, .model and .kepler.
// result.source identifies the selected provider; skippedSources records reasons.
```

Both file sources precede the analytic/fitted and Kepler tiers. The new overload
permits only `SpiceSPKResolutionError.missingLink` and `.coverageGap` from the SPK
slot to fall through, and records the missing body ID. Use `.declaredCoverageOnly`
to propagate those failures too. No other error permits retry, including unsupported
SPK types, malformed coefficients, cycles, work limits, missing frame bindings,
time-model failures or returned state-contract violations. If no fallback succeeds,
the last permitted SPK availability error is rethrown. Generic `outsideCoverage`
errors from an otherwise eligible provider are not interpreted as SPK gaps.

Assets remain caller-owned and optional: omit absent providers; loading or constructing
an invalid provider must fail before composing the chain. SPK construction still validates
coverage endpoints, so a wholly unavailable dataset is not an implicit optional provider.
Interior SPK gaps can be handled by the composed chain. Every supplied provider must
already share exact target, origin and frame identities. Retained SPK snapshots preserve
selection after pool unload/reload. The existing acceleration normalization and coverage
union checks still apply; source metadata records provider-level provenance, and
`SpiceSPKStateProvider.evaluation` remains available for per-segment chain details.

`DefaultEphemerisModel.makeProvider(for:coverage:tdbModel:frames:velocityStep:accuracy:)`
constructs the native-frame default where one exists:

| Body identity | Built-in model |
| --- | --- |
| Sun, Mercury, Venus, Earth | VSOP87E |
| Mars, Jupiter, Saturn, Uranus, Neptune system barycentres | VSOP87E |
| Moon | Existing legacy ELP82-labelled truncated Meeus lunar series |
| Pluto system barycentre | Goffin |
| Other identities | None; use a supplied fitted model or Kepler |

The current lunar implementation is **not full ELP82** and retains the honest
`meeusMoon` name. The factory returns nil for bodies without a default, rather
than substituting a different body identity. A caller-supplied fitted provider can
occupy the model tier for any body.

Compose providers only after converting them to the same origin and frame using
the frame system and explicit metadata adapters. DE441 is usually barycentric
ICRF; VSOP87E needs an axes conversion, the Moon needs Earth translation, and
Goffin needs Sun translation. Retain the tree used by converted providers, and
restrict coverage to all transform dependencies. The chain rejects different
body identities, origins, frames, and gaps in its combined coverage.

`KeplerEphemeris(body:propagator:coverage:provenance:)` wraps a `TwoBodyPropagator`
with an explicit fixed reference state and central GM. The seed state must be
relative to the central body in a suitable inertial frame. Convert this provider's
origin/frame before combining it with barycentric providers. Supply a sensible
usable coverage interval; there is no automatic orbit inferred from a body ID or
seeding from the last query. Unknown bodies can therefore use a Kepler-only chain.

```swift
let ephemeris = try FallbackEphemeris(
    de441: optionalDEProvider, // already converted to common coordinates
    model: analyticOrFittedProvider,
    kepler: optionalKeplerProvider)
let answer = try ephemeris.selection(at: epoch)
// answer.state is SI; answer.tier and answer.source report the selected model.
```

`selection(at:)` reports the selected source metadata, including its provenance
and accuracy. Chain metadata lists all participating time scales. Acceleration is
returned only if every tier supplies it, so its availability stays consistent
when switching sources. Selection is deterministic and independent of query
order. Switching models can introduce a jump; this policy does not blend states
or assert that the tiers agree at their boundaries.


## Native SPK geometric providers

`SpiceSPKStateProvider` adapts a retained `SpiceKernelSnapshot` to the common SI
interface using exact NAIF target/origin identifiers, an explicit `TDBTimeModel`
and provenance-backed `SpiceFrameBindings`. Its frame tree converts each relative
link's axes and rotation rates; origin translation is already resolved by the SPK
chain. It evaluates full type 2/3 coefficients, returns atomic-second velocity and
leaves acceleration unavailable. Non-SSB origins require an explicitly centred
output frame. Endpoint validation does not certify interior chain coverage; gaps
and corrupt/unsupported selected segments propagate without automatic fallback.
See [SPICE_KERNELS.md](SPICE_KERNELS.md) for construction, dataset-specific J2000/ICRF
interpretation, retained snapshot ownership, frame-tree lifetime and provenance.

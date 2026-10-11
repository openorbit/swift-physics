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

There is no shared cache or dependence on query order. Collections do not make
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

# Orbital mechanics

## Kepler equations

`KeplerSolver.elliptic(meanAnomaly:eccentricity:options:)` solves
`M = E - e sin(E)` for `0 <= e < 1`. Anomalies are radians. The result
preserves the input mean anomaly's revolution offset, including negative
anomalies. At very large magnitudes, a `Double` cannot preserve small orbital
phase changes; use a separately tracked revolution count when that matters.

`KeplerSolver.hyperbolic(meanAnomaly:eccentricity:options:)` solves
`M = e sinh(H) - H` for `e > 1`. Hyperbolic mean anomaly is dimensionless and
is not reduced modulo a revolution. Eccentricity exactly one is rejected:
parabolic motion requires a different equation and time parameterization.

```swift
let options = try KeplerSolverOptions(
    absoluteTolerance: 1e-13,
    relativeTolerance: 1e-14,
    maximumIterations: 128)
let solution = try KeplerSolver.elliptic(
    meanAnomaly: 0.1, eccentricity: 0.9, options: options)
let eccentricAnomaly = solution.anomaly
```

Both solvers use a monotonic root bracket and safeguarded Newton steps, with
bisection when Newton would make insufficient bracket progress. Small-anomaly
series avoid subtracting nearly equal values near eccentricity one. This
approach follows the combination of Newton iteration and bisection described
by [JPL's KPSOLV documentation](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/kpsolv.html);
this library implements the scalar equations directly in Swift.

Options specify absolute and relative **anomaly** tolerances, rather than a
residual tolerance that can hide large root errors near a flat derivative.
A bracket half-width within `absoluteTolerance + relativeTolerance * |anomaly|`,
an exactly zero computed residual, or adjacent representable bracket values
ends iteration. Floating-point rounding limits attainable accuracy; these
settings are not formal interval-arithmetic error guarantees.
`KeplerSolution.iterations` reports the iteration count (zero for an exact
endpoint or the circular special case).

`KeplerSolverError` distinguishes invalid eccentricity, nonfinite mean anomaly,
invalid options, arithmetic range limits, and exhausted iterations. Exhaustion
includes the final positive-anomaly bracket and iteration count; for negative
mean anomalies, negate and reverse the bracket to interpret the signed root.
Extreme hyperbolic inputs may overflow intermediate values even if the
mathematical root is finite; the solver reports `numericalRangeExceeded`.

Tests include circular and near-parabolic eccentricities, signed anomalies,
multiple elliptic revolutions, invalid inputs, and iteration exhaustion.
Independent reference cases use 80-digit Decimal arithmetic and bisection;
regenerate them with `python3 Tools/generate_kepler_fixtures.py`.

## Cartesian states and osculating elements

`OsculatingOrbitalElements` carries epoch, frame, central-body identity and the
full sourced `GravitationalParameter`, including its coordinate-time convention.
It represents the instantaneous two-body conic associated with a position and
velocity, rather than the legacy Meeus `PlanetaryOrbitalElements` mean series.

Supply a state already relative to the central body, in right-handed,
pseudo-inertial axes centred on that body. The caller must choose a GM
compatible with the state's coordinate-time convention. Frame metadata cannot
verify the body's identity or time scaling; the converter checks the
`pseudoInertial` flag and does not infer an origin or transform rotating states.
Use the frame system to convert states before extracting elements.

```swift
let elements = try OsculatingOrbitalElements(
    state: centralBodyRelativeState,
    centralBody: .earth,
    gravitationalParameter: earthGM)
let reconstructed = try elements.cartesianState()
```

The six geometric fields are semi-latus rectum `p` in metres, eccentricity `e`,
inclination `i`, ascending-node longitude `Omega`, periapsis argument `omega`,
and true anomaly `nu`. Angles are radians. Node/periapsis angles are normalized
to `[0, 2pi)` and true anomaly to `[-pi, pi]`. `p = h² / mu` avoids dividing
by `1 - e²` during conversion and supports all three conics. The classical
state-vector geometry is described in
[NASA/TM-2004-213230, Appendix A](https://ntrs.nasa.gov/api/citations/20040084254/downloads/20040084254.pdf).
`semiMajorAxisMetres` is derived, negative for hyperbolas, and nil for an exact
parabola; near eccentricity one its value is ill-conditioned and may overflow.
Use `p` and `periapsisDistanceMetres` for more stable size information.

The direct element initializer supports elliptical, exactly parabolic, and
hyperbolic cases. It rejects anomalies with `1 + e cos(nu) <= 0`, outside the
finite physical branch. Conversion back uses perifocal position/velocity and
an explicit right-handed rotation. It preserves epoch and frame and supplies
no acceleration: osculating elements alone do not assert the actual forces.

`angleConvention` makes singularities explicit:

- Classical: all three orientation angles have their usual meanings.
- Circular inclined: periapsis argument is zero; true anomaly holds argument of latitude.
- Eccentric equatorial: node is zero; periapsis is measured from +X along orbital motion.
- Circular equatorial: node and periapsis are zero; true anomaly holds longitude along orbital motion.

For equatorial retrograde motion, angles along orbital motion run clockwise
when viewed from +Z. The inclination remains pi, retaining handedness.
The direct initializer canonicalizes exact singularities while preserving
combined phase. State extraction defaults to circular eccentricity and
absolute sine-inclination thresholds of `1e-12`. Values inside these thresholds
are snapped to circular/equatorial geometry; reconstruction is therefore
slightly approximate. `OrbitalElementConversionOptions` exposes both thresholds,
including zero to disable tolerance-based snapping. This does not recover
angles that floating-point precision cannot resolve.

`OrbitalElementsError` reports invalid elements/options, noninertial axes,
zero position, radial/zero-velocity trajectories, arithmetic range limits,
and unreachable conic anomalies. Radial trajectories do not define an orbital
plane and require a different representation. Very small nonzero angular
momentum and near-asymptote states can remain poorly conditioned; a successful
conversion does not guarantee well-determined classical angles. Conic type
uses the computed eccentricity without snapping to one, so roundoff can place
an extracted exact parabolic state slightly on either side of unity.

Tests include an independently derived polar state, prograde/retrograde
circular Cartesian fixtures, and a reconstruction grid across inclined,
equatorial, circular, nearly parabolic, parabolic and hyperbolic states.

## Dynamical anomalies and conic boundaries

`ConicAnomaly` distinguishes `.eccentric(E)`, `.hyperbolic(H)` and `.barker(D)`.
Convert from a true anomaly with `ConicAnomaly.from(trueAnomaly:eccentricity:)`
or obtain it from `elements.conicAnomaly()`. The tagged value's
`trueAnomaly(eccentricity:)` and `meanAnomaly(eccentricity:)` methods reject a
conic/eccentricity mismatch and nonfinite parameters.

`ConicAnomaly.from(meanAnomaly:eccentricity:options:)` uses the bounded Kepler
solvers for elliptic and hyperbolic cases. At exactly `e == 1`, its input is
Barker's `B = D + D³/3`, where `D = tan(nu/2)`, and it uses a closed-form
inverse. These conventions match the parabolic time-of-flight equation in
[JPL's PARTOF documentation](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/partof.html).
Barker B is dimensionless; it is not interchangeable with either conic's M.
Solver options apply to E/H iteration; the closed-form Barker inverse does not
iterate. Near zero it uses the linear limit to retain subnormal inputs.

True-anomaly conversion yields a principal phase. Elliptic anomalies supplied
through the mean-anomaly solver retain revolution offsets, and converting
those E values back to M retains those offsets. Converting to true anomaly
necessarily removes them. For circular states, E/M follow the canonical phase
convention described above rather than an undefined physical periapsis.

Mean anomaly evaluation uses cancellation-resistant series near E/H zero.
True-anomaly conversions use half angles and separate square roots of `1-e`
and `1+e`, without snapping a near-parabolic eccentricity to one. As e approaches
one, the mapping from E/H to true anomaly becomes very sensitive: a tolerance
on E/H does not guarantee the same tolerance on true anomaly. Choose tighter
solver options when resolving this phase is necessary. Future propagation
will use a common formulation instead of switching between ill-conditioned
classical time parameters near unity.

The finite parabolic branch requires `|nu| < pi`; a hyperbola requires
`|nu| < acos(-1/e)`. Boundaries and angles beyond them are rejected explicitly.
For extremely large D/H, converting to a `Double` true anomaly can round to
an asymptote and is also rejected. Near the boundary, state reconstruction
uses half-angle expressions for radius and transverse velocity to retain
small terms that `1 + e*cos(nu)` or `e + cos(nu)` can lose. These formulas
improve numerical conditioning but do not recover precision lost in an input
angle already rounded to the boundary.

## Two-body propagation

`Propagator` extends `StateProvider` with `referenceEpoch` and `centralBody`.
`TwoBodyPropagator` takes a `CartesianState`, explicit central-body identity,
sourced GM and validated `TwoBodyPropagationOptions`. It solves Newtonian
point-mass motion directly in Cartesian coordinates using scaled universal
variables, Stumpff functions and Lagrange f/g coefficients. This common
formulation supports circular, elliptic, parabolic and hyperbolic motion
without a conic-specific switch near eccentricity one. See the universal
variable approach described by
[JPL's PROP2B documentation](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/prop2b.html)
and [Wisdom and Hernandez (2015)](https://arxiv.org/abs/1508.02699).
The implementation here uses its own Swift Stumpff evaluation and a
safeguarded Newton/bisection solver, rather than the latter paper's algorithm.

```swift
let orbit = try TwoBodyPropagator(
    initialState: centralBodyRelativeState,
    centralBody: .earth, gravitationalParameter: earthNewtonianGM)
let future = try initialState.epoch.advanced(by: Duration(seconds: 600))
let predicted = try orbit.state(at: future)
let diagnostics = try orbit.propagation(at: future)
```

The initial position and velocity must already be relative to the central body
in a pseudo-inertial frame centred on it. They define a massless test-body
trajectory under `acceleration = -mu*r/|r|³`; incoming acceleration is not an
initial condition. Every result supplies the model acceleration, including at
the reference epoch. The central body's motion belongs in the frame tree.
No planetary perturbations, surface impacts or finite thrust are included.

Elapsed time and returned derivatives use continuous atomic/TT seconds.
GM must be `.newtonian` or `.ttCompatible`. A `.tdbCompatible` or
`.tcbCompatible` catalogue entry is rejected, because using it with atomic
elapsed time would silently mix coordinate-time conventions. A game may
explicitly define its own Newtonian GM with provenance; the propagator does
not relabel a catalogue value or infer a relativistic scale conversion.

Every query starts from the immutable initial state, so pause, time warp,
reverse time, and query order have no persistent effects. Elliptic queries
reduce elapsed time by the model period before solving, avoiding enormous
trigonometric arguments after many revolutions; the result retains the
requested absolute epoch. Converting a very large split duration to the
solver's `Double` time argument can still lose fractional phase precision.
No coverage interval is imposed, but arithmetic or convergence can fail for
extreme queries. An ephemeris adapter can impose an application coverage policy.

Options expose absolute time-equation residual tolerance (seconds), relative
residual tolerance, maximum root iterations and maximum bracket expansions.
These are numerical solver controls, not physical error guarantees. Near-zero
Stumpff functions use series; other arguments use trigonometric/hyperbolic
forms with a half-angle C expression. The root is bracketed with bounded
expansion and Newton steps are accepted only when they contract the bracket.
Iteration ends when the residual is within tolerance or the bracket reaches
adjacent representable numbers. `propagation(at:)` returns the state and
iteration/expansion counts.

`PropagationError` reports invalid options, incompatible GM convention,
noninertial frame, zero position, radial trajectory, arithmetic range limits,
bracketing failure and iteration exhaustion. Radial motion is excluded to
avoid propagation through a point-mass collision. Very small angular momentum,
extreme energy and long unbound flight can be ill-conditioned even before an
arithmetic limit is reached. Physical-model fidelity depends on the chosen GM
and the validity of the two-body approximation independently of solver tolerance.

## Impulsive maneuvers

`ImpulsiveManeuver` specifies epoch, frame identity and finite delta-v in m/s.
`ManeuveredTwoBodyPropagator` builds immutable coast segments from a base
`TwoBodyPropagator`, sorting burns by epoch. Burns before the reference epoch
and duplicate burn epochs are rejected. Combine simultaneous delta-v vectors
in one explicit frame before construction if that is the intended operation.

```swift
let burn = try ImpulsiveManeuver(
    epoch: future, frame: orbit.frame,
    deltaVelocityMetresPerSecond: SIMD3(0, 20, 0))
let trajectory = try ManeuveredTwoBodyPropagator(base: orbit, maneuvers: [burn])
let before = try trajectory.state(at: future, side: .before)
let after = try trajectory.state(at: future) // Post-burn at the exact epoch.
```

With a different burn frame, pass `using: frameTree`. Delta-v is a free vector:
only the rotation from burn axes to propagation axes is applied. Frame origin
translations, origin velocities and angular-rate terms do not apply to an
instantaneous velocity *difference at the same position*. Position is
continuous, velocity changes by the resolved delta-v, and acceleration remains
the central-gravity value on each side. Impulses are prescribed changes;
spacecraft mass, thrust duration and propellant are not inferred.

Frame transforms are evaluated at burn epochs during construction, and
`resolvedManeuvers` exposes the rotated vectors plus pre/post-burn states.
The frame tree is not retained. If source frame providers or the burn plan
change, construct a new trajectory. Any spacecraft-relative frame supplied by
a caller must describe the intended pre-burn trajectory, including earlier
burns; construction does not retarget arbitrary frame providers automatically.
At a burn epoch the trajectory is right-continuous by default; `.before` selects
the preceding segment. Queries earlier than the initial epoch follow the
unmodified base orbit. Queries in any order select the same fixed segment.

## Sampling and the game integration checkpoint

`TrajectorySampler.sample(_:at:in:using:maximumSamples:)` samples any
`StateProvider`, preserving epoch order and duplicates. A uniform overload
accepts `from:through:step:`: a positive step works in either time direction,
both endpoints are included, and the final interval may be shorter. Fractional
steps are anchored to the start epoch to avoid accumulated schedule drift.
Schedules are bounded (default 100,000 samples) and validated before querying.
A step below the epoch representation's resolution is rejected if it cannot
advance the schedule. Provider/frame errors propagate without a partial result.

```swift
let points = try TrajectorySampler.sample(
    trajectory, from: initialState.epoch, through: future,
    step: Duration(seconds: 10),
    in: bodyFixedFrame, using: frameTree)
```

Each result carries its exact requested epoch, SI state and destination frame.
Provider epoch/frame contracts are checked before conversion. Sampling does not
interpolate, integrate incrementally, or round epochs. To display velocities
on each side of an impulse, query its two sides explicitly; a uniform schedule
does not automatically insert burn epochs. Rendering interpolation and `Float`
conversion belong at the application boundary.

Run the complete checkpoint with:

```sh
swift run --disable-sandbox TwoBodyGameDemo
```

[Examples/TwoBodyGameDemo/main.swift](Examples/TwoBodyGameDemo/main.swift) queries
VSOP87 Earth, creates an Earth-centred inertial frame, propagates a spacecraft,
applies a scheduled impulse, and emits JSON samples in a prescribed rotating
Earth frame. The spin is a game model rather than observed Earth orientation.
The explicit `.assumeZero` origin-acceleration policy is used because VSOP87's
state provider supplies velocity but no acceleration. This affects transformed
accelerations, not the propagated Earth-relative coast solution.

Validation includes independent 80-digit Decimal conic-reference fixtures
(regenerate with `python3 Tools/generate_two_body_fixtures.py`), a non-periapsis
analytic ellipse, two-body energy/angular-momentum conservation, forward/backward
recovery, multiple revolutions, maneuver discontinuities and rotating-frame
burn vectors. The end-to-end test queries an ephemeris-driven Earth alongside
a spacecraft under pause/time warp and verifies rotating-frame conversion.

## Hohmann transfers

`HohmannTransfer` calculates an impulsive Newtonian transfer between two coplanar
circular orbits around the same central body. Supply centre-relative radii in
metres and a sourced Newtonian or TT-compatible GM. Surface altitudes must be
converted to radii by the caller. The result exposes transfer semi-major axis,
eccentricity, circular endpoint speeds, signed departure/arrival delta-v, total
absolute delta-v and the half-ellipse flight duration in atomic seconds.

Delta-v is signed along **local orbital motion**: both burns are positive for an
outward transfer and negative for an inward transfer. At arrival the spacecraft
is on the opposite side of the central body, so its tangential axes are reversed
relative to departure. Equal radii yield zero burns and a half-period circular
coast, rather than a zero-duration plan.

```swift
import CelMek
import Math

let gm = try GravitationalParameter(
    cubicMetresPerSecondSquared: 3.986004418e14,
    convention: .newtonian,
    provenance: PhysicalDataProvenance(source: "Application Earth GM", version: "1"))
let transfer = try HohmannTransfer(
    departureRadiusMetres: 7_000_000,
    arrivalRadiusMetres: 42_164_000,
    gravitationalParameter: gm)
let frame = Frame(name: "Earth-centred inertial", pseudoInertial: true)
let plan = try transfer.plan(
    departureEpoch: .j2000, frame: frame, centralBody: .earth,
    radialDirection: SIMD3(1, 0, 0),
    progradeDirection: SIMD3(0, 1, 0))
let arrival = try plan.trajectory.state(at: plan.arrivalEpoch)
```

`plan` normalizes the supplied finite, nonzero radial and prograde directions and
requires perpendicularity within 1e-12 in their unit-vector dot product. They
define the departure position, orbital plane and sense of motion; no global XY
plane or particular spin direction is assumed. The frame must be pseudo-inertial
and actually centred on the supplied central body. The planner cannot establish
that physical origin from a frame name. It builds an initial circular state, two
`ImpulsiveManeuver`s and a `ManeuveredTwoBodyPropagator`; query `.before` to inspect
pre-burn states, or use the existing sampler for rendering and prediction.
The trajectory continues on the departure orbit before the first burn and on
the destination orbit after the second. Queries at burn epochs return post-burn
states by default, following the existing maneuver API.

This computes arrival at an **orbit**, not rendezvous with a particular target.
Target phase, launch windows, plane changes, finite burns, obstacles/body surface
clearance and perturbations are separate constraints. It makes no claim of
minimum delta-v across other transfer families. Invalid radii, directions,
incompatible GM conventions, non-inertial frames and unrepresentable numerical
ranges fail explicitly; propagation/epoch errors retain their existing types.
Scaled radius ratios and rationalized speed differences preserve small changes
near equal radii. The executable plan remains subject to the two-body propagator's
numerical range and tolerances.

Verification uses independent 80-digit Decimal energy/vis-viva and Kepler-period
calculations (`Tools/generate_hohmann_fixtures.py`,
`Tests/CelMekTests/Fixtures/hohmann-reference.json`). Cases include inward/outward,
equal and nearly equal radii and a large radius ratio. Propagation tests exercise
an inclined plane, reversed sense of motion, both burn sides, arrival
circularization, later radius conservation and sampling/query-order independence.
For the classical transfer geometry and half-ellipse interpretation, see
[NASA/JPL's orbital mechanics chapter](https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/MPG%20Book/Release/Chapter7-OrbitalMechanics.pdf)
and [JPL's launch-window exercise](https://www.jpl.nasa.gov/edu/resources/lesson-plan/lets-go-to-mars-calculating-launch-windows/).

## Ideal brachistochrone planning

`IdealBrachistochroneTrajectory` plans minimum-time Newtonian motion between
stationary endpoints in one inertial frame, with an explicit positive maximum
acceleration magnitude. It accelerates along the straight-line displacement at
that limit, then reverses acceleration at the midpoint to stop at the destination.
The model treats both endpoint velocities as zero in the supplied frame.

For distance `d` and acceleration bound `a`, the half-duration is `sqrt(d/a)`,
total duration is `2*sqrt(d/a)`, peak speed is `sqrt(a*d)` and integrated
acceleration magnitude (total delta-v) is twice the peak speed. These quantities
use metres and atomic SI seconds. The minimum-time claim applies to these
rest-to-rest, gravity-free assumptions with an instantaneous acceleration reversal;
it does not optimize a gravity-driven spacecraft trajectory.

```swift
import CelMek
import Math

let frame = Frame(name: "Planning inertial frame", pseudoInertial: true)
let trajectory = try IdealBrachistochroneTrajectory(
    departureEpoch: .j2000,
    frame: frame,
    departurePositionMetres: SIMD3(0, 0, 0),
    arrivalPositionMetres: SIMD3(1_000_000, 0, 0),
    maximumAccelerationMetresPerSecondSquared: 1)
// 2,000 seconds total, switch at 1,000 seconds, peak speed 1,000 m/s.
let middle = try trajectory.state(at: trajectory.switchEpoch)
let samples = try TrajectorySampler.sample(
    trajectory, from: trajectory.departureEpoch, through: trajectory.arrivalEpoch,
    step: Duration(seconds: 20))
```

The immutable provider exposes departure/switch/arrival epochs, coverage, endpoint
positions, distance, duration, peak speed, total delta-v, acceleration direction
and braking direction. `phase(at:)` reports accelerating, braking or arrived.
`state(at:)` supplies position, velocity and acceleration at the requested epoch.
The switch is right-continuous with braking acceleration; arrival returns the
exact destination, zero velocity and zero acceleration. Position and velocity
remain continuous at the switch. Queries outside departure-through-arrival
coverage fail explicitly. Coincident endpoints produce a zero-duration stationary
state at the departure epoch, with zero direction, speed and acceleration.

All positions share the supplied frame and origin; the caller is responsible for
that geometry and for stationary endpoints. Construction requires a pseudo-inertial
frame. Existing frame conversion and sampling APIs can display the trajectory in
rotating frames without changing its native planning assumptions. No force, GM,
central body, spacecraft mass, exhaust speed or propellant budget is inferred.
Constant acceleration is prescribed kinematics; finite thrust with changing mass,
gravity, moving targets, nonzero endpoint velocities, attitude slew limits and
numerical trajectory optimization remain subsequent dynamics/planning work.

Invalid acceleration, nonfinite endpoints, non-inertial frames and numerical-range
failures have explicit `BrachistochroneError` cases. Epoch/duration arithmetic
errors propagate. Square-root scaling avoids forming an overflowing `d/a`; braking
is evaluated from remaining time, and positions use bounded endpoint weights.
Binary64 coordinate resolution still limits small displacements near large origins.
A nonzero trip whose switch/arrival cannot be represented as distinct epochs is
rejected. Queries are independent of order and do not integrate hidden steps.

Independent 80-digit Decimal fixtures integrate the constant-acceleration phases
(`Tools/generate_brachistochrone_fixtures.py`,
`Tests/CelMekTests/Fixtures/brachistochrone-reference.json`). Swift Testing checks
3D endpoint geometry, large coordinates, coincident endpoints, continuous
position/velocity, acceleration jumps, derivative consistency, coverage failures,
reverse sampling and rotating-frame derivatives. The bounded-acceleration
minimum-time control background is described in
[MIT's double-integrator discussion](https://underactuated.csail.mit.edu/dp.html)
and [Liberzon's double-integrator example](https://liberzon.csl.illinois.edu/teaching/cvoc/node85.html).

## Fly-through waypoint routes

`WaypointBrachistochroneRoute` plans gravity-free motion from A to B through an
ordered list of fixed inertial waypoints. It passes through intermediate points
with continuous **nonzero velocity**; waypoints are neither stops nor impulses.
Each leg has its own positive acceleration limit, allowing prescribed lighter
burns for cooldown or other operational constraints. No engine temperature model
is inferred. Route endpoints are at rest in this first version. Optional dwell
is accepted only at the initial/final endpoint; intermediate dwell is rejected.

```swift
let route = try WaypointBrachistochroneRoute(
    departureEpoch: departure, frame: frame, origin: .sun,
    waypoints: [
        BrachistochroneWaypoint(positionMetres: start),
        BrachistochroneWaypoint(positionMetres: avoidanceWaypoint),
        BrachistochroneWaypoint(positionMetres: destination)
    ],
    legAccelerationsMetresPerSecondSquared: [1, 0.25],
    obstacles: [SphericalExclusionRegion(
        identifier: "Sun", frame: frame, origin: .sun, centreMetres: .zero,
        radiusMetres: suppliedSolarRadius, safetyMarginMetres: clearanceMargin)]
)
let crossing = try route.state(at: route.waypointArrivalEpochs[1])
// crossing.velocity is nonzero; the next leg carries that exact velocity.
```

This is a **feasible smooth approximation**, not a minimum-time brachistochrone
optimizer. Each leg is a time-parameterized cubic Hermite curve, represented as a
cubic Bézier segment. Initial nominal durations use the ideal rest-to-rest model
with each leg's acceleration limit. Internal velocity follows the chord through
neighboring waypoints divided by the adjacent nominal durations. A common time
scale then reduces acceleration until every leg satisfies its own bound, with a
1% scale margin against floating-point rounding. Scaling one leg's requested
limit can also change neighboring passage speeds and route-wide timing.
The final controls and acceleration bounds are checked before returning a route.

For leg duration T and Bézier controls P0...P3, velocity is quadratic and
acceleration is affine in normalized time u. Its maximum magnitude is bounded
by the larger endpoint acceleration magnitude, giving a whole-leg acceleration
check rather than a sampled approximation. Position and velocity match at each
waypoint; acceleration may jump when switching legs. This assumes instantaneous
thrust-direction/magnitude changes with finite acceleration; jerk/attitude slew
constraints belong to later refinement. Consecutive coincident points and
neighbor geometry that yields zero passage velocity are rejected, including an
exact retracing turn whose neighboring endpoints coincide. More geometric points
are required to describe a nonstopping turn in that case.

`SphericalExclusionRegion` supplies an explicit frame/origin, fixed centre, body
radius and safety margin. Clearance checks follow the **curves**, not the waypoint
polyline or output samples. The [Bézier convex-hull property](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node12.html)
allows control-point bounding boxes to certify exclusion; uncertain boxes are
split with [de Casteljau subdivision](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node13.html).
Sampled points can demonstrate obstruction, but clearance is accepted only when
all boxes have positive distance from the exclusion sphere, including a coordinate
scale-dependent floating-point margin. Reports contain conservative lower bounds
on clearance, not the exact closest distance. Tangency/near-contact may produce
an obstruction or explicit unresolved-clearance error; unresolved cases never
silently pass. This is conservative floating-point geometry, not formal interval
arithmetic certification.

Subdivision depth is capped at 24, total clearance-node checks at 1,000,000,
and route legs/obstacles at 10,000 each, with at most 1,000,000 leg-obstacle pairs.
Invalid inputs, incompatible obstacle frames/origins, obstructed legs, unresolved
clearance and exhausted budgets throw explicit errors. Callers supply the route;
no automatic avoidance-waypoint generation is performed. Obstacles are fixed
throughout the route, so a Sun-centred scenario uses an explicit Sun-centred frame
and origin. Moving bodies need ephemeris-aware clearance refinement. The planner
has no gravity, mass, propellant or thermal evolution, and supplies no globally
optimal travel-time claim.

Arrival is the final waypoint crossing; completion additionally includes any final
endpoint dwell. Coverage is departure through completion. Querying a shared leg
boundary selects the outgoing leg, with unchanged position/velocity and its new
acceleration. Completion returns zero velocity/acceleration. The route and each
leg implement `StateProvider`, so existing sampling and explicit frame conversion
remain available without hidden clocks or mutable query history.

Swift Testing checks fly-through position/velocity continuity, nonzero passage
speed, unequal per-leg acceleration limits, independent cubic state values,
endpoint dwell, a solar-avoidance route, an interior curve obstruction with clear
endpoints, rejected intermediate dwell and degenerate retracing geometry.

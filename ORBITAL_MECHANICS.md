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
is inferred. Route endpoints default to rest and support explicit velocities. Optional dwell
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
acceleration. Completion preserves the specified final velocity (zero by default)
and returns zero geometric acceleration. This endpoint value does not prescribe
motion beyond coverage. The route and each
leg implement `StateProvider`, so existing sampling and explicit frame conversion
remain available without hidden clocks or mutable query history.

Swift Testing checks fly-through position/velocity continuity, nonzero passage
speed, unequal per-leg acceleration limits, independent cubic state values,
endpoint dwell, a solar-avoidance route, an interior curve obstruction with clear
endpoints, rejected intermediate dwell and degenerate retracing geometry.

### Numerical execution of a prescribed waypoint route

`WaypointRoutePropagator` executes a `WaypointBrachistochroneRoute` with
`SpacecraftMassProperties`, explicit `backgroundDynamics`, and one
`WaypointRouteEngineLimits` per leg. Each engine limit supplies a maximum thrust
acceleration, maximum force in newtons, and effective exhaust velocity. Supply
gravity through the existing force models; the route and dynamics must share
the exact frame and physical origin. Background models must have zero mass flow.

The engine commands the reference cubic acceleration minus the background
acceleration evaluated at the integrated spacecraft state. Required thrust is
the norm of this difference times current total mass; propellant consumption is
thrust divided by exhaust velocity. This assumes ideal instantaneous thrust
steering and perfect knowledge of the environment. Geometry bounds net
acceleration, while the engine limits bound thrust acceleration including gravity
compensation: those limits can differ substantially near a massive body.

Propagation starts at the route's departure state and replays forward for each
query. Integration splits at every waypoint, carrying the actual position,
velocity and remaining propellant into the next leg without resetting them to
the reference. Intermediate velocities remain nonzero. Endpoint dwells are
currently rejected. Completed-leg checkpoints report position and velocity
residuals, spacecraft mass, and integrator diagnostics. Explicit positive
endpoint tolerances reject excessive residuals; local RK tolerances alone are
not a global endpoint guarantee. Shared boundaries expose the outgoing leg's
acceleration; the final arrival exposes the final leg's left-limit acceleration.
No continuation beyond arrival is prescribed.

Engine constraints are checked at every derivative evaluation, including trial
stages. Violations, dry firing, unsupported background mass flow and integrator
failures throw; a feasible prefix query does not establish that the whole route
is feasible. These checks can conservatively reject a trial stage and are not a
continuous bound between evaluations. The geometric route's obstacle clearance
continues to describe its reference curve; endpoint residuals do not certify
whole-path clearance of the numerical execution.

This executes an existing feasible geometric plan under a chosen ideal engine;
it does not optimize waypoint locations or times, solve a minimum-time transfer,
or model attitude, thermal limits, moving obstacles or automatic avoidance.
The broader numerical refinement backlog remains open.

### Timed engine-limit phases

`WaypointRouteEnginePhase` specifies engine limits over `[startEpoch, endEpoch)`.
Its alternative initializer accepts a departure epoch, nonnegative start offset
and positive duration. These intervals are independent of the spatial waypoints:
a phase may begin within a leg, span several waypoint passages, or end exactly
at a passage. Input order is arbitrary. `WaypointRoutePropagator` sorts phases,
rejects overlaps (including duplicate intervals), and rejects intervals outside
the route's coverage. Adjacent intervals are allowed. At most 10,000 phases are
accepted.

Pass phases using `timedEnginePhases:` alongside the existing per-leg engine
limits. An active phase **replaces** all three per-leg engine limits; gaps use
the per-leg limits. These are upper limits on the engine, not commands to burn
at that magnitude. The inverse-dynamics controller still supplies the thrust
needed by the reference curve. To express 0.3g, use
`0.3 * WaypointRouteEngineLimits.standardGravityMetresPerSecondSquared`
(standard gravity is 9.80665 m/s²).

For example, given a route lasting at least three hours:

```swift
let full = try WaypointRouteEngineLimits(
    maximumAccelerationMetresPerSecondSquared: 9.80665,
    maximumThrustNewtons: 1_000_000,
    effectiveExhaustVelocityMetresPerSecond: 30_000)
let reduced = try WaypointRouteEngineLimits(
    maximumAccelerationMetresPerSecondSquared: 0.3 * 9.80665,
    maximumThrustNewtons: 1_000_000,
    effectiveExhaustVelocityMetresPerSecond: 30_000)
let phases = try [
    WaypointRouteEnginePhase(departureEpoch: route.departureEpoch,
        startOffset: .zero, duration: Duration(seconds: 3600), limits: full),
    WaypointRouteEnginePhase(departureEpoch: route.departureEpoch,
        startOffset: Duration(seconds: 3600), duration: Duration(seconds: 7200), limits: reduced)
]
```

Integration splits exactly at phase starts, phase ends and waypoint epochs,
using the incoming limits for the final RK stage of an interval and the outgoing
limits at a boundary query. Position, velocity and fuel are carried unchanged
across each split. `propagation(at:)` returns execution `segments` with initial
and final spacecraft states, the spatial leg index, the active sorted phase
index (or nil for per-leg defaults), and numerical diagnostics. Waypoint
checkpoint diagnostics aggregate all segments in that leg. Coincident temporal
and spatial boundaries cause a single split, without a stop or duplicate
positive-duration segment. Final arrival retains the incoming left limit.

This checks the existing spatial route against a temporal engine schedule. It
does not change the route timing to fit that schedule. The geometric planner
supports arrival windows through uniform retiming as described below. Solving
the relative order of untimed passages and temporal constraints, and
commanded burns with explicit directions remain planning extensions. Engine
checks retain the sampled-feasibility limitations described above.

### Waypoint arrival epochs and windows

`BrachistochroneWaypoint` accepts an optional inclusive `arrivalWindow`, or an
`arrivalEpoch` for an exact passage time. The first waypoint's window must
contain the route departure epoch. Arrival means reaching the waypoint, before
any permitted endpoint dwell. Intermediate waypoints still have no dwell.

The geometric planner first computes its existing acceleration-bounded route.
If passage constraints are supplied, it uniformly slows all flight legs to fit
them. For each constrained waypoint, the allowed time relative to flight start
divided by the nominal cumulative flight duration defines an interval of
allowed scale factors. The planner intersects these intervals with scales at
least one and selects the smallest common factor. Initial/final dwells are not
scaled. The route then rebuilds epochs and velocities, validates acceleration
bounds and checks obstacle clearance on the resulting curves. Exact arrival
epochs are preserved as `Epoch` values; fractional-second rounding at window
boundaries is resolved before the final checks.

```swift
let passage = try BrachistochroneWaypoint(positionMetres: waypointPosition,
    arrivalWindow: earliestPassage...latestPassage)
let destination = try BrachistochroneWaypoint(positionMetres: destinationPosition,
    arrivalEpoch: requiredArrival)
```

Slower timing reduces passage velocities without introducing a stop or waiting
at intermediate waypoints. Independently defined engine phases can then be
supplied to `WaypointRoutePropagator`; their epochs do not move when geometry
is retimed. Engine, gravity and fuel feasibility remain execution checks.

`uniformRetimingUnavailable` means these arrival constraints cannot be met by
uniformly slowing this nominal route. A different distribution of leg times,
curve or thrust policy may still work. This method does not search those
alternatives, optimize times under engine phases, or claim minimum travel time.
An acceleration failure after reconstructing precise epochs is reported as
`timingAccelerationLimitExceeded(leg:)`; geometry and clearance errors retain
their existing explicit diagnostics. Local acceleration bounds on the reference
curve and numerical endpoint tolerances retain their distinct meanings.

### Searching independent leg timings

`WaypointBrachistochroneRoute` optionally accepts `legDurations: [Duration]`.
These positive durations replace automatic timing; acceleration bounds, arrival
windows, nonzero passage velocities and whole-curve geometric clearance are
still checked on the resulting route. The planner never stretches an explicit
duration to repair a failed constraint. This allows schedules that uniform
retiming cannot satisfy.

`WaypointRouteTimingPlanner.search` accepts one array of candidate arrival
epochs for every waypoint after departure. It sorts and deduplicates each
array, then considers their Cartesian product in lexicographic order, prioritizing
earlier intermediate passage times. Spatial waypoints retain their input order;
candidate epochs determine their relative order with temporal phase boundaries.
Each ordered candidate builds the route with exact epoch-difference durations
and executes it using the supplied background dynamics, mass, engine limits,
timed phases and numerical options. This includes carried fuel, sampled thrust
feasibility, endpoint residuals and reference-curve obstacle checks.

The result reports `feasible` for the first candidate that passes all checks,
`candidateSetExhausted` when every combination fails, or
`evaluationLimitReached` when the budget ends with candidates remaining. A
feasible result includes the route and full propagation result. `failedAttempts`
retain arrival epochs and distinguish ordering, geometric, engine and numerical
integration failures. Integrator budget/tolerance failures are not proof of
physical infeasibility. Unexpected force/provider errors propagate to the caller;
they are not hidden as rejected timings.

The default budget is 1,000 considered combinations; callers may specify
1...10,000. Invalid ordering and rejected constraints count against that budget.
Input candidates are limited to 100,000 entries in total. Existing numerical
work limits apply separately to each propagation segment; there is no aggregate
derivative budget across candidates. Endpoint dwells are rejected by numerical
execution. Temporal phases must start at or after departure; a candidate ending
before a phase finishes is rejected and a later candidate may succeed.

For example, a physical passage constrained to 10...20 seconds followed by an
exact arrival at 50 seconds cannot fit equal nominal leg durations by uniform
scaling. Candidate passage epochs of 10 and 20 seconds produce independently
timed legs of 10/40 and 20/30 seconds. If a low-acceleration phase occupies the
first 15 seconds, the former candidate can fail while the latter passes,
moving the physical passage across the temporal boundary without stopping.
Swift Testing verifies this case, mass continuity, candidate normalization and
the distinction between candidate exhaustion and budget exhaustion.

This is a finite feasibility search, not continuous timing optimization. The
first feasible candidate need not have the earliest final arrival, and exhaustion
only describes the supplied set. Commanded-burn planning and continuous
optimization remain separate extensions. Numerical
execution does not certify whole-path clearance or continuous engine feasibility.

### Generating timing candidates from windows

`WaypointRouteTimingPlanner.arrivalCandidates(in:spacing:including:maximumTotalCandidates:)`
generates the per-waypoint arrays used by `search`. Supply a finite arrival
window for each spatial waypoint after departure and a positive `Duration`
spacing. Each array includes both inclusive window endpoints, a regular grid
anchored at the lower endpoint, and all supplied engine-phase starts/ends within
that window. Arrays are sorted and deduplicated; exact windows yield one epoch.
Off-grid phase boundaries are included exactly, preserving fractional epochs.

```swift
let candidates = try WaypointRouteTimingPlanner.arrivalCandidates(
    in: arrivalWindows,
    spacing: Duration(seconds: 3600),
    including: enginePhases,
    maximumTotalCandidates: 10_000)
// Pass candidates as candidateArrivalEpochs: to search with the same phases.
```

Windows are independent of phase order. A physical waypoint can have candidate
passages before, exactly at, or after a burn-phase boundary. Spatial passage
order is checked by the search. For an untimed waypoint, the caller supplies a
finite search window; no horizon or arrival deadline is inferred. A window
that admits one candidate need not admit any physically feasible transfer.

Candidate generation and candidate evaluation have separate budgets. The
generation limit defaults to 100,000 total entries across all arrays and accepts
1...100,000. More than 10,000 windows or phases is rejected. Generation stops
with `candidateGenerationLimitExceeded` instead of truncating a window; zero or
negative spacing is rejected. Precise epoch advancement must increase time.
The Cartesian product is enumerated lazily by the existing search, subject to
its own evaluation budget.

This supplies an automatic grid and phase-boundary policy. Adaptive refinement
and local timing improvement are available through the APIs below; global
optimization remains open.
Smaller spacing can find additional candidates at
greater cost, and candidate-set exhaustion does not establish global infeasibility.

### Adaptive timing feasibility search

`WaypointRouteTimingPlanner.adaptiveSearch` combines candidate generation and
execution over successive grid resolutions. Supply finite `arrivalWindows`
for spatial passages after departure and `WaypointTimingRefinementOptions`
with initial and minimum spacing in SI seconds. The search intersects these
windows with the waypoint's own arrival constraints; disjoint windows are an
invalid search configuration. Exact windows remain exact at every level.

At each level the search includes window endpoints and timed engine boundaries,
then considers the candidate product in the existing lexicographic order. It
stops at the first feasible route. If that grid is exhausted, spacing is halved,
bounded below by the requested minimum, and another level is searched.
Previously rejected exact arrival combinations are skipped across levels.
Providers and forces must remain pure and deterministic throughout the search;
inputs, force models and numerical tolerances are identical for every level.
Refinement is global grid refinement, not a local optimizer guided by a smooth
constraint residual. No monotonic feasibility assumption or bisection across
feasible/infeasible candidates is made.

The result includes aggregate candidate evaluations, failed attempts, and
per-level spacing, generated entries, new evaluations and discrete-search
status. Successful results include the selected route and its full numerical
execution. Termination distinguishes:

- `feasible`: a sampled candidate passed the route and execution checks.
- `resolutionLimitReached`: the minimum-spacing grid is exhausted, or all
  windows are exact and their one combination fails.
- `evaluationLimitReached`: the total candidate evaluation budget is consumed.
- `refinementLimitReached`: all permitted refinement levels are exhausted.
- `candidateGenerationLimitReached`: the next grid exceeds the generation budget.

Options allow 1...32 levels, 1...10,000 total new candidate evaluations and
1...100,000 generated entries per level. Defaults are eight levels, 1,000
evaluations and 100,000 generated entries. Generation and numerical budgets
retain their separate meanings; RK work is bounded per propagation segment.
Precise-epoch resolution errors and unexpected force/provider errors still throw.

The finest sampled grid can miss a narrow feasible interval. These termination
diagnostics do not prove global infeasibility or continuous convergence, and
the first feasible result does not minimize arrival time. Reference-clearance
and sampled-engine limitations also remain. Tests include a route whose coarse
arrival endpoints both fail engine limits but whose midpoint succeeds, verifying
that refinement finds it without evaluating the failed endpoints again.

### Local timing improvement from a feasible seed

`WaypointRouteTimingPlanner.optimizeTiming` starts with `seedArrivalEpochs`
and finite `arrivalWindows`, then tries to reduce final arrival time. A route
returned by discrete/adaptive search can supply its arrival epochs after dropping
the departure entry. `WaypointTimingOptimizationOptions` supplies positive initial
and minimum steps in SI seconds, an iteration limit and a total candidate
evaluation limit. Waypoint windows intersect the supplied search windows; seed
epochs must lie inside those intersections. The seed itself is fully executed
and counts against the evaluation budget.

At each iteration the optimizer tries an earlier final arrival with unchanged
intermediate epochs, proportional compression of all passage times, and positive
or negative changes to each individual intermediate passage paired with that
earlier arrival. Proposals are clipped to their arrival windows and generated
lazily. Already tested exact combinations are skipped. Only feasible proposals
with strictly earlier final arrival are accepted, and every candidate uses the
same geometry, obstacle checks, physical forces, fuel model and numerical options
as the existing search. If no pattern improves the incumbent, the step is halved
until the requested minimum. Temporal phase epochs stay fixed.

Results retain the best feasible route and full spacecraft execution even when
an evaluation or iteration limit interrupts improvement. `seedRejected` returns
no route and the seed's explicit failure. `stepLimitReached` means the tested
patterns found no improvement at the requested minimum step; it is not a proof
of stationarity or global optimality. `evaluationLimitReached` and
`iterationLimitReached` identify the other termination conditions. Diagnostics
include aggregate evaluations, candidate failures and each iteration's step,
new evaluations, improvement flag and best arrival epoch.

Defaults permit 256 iterations and 1,000 candidate evaluations; explicit limits
accept 1...1,000 iterations and 1...10,000 evaluations. Numerical limits remain
per propagation segment. Exact intermediate arrival constraints remain exact;
intermediate waypoints continue to be fly-through passages. If final arrival is
fixed by its window, this objective cannot improve it. No fuel-only objective,
equal-duration rearrangement or temporary acceptance of worse schedules is used.

This is a local direct-search heuristic in continuous epoch variables. Its
limited directions and strict feasible-improvement policy can miss routes that
require coordinated changes to several passages, a temporary worse arrival or
crossing an infeasible region. Multi-start/global optimization, commanded burns
and continuous constraint certification remain open. It inherits the sampled
engine and reference-clearance limitations of route execution.

Swift Testing checks an independent rest-to-rest cubic threshold:
`max |a| = 6*d/T²`, so a four-metre cubic at 1 m/s² requires at least
`sqrt(24)` seconds. The optimized result approaches this threshold within its
timing step; this does not make it the faster ideal bang-bang brachistochrone.
Tests also cover monotonic improvement, determinism, budget-retained feasible
results, rejected seeds and exact fly-through passage through a timed engine
schedule with mass continuity.

### Nonzero endpoint and prescribed passage velocities

`BrachistochroneWaypoint` accepts optional `velocityMetresPerSecond`, expressed
in the route's inertial axes relative to its physical origin. The value is a
constraint, not a timing-dependent hint. Unspecified endpoints remain at rest;
unspecified intermediate velocities use the existing neighboring-chord rule.
Explicit intermediate velocities must be nonzero, preserving fly-through
behavior. An explicit velocity can also specify a turn that the neighboring
chord rule would otherwise classify as degenerate; acceleration and clearance
checks still apply to the full resulting curve.

```swift
let departure = try BrachistochroneWaypoint(positionMetres: departurePosition,
    velocityMetresPerSecond: currentVelocity)
let arrival = try BrachistochroneWaypoint(positionMetres: arrivalPosition,
    arrivalEpoch: requiredArrival, velocityMetresPerSecond: requiredVelocity)
```

All route timing search, refinement and local optimization APIs carry these
constraints through the supplied waypoint values. Numerical execution starts
with the specified departure velocity, carries actual velocity/mass across
legs, and checks the prescribed arrival velocity with its endpoint tolerance.
The final geometric state preserves the specified velocity. Numerical execution
retains its existing final incoming acceleration convention. Neither provider
prescribes continuation after its final epoch.

Fixed velocities do not decrease when durations increase. Automatic geometric
timing therefore re-evaluates each cubic's actual acceleration bound, doubling
all durations as necessary, with a 64-check limit and a one-percent acceleration
margin. Failure reports `velocityTimingLimitReached` or an explicit numerical
range error. This is conservative feasible timing, not minimum time. Arrival
window retiming rebuilds and rechecks the resulting acceleration bounds; with
fixed velocities, slowing a curve is not assumed to reduce its acceleration
monotonically. A failed automatic schedule may still admit an explicit-duration
schedule or a different timing candidate.

Explicit leg durations are never extended to accommodate velocity constraints.
They must satisfy the acceleration and arrival-window checks as given. Nonzero
velocity with endpoint dwell is rejected as `velocityDwellConflict(index:)`,
since a stationary dwell would require an unmodelled impulse. Nonfinite velocity
components are rejected. Large velocities and long flight times can bend a
cubic beyond its waypoint chords, so whole-curve clearance is checked on the
actual controls. There is no implicit turn smoothing, jerk/attitude constraint,
implicit moving-target lookup. Explicit moving endpoints are available below;
there is no general orbital rendezvous or optimal-control solver.

Swift Testing uses independent constant-velocity and constant-deceleration
references, analytic exhaust-velocity fuel consumption, constrained timing
optimization, intermediate velocity continuity and invalid velocity/dwell cases.

### Moving-target rendezvous endpoints

`MovingWaypointRouteTarget` adapts an `EphemerisProvider` or another
`StateProvider`. Its ephemeris initializer uses the provider's physical origin,
frame and declared coverage and invokes `checkedState(at:)`. For arbitrary
state providers the caller supplies and guarantees origin and coverage. Every
query checks coverage before evaluation and verifies the returned epoch and
exact frame identity. Providers must be pure and deterministic during planning.

Pass this adapter as `movingTarget:` to discrete `search`, `adaptiveSearch` or
`optimizeTiming`. For each candidate final epoch, the planner queries the target
anew and replaces the final waypoint's nominal position with the target position
and its velocity with the target velocity. The final waypoint still supplies its
arrival window. Its nominal position is ignored, and it must leave its optional
velocity unspecified; an explicit static velocity would conflict with matching
the moving target and is rejected. Intermediate waypoints default to fixed spatial
constraints and can use the moving-passage API below. Endpoint dwells remain
unsupported by numerical execution.

```swift
let target = try MovingWaypointRouteTarget(ephemeris: targetEphemeris)
// Pass movingTarget: target to the timing planner. The final waypoint supplies
// the arrival window; its position is a nominal placeholder and its velocity nil.
```

For a single chosen arrival epoch, `target.waypoint(at:)` returns a waypoint with
that exact epoch and sampled target position/velocity, usable with explicit leg
durations. This fixed-time helper does not track later timing changes; use the
planner's `movingTarget:` parameter when searching or optimizing arrival time.
Keep the original waypoint template for those searches; a returned route's final
waypoint already contains a resolved velocity at its selected arrival epoch.

The target must already use the route's inertial frame and physical origin.
Use explicit frame/origin conversion through the existing frame system before
constructing the target adapter; no implicit transformation or origin translation
is applied here. Matching uses geometric SI position/velocity at the same epoch,
without light-time or apparent-position corrections. Target acceleration is not
an arrival constraint, and the result prescribes no continuation beyond arrival.

Discrete out-of-coverage arrival candidates are recorded as
`targetOutsideCoverage`, permitting later candidates to be tested. Adaptive search
and local optimization intersect the final search window with target coverage;
a disjoint window is an invalid search configuration, and an optimizer seed must
lie inside the intersection. Frame/origin mismatches, malformed returned states
and unexpected provider errors propagate instead of being hidden as rejected
timings. Coverage is never extrapolated. One target lookup is made per ordered,
newly evaluated candidate; force-model provider calls retain their own numerical
work limits.

The resulting route and numerical execution match the candidate target state
within the existing endpoint tolerances, with carried propellant and sampled
engine feasibility checks. Timing improvement remains local to the cubic route
family, with no global minimum-time claim or guarantee of continuous clearance.
Targets may be ephemerides or spacecraft propagators. Moving intermediate
waypoints are supported below; moving obstacle avoidance remains separate.

Swift Testing uses an independent linear target `p(T)=4+0.5*T`, `v(T)=0.5`.
The rest-to-rendezvous cubic's departure acceleration is `24/T²+2/T`, giving a
six-second threshold at 1 m/s². Tests check candidate resampling, matching both
position and velocity, timing improvement against this threshold, ephemeris
adaptation, coverage handling and explicit contract failures.

### Moving intermediate waypoint passages

`MovingWaypointPassage` supplies a moving spatial constraint for an intermediate
waypoint using the same checked provider adapter as rendezvous targets. Pass
`movingWaypoints: [waypointIndex: passage]` to discrete `search`, `adaptiveSearch`
or `optimizeTiming`. Indices refer to the complete ordered waypoint array:
departure is zero, intermediate indices are `1..<waypoints.count-1`, and final
arrival uses the separate `movingTarget:` rendezvous option. Other indices are
rejected. Multiple moving passages can be combined with a moving final target.

```swift
let body = try MovingWaypointRouteTarget(ephemeris: bodyEphemeris)
let flyby = MovingWaypointPassage(provider: body)
// Pass movingWaypoints: [1: flyby] to a timing planner with at least three points.
```

Every new, ordered timing candidate resolves each moving passage at its own
candidate passage epoch, in waypoint-index order, before constructing the cubic
legs. The resolved position replaces the template's nominal position; its
arrival window and optional explicit route velocity remain constraints. The
default `.preserveRouteVelocity` mode preserves an explicit waypoint velocity
or generates the existing neighboring-chord fly-through velocity. It does not
match the provider's velocity and therefore does not stop at a stationary body.

`.matchProviderVelocity` instead constrains spacecraft passage velocity to the
provider's velocity. A separate explicit template velocity conflicts with that
mode and is rejected as `conflictingWaypointVelocity(index:)`. Matching a zero
provider velocity at an intermediate point is rejected by the route's existing
nonzero-passage rule. Both modes preserve position/velocity/mass continuity
across legs and undergo the same acceleration, engine, fuel and geometric
clearance checks as fixed spatial waypoints. Provider acceleration is not used
as a passage constraint.

The provider must share the exact route frame and physical origin; declared
coverage and returned state contracts use the existing target checks. Discrete
coverage misses are recorded as `movingWaypointOutsideCoverage(index:)`, while
unexpected provider errors propagate. Adaptive search and optimization intersect
each moving waypoint's search window with that provider's coverage; disjoint
windows are invalid configuration, and optimized seeds must lie inside the
intersection. Timing changes trigger fresh provider queries. No extrapolation,
automatic origin conversion, departure-point motion or reordering of spatial
constraints is performed.

`passage.waypoint(at:template:)` can sample a fixed-time passage for manual
construction; use explicit durations consistent with that sampled epoch. Keep
the original template and `movingWaypoints:` when searching different times.
Search work includes at most one query per moving passage and one final-target
query for each new ordered candidate, stopping on the first rejected lookup.
Existing candidate, waypoint and numerical work limits still apply.

A moving passage is a position constraint at one epoch, not a moving exclusion
region or whole-path avoidance guarantee. Moving obstacle bounds, target-relative
offset geometry, light-time corrections and global optimal control remain open.
Tests cover time-dependent positions, coverage rejection, independent velocity
choices, stationary-body fly-through, zero-velocity matching rejection,
provider-coverage-limited refinement/optimization and carried mass continuity.

### Moving spherical exclusions

Routes and all three timing planners accept `movingObstacles: [MovingSphericalExclusionRegion]`.
Each region uses a checked `MovingWaypointRouteTarget` centre provider, a radius,
a safety margin, and an explicit `maximumCentreSpeedMetresPerSecond`. The caller
must guarantee that speed bound throughout the provider's declared coverage in
the route's native frame and physical origin. An instantaneous ephemeris velocity
does not establish this bound.

For each subdivided cubic interval, the centre at its midpoint and the speed
bound enclose the obstacle's motion in a sphere of radius
`maximumCentreSpeed * halfInterval`. The distance from that centre to the cubic
control polygon's bounding box, minus the motion radius and exclusion radius,
gives a conservative clearance bound. Positive bounds, reduced by a floating
point allowance, accept intervals; endpoint/midpoint samples can demonstrate
collision but cannot establish safety. Ambiguous intervals subdivide up to depth
24, with a shared static/moving budget of one million interval checks per route.
Unresolved clearance and exhausted work budgets are explicit failures.

Reports in `route.clearances` cover each entire geometric leg. Obstacle coverage
must include both leg endpoints; coverage misses appear as indexed geometry
failures in timing search. Every candidate checks the obstacle at its own epochs,
including adaptive refinement and local optimization. Provider errors propagate.
Identifiers must be unique across fixed and moving regions. Routes with endpoint
dwells currently reject moving exclusions explicitly rather than leaving those
stationary intervals unchecked.

The guarantee depends on the declared motion bound and uses floating point
geometry, not formal interval arithmetic. It applies to the prescribed reference
curve; numerical integration residuals do not certify executed whole-path
clearance. Automatic avoidance waypoint generation and relative waypoint offsets
remain separate extensions.

### Relative passages and rendezvous

A checked moving target can be shifted with
`target.offset(positionMetres: offset)`. The offset is constant in its native
axes, so the target velocity and acceleration are retained. To express the
displacement in body-fixed or other rotating axes, use
`target.offset(positionMetres: offset, in: axesFrame, using: frameTree)`.
The frame system rotates the displacement and contributes its angular velocity
and angular acceleration terms. Frame translations are deliberately excluded:
the target provider supplies the physical centre. Missing target acceleration
remains missing. Both forms preserve physical origin and declared coverage and
can be used as a moving passage or final rendezvous target in every timing planner.
Providers and the frame tree must remain deterministic during a search.

### Exploration from multiple timing seeds

`optimizeTimingFromSeeds` runs bounded local improvements from caller-selected
arrival vectors and retains the earliest feasible result. Seeds can span separate
timing basins within the arrival windows; identical seeds are visited once, and
rejected seeds do not discard another seed's feasible result. Per-start diagnostics
are retained, including local termination and rejected proposals, with an aggregate
candidate evaluation limit. Local iteration limits apply separately to each start.
The local patterns now include joint intermediate shifts and shifts of each
prefix and suffix, paired with earlier final arrival, in addition to final-only,
proportional and individual-passage changes.

This is finite exploration plus local improvement. Exhausting the seeds or the
requested step resolution proves neither global infeasibility nor stationarity
or global minimum time. All proposals re-evaluate moving targets, exclusions,
engine schedules and physical constraints.

### Prescribed attitude and thermal constraints

`WaypointRoutePhysicalConstraints` is accepted by numerical route execution and
all timing planners. A `WaypointThrustCone` supplies a pure epoch-dependent axis,
native frame and half-angle in radians. Execution checks the *required engine
acceleration*, including compensation for background gravity, against that cone
at derivative stages. No thrust means no direction restriction. Cone checks are
sampled, like the engine limits, and do not certify between-stage compliance.
These are prescribed attitude constraints; rigid-body attitude and slew dynamics
belong to Milestone 8.

A `WaypointThermalEnvelope` explicitly defines stored heat in joules, a maximum,
heating watts per newton of thrust, and a linear cooling rate per second:
`dH/dt = heatingWattsPerNewton * thrust - coolingRatePerSecond * H`.
The engine ceiling for each spatial/temporal execution segment bounds heating.
Its constant-power solution provides a conservative carried envelope, including
a small floating arithmetic allowance; both endpoints bound that solution's
maximum. Heat is carried across waypoints and engine-phase boundaries, and
segments report their initial and final heat bounds. Zero acceleration/thrust
ceilings permit explicit coast/cooldown phases; exhaust velocity remains positive.

The envelope is conditional on the engine ceiling being respected continuously.
It is not a reconstructed actual temperature history, and it can reject physically
cooler trajectories. `thermalBoundExceeded` means this conservative model could
not establish thermal feasibility; it is not a proof of physical infeasibility.
No heat model is inferred merely because a lighter-burn phase was requested.

### Commanded-burn shooting

`CommandedBurnPlanner.search` evaluates a bounded list of
`CommandedBurnCandidate(arrivalEpoch:schedule:)` values. Schedules contain explicit
inertial `FiniteBurn` directions, thrust and exhaust velocity. The planner uses
`ScheduledThrustPropagator` under caller-supplied background forces; it does not
replace commands with inverse-dynamics thrust. Candidates are ordered by arrival
epoch (stable input order for ties), and the earliest supplied feasible one wins.

The final checked moving target requires position and velocity matching within
explicit tolerances. Independent `CommandedBurnPassage` values require position
matching and optionally velocity matching, at strictly increasing physical
passage epochs. A passage may occur inside a burn or coast, with nonzero velocity.
Numerical execution carries position, velocity and mass through these passages;
ignition/cutoff boundaries remain separate. Burns must lie within the candidate
flight, and target coverage, residuals, fuel reserve at checkpoints and numerical
failures are diagnosed. Failed attempts retain available residual checkpoints.
Provider/configuration errors propagate rather than being classified as
trajectory infeasibility. Dry-mass checks occur inside numerical propagation.

Candidate count/evaluation budgets are bounded; integrator limits apply to each
propagation segment. Candidate-set exhaustion is not global infeasibility.
This planner does not automatically synthesize burn directions, enforce waypoint
route cones/thermal envelopes, or certify between-checkpoint obstacle clearance.
Its commands are chosen explicitly by the caller; obstacle and physical
constraint planning for prescribed cubic routes remains the separate route API.

Milestone 6 provides bounded planning and numerical simulation with explicit
limits. General optimal control, automatic avoidance geometry and formal
continuous-constraint certification are future extensions.

### Circular-target Hohmann rendezvous

`CircularHohmannOrbit` supplies an explicit circular two-body model with its own
reference epoch, radius, radial/prograde directions, central body, native frame
and GM. Directions are normalized; a perpendicularity tolerance of 1e-12 is
checked and tiny residual nonorthogonality is removed. The model provides
analytic position, velocity and acceleration. Its direction-based initializer
does not infer orbital shape; checked state/provider adapters are described below.

`HohmannRendezvousPlanner.search` accepts departure and target circular models
and an inclusive departure window. It requires exact frame identity, the same
central body, matching GM value/convention and coplanar orbital normals agreeing
within 1e-12. Opposite orbital directions are rejected; global clockwise or
counterclockwise motion is otherwise supported by the supplied axes. The caller
guarantees the physical frame origin.

The target's required departure phase is
`phi = pi - nTarget * transferDuration`, reduced to [0, 2*pi), positive in the
departure orbit's direction of motion. Relative phase evolves at
`nTarget - nDeparture`; the planner computes the first analytic opportunity in
the window and subsequent synodic opportunities, rather than sampling an epoch
grid. A stable mean-motion difference is used for nearly equal radii. This is
the circular interplanetary phase construction described in the
[FreeFlyer University Guide](https://ai-solutions.com/wp-content/uploads/2023/04/FreeFlyerUniversityGuide.pdf).

Every evaluated opportunity builds the existing two-impulse propagated plan and
compares its post-arrival-burn position and velocity against the target at the
same exact arrival epoch. Results retain departure phase error and actual
position/velocity residuals. The earliest verified opportunity is returned;
arrival may lie after the departure window. Work is bounded by an explicit
opportunity evaluation limit, with distinct window-exhausted and evaluation-limit
statuses. Propagation, epoch and unrepresentable-motion errors remain explicit.

Equal radii have no relative phase drift in this model. Coincident position and
velocity within the requested tolerances return `alreadyRendezvoused` at the
window's lower bound with no artificial transfer plan. Separated co-orbital
targets return `sameOrbitPhasingRequired`; waiting or the existing zero-burn
half-period coast cannot resolve their separation.

```swift
let departure = try CircularHohmannOrbit(
    referenceEpoch: .j2000, frame: frame, centralBody: .earth,
    radiusMetres: 7_000_000,
    radialDirection: SIMD3(1, 0, 0), progradeDirection: SIMD3(0, 1, 0),
    gravitationalParameter: gm)
let target = try CircularHohmannOrbit(
    referenceEpoch: .j2000, frame: frame, centralBody: .earth,
    radiusMetres: 42_164_000,
    radialDirection: SIMD3(1, 0, 0), progradeDirection: SIMD3(0, 1, 0),
    gravitationalParameter: gm)
let result = try HohmannRendezvousPlanner.search(
    departureOrbit: departure, targetOrbit: target,
    departureWindow: Epoch.j2000...Epoch.j2000.advanced(by: Duration(seconds: 100_000)))
if let plan = result.plan {
    let rendezvous = try plan.trajectory.state(at: plan.arrivalEpoch)
}
```

This is rendezvous for explicit circular two-body models, with caller-selected
position/velocity tolerances; it is not a planetary capture plan or numerical
multi-body targeting solution. Large epoch offsets and poorly conditioned
near-equal-radius opportunities can lose floating phase precision. Residual
checks protect against accepting a model opportunity that fails the requested
tolerances, but do not certify fidelity to real ephemerides. Checked circular
state/provider adapters are available below; parking-orbit escape/capture and
general transfers remain separate backlog items. The original `HohmannTransfer.plan` continues to
provide an orbit transfer without requiring a target.

### Checked circular states and providers

`CircularHohmannOrbit(referenceState:centralBody:gravitationalParameter:validationOptions:)`
validates a local circular approximation. `CircularHohmannValidationOptions`
specifies relative circular-speed error, radial-speed fraction of circular speed
and the permitted difference between unit plane normals. Defaults are 1e-8,
1e-8 and 1e-10 respectively; all must be finite, nonnegative and below one.
Circular speed is evaluated from the supplied central-body GM and sampled radius.
Accepted radial velocity and circular-speed error are projected out of the
model explicitly; eccentric or inclined input beyond the selected tolerances
is rejected rather than automatically redirected or silently approximated.

The provider initializer samples a checked `MovingWaypointRouteTarget` at an
explicit epoch and requires its physical origin to equal the central body.
The existing ephemeris initializer on that wrapper supplies native metadata and
coverage checks. No frame conversion or origin inference is performed.

`CheckedHohmannTransferPlanner.plan` takes departure and target states at the
same epoch, with exact frame identity and an explicit central body/GM. It checks
circularity, coplanarity and co-directional motion before constructing a transfer.
Its trajectory starts from the **actual sampled departure position and velocity**.
The first impulse removes the accepted radial/circular-speed residual and inserts
the spacecraft into the circular-fit Hohmann transfer; the arrival impulse
circularizes in the transfer plane. It does not substitute a remote target's
velocity when its physical position has been missed.

The result exposes each velocity projection error, the actual burn vectors and
their summed norms, plus target position/velocity residuals and
`rendezvousSatisfied`. Scalar delta-v on the underlying `HohmannTransfer`
describes the ideal circular model; the checked plan's impulse norms include
the accepted departure-state corrections. Approximation tolerances and
rendezvous tolerances are independent. A relaxed input tolerance does not make
a missed target feasible. With state-only inputs, the arrival target is the
fitted circular model; this is explicitly a model prediction.

`CheckedHohmannTransferPlanner.search` instead accepts checked departure/target
providers and a bounded list of candidate departure epochs. It sorts and
deduplicates dates, samples both providers afresh at each date and queries the
target again at the exact transfer arrival. The arrival state must also meet
circularity and plane checks, and actual arrival residuals decide feasibility.
The earliest feasible supplied date wins. Failed attempts retain available plans
and residuals, with distinct coverage, unsupported-geometry and rendezvous failures.
Provider contract/data errors and numerical propagation errors propagate.

For equal sampled radii, already-coincident states return
`alreadyRendezvoused` with no transfer and no requirement for future provider
coverage. Separated states record `sameOrbitPhasingRequired` for that date;
the search may continue because a provider is not assumed to obey one circular
model throughout its coverage. Candidate-set and evaluation-limit exhaustion
are distinct and do not prove global infeasibility.

Use the analytic circular rendezvous planner to obtain launch-date initial guesses
from checked circular models, then pass selected dates to the provider-backed
search for verification. This does not continuously refine provider-backed
launch dates or validate their full trajectories as circular. Lambert targeting,
noncoplanar transfers and numerical refinement remain separate extensions.

### Parking-orbit escape and capture

`ParkingOrbitTransfer` computes a local periapsis impulse between a circular
parking orbit and a body-relative hyperbola, using explicit parking radius,
body radius, minimum periapsis altitude, excess speed and body GM. Radii are
centre-relative metres; the parking radius must exceed the body radius and
meet the inclusive minimum-altitude bound. A positive body radius and a
Newtonian or TT-compatible GM are required.

For circular speed `vc = sqrt(mu/rp)`, escape speed `vesc = sqrt(2)*vc`
and specified excess speed `vinf`, the hyperbolic periapsis speed is
`vp = hypot(vinf, vesc)`. The impulse magnitude is `vp-vc` for either
injection or capture; its direction is opposite for the two maneuvers.
The result also exposes characteristic energy `C3 = vinf²`, eccentricity
`e = 1 + rp*vinf²/mu` and the acute asymptote/apse-line angle
`beta = acos(1/e)`, evaluated without cancellation near the parabolic limit.
These are the local escape/capture energy relations in the
[FreeFlyer patched-conics tutorial](https://freeflyer.com/_freeflyeruniversityguide/patched_conics_transfer.htm).

`plan(kind:periapsisEpoch:frame:centralBody:hyperbolicExcessVelocityMetresPerSecond:orbitalNormal:)`
constructs an explicit impulse and sampleable body-relative trajectory.
The excess vector is the **velocity direction**, incoming for capture and
outgoing for escape; it is not a position vector or a heliocentric velocity.
Its magnitude must match the scalar excess speed within the selected relative
tolerance (default 1e-12). The normal supplies the right-handed orbital plane
and sense of motion, with perpendicularity checked within 1e-12.

For unit periapsis radial/prograde directions `rhat, that`, the asymptotic
velocity directions are `-rhat/e + sqrt(1-1/e²)*that` on the outgoing
branch and `+rhat/e + sqrt(1-1/e²)*that` on the incoming branch.
The planner solves for those periapsis axes from the supplied vector and normal,
so the local burn is generally not parallel to the excess vector.

Escape accelerates from the circular parking state to the hyperbolic periapsis
state; capture brakes from the hyperbolic periapsis state to the circular state.
The result contains pre/post-burn states, periapsis axes, an `ImpulsiveManeuver`
and a `ManeuveredTwoBodyPropagator`. Its queries before and after the impulse
follow the corresponding circular/hyperbolic branches, including backward
queries into an incoming capture approach. The frame must be pseudo-inertial
and physically centred on the explicitly supplied body; names do not establish
that origin.

A zero excess speed is supported by scalar calculations as the parabolic escape
limit, but has no nonzero velocity direction at infinity. Vector-defined plans
reject it explicitly. Unrepresentable ranges, inconsistent vector magnitude,
unsafe radius and invalid plane/direction inputs fail explicitly.

These costs belong to the local body-relative maneuver budget. They must not
be substituted for, or blindly added to, heliocentric transfer impulses whose
role is already represented by the departure/arrival excess velocities.
The plan assumes a spherical clearance boundary, an impulsive tangential burn
and a circular parking orbit in the selected encounter plane. It does not solve
parking-orbit phasing, finite-burn losses, atmospheric clearance, encounter
boundary handovers or a full interplanetary trajectory.

### Plane changes and combined impulses

`PlaneChangeImpulse` specifies before/after tangential speeds and a velocity-plane
rotation angle in [0, pi] radians at a common node. Its cost is
`hypot(vAfter-vBefore, 2*sqrt(vBefore*vAfter)*sin(angle/2))`, evaluated with
scaled factors to avoid unnecessary overflow and cancellation. Equal speeds
give the pure plane-change result; unequal speeds combine speed and plane
changes in one impulse. This is the vector relation described in
[ASU's orbital-maneuver lecture](https://control.asu.edu/Classes/MAE462/462Lecture09.pdf).

The `maneuver` method returns an explicit-epoch/native-frame impulse from an
outward radial direction and initial prograde direction. Positive sense uses
the right-hand rule about that outward radius; negative sense reverses the
rotation. The directions must be perpendicular within 1e-12 and the frame
pseudo-inertial. The method subtracts the actual before/after velocity vectors;
costs of separate burns are not added to approximate a combined burn.
Both velocities are tangential. General radial-velocity changes and choosing
the appropriate intersection/node epoch remain caller responsibilities.

### Bi-elliptic transfers

`BiEllipticTransfer` joins two half-ellipses through a finite, explicitly selected
intermediate apoapsis strictly above both circular endpoint radii. The first
impulse raises apoapsis, the second changes periapsis at that apoapsis and the
third circularizes at the destination. Inward and outward endpoint transfers
are supported; the signed intermediate speed change reverses with their order.
The model exposes each half-ellipse's semi-major axis and exact `Duration`,
signed endpoint speed changes, the middle impulse and total absolute delta-v.
See the [FreeFlyer bi-elliptic tutorial](https://freeflyer.com/_freeflyeruniversityguide/bi_elliptic_transfer.htm)
for the vis-viva construction.

An optional `planeChangeRadians` combines a plane rotation with the intermediate
speed change. `plan` chooses its sense and emits three explicit impulses plus a
`ManeuveredTwoBodyPropagator`, intermediate/arrival epochs and expected arrival
state. The plane rotates about the **outward apoapsis radius**, opposite the
initial departure radius. Position is continuous at the node, and the final
circular orbit lies in the rotated plane. No destination body rendezvous or
parking-orbit capture is implied.

A rationalized difference preserves the small signed middle speed change for
nearly equal endpoint radii. `intermediateDeltaVelocityMetresPerSecond` supplies
the middle cost used in total accounting; for zero plane rotation it uses that
stable difference. The middle helper's separately rounded before/after speeds
may lose such a tiny difference when inspected independently. Equal endpoint
radii describe a raised-apoapsis excursion, not a zero-duration identity.
Invalid intermediate radii, non-inertial frames, unsupported GM conventions,
unrepresentable ranges and epoch/propagation failures remain explicit.

### Comparing circular transfer strategies

`CircularTransferComparison` ranks the existing coplanar Hohmann transfer and
bi-elliptic transfers at caller-supplied intermediate radii. It deduplicates and
sorts those radii, accepts at most 10,000 entries and supports
`minimumDeltaVelocity` or `minimumDuration`. Optional inclusive duration and
delta-v bounds filter candidates; results retain ranked and rejected candidates.
`best` is nil when none meets the bounds. Ties use the other objective, then
deterministic input strategy order. Invalid model inputs throw rather than
being mistaken for a budget rejection.

This comparison uses ideal two-body, coplanar circular costs with no plane
change. Compare combined noncoplanar maneuver costs explicitly using the helpers;
a coplanar Hohmann plan is not silently promoted to a plane-changing strategy.
The comparison does not certify target phasing, body clearance, perturbations,
finite-burn losses or propellant feasibility. Its best result is best among the
supplied radii and strategies, not a global optimum over all transfer families.

### Finite-burn transfer realization

`FiniteTransferPlanner.initialGuess` accepts the impulses from Hohmann,
bi-elliptic, parking-orbit or other plans, together with an actual spacecraft
state, arrival epoch, thrust and effective exhaust velocity. It computes each
constant-thrust duration from `mConsumed = mBefore * -expm1(-deltaV / ve)` and
`duration = mConsumed * ve / thrust`, carrying mass between burns. Zero impulses
are omitted. Burns are centered on their impulse epochs, shifted inside the
flight at its endpoints; overlapping burns, inadequate fuel and burns that do
not fit are rejected. Frames must match exactly. The input state is the physical
pre-burn state; the context origin must describe the same physical origin as the
transfer and target.

`FiniteTransferPlanner.refine` performs bounded local coordinate shooting against
an explicit Cartesian arrival state at a fixed epoch. It adjusts ignition,
duration and inertial direction while keeping each engine's thrust/exhaust
velocity fixed, preserving burn order and non-overlap. Supply `CentralGravity`
and any other background forces explicitly. Each trial numerically propagates
actual commands and evolving mass, accounting for finite-burn motion under those
forces. Position and velocity errors are scaled by caller tolerances; convergence
requires both tolerances separately. Results retain the initial and best schedules,
actual arrival states, endpoint errors, fuel consumed, evaluation count and
rejected trial count, even when the evaluation or step limit is reached. Remaining
fuel is available in each arrival state's mass; a caller reserve is checked.
Initial propagation failures throw; numerical/fuel failures in trial schedules
are rejected, while other force-model errors propagate to the caller.

This is a local refinement of a seed, with no global optimum or feasibility
certificate for collision clearance, attitude, thermal limits or intervening
waypoints. Evaluation limits bound propagation calls; integrator budgets apply
per scheduled segment. A positive fuel margin is required to avoid firing exactly
through depletion. Rocket-equation seeds do not imply impulsive endpoint
equivalence: tests compare actual Hohmann realizations under central gravity and
verify that progressively shorter burns approach the impulsive result.

### Zero-revolution Lambert arcs

`LambertSolver.solve` takes two central-body-relative positions in metres,
explicit departure/arrival epochs, an inertial frame, central body and compatible
GM. It returns both endpoint velocities in metres per second and Cartesian states,
plus a `trajectory()` convenience method for independent two-body propagation.
Input states can supply their positions and epochs; their velocities do not
participate in the boundary-value solution. The caller asserts the physical origin
and must express both positions in the same supplied frame.

A nonzero `referencePole` and `LambertDirection` are mandatory. Positive means
angular momentum has positive projection on that pole; negative means negative
projection. The pole need not be the orbit normal. Geometry and that direction
select the short-way (angle below pi) or long-way (above pi) zero-revolution arc.
The result reports the branch, angle and actual unit orbital normal. It does not
silently choose the shortest arc. This API does not match spacecraft or body
velocities or generate rendezvous/capture burns.

The native Swift solver scales radius, velocity and time, brackets the universal
time equation below its first positive singularity, and bisects within caller
work limits. The universal-variable f/g equations follow the formulation described
in [Vallado's Lambert routine](https://github.com/CelesTrak/fundamentals-of-astrodynamics/blob/main/software/matlab/lambertu.m).
Small-argument Stumpff series and half-angle differences reduce cancellation.
The absolute/relative tolerances bound the time-equation residual, not physical
endpoint errors; tighter reference comparisons request tighter time tolerances.
Results include iteration count, bracket evaluation count and signed time residual.

Collinear/coincident/opposite endpoints, and near-collinear geometry with normalized
cross-product magnitude at or below `1e-8`, are explicitly rejected. A pole whose
projection on the geometry's unit normal is at or below `1e-12` cannot define a
direction and is rejected. Revolution counts outside `0...32` are unsupported. Invalid
inputs, incompatible GM conventions, rotating frames, numerical range limits,
bracket limits and convergence limits report errors. Extreme flight times can
fail within the bounded numerical domain; the solver does not change branches
or relax tolerances to return an answer.

Tests compare velocities against independent circular and Kepler-anomaly elliptic,
hyperbolic and Barker-equation parabolic references, then propagate the solutions
to their arrival positions and velocities. Both directions, long-way motion,
inclined planes, time reversal, near-collinear accepted cases and bounded failures
are covered. Multi-revolution support is described below; ephemeris-driven transfer windows
remain a separate backlog item.


### Multi-revolution Lambert branches

`LambertSolver.solve` supports `revolutions` from 1 through 32 when an explicit
`multiRevolutionBranch` is supplied. `lowerUniversalParameter` and
`upperUniversalParameter` select roots on the decreasing and increasing sides of
the time-of-flight curve's minimum. These names describe numerical root identity;
they do not rank energy, delta-v or focus position. Short/long-way geometry and
orbital direction retain the existing reference-pole conventions independently
of this choice. Zero-revolution solves reject a multi-revolution branch argument.

`solveBranches` enumerates both roots for a caller-selected positive revolution
count in lower/upper order. It returns one solution when the roots coalesce within
the caller's time tolerance. It throws when a solve fails instead of silently
omitting a branch. No branch switching is performed. The reported transfer angle
excludes completed revolutions: total angular sweep is `2*pi*revolutions + angle`.

The solver searches the minimum between consecutive universal-variable
singularities using bounded golden-section minimization, then brackets and bisects
the requested side. The minimum search stops at relative universal-parameter
interval width `1e-12`; the reported minimum time is a numerical estimate.
`maximumIterations` bounds the minimum search and root solve separately;
`maximumBracketExpansions` bounds the selected root's bracketing. Enumeration
applies those budgets per solve. The solution reports revolution count, selected
root, universal parameter, minimum time, minimum-search iterations, a coalescence
flag and the existing root/time-residual diagnostics. Below-minimum flights beyond
the requested time tolerance produce `belowMinimumFlightTime` with the estimated
minimum; within that tolerance the minimum orbit is returned and explicitly
flagged as coalesced. Minimum-search and bracketing failures remain explicit errors.

The interval/branch structure follows the universal-variable formulation in
[Vallado's routine](https://github.com/CelesTrak/fundamentals-of-astrodynamics/blob/main/software/matlab/lambertu.m);
[Izzo's discussion](https://arxiv.org/abs/1403.2705) describes the multiple-root
and minimum-time structure. This implementation uses local anomaly reduction
between singularities so completed-turn parity does not alter the geometry.
Tests cover both roots, both directions, short/long-way motion and revolution
counts 1, 2, 5 and 32. Independent circular and eccentric Kepler solutions verify
velocities; propagated endpoint states and sampled angular sweeps verify every
returned root and its revolution count. Minimum-time coalescence, nearby distinct
roots, infeasible flight times and work limits are also checked.

### Moving-body Lambert transfer plans

`LambertTransferPlanner.plan` accepts explicit endpoint states, checked
`MovingWaypointRouteTarget` sources or `EphemerisProvider` sources. The moving
sources are sampled at departure and arrival respectively, so the destination is
resolved at arrival rather than frozen at departure. Coverage, exact state epochs,
frame identity and physical origin are checked; provider errors propagate.
Both providers must be relative to the selected central body in the same inertial
frame. For the direct-state overload the caller guarantees those physical origins.
No automatic rotation, origin shift or apparent-position correction is applied.

The caller explicitly chooses `intercept`, `flyby` or
`velocityMatchedRendezvous`, together with the Lambert direction, revolution count
and root where required. The departure impulse is the Lambert departure velocity
minus the sampled source velocity. The coast is independently propagated to
arrival. Intercept and flyby retain that incoming velocity and have no arrival
impulse. Rendezvous uses target velocity minus the actual propagated incoming
velocity, preserving actual position. The plan reports source/target states,
pre/post-arrival states, relative velocity vectors, actual impulse norms and their
sum, endpoint errors and whether the selected endpoint condition meets the caller's
tolerances. Intercept/flyby require position matching; rendezvous additionally
requires velocity matching. A residual position miss is retained and cannot be
made feasible by an arrival velocity impulse.

`LambertTransferTrajectory` covers only departure through arrival. Its default
queries are right-continuous at impulses, with `.before` available at each
endpoint. It retains the supplied physical pre-departure state, the transfer coast
and the post-arrival impulse state; queries outside the flight throw. No preceding
spacecraft history or following encounter orbit is inferred. This also supports
stationary or radial target velocities at rendezvous without inventing an orbit.
The plan's `maneuvers` can seed `FiniteTransferPlanner` using an actual spacecraft
state, engine and mass.

When departure/arrival sources represent body centres, relative velocity vectors
are ideal patched-conic excess-velocity inputs. Their norms are central-frame
velocity-difference costs, not local parking-orbit injection/capture costs. Use the
body's own GM and parking radius with `ParkingOrbitTransfer` for those local energy
changes. Matching a body-centre velocity does not establish a body-bound parking
orbit. A flyby here describes the incoming encounter endpoint; deflection, periapsis
clearance, B-plane targeting and outgoing encounter arcs require the separate
flyby/patch models in the backlog. Offset providers may instead specify explicit
physical passage or rendezvous locations.

Swift Testing cases use independent circular endpoint velocities and moving linear
sources to verify arrival sampling, actual impulses, mode-dependent costs and
velocity retention. They also check ephemeris metadata integration, both
multi-revolution roots, stationary rendezvous, residual preservation, bounded
trajectory queries and coverage/frame/origin errors.

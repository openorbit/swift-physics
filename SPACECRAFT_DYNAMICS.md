# Spacecraft dynamics foundations

The initial interfaces describe inertial translational dynamics in SI units,
using explicit epochs and existing Cartesian/frame identities. They evaluate
an ordinary differential equation's right-hand side, used by the adaptive
numerical propagator described below.

`SpacecraftMassProperties` holds positive finite dry mass and nonnegative finite
propellant mass in kilograms, with a finite total. Dry mass includes everything
outside the one represented propellant reservoir and remains constant in this
initial model. Inertia tensors, centres of mass and rotational evolution belong
to the attitude/rotational-dynamics milestone.

`SpacecraftState` combines a `CartesianState` with instantaneous mass properties.
Position and velocity retain their epoch, frame and SI conventions. A Cartesian
acceleration attached to that state is optional observational data: dynamics
computes acceleration from force models and does not add that field again.

`ForceEvaluationContext` names the frame's physical origin and optionally carries
an explicitly selected `EphemerisCollection`. The caller must ensure the frame
actually represents that origin. Environmental body queries must be converted
explicitly into compatible frames/origins; the context does not relabel or
transform ephemeris states. Coverage and data errors propagate from providers.
No global clock, body catalogue, implicit GM or ephemeris fallback is selected.

A `ForceModel` evaluates `contribution(for:in:)` at a supplied spacecraft state,
returning a `ForceContribution` carrying that same epoch and frame. Its fields
are acceleration in m/s² and propellant mass rate in kg/s, per atomic SI second.
Gravity models can supply acceleration directly; thrust models divide force by
the current total spacecraft mass. This avoids forcing gravitational acceleration
to be represented as a spacecraft-mass-dependent force in newtons.

Negative mass rates represent propellant consumption; positive rates represent
explicitly modelled inflow. The derivative contract supports both. Evaluation
does not mutate mass, clip consumption, or silently disable thrust. Models must
handle their operating constraints explicitly. The integrator rejects negative
propellant stages; scheduled constant-flow burns additionally support the analytic
cutoff policy described below. This single-reservoir interface does not yet model
multiple tanks, staging or changing dry mass.

```swift
import CelMek
import Math

struct ExampleThrust: ForceModel {
    func contribution(for state: SpacecraftState,
                      in context: ForceEvaluationContext) throws -> ForceContribution {
        // Application-defined thrust along the state's X axis, not spacecraft axes.
        let enabled = state.mass.propellantMassKilograms > 0
        return try ForceContribution(
            epoch: state.epoch, frame: state.frame,
            accelerationMetresPerSecondSquared: enabled
                ? SIMD3(100, 0, 0) / state.mass.totalMassKilograms : .zero,
            propellantMassRateKilogramsPerSecond: enabled ? -0.1 : 0)
    }
}

let frame = Frame(name: "Earth-centred inertial", pseudoInertial: true)
let state = SpacecraftState(
    cartesian: try CartesianState(epoch: .j2000, frame: frame,
        position: SIMD3(7_000_000, 0, 0), velocity: SIMD3(0, 7_500, 0)),
    mass: try SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50))
let dynamics = try SpacecraftDynamics(
    models: [ExampleThrust()],
    context: ForceEvaluationContext(frame: frame, origin: .earth))
let derivative = try dynamics.derivative(at: state)
// Position rate (0, 7500, 0) m/s, velocity rate (1, 0, 0) m/s²,
// propellant rate -0.1 kg/s. This example deliberately has no gravity model.
```

`SpacecraftDynamics` evaluates the models in declared order, validates each
contribution's epoch/frame, and sums acceleration and propellant rate. Its
`SpacecraftDerivative` carries position rate (the state's velocity), velocity
rate (summed acceleration) and propellant mass rate. Empty models mean unforced
motion with constant mass. Finite contribution validation and overflow checks
prevent an invalid derivative from escaping. Model errors remain visible; no
partial derivative or hidden fallback is returned.

Construction requires a pseudo-inertial context frame; state queries must use
that exact frame identity. Rotating-frame integration and its inertial terms are
later work. Existing frame transforms can still be used to view states in other
frames. Model evaluators must be deterministic for fixed inputs and tolerate
repeated/out-of-order epochs, as an adaptive solver may revisit or reject a step.
They must not advance an internal clock, consume propellant, or otherwise commit
simulation state during evaluation. Mutable caches are permissible only if they
do not change results; thread safety remains each provider/model's responsibility.

Swift Testing covers state-dependent thrust acceleration, additive gravity-like
acceleration, absence of double-counting, unchanged mass inventories, unforced
motion, invalid mass/contributions, frame/epoch contract failures, model errors
and overflow. Point-mass gravity models are documented below; finite thrust,
additional force settings and propellant-event handling are follow-up work.

## Adaptive numerical propagation

`NumericalSpacecraftPropagator` integrates Cartesian position/velocity and
propellant mass with the native Swift Dormand–Prince embedded 5(4) method.
Dry mass stays constant. It implements `StateProvider` for the existing frame
conversion/sampling APIs; `spacecraftState(at:)` also returns mass, and
`propagation(at:)` returns that state with diagnostics.

```swift
let options = try NumericalPropagationOptions(
    absolutePositionToleranceMetres: 1e-3,
    absoluteVelocityToleranceMetresPerSecond: 1e-6,
    absolutePropellantToleranceKilograms: 1e-9,
    relativeTolerance: 1e-9,
    initialStepSeconds: 10,
    minimumStepSeconds: 1e-6,
    maximumStepSeconds: 300,
    maximumStepAttempts: 100_000,
    maximumDerivativeEvaluations: 700_001)
let propagator = try NumericalSpacecraftPropagator(
    initialState: state, dynamics: dynamics, options: options)
let result = try propagator.propagation(
    at: state.epoch.advanced(by: Duration(seconds: 60)))
// result.state includes mass and freshly evaluated acceleration.
// result.diagnostics reports the work and accepted-step error estimates.
```

Each query restarts from the same immutable seed. Backward queries integrate
with negative time steps, including the signed mass derivative; this is a replay
of the prescribed equations, not a reversal of committed game events. Repeating
or reordering queries does not change the trajectory. Sampling many output epochs
currently performs a separate integration for each epoch, trading predictable
semantics for repeated work. There is no dense interpolation or trajectory cache.
Rendering cadence neither advances the provider nor fixes its adaptive steps.

The fifth-order solution is accepted when every state component's normalized
embedded error is at most one. Each scale is its positive absolute tolerance
plus `relativeTolerance * max(abs(previous), abs(candidate))`; position,
velocity and propellant have separate SI tolerances. The relative tolerance may
be zero. These are **local numerical error controls**, not physical model accuracy
or a guaranteed final-position error. Absolute position is measured in the
selected frame; at a large barycentric origin a relative tolerance can permit
large position errors, so select a suitable origin and component tolerances.
The embedded error is formed directly from weighted stage derivatives rather
than subtracting two large endpoint states.

Initial/minimum/maximum step sizes are positive finite seconds with
`minimum <= initial <= maximum`. The final step is shortened to land on the exact
requested epoch and may be smaller than the minimum. Rejected trials reduce the
step; unmet tolerances at the minimum fail explicitly. Accepted/rejected attempts
and derivative evaluations each have bounded work limits. Every RK stage uses an
explicit epoch, and a step whose interior epochs cannot be represented distinctly
fails rather than looping or querying collapsed stages.

Diagnostics include accepted/rejected step counts, derivative evaluations,
smallest/largest accepted step sizes (nil for a zero-duration query), and maximum
accepted normalized error. Seven derivative evaluations are used per complete
trial plus one final evaluation; a rejected negative-mass trial can terminate
earlier. A reference-epoch query performs one derivative evaluation, no steps,
and returns seed position/velocity/mass with computed acceleration rather than
trusting the seed's optional acceleration field.

Trial states with negative propellant are rejected before force evaluation.
Step reduction may recover a valid trial; continued depletion at the minimum
fails with `propellantDepleted`. No propellant is clipped, no engine is implicitly
disabled by this generic propagator, and no partial final state is returned. General depletion localization,
staging, impacts, impulsive burns and other discontinuity events remain future
work. Segment integration explicitly at known force discontinuities; the adaptive
error estimate alone does not guarantee that a discontinuity is discovered.
Nonfinite arithmetic, exhausted work limits and epoch-resolution failures are
explicit `NumericalPropagationError` cases. Model, frame, environmental-coverage,
mass-validation and time-arithmetic errors propagate without provider fallback.

Tests compare polynomial motion/linear mass change, an exponential velocity
solution and central-gravity motion against the existing independently validated
two-body propagator. They check tighter-tolerance convergence, forward/backward
integration, rejection behavior, exact requested epochs, final-step shortening,
query-order independence, work limits and depletion/provider failures.
The Dormand–Prince coefficient table is also documented in
[SciPy's RK45 source](https://github.com/scipy/scipy/blob/v1.16.0/scipy/integrate/_ivp/rk.py);
no external numerical library is required at runtime.

Shared accuracy hints can select local-error presets via
`NumericalSpacecraftPropagator(initialState:dynamics:accuracy:options:)`.
Explicit options override the preset; the original direct-options initializer
remains available. The selected settings and requested/preset metadata are
inspectable. See [accuracy policy](ACCURACY.md) for the versioned preset table
and its interpretation; force/environment fidelity stays independently selected.

## Point-mass gravity

`CentralGravity(body:gravitationalParameter:)` supplies `-μ r/|r|³` for a
central body located at the context's declared origin. Its body identifier must
match `context.origin`. It does not query an ephemeris for an origin already at
zero, infer GM from mass, or apply surface-radius softening. Its sourced GM must
be Newtonian or TT-compatible, matching the atomic-second dynamics convention;
raw TDB/TCB-compatible GM values require an explicit justified conversion first.

`ThirdBodyGravity(body:centralBody:gravitationalParameter:)` supplies a relative
perturbation in nonrotating axes centred on the named central body:

```text
a = μ [(R-r)/|R-r|³ - R/|R|³]
```

Here `r` is spacecraft position relative to the central body, and `R` is the
perturber's position relative to that same body. The second, indirect term
removes the perturber's acceleration of the origin. This is not an absolute
point-mass attraction model for arbitrary barycentric propagation, nor does it
infer the actual acceleration of an independently prescribed origin. The caller
must choose a physically consistent central-body relative model. A barycentre
cannot simply replace that central body while retaining the same indirect term.

```swift
let central = try CentralGravity(body: .earth, gravitationalParameter: earthGM)
let lunar = try ThirdBodyGravity(
    body: .moon, centralBody: .earth, gravitationalParameter: moonGM)
let dynamics = try SpacecraftDynamics(
    models: [central, lunar],
    context: ForceEvaluationContext(
        frame: earthCentredInertialFrame, origin: .earth,
        ephemerides: earthRelativeProviders))
// Integrate with NumericalSpacecraftPropagator as above.
```

The third-body model resolves its target through the injected collection on each
query, then uses the provider's checked state at the spacecraft epoch. Provider
origin must equal the central body and frame identity must equal the propagation
frame. Convert native DE/VSOP/lunar states and their metadata explicitly first;
there is no implicit translation or axes relabeling. Collections can contain
explicit fallback/cache adapters, but gravity adds no separate fallback policy.
Missing collections, unsupported bodies, bad frame/origin metadata and uncovered
queries remain explicit errors. Target and central body cannot be identical.

Both models return acceleration independent of spacecraft mass, with zero
propellant rate, and enforce frame identity and pseudo-inertial axes even when
queried directly outside `SpacecraftDynamics`. Body identifiers preserve planet
centre versus system-barycentre distinctions; supply a GM appropriate to the
represented mass and geometry. No collision radius, gravity harmonics, tidal
response, relativistic correction or spacecraft self-gravity is added.

When `|r|/|R| <= 0.1`, third-body evaluation avoids subtracting two nearly equal
large accelerations. With `u=r/|R|`, `n=R/|R|`, `s=|n-u|` and
`q=u·(u-2n)`, it rationalizes `s³-1` as
`f=q*(s²+s+1)/(s+1)` and evaluates `-(μ/|R|²)*(u+f*n)/s³`.
Larger relative separations use direct acceleration subtraction, preserving
precision near the perturber. A spacecraft at the origin gives exactly zero
third-body perturbation after the provider's coverage/contract checks. The
perturber at the origin or the spacecraft coinciding with an attracting point
mass is singular and fails explicitly. Finite input vectors can still exceed
supported floating arithmetic range; nonfinite intermediate/results fail rather
than being softened or extrapolated.

Independent 80-digit Decimal fixtures evaluate direct Newtonian terms and their
difference (`Tools/generate_gravity_fixtures.py`, `Fixtures/gravity-reference.json`).
They include a tiny distant-body tide that would lose precision under ordinary
subtraction, as well as central, lunar-distance and general 3D cases. The numerical
orbit/convergence tests now exercise production `CentralGravity` against the
analytic two-body propagator. Further tests check moving provider epochs,
force composition, mass independence, origin/frame/GM errors and singularities.
Combined central/third-body propagation is also tested with a moving injected
provider, unchanged mass and acceleration consistent with the summed derivative.
The direct/indirect distinction is also described in
[Orekit's third-body gravity API](https://www.orekit.org/site-orekit-7.1/apidocs/org/orekit/forces/gravity/ThirdBodyAttraction.html).

## J2 gravity

`J2Gravity` adds only the axisymmetric degree-two correction. Compose it with
`CentralGravity` for central attraction plus oblateness; do not add it to another
model already including degree two/order zero. The signed dimensionless coefficient
is J2 = −C20 for **unnormalized** C20, not a normalized harmonic coefficient.
This convention follows [Orekit's J2-only model](https://javadoc.io/static/org.orekit/orekit/12.2.1/org/orekit/forces/gravity/J2OnlyPerturbation.html).

Supply the body identifier, Newtonian or TT-compatible GM, a sourced equatorial
radius in metres, a finite J2 with its provenance, and an `OrientationProvider`.
No body constants are inferred. Zero and negative J2 are supported. The context
origin must be that body, and propagation requires a pseudo-inertial frame.
The orientation provider and its result must use the exact propagation frame
and requested epoch; provider failures propagate. An ICRF orientation cannot
be silently relabeled to a different frame.

The provider maps reference axes into body axes with the pole along body z.
For unit body position n, distance r, and radius R, the correction is
`1.5 J2 (GM/r²) (R/r)² [nx(5nz²−1), ny(5nz²−1), nz(5nz²−3)]`.
Rotate this acceleration vector back into propagation axes without kinematic
rotation terms. Spin around a fixed pole leaves this axisymmetric force unchanged;
a changing pole changes the force. Mass and propellant remain unchanged.

This is an exterior gravity expansion, with no surface collision model, radius
cutoff, higher harmonics, tides, or relativistic corrections. The caller must
select a physically valid domain. The origin is singular, even for zero J2;
nonfinite arithmetic fails explicitly.

Swift Testing checks equatorial/polar values, independent numerical gradients of
the potential, signed/zero coefficients, tilted poles and spin invariance, and
central-plus-J2 integration. A ten-orbit weak-J2 case agrees within 1% with the
[first-order secular nodal regression](https://ai-solutions.com/_freeflyeruniversityguide/j2_perturbation.htm)
and checks conservation of angular momentum about the pole. This validates the
force and integrator; it does not validate a particular body's chosen constants.

## Finite thrust and scheduled burns

`FiniteThrust` prescribes constant thrust in newtons, positive effective exhaust
velocity in m/s, and a normalized direction in an explicit pseudo-inertial frame.
The effective exhaust velocity incorporates pressure effects in the selected
engine operating condition. Acceleration is `direction * thrust / totalMass`,
and propellant rate is `-thrust / effectiveExhaustVelocity`. These follow the
[NASA thrust and specific impulse convention](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/).
As fuel burns, constant force produces increasing acceleration. This model has
no implicit gravity, attitude steering, atmosphere correction, throttle history,
staging, or multiple reservoirs. A zero-thrust engine consumes no fuel; an active
engine with an empty reservoir throws `propellantDepleted`.

`FiniteBurn` specifies a positive-duration interval. `ScheduledThrust` activates
burns on `[start, end)`. Adjacent burns switch at their shared epoch; overlapping
burns represent independent engines and add both thrust acceleration and fuel
consumption. Each burn can supply different engine parameters and direction.
The schedule is immutable and evaluates arbitrary epochs without consuming fuel
or changing hidden state. All engines must share the exact propagation frame.

Use `ScheduledThrustPropagator` for numerical schedule execution. It divides
forward or backward propagation at every ignition/cutoff boundary and integrates
with the segment's active engines held through its RK endpoint evaluations.
It then reports acceleration according to the right-continuous schedule at the
requested epoch. Position, velocity, dry mass and propellant carry continuously
between segments. Short burns cannot disappear between adaptive steps, even if
both their endpoints lie between requested samples. Every query restarts from
the immutable seed; backward replay adds back the prescribed consumed fuel and
reconstructs earlier motion, rather than reversing engine direction.

```swift
let engine = try FiniteThrust(
    frame: initialState.frame,
    thrustNewtons: 20,
    effectiveExhaustVelocityMetresPerSecond: 100,
    direction: SIMD3(1, 0, 0))
let ignition = try initialState.epoch.advanced(by: Duration(seconds: 10))
let cutoff = try ignition.advanced(by: Duration(seconds: 20))
let schedule = try ScheduledThrust(frame: initialState.frame, burns: [
    FiniteBurn(start: ignition, end: cutoff, thrust: engine)
])
let trajectory = try ScheduledThrustPropagator(
    initialState: initialState, schedule: schedule,
    backgroundModels: [centralGravity], context: context,
    options: .preset(for: .high))
let result = try trajectory.spacecraftState(at: cutoff)
```

The schedule accepts at most 10,000 burns. Numerical options and work budgets
apply separately to each segment; there is no aggregate step-diagnostics API yet.
Background models must be continuous between scheduled boundaries. Other force
events still require explicit handling. Passing `ScheduledThrust` directly to
the generic numerical propagator does **not** provide boundary splitting.

The default `FuelDepletionPolicy.fail` preserves explicit exhaustion failures,
including potentially at exact scheduled cutoffs. Select `.cutoff` to localize
constant-flow scheduled depletion and continue under background forces, as
specified below. The standalone `FiniteThrust` and generic numerical propagator
retain their failure semantics.

Swift Testing checks an independently integrated constant-thrust rocket solution
for velocity, position and fuel, short burns followed by coast, exact boundary
acceleration, overlapping/adjacent engines, reverse recovery, query determinism,
invalid inputs, frame mismatches and depletion failures.

## Drag and solar radiation pressure

`AtmosphericDrag` uses the isotropic [drag equation](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-equation/):
`a = -rho * Cd * A * |vRelative| * vRelative / (2 * mass)`.
Area is in m²; Cd is dimensionless. Both are finite and nonnegative; zero disables
the force. `AtmosphereProvider` supplies a local `AtmosphereSample` containing
nonnegative density in kg/m³ and medium velocity in m/s at the requested epoch.
The sample must have the exact spacecraft frame and context origin. Winds and
atmospheric rotation belong in that velocity, expressed in propagation axes;
zero medium velocity explicitly describes a stationary atmosphere. Provider
coverage/errors propagate. `UniformAtmosphere` supports controlled scenarios;
it is not an Earth atmosphere or altitude/density approximation. No empirical
planetary density model is selected automatically.

`SolarRadiationPressure` uses a point source and an isotropic constant area with
a dimensionless radiation coefficient, following the single-coefficient approach
used by [Orekit](https://www.orekit.org/site-orekit-12.2.1/apidocs/org/orekit/forces/radiation/IsotropicRadiationSingleCoefficient.html).
Supply a sourced reference pressure in pascals, sourced reference distance in
metres, area in m², and coefficient. Acceleration points away from the source,
with magnitude `pressure * (referenceDistance/distance)² * coefficient * area /
mass * illuminatedFraction`. If the source is the context origin its position
is zero; otherwise its checked ephemeris must have the same physical origin and
exact frame. No automatic origin conversion or fallback is selected. The source
position is evaluated at the spacecraft epoch without light-time correction.
Coincidence with the source is singular, even for zero illumination.

An explicit `RadiationIlluminationProvider` supplies a finite fraction in 0...1
at each state. `ConstantRadiationIllumination(illuminatedFraction: 1)` means full
illumination; zero means shadow and intermediate values attenuate incident flux.
Injected providers can use ephemerides and geometry for eclipses. No built-in
umbra/penumbra geometry or occultation detection is included, and discontinuous
shadow transitions require event handling or known-boundary segmentation. A
constant fraction alone makes no claim that a path is free from eclipses.

Both forces use instantaneous total mass, consume no propellant, and require
pseudo-inertial propagation. They do not model attitude-dependent projected
area, lift, aerodynamic heating, optical panels, albedo, infrared radiation,
or relativistic effects. The caller selects coefficients and valid domains.

```swift
let medium = try UniformAtmosphere(
    frame: initialState.frame, origin: context.origin,
    densityKilogramsPerCubicMetre: 1e-12,
    velocityMetresPerSecond: SIMD3(0, 0, 0))
let drag = try AtmosphericDrag(
    referenceAreaSquareMetres: 10, dragCoefficient: 2.2, atmosphere: medium)
let radiation = try SolarRadiationPressure(
    referencePressurePascals: pressureAtReferenceDistance,
    referenceDistanceMetres: referenceDistance,
    referenceAreaSquareMetres: 10, radiationCoefficient: 1.3,
    illumination: ConstantRadiationIllumination(illuminatedFraction: 1))
let dynamics = try SpacecraftDynamics(
    models: [centralGravity, drag, radiation], context: context)
```

Swift Testing checks relative-flow drag vectors and an independent analytic
quadratic-drag trajectory, radiation distance/mass/illumination scaling and
source direction, explicit ephemeris origins, invalid environments and provider
contracts. Synthetic parameters test the equations without asserting planetary
environment fidelity.

## Event detection

`SpacecraftStateProvider` extends `StateProvider` with mass-aware queries; both
numerical propagators conform. `EventDetectedPropagation` searches any such
trajectory from an explicit start epoch to a target, forwarding provider errors
and checking returned epoch/frame contracts. It reports occurrences and a final
state, whether the target was reached, and the number of state evaluations.
It is deliberately not a `StateProvider`: a terminal event can end before the
requested epoch, and that result must not masquerade as a target-epoch state.

A `PropagationEventDetector` supplies a unique nonblank identifier, a continuous
finite scalar value, crossing direction, and `.record` or `.stop` action.
Directions describe increasing physical time, even during backward searches.
`EpochEvent` localizes an explicit epoch exactly, useful for maneuver checkpoints.
`SphericalImpactEvent` uses distance from an explicit fixed centre minus radius,
with decreasing crossings selected by default. Its frame must match exactly;
the caller assigns the physical origin and radius. Moving bodies, planes,
clearance margins and illumination thresholds can use custom detectors.

```swift
let impact = try SphericalImpactEvent(
    identifier: "surface", frame: trajectory.frame, radiusMetres: bodyRadius)
let search = try EventDetectedPropagation(
    trajectory: trajectory, detectors: [impact],
    options: PropagationEventOptions(
        maximumCheckIntervalSeconds: 1, timeToleranceSeconds: 1e-6))
let result = try search.propagation(from: initialState.epoch, to: requestedEpoch)
// result.state may be at an impact, before requestedEpoch.
```

Search samples at bounded time intervals and bisects detected sign changes,
querying the actual trajectory at every trial epoch. There is no linear state
interpolation or dense-output assumption. Occurrences are sorted in traversal
order; the first terminal occurrence stops the result. Equal localized epochs
use detector input order. Initial exact zeros are reported by default regardless
of crossing direction; `includeInitialEvents: false` suppresses these. Left-end
zeros are not emitted again on the next interval. Only events at or before a
terminal occurrence in traversal order are returned.

`PropagationEventOptions` bounds check interval, time bracket tolerance,
localization iterations, state queries and occurrence count. Localization reports
the final bracket width, not a physical state-error guarantee. Provider numerical
tolerances remain independent. Nonfinite values, exhausted budgets, insufficient
epoch resolution and failed localization throw explicit errors; no partial result
is returned. All detector evaluations must be pure and may occur out of order.

This is sign-crossing detection: an even number of crossings within one check
interval, a tangency, an initially interior state, or a function identically zero
requires a caller-specific policy. Choose check spacing to separate crossings;
this is not a certified continuous collision-clearance algorithm. The provider
and every detector must remain valid throughout sampled/localized brackets.
A terminal search may evaluate past the root while forming its bracket, so it
cannot prevent a singular force or depleted-fuel provider error there. A positive
fuel-reserve detector can be used where the trajectory is still valid, but exact
zero-fuel localization for arbitrary mass-flow models requires integration-level
handling. Constant-flow scheduled burns support the analytic cutoff policy below.

Event actions report or terminate; they do not mutate forces or automatically
apply maneuvers. Restart from the reported state with explicit changed inputs to
model a discontinuity. Scheduled ignition/cutoff splitting remains handled by
`ScheduledThrustPropagator`. Underlying queries restart from their immutable seed,
so this general search is intended as a correctness-first reference implementation,
not an optimized accepted-step observer. Numerical work limits apply per provider
query in addition to the search's state-query budget.

Swift Testing validates analytic straight-line sphere crossings, backward crossing
direction, record/terminal ordering, initial-root suppression, exact epoch events,
scheduled-burn mass state, invalid detector contracts and bounded search failure.

## Automatic scheduled fuel cutoff

Set `fuelDepletionPolicy: .cutoff` on `ScheduledThrustPropagator`. The default is
`.fail` for compatibility. For every segment, sum the active engines' positive
constant consumption rates q. If the forward segment would exhaust available
propellant P, split at `segmentStart + P/q`, integrate thrust up to that epoch,
set the known boundary fuel to exactly zero, and continue with background forces.
All scheduled engines draw from this one reservoir and shut down together;
later ignitions also remain off. `propagation(at:)` returns the state and optional
`depletionEpoch`, including depletion coincident with a scheduled cutoff.

Fuel evolves analytically within each constant-flow segment. Cartesian integration
supplies that exact inventory to every background force evaluation, so drag,
radiation pressure and other mass-dependent forces see evolving total mass.
The RK reservoir coordinate is held fixed internally rather than generating
negative fuel stages. Endpoint stage evaluations use the smooth thrust left
limit at dry mass; the public returned state uses zero fuel and thrust-free
right-limit acceleration. This is a localized boundary condition, not clipping
arbitrary negative inventories. Numerical position/velocity tolerances still
apply; the exhaustion time inherits floating-point P/q precision and epoch
resolution. An unrepresentable positive-duration depletion boundary fails.

```swift
let trajectory = try ScheduledThrustPropagator(
    initialState: initialState, schedule: schedule,
    backgroundModels: [centralGravity], context: context,
    fuelDepletionPolicy: .cutoff)
let result = try trajectory.propagation(at: requestedEpoch)
// result.depletionEpoch reports reservoir exhaustion; result.state may be a coast state.
```

Cutoff mode requires background models to return **zero propellant mass rate**.
Unexpected background consumption or inflow throws `backgroundMassFlowUnsupported`;
it cannot silently invalidate the known depletion time. Variable throttle,
refueling, staging and multiple reservoirs remain unsupported; continuous variable
consumption uses `FuelLimitedNumericalPropagator` below.
An initially empty reservoir coasts forward with scheduled engines disabled.
Backward queries from the immutable wet seed integrate the prescribed mass flow
in reverse. A new empty seed cannot infer historical cutoff/engine activity, so
backward traversal through a positive-flow burn fails explicitly. Querying from
the original wet seed after depletion remains deterministic.

Swift Testing compares exhaustion velocity/position and subsequent coast against
the rocket equation, verifies exact scheduled cutoff, overlapping engines,
empty seeds, later suppressed ignitions, evolving mass in background forces,
unsupported background mass flow, and missing backward-history failures.

## Accepted-step events and explicit transitions

`EventDrivenNumericalPropagator` checks detector values at error-controlled
accepted step endpoints. It reuses those states; only a candidate crossing
requires localization queries. Maximum step size is capped by the event check
interval independently of the numerical tolerance preset. Localization restarts
from the crossing step's left state, avoiding full-seed replay for every root.
The accepted endpoint past a root is discarded; integration restarts at the
localized state with the transition's force composition.

Each detector is **one-shot per query by default**. Optional rearming thresholds
are described below. Initial exact roots are included; fired
identifiers are disabled for the rest of that query, preventing repeated resets
at the same root. Recurring events opt into a caller-selected nonzero rearming band.
Equal exact epochs use detector input order; remaining detectors are rechecked
at every reset epoch. `.stop` returns immediately without applying a transition.
A `.record` occurrence continues and optionally invokes the transition closure.
Returned occurrences contain the pre-transition state. A reset at the target
is applied and its resulting acceleration recomputed before returning.

The closure receives `(occurrence, forward, currentDynamics)` and returns an
`EventStateTransition(state:dynamics:)`. It can alter position/velocity/mass
explicitly or replace the composed forces. State epoch/frame and dynamics
frame/origin must match the event's propagation context. The caller supplies
inverse reset behavior for backward replay; arbitrary changes are not implicitly
inverted. Closures and detectors must be deterministic without external mutation.
Each query starts from the original seed and original dynamics.

```swift
let checkpoint = try EpochEvent(
    identifier: "impulse", epoch: burnEpoch, action: .record)
let trajectory = try EventDrivenNumericalPropagator(
    initialState: initialState, dynamics: dynamics, detectors: [checkpoint]
) { event, forward, currentDynamics in
    let state = event.state
    let velocity = state.cartesian.velocity + (forward ? deltaV : -deltaV)
    let reset = SpacecraftState(
        cartesian: try CartesianState(epoch: state.epoch, frame: state.frame,
            position: state.cartesian.position, velocity: velocity), mass: state.mass)
    return EventStateTransition(state: reset, dynamics: currentDynamics)
}
let result = try trajectory.propagation(at: requestedEpoch)
```

Known discontinuities must be changed by the reset closure, rather than hidden
inside a force's epoch-dependent branch. For example, use a continuous pre-event
force set and return the post-event force set at an `EpochEvent`. Prescribed
constant-flow engine schedules still use `ScheduledThrustPropagator` and its
specialized fuel-cutoff handling.

The same sign-crossing limits apply: tangencies and multiple crossings inside
one check interval are not guaranteed. All forces must remain valid through the
trial step and crossing bracket, so arbitrary zero-fuel events cannot yet recover
negative fuel stages or singular-force evaluations. Variable-flow depletion
uses the dedicated consuming-force boundary contract described below.
Numerical work limits apply per restart/localization query; event query, accepted
step-check and occurrence counts are bounded across the complete event run.
`stateEvaluations` reports localization/search provider queries, not RK derivative
evaluations. General accepted-step events have no aggregated RK diagnostics yet.

Swift Testing verifies a terminal plane crossing, an impulse with inverse replay,
two force-switch events against analytic motion, deterministic replay, initial
one-shot roots, and rejected resets with incompatible epochs. Existing numerical
propagation tests verify unchanged behavior when no observer is installed.

## Variable-flow depletion and recurring events

`VariableThrust` accepts a pure command closure `(state, context) -> FiniteThrust`.
A command can change throttle, effective exhaust velocity and inertial direction
continuously with epoch/state. Engine validation and frame checks apply at every
evaluation. Commands must be defined at dry mass for depletion localization and
must not mutate hidden state. Known discontinuities still use explicit resets.

`FuelConsumingForceModel` marks a model as drawing from the shared reservoir and
provides `depletionBoundaryContribution`. For positive fuel this must equal normal
force evaluation; at zero fuel it provides the continuous dry-mass limit for root
localization. This method must return a nonpositive mass rate. `FiniteThrust` and
`VariableThrust` conform. Normal evaluation still rejects firing with zero fuel.
Unmarked background forces must return exactly zero mass rate; incompatible
consumption/inflow throws rather than silently changing the depletion model.

`FuelLimitedNumericalPropagator` integrates the coupled position, velocity and
variable fuel inventory. During controlled boundary trials, the private reservoir
coordinate may go negative while force models see zero propellant and dry total
mass. This extension exists only to bracket depletion; negative inventory is
never exposed in a returned spacecraft state. An accepted crossing triggers
bounded bisection through reintegration from the valid side. The reported depletion
epoch is the last valid side of a time bracket no wider than `timeToleranceSeconds`;
its bounded residual inventory is replaced with the known zero-fuel boundary.
`FuelLimitedPropagationResult` reports the state, optional depletion epoch and
final bracket width. All consuming models then switch off together and the
background forces continue.
No arbitrary invalid stage is accepted or silently clipped in the generic solver.

```swift
let engine = try VariableThrust(frame: initialState.frame) { state, context in
    let elapsed = try state.epoch.duration(since: initialState.epoch).seconds
    return try FiniteThrust(frame: state.frame,
        thrustNewtons: 10 * (1 + elapsed),
        effectiveExhaustVelocityMetresPerSecond: 100,
        direction: SIMD3(1, 0, 0))
}
let powered = try SpacecraftDynamics(models: [centralGravity, engine], context: context)
let trajectory = try FuelLimitedNumericalPropagator(
    initialState: initialState, dynamics: powered,
    eventOptions: PropagationEventOptions(timeToleranceSeconds: 1e-8))
let result = try trajectory.propagation(at: requestedEpoch)
```

The boundary time has numerical integration error in addition to bracket width;
a time tolerance does not certify a fuel/state error bound. Use appropriate
mass/position/velocity tolerances. Event options bound localization iterations
and all integration queries (initial search, trials, coast); RK work limits apply
per query. The check-interval and occurrence-count event options are not used by
this dedicated single, monotonic depletion search. Continuous consumption with a
well-defined dry-mass limit is required. Refueling, staging, multiple tanks and
other simultaneous state resets are outside this propagator. Backward queries
use the original ordinary dynamics and require an evaluable wet seed; a dry seed
cannot reconstruct burn history. The generic solver's depletion failure behavior
is unchanged.

For recurring general events, supply `rearmingThresholds: [identifier: threshold]`
to `EventDrivenNumericalPropagator`. Thresholds must be finite and positive and
identify registered detectors. They use each detector's scalar units, so the
caller chooses a meaningful band (metres for clearance, for example). After firing,
a detector stays disabled until an accepted step starts outside that band, then
can trigger a later crossing. Reset states are not immediately rearmed, and the
same root cannot repeatedly fire while remaining in the band. Check intervals
must still separate roots and resolve the excursion outside the band. Detectors
without a threshold retain one-shot behavior. Occurrence and evaluation limits
bound repeated events, and backward searches use the same physical-time direction.

Swift Testing validates a variable-thrust exhaustion against analytic fuel/time,
rocket-equation velocity and independent Simpson quadrature of position, post-cutoff
coast, unchanged pre-depletion states, deterministic replay, recurring crossings
forward/backward, invalid bands and unsupported background mass flow.

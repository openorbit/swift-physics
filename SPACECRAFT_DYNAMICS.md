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
handle their operating constraints explicitly; the future integrator must enforce
nonnegative propellant and handle depletion as an event rather than integrating
through a negative reservoir. This single-reservoir interface does not yet model
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
disabled, and no partial final state is returned. Exact depletion localization,
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

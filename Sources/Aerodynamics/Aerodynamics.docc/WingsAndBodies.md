# Segmented wings and whole-body loads

Evaluate aerodynamic loads for supplied wing or capsule geometry and coefficient data.

These evaluators remain pure aerodynamic physics. They do not integrate motion,
assemble vehicles, manage control systems or calculate reentry heating. All body
vectors use the same caller-selected origin and X-forward/Y-right/Z-down axes.

## Wing strips

``SegmentedWing`` aggregates ``WingSegment`` loads about a supplied moment reference.
A rectangular helper defaults to eight strips; other counts and arbitrary segment
geometry are supported. ``WingSectionAxes`` defines the segment's chord/span/down
basis, allowing twist, dihedral, sweep or a vertical tail orientation through explicit
orthonormal directions. Positive incidence raises the leading edge.

Each strip samples body-relative-to-air velocity at its own reference location,
projects that velocity into the chord plane and evaluates its section coefficients.
This independent-strip approximation discards spanwise velocity: dynamic pressure,
Mach and Reynolds number use the chord-plane speed. It is not a complete swept-wing
or finite-wing model.

Positive lift is perpendicular to chord-plane velocity, toward section-up; drag
opposes that velocity. Positive CM rotates about section +span. Force is `q S C`;
section moment is `q S chord CM`. The segment position must match the polar's moment
reference, such as quarter chord. Summation adds the lever-arm moment exactly once.
See [NASA aerodynamic-center conventions](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/aerodynamic-center/).

```swift
import Aerodynamics
let ambient = try ISAAtmosphere().properties(at: 0)
let flow = try FluidFlow(relativeVelocityMetresPerSecond: SIMD3(40, 0, 0), properties: ambient)
// `airfoil` is a caller-supplied AirfoilCoefficientModel with appropriate coverage.
let wing = try SegmentedWing.rectangular(spanMetres: 8, chordMetres: 2,
    incidenceRadians: 0.08, model: AirfoilWingSectionModel(airfoil: airfoil))
let result = try wing.evaluate(flow: flow,
    angularVelocityRadiansPerSecond: SIMD3(0.1, 0, 0))
```

The uniform-flow convenience adds `ω×r`, with offsets from an explicit velocity
reference point. A separate sampler overload accepts local flow for each segment,
allowing atmosphere gradients, winds and prescribed induced-flow corrections.
Results include each strip's section flow, coefficients and physical loads as well
as the total. Zero chord-plane dynamic pressure yields zero loads and nil coefficients
without querying an undefined incidence/Reynolds number.

``WingSectionCoefficientModel`` receives flow, chord and caller-supplied control
surface deflection. A calibrated provider may use deflection to select/interpolate
its aerodynamic data. ``AirfoilWingSectionModel`` adapts existing checked polars;
it rejects nonzero deflection and missing Mach/Reynolds inputs instead of ignoring
them. No generic flap effectiveness, stall or wake field is invented. The optional
<doc:FiniteWing> correction supplies equivalent uniform-wing induced angle/drag
under explicitly bounded conditions; it does not couple arbitrary strips.
Use a provider and dataset appropriate to the flight regime; this API does not
extend low-speed airfoil data to hypersonic entry automatically.

## Capsule and fuselage models

``WholeBodyAerodynamics`` consumes ``BodyAerodynamicCoefficientModel`` and reference
geometry. It uses full freestream dynamic pressure, scales all six body-axis
coefficients and transfers the resulting moment to the requested point if needed.
It never adds a second centre-of-pressure moment to a complete coefficient moment.

``BodyAerodynamicQuery`` carries flow, Reynolds reference length, moment reference
lengths and body angular velocity. The normalized angular rates are
`(p span, q chord, r span)/(2 V)` when speed is nonzero. A custom provider can use
these for calibrated damping. No damping is applied by the static table provider.

``TabulatedBodyCoefficients`` linearly interpolates a complete rectangular grid of
angle of attack, sideslip, Mach and Reynolds conditions. Up to 16 corner coefficients
contribute to a query. All axes default to rejecting out-of-coverage requests;
explicit clamping uses their nearest boundary values. A singleton axis covers its
single value rather than implying independence. Empty, duplicate, malformed or
incomplete grids are rejected, including during Codable decoding.

``BodyCoefficientTableLoader`` reads the Codable JSON table format from supplied data
or files; I/O is separate from evaluation. The JSON fields are `samples`,
`boundaryPolicy`, `source` and `version`. Samples contain `angleOfAttackRadians`,
`sideslipRadians`, `machNumber`, `reynoldsNumber` and `coefficients`; coefficients
contain three-element `force` and `moment` arrays. Source/version identify supplied
data. Geometry dimensions and moment references must match that dataset; the
caller is responsible for this compatibility.

``WindForceCoefficients`` converts positive drag, signed side force and positive
lift into body-axis force coefficients using incidence and sideslip. Wind X is along
the body's relative velocity, wind Y is right and wind Z is down; force components
in wind axes are `(-CD, CY, -CL)`. Moments are supplied separately in body axes,
with body-axis moment reference lengths. The conversion convention is consistent
with [NASA's component/velocity-axis definitions](https://rotorcraft.arc.nasa.gov/Publications/files/NDARCTheory_v1_6_938.pdf).

``ConstantBodyCoefficients`` provides explicitly fixed body-axis coefficients for
controlled scenarios. It does not represent wind-aligned drag at arbitrary incidence.
For zero dynamic pressure the whole-body evaluator returns zero loads without
querying coefficients. A table query with undefined incidence (pure spanwise flow)
or missing Mach/Reynolds properties fails explicitly.

No capsule-specific flight database is bundled. A lifting or ballistic capsule
must supply measured, computed or explicitly idealized data and geometry. Heat flux,
aerodynamic damping, rarefied-flow transitions and hypersonic surface-pressure
models remain separate workstreams.

## Runnable example and verification

```sh
swift run FluidPhysicsDemo
```

The example evaluates an eight-strip wing using Earth sea-level ISA and a lifting
capsule using ISA at 50 km. Its coefficients are synthetic demonstrations of the
interfaces, not validated flight or entry predictions.

Tests compare uniform strip loads with analytic `q S C` totals, validate moment
reference transfers, examine rotational local flow and check midpoint convergence
against an independently integrated spanwise drag distribution. Body tests use an
exact affine four-dimensional coefficient field, checked JSON round trips and
wind-axis transformations. These checks validate implementation conventions;
vehicle-specific physical accuracy requires independent aerodynamic datasets.

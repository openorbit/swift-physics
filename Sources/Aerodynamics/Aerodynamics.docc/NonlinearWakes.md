# Nonlinear sections and coupled wakes

Solve a bounded fixed-wake lifting-line problem across supplied aerodynamic strips.

``NonlinearLiftingLine`` accepts ``WingSegment`` geometry and the existing
``WingSectionCoefficientModel`` contract, including tabulated nonlinear/control
polars. Strips from separate wings or tails may share a solve. The caller samples
incident flow at each strip, including local rotation, atmosphere gradients,
external body interference or prescribed winds. No vehicle assembly is owned.

Each strip is a constant-circulation horseshoe: a bound vortex along the supplied
section span through its moment-reference point, with width `area/chord`, and two
straight finite trailing legs. Use quarter-chord geometry and a matching section
moment convention. Wake direction is physical air motion (normally body −X).
It must be perpendicular to each strip's span; the model supports unswept lifting
lines with incidence/dihedral and separated surfaces, rather than arbitrary swept
or intersecting lifting lines. Distinct control points alone do not establish a
well-resolved geometry. Avoid other bound/wake vortex cores and endpoints.

The solver sums regularized Biot–Savart air velocity, subtracts it from each
incident body-relative velocity, queries the section model and relaxes toward
`Gamma = 0.5 V_section chord CL`. It reports circulation, induced air velocities,
section results, loads, iteration count and maximum circulation residual normalized
by incident `0.5 V chord`. A failed solve throws `nonConvergence`; it never returns
an unconverged best guess. Multiple nonlinear roots or convergence from a zero
initial state do not establish uniqueness or post-stall physical accuracy.

```swift
let sections = try SegmentedWing.rectangular(spanMetres: 8, chordMetres: 1,
    model: sectionModel).segments
let solver = try NonlinearLiftingLine(segments: sections,
    wakeLengthMetres: 800, vortexCoreRadiusMetres: 0.0001)
let result = try solver.evaluate { segment in
    try incidentFlow.atLocalPoint(offsetMetres: segment.positionMetres,
        angularVelocityRadiansPerSecond: suppliedBodyRate)
}
let airWake = try solver.inducedAirVelocity(at: downstreamBodyPoint,
    circulation: result.circulationSquareMetresPerSecond)
let downstreamFlow = try incidentFlow.replacingRelativeVelocity(
    incidentFlow.relativeVelocityMetresPerSecond - airWake)
```

Both vectors must share axes, units and frame. `FluidFlow/replacingRelativeVelocity(_:)`
also works with body, panel and parachute providers. For mutually interacting
lifting surfaces, put all their strips in one solver; successive isolated solves
provide only one-way interaction. Body shielding, wake deficits and arbitrary
interference may be supplied through the flow sampler; they are not inferred.

Loads use the converged effective flow, so its tilted lift already produces induced
drag. Supply **section profile drag**; do not wrap the polar in ``FiniteWingAirfoil``
or add another induced-drag force. Section moments and force lever arms are summed
once about the requested reference point.

Limits are Mach ≤ 0.3, incident speed ≥ 1 m/s, incident alignment within 0.35 rad
of the upstream wake axis, incident spanwise speed ≤ 10% of speed, section incidence
within ±0.35 rad and |CL| ≤ 2. Provider coverage also applies at every iterate.
No 3-D separation, transonic shocks, dynamic stall, free-wake rollup or ground
images are inferred. A supplied post-stall polar remains an approximation of local
2-D behavior. Zero density/rest are outside this steady circulation solver's coverage.

At most 256 strips are accepted. Setup/memory cost O(N²), evaluation O(I N²).
Relaxation, residual tolerance, iteration cap, wake length and core radius are
explicit. Core radius is limited to 2% of the narrowest strip; this is a numerical
restriction, not a vortex-physics model. Refine strip count, increase wake length
and decrease core radius independently to assess numerical sensitivity.

Tests compare a single horseshoe against an independent semi-infinite-vortex/scalar
root solution, check 8/16/32/64-strip convergence against Fourier lifting line, and
exercise nonlinear polar response, two-way wing/tail coupling, nonuniform flow,
reference-point transfer and explicit failures.

Background: [CN-AeroModels nonlinear lifting-line formulation](https://lheea.gitlab.io/CN-AeroModels/nonlinear_lifting_line.html).

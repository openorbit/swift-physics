# Coupled spanwise lifting line

``LiftingLineWing`` solves a straight, unswept planar wing in uniform incident
flow. Unlike independent strips or an equivalent whole-wing correction, its
circulation and induced incidence are coupled across the span. It accepts supplied
linear section polars with spanwise chord, incidence and zero-lift-angle variation.

Use it for attached, low-speed flow: aspect ratio ≥ 4, Mach ≤ 0.3, zero sideslip
within numerical tolerance, speed ≥ 1 m/s and small angles. Every station must
satisfy its declared effective-incidence and Reynolds coverage. Lift slopes are
positive and bounded to 0.01…20 per radian; geometric/zero-lift angles and declared
effective-incidence bounds stay within ±0.35 rad. These implementation limits do
not constitute validated accuracy bounds for a particular wing.

## Geometry and section inputs

The reference line lies along body Y, with a supplied centre in body axes. All
section moments refer to that same straight line; positions of swept or staggered
quarter-chord lines are not represented. The constructor samples
`sectionAtSpanFraction` at cosine-spaced stations, from left (−1) to right (+1).
``LiftingLineSection`` supplies chord, lift slope, zero-lift angle, incidence,
constant profile drag, pitching moment and coverage. Profile drag excludes induced
drag. Supplied control/camber changes may alter these parameters externally; no
control-effectiveness factor or actuator state is inferred.

Reference area sets coefficient normalization and must match the geometry/dataset.
Rectangular area is span times chord; elliptical area is π span rootChord / 4.
Stations do not include the tips, where circulation vanishes. Chord/span ratios at
sampled stations must lie within 1e-6…1. Station count is configurable from 2 to 64,
with eight as the initial default rather than an accuracy claim.

```swift
let wing = try LiftingLineWing(spanMetres: 8,
    referenceAreaSquareMetres: 4 * .pi, stationCount: 16,
    source: "Ideal elliptical wing", version: "1") { fraction in
    try LiftingLineSection(chordMetres: 2 * sqrt(1 - fraction * fraction),
        liftSlopePerRadian: 2 * .pi, profileDragCoefficient: 0.02,
        effectiveAngleRangeRadians: -0.2...0.2, reynoldsRange: 1e3...1e8)
}
let result = try wing.evaluate(flow: flow)
```

The example parameters are idealized. Flow supplies sound speed and viscosity for
coverage checks, and is already relative to the atmosphere in the stated body axes.
This model does not sample local rotation, varying incident speed or wind gradients.
Those require a more general lifting-line/wake formulation; the independent-strip
sampler remains available for externally prescribed local flow.

## Solution and loads

The sine-series collocation system is solved with a pivoted factorization cached
at construction. Setup costs O(N³) and each query O(N²). ``LiftingLineResult`` returns
Fourier coefficients, per-station circulation/effective incidence, lift, profile
and induced drag coefficients, span efficiency and maximum equation residual.
Numerically zero total lift has undefined span efficiency, even when asymmetric
circulation produces nonzero induced drag.

Lift, induced drag and the spanwise lift moment use Fourier integrals. Profile drag,
section pitching moments and the spanwise drag moment use numerical quadrature;
increase station count to check their convergence. Local induced-drag allocation
uses the small-angle lifting-line approximation. The returned force/moment are
about an explicit point and can feed CelMek's aerodynamic wrench adapter.

Do not add another finite-wing correction or induced-drag load for the same wing.
The linear section data does not supply nonlinear stall, dynamic stall, compressibility,
surface interference, multiple-wing wakes or a validated control database.
`FiniteWingAirfoil` remains the cheaper equivalent uniform-wing option.

Tests reproduce elliptical lift slope and minimum induced drag at several station
counts, an analytic second-harmonic twist/rolling moment, rectangular convergence,
constant profile drag, reference-point transfer and unsupported coverage.

References: [classical lifting-line Fourier formulation](https://www.mdpi.com/2076-3417/15/6/3383)
and [NASA induced-drag relationship](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/induced-drag-coefficient/).

# Hypersonic surface pressure

`SurfacePressureIntegration` integrates gauge pressure over supplied flat
`FluidSurfacePanel` elements. Each panel has area in m², an outward normal,
a centre in metres and a caller-supplied exposure flag. Geometry and flow use the
same body axes. Subdivision is caller-selected; evaluation costs O(panel count).
It returns per-panel pressure coefficients and gauge pressures plus force and
moment about an explicit point, ready for CelMek's aerodynamic wrench adapter.

`SurfacePressureCoefficientModel` permits replacement pressure laws. The initial
``ModifiedNewtonianPressure`` approximation uses a calorically perfect gas with
explicit specific-heat ratio γ. It computes downstream normal-shock Mach and static
pressure, then isentropically decelerates to obtain stagnation pressure. The law is:

`CpMax = (pStagnationBehindShock/pFreeStream − 1)/(γ Mach²/2)`

`Cp = CpMax max(0, outwardNormal · vehicleRelativeVelocityDirection)²`

Vehicle-relative velocity points into the windward outward normal. The pressure
force is inward: `F = −normal Cp q area`. Panel moments use the force lever arm;
ambient pressure is subtracted. For an open surface, these are aerodynamic gauge
loads, not total structural pressure loads. Leeward, tangent and explicitly
unexposed panels have zero pressure coefficient in this approximation.

```swift
let panel = try FluidSurfacePanel(areaSquareMetres: 1,
    centreMetres: SIMD3(2, 0, 0), outwardNormal: SIMD3(1, 0, 0))
let pressureModel = try ModifiedNewtonianPressure(specificHeatRatio: 1.4)
let surface = try SurfacePressureIntegration(panels: [panel],
    referenceLengthMetres: 1, model: pressureModel)
let flow = try FluidFlow(relativeVelocityMetresPerSecond: SIMD3(1000, 0, 0),
    densityKilogramsPerCubicMetre: 0.02, speedOfSoundMetresPerSecond: 100,
    meanFreePathMetres: 0.001)
let result = try surface.evaluate(flow: flow)
```

This synthetic example supplies γ and sound speed consistently with the chosen gas
approximation. Thermodynamic profiles do not yet supply mean free path; it must be
provided independently. The model requires known Mach and Knudsen numbers for
nonzero speed, defaults to Mach 5…25 and Kn ≤ 0.01, and allows narrower caller
coverage within Mach 5…50. Reference length must be positive and should represent
the relevant body or local geometry scale. These are implementation coverage
bounds, not guarantees of accuracy or proof that perfect-gas behavior applies.
Zero speed returns zero aerodynamic pressure load without requiring gas properties.
Unsupported coverage and missing properties throw even for unexposed panels.

Use the panel-flow sampling overload to account for vehicle and atmospheric
rotation, winds or gradients at each panel centre. Uniform-flow evaluation assumes
the same flow everywhere. No geometry shadowing or shock interaction is inferred.
Normals are normalized during construction; outward orientation and correct areas
remain the caller's responsibility.

The approximation is intended for idealized blunt-body windward pressure. It does
not model viscous skin friction, base suction, expansion pressures, heating,
dissociation or variable heat capacity. Add separately supported contributions
without counting pressure drag twice; do not silently extrapolate it into another
flight regime or treat constant γ as a real-gas entry model.

Tests check a normal-shock stagnation reference, inclined-panel force and lever-arm
moment, coverage and rest behavior, and convergence to the analytic unit-sphere
result `CD = CpMax/2` based on projected area. This validates the implementation of
the approximation, not a particular vehicle's flight accuracy.

The formulation follows [NASA normal-shock relations](https://www.grc.nasa.gov/WWW/BGH/normal.html)
and the modified Newtonian coefficient discussed in
[NASA TR R-127](https://ntrs.nasa.gov/citations/19630006549).

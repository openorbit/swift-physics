# Local surface convection

Evaluate bounded laminar or fully turbulent flat-plate heat transfer separately from capsule stagnation heating.

``SurfaceConvectiveHeatingModel`` accepts local boundary-layer-edge flow, distance
along the streamline from the leading edge, a supplied wall temperature and
``AerothermalTransportProperties``. The latter carries gas identity/provenance and
film-temperature density, viscosity, conductivity and specific heat. No body-specific
transport defaults or chemistry are inferred. Film temperature is `(Te + Tw)/2`;
the evaluator checks the supplied property temperature.

``LaminarFlatPlateHeating`` implements the isothermal, zero-pressure-gradient,
attached parallel-flow approximation `Nu_x = 0.332 sqrt(Re_x) Pr^(1/3)`, with
`Re_x = rho_f V_edge x / mu_f`, `Pr = mu_f cp_f / k_f` and `h_x = Nu_x k_f / x`.
The edge-flow density/viscosity are not used in place of film properties.
``SurfaceConvectionResult`` retains signed `q = h_x (Te - Tw)`: positive heats the
wall and negative cools it. This differs deliberately from the incoming-only
``AerodynamicHeatFlux`` used by stagnation correlations. Do not clamp cooling to
zero or add this flux to a stagnation correlation at the same point.

Implementation bounds are Mach ≤ 0.3, local Reynolds 1,000–500,000, Prandtl
0.6–15, local Knudsen ≤ 0.01 and wall/edge temperature ratio 0.8–1.25.
These conservative limits are restrictions, not an accuracy guarantee or a
transition detector. The caller must establish an attached laminar boundary layer,
parallel edge flow and approximately uniform wall temperature and properties.
Shock interactions, pressure gradients, sweep, separation, automatic turbulent transition,
compressible recovery temperature and high-enthalpy chemistry need other models.
Zero speed is outside forced-convection coverage; no natural-convection estimate
is substituted. The singular leading edge at x=0 is rejected.

```swift
let transport = try AerothermalTransportProperties(
    gasIdentifier: "supplied-air", source: "caller property dataset",
    temperatureKelvin: 315, densityKilogramsPerCubicMetre: 1.12,
    dynamicViscosityPascalSeconds: 1.92e-5,
    thermalConductivityWattsPerMetreKelvin: 0.0274,
    specificHeatJoulesPerKilogramKelvin: 1007)
let edge = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(10, 0, 0),
    densityKilogramsPerCubicMetre: 1.07, temperatureKelvin: 330,
    speedOfSoundMetresPerSecond: 364, meanFreePathMetres: 7e-8)
let result = try LaminarFlatPlateHeating().evaluate(SurfaceHeatingQuery(
    edgeFlow: edge, distanceFromLeadingEdgeMetres: 0.1,
    surfacePositionMetres: SIMD3(0.1, 0, 0),
    wall: AerothermalWallState(temperatureKelvin: 300), transport: transport))
```

Evaluation costs O(1) per point. Callers sample a surface distribution and integrate
flux over area externally; no surface mesh, thermal evolution or assembly is owned
by the model. Position is an output label, not a substitute for streamline distance.
This establishes a surface-convection extension point and a low-speed reference
model. Use <doc:EntryHeating> for supplied high-enthalpy surface data and a bounded
swept-cylinder stagnation-line approximation.

The correlation and independent numerical similarity coefficients are given in
[NPTEL External Flows, lecture 8, table 2.2](https://archive.nptel.ac.in/content/storage2/courses/112104159/module2/Lec8.pdf).
Tests compare those coefficients, integrate the local distribution over a bounded
plate interval and exercise heating/cooling, coverage and property requirements.


``TurbulentFlatPlateHeating`` explicitly selects a fully turbulent boundary layer
from the leading edge, with `Nu_x = 0.0296 Re_x^0.8 Pr^(1/3)` and Reynolds
500,000–10,000,000. It retains the same Mach, Kn, Pr and temperature-ratio limits,
film-property requirements and signed flux. It does not locate transition or model
an initially laminar region. The different spatial dependence is tested alongside
heating/cooling and out-of-range rejection. Source:
[NPTEL External Flows, lecture 12](https://archive.nptel.ac.in/content/storage2/courses/112104159/module2/Lec12.pdf).

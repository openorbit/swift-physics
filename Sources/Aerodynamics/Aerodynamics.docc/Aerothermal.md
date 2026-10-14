# Stagnation-point heating

``ConvectiveHeatingModel`` and ``RadiativeHeatingModel`` accept a
``StagnationHeatingQuery`` and return incoming heat flux in W/m² at an explicit
surface position. They are separate contributions: no radiative model means
radiation is unavailable, not zero. There is no built-in radiative correlation yet.
Wall radiation emission, temperature evolution, conduction and ablation belong to
external thermal models. These interfaces do not advance a clock or mutate vehicle
state. Richer gas chemistry or catalytic-surface data can be supplied through a
replacement model's provider; these initial contracts are for stagnation heating.

``SuttonGravesHeating`` supplies a laminar cold-wall convective estimate for an
axisymmetric blunt body:

`qConvective = k sqrt(ρ / Rnose) V³`

Density is kg/m³, nose radius is metres, relative speed is m/s and output is W/m².
The coefficient `k` has SI dimensions kg^0.5/m. Its numerical value depends on
the gas and the coefficient's unit convention; a W/cm² coefficient cannot be used
unchanged with the W/m² interface.

``SuttonGravesParameters`` requires a gas identifier, coefficient, source and
version. ``SuttonGravesParameters/earthAirColdWall`` explicitly selects the Earth-air
value `1.7415e-4` used in the NASA TFAWS 2012 course. No planet or atmosphere lookup
selects it automatically. The query's gas identifier must match the parameters;
this is a dataset identity check, not composition analysis. Custom parameter sets
allow independently justified approximations for other gas mixtures.

```swift
let coverage = try StagnationHeatingCoverage(
    speedMetresPerSecond: 1000...8000,
    densityKilogramsPerCubicMetre: 1e-6...0.01,
    noseRadiusMetres: 0.1...5, machNumber: 5...40)
let model = SuttonGravesHeating(parameters: .earthAirColdWall, coverage: coverage)
let flow = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(3535, 0, 0),
    densityKilogramsPerCubicMetre: 3.1459e-4, speedOfSoundMetresPerSecond: 300,
    meanFreePathMetres: 1e-4)
let query = try StagnationHeatingQuery(flow: flow, gasIdentifier: "earth-air",
    noseRadiusMetres: 1, surfacePositionMetres: SIMD3(1, 0, 0))
let heat = try model.heatFlux(for: query)
```

These example bounds are declared simulation coverage, not experimentally validated
limits. The caller must establish that the correlation suits the gas, geometry and
entry conditions. The model restricts declared Mach coverage to hypersonic flow
(Mach ≥ 5) and continuum coverage to Kn ≤ 0.01, using nose radius as the reference
length. Missing sound speed or mean free path throws for nonzero speed. Speed,
density, radius, Mach and Knudsen coverage are checked before evaluating the
correlation. Zero speed produces zero aerodynamic convective flux after checking
gas identity, wall assumptions and radius, without requiring Mach or mean free path.
Zero density at nonzero speed is outside the declared positive-density coverage.
Low density alone does not establish continuum applicability.

The stagnation point and effective nose radius are caller-supplied. The model does
not locate a stagnation point or extrapolate a capsule value over a wing. Local flow
must already account for atmospheric co-rotation and winds, as described in
<doc:Foundations>. Use consistent thermodynamic and transport inputs; current
atmosphere profiles do not supply mean free path automatically.

This cold-wall approximation does not account for finite wall temperature. A query
with ``AerothermalWallState`` therefore throws `unsupportedQuery`; a different
correlation must use that state. No turbulent multiplier, catalytic correction,
radiative contribution or species-equilibrium solver is inferred. The estimate
may be integrated externally as heat load in J/m², but does not itself predict
wall temperature or protection-system survival.

See <doc:PlanetaryHeating> for a Mars preset, explicit Venus mixture parameters
and a separate hydrogen/helium cold-wall correlation with different wall assumptions.

Tests reproduce the NASA course's rounded shuttle-like reference of 13.6 W/cm²,
check cubic velocity, square-root density and inverse-square-root radius scaling,
and exercise gas mismatch, custom coefficients, rest and coverage boundaries.

Sources: [Sutton and Graves, NASA TR R-376](https://ntrs.nasa.gov/citations/19720003329)
and [NASA TFAWS 2012 Aerothermodynamics Course](https://tfaws.nasa.gov/TFAWS12/Proceedings/Aerothermodynamics%20Course.pdf).

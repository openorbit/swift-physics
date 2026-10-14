# Planetary heating parameters

Heating parameters are explicit gas datasets, independent of the atmosphere
catalogue. Density and temperature profiles alone do not identify a heating model.
Use a query gas identifier matching the chosen parameter set.

## Mars

``SuttonGravesParameters/marsCarbonDioxideNitrogenColdWall`` supplies
`k = 1.83e-4` for an idealized 97% CO₂ / 3% N₂ mixture **by mass**. Use it with
``SuttonGravesHeating`` and caller-declared coverage. It is the simple cold-wall
estimate, not the newer CFD correlation in the same paper.

Source: [West and Brandis, 2020, equation 2](https://ntrs.nasa.gov/api/citations/20200002354/downloads/20200002354.pdf).

## Venus and explicit mixtures

``SuttonGravesMixtureParameters`` accepts mass or mole fractions of nine base
gases. It derives the enthalpy-form coefficient using equation 42 and the rounded
transport products in table I. Mole fractions are converted to mass fractions.
``SuttonGravesMixtureParameters/venusCarbonDioxideNitrogen`` assumes 96.5% CO₂ /
3.5% N₂ by mole, ignoring trace species. This is a derived mixture approximation,
not a fitted Venus flight dataset.

``SuttonGravesMixtureHeating`` requires ``StagnationGasState`` in the query:
stagnation pressure in Pa and boundary-layer-edge specific enthalpy in J/kg.
It uses `q = K sqrt((pStag/101325)/Rnose) hEdge`, neglecting wall enthalpy.
The coefficient's source units contain atm and differ from the velocity-cubed
coefficient; the evaluator performs the pressure conversion. Wall state must be
absent for this cold-wall estimate. Flow and continuum coverage are also checked.
Pressure and edge enthalpy are limited to the original study's 0.001…100 atm and
2.3…116.2 MJ/kg ranges. These bounds do not establish the cold-wall approximation's
accuracy for every state. The model does not solve equilibrium chemistry or infer
shock-layer properties from an ISA-like profile.

Sources: [Sutton–Graves report](https://ntrs.nasa.gov/citations/19720003329) and
[NASA Venus composition](https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html).

## Hydrogen and helium

``HydrogenHeliumColdWallHeating`` implements the separate nonlinear cold-wall
fit. Its baseline is 81.55% H₂ / 18.45% He by mole, with a fully catalytic,
non-ablating spherical wall at 300 K and freestream temperature 141.2 K. A missing
wall adopts 300 K; other wall temperatures and inconsistent freestream temperature
are rejected. It requires Mach ≥ 5, Kn ≤ 0.003, radius 0.2…1 m and density
1e-5…0.1 kg/m³. Speed coverage follows a log-density interpolation of table-2
boundary points rather than a full 7…45 km/s rectangle.

``HydrogenHeliumHeatingMixture`` provides the Neptune baseline and Saturn/Uranus
approximations. Their factors 0.952/0.980 originate in hot-wall composition studies
and were also applied by the authors to cold-wall comparisons. They are not
independent cold-wall fits. The model preserves the published branch switch at
29 km/s, including an approximately 5–9% downward step for supported radii near
that boundary; it does not smooth or extrapolate it.
Diagnostics expose fitted pressure, enthalpy difference, K₂ and selected branch.
Tests compare both flux and K₂ with all 180 Appendix-B CFD cases within 20%; this
checks implementation against the dataset, not arbitrary trajectories.

Source: [Carroll and Brandis, 2022](https://ntrs.nasa.gov/api/citations/20220018610/downloads/AIAA_Aviation_2022_Paper-6.pdf).

## Example

```swift
let mixture = HydrogenHeliumHeatingMixture.neptuneReference
let model = HydrogenHeliumColdWallHeating(mixture: mixture)
let flow = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(19000, 0, 0),
    densityKilogramsPerCubicMetre: 1e-4, temperatureKelvin: 141.2,
    speedOfSoundMetresPerSecond: 1000, meanFreePathMetres: 1e-5)
let query = try StagnationHeatingQuery(flow: flow, gasIdentifier: mixture.gasIdentifier,
    noseRadiusMetres: 0.4, surfacePositionMetres: SIMD3(1, 0, 0),
    wall: AerothermalWallState(temperatureKelvin: 300))
let result = try model.evaluate(query)
```

Example sound speed and mean free path are supplied assumptions. No radiative flux
is included. Jupiter, Titan and tenuous-moon heating presets remain unsupported.
Current nominal atmosphere profiles also do not establish upper-atmosphere or
transport coverage for these entry calculations. See <doc:Aerothermal> for the
external thermal-integration boundary.

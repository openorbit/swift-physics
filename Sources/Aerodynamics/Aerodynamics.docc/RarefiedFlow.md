# Free-molecular surface loads

`SurfaceTractionModel` returns signed gauge pressure and tangential stress in Pa.
`SurfaceTractionIntegration` integrates these over supplied surface panels,
including force lever arms and moment-reference transfer. It validates that shear
is tangential to each panel. `PressureOnlyTraction` adapts existing continuum
pressure laws into the same contract with zero shear; it is an alternative evaluator,
not an additional force to count alongside pressure-only integration.

## Gas and surface inputs

``MaxwellianFreeMolecularTraction`` describes independent impacts from an ideal
single-species Maxwellian gas. ``FreeMolecularGas`` supplies specific gas constant
in J/(kg K), identity and provenance. Flow supplies density, temperature, local
vehicle-relative-to-air velocity and mean free path. Sound speed and viscosity are
not required. Ambient pressure is computed as `ρ R T`; any supplied flow pressure
must be consistent with those inputs. A mixture-average R is an additional
approximation, not a species-resolved mixture treatment.

``GasSurfaceReflection`` sets the specular fraction and wall temperature. Specular
reflection transfers normal momentum without shear. The remaining fraction is fully
accommodated diffuse re-emission at wall temperature, with outgoing molecular flux
balancing incoming flux. This reflection fraction is not a general independent
normal/tangential accommodation law, and is not inferred from material name.

For outward normal `n`, let `U = vehicleVelocity · n`, `σ = sqrt(R T)` and `a = U/σ`.
Incoming mass and normal momentum fluxes follow Gaussian half-space moments:

`J = ρ σ [φ(a) + a Φ(a)]`

`Pincident = ρ σ² [(a² + 1) Φ(a) + a φ(a)]`

The diffuse outgoing pressure is `J sqrt(π R Twall / 2)`. The model combines that
pressure with specular reflection, subtracts ambient pressure, and supplies diffuse
shear opposite the tangential vehicle velocity. Thermal molecules can hit leeward
faces, so they are not discarded by an incidence sign test. Deep leeward Gaussian
tails use an asymptotic expansion to avoid subtraction cancellation.

## Coverage and geometry

The evaluator requires Kn ≥ 10 using the supplied reference length and mean free
path. Missing temperature or mean free path throws, including at rest; low density
does not imply vacuum or select a regime. A separate numerical bound requires
`abs(U/σ) ≤ 1e6`. Zero speed does not necessarily mean zero pressure load when wall
and gas temperatures differ. At equal temperatures, stationary gas produces zero
gauge pressure and shear. A uniform closed body at rest has no net pressure force.

Normals point outward. Gauge pressure may be negative on leeward panels; absolute
surface pressure is obtained by adding `ρ R T`. Gauge subtraction removes the
local ambient contribution, so these loads do not supply a separate hydrostatic
buoyancy force. Sample consistent flow and any needed buoyancy separately.

The exposure flag is a caller-supplied inclusion mask. A false flag skips the panel's
traction; it does not solve molecular occlusion, reflected-particle reimpact or
thermal-view geometry. Independent panels suit externally exposed geometry;
concavities and nearby interacting surfaces need a more complete model. All flow,
stress and location vectors use the same axes, with atmospheric co-rotation already
included in relative flow. The panel sampling overload supports local rotation and
external gradients.

## Example

```swift
let gas = try FreeMolecularGas(identifier: "Supplied single-species gas",
    specificGasConstantJoulesPerKilogramKelvin: 300, source: "Caller gas dataset")
let surface = try GasSurfaceReflection(specularFraction: 0.2, wallTemperatureKelvin: 300)
let model = try MaxwellianFreeMolecularTraction(gas: gas, surface: surface)
let panel = try FluidSurfacePanel(areaSquareMetres: 1,
    centreMetres: SIMD3(1, 0, 0), outwardNormal: SIMD3(1, 0, 0))
let evaluator = try SurfaceTractionIntegration(panels: [panel],
    referenceLengthMetres: 1, model: model)
let flow = try FluidFlow(relativeVelocityMetresPerSecond: SIMD3(7000, 0, 0),
    densityKilogramsPerCubicMetre: 1e-8, temperatureKelvin: 200, meanFreePathMetres: 20)
let loads = try evaluator.evaluate(flow: flow).loads
```

These are demonstration inputs, not an atmosphere or material preset. No chemical
reactions, ablation, heat flux or gas inventory is evolved. Do not extrapolate either endpoint model into slip/transitional flow.
``TabulatedRarefiedTraction`` accepts a supplied consistent transition dataset
through the same integration contract; no universal endpoint blend is inferred.

Tests compare molecular fluxes with independent Simpson quadrature, equal-temperature
rest and hot-wall limits, thermal leeward impacts, specular high-speed pressure,
force/moment transfer and equivalence of the continuum pressure adapter.

References: [NASA free-molecular pressure/shear formulation](https://ntrs.nasa.gov/api/citations/19890015753/downloads/19890015753.pdf)
and [NASA rarefied-flow regime definitions](https://www.nas.nasa.gov/assets/nas/pdf/ams/2022/AMS_20220802_Singh.pdf).


## Supplied transitional traction

``TabulatedRarefiedTraction`` uses matched pressure/shear `FluidScalarTable`
axes, ordered `log10Kn`, `incidenceCosine`, `speedRatio`, `wallTemperatureRatio`,
all dimensionless. Incidence cosine is outward normal dotted with body-relative
velocity direction; speed ratio is `V / sqrt(R T)`, and temperature ratio is
`Tw / T`. Pressure values are gauge pressure divided by `rho R T`; nonnegative
shear magnitudes use the same normalization and oppose tangential relative velocity.
This normalization differs from conventional dynamic-pressure coefficients.
Normals/rest have no tangential direction and must have zero interpolated shear;
inconsistent data throws. Nonzero molecular thermal pressure at rest can remain.

Gas R and wall temperature are explicit; source/version/assumptions must identify
the gas composition, surface accommodation, geometry and data method (for example
DSMC or a calibrated bridge). A complete rectangular dataset may include continuum,
slip, transition and free-molecular conditions, with continuous multilinear output
inside its stated coverage. The model never changes mean free path to trick the
asymptotic evaluators into accepting unsupported states. Missing temperature/path
or an out-of-range coordinate throws; no density-based vacuum fallback exists.
Caller exposure flags remain the only shielding input.

Tests use a synthetic affine dataset to verify log-Kn interpolation/continuity,
normalization, shear direction and integrated moments. These tests do not validate
a physical transition law; no calibrated transition dataset ships by default.

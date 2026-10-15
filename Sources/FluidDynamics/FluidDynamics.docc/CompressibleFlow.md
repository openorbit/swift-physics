# Conservative compressible flow

Capture shocks and expansion waves with a separate energy-conserving Euler formulation.

``CompressibleFlowGrid`` evolves density, three momentum-density components and
total energy density on ``CartesianFluidGrid`` geometry. It uses an unsplit,
first-order finite-volume method with Rusanov/local Lax–Friedrichs face fluxes.
The same interface flux updates both adjacent fluid cells. Its numerical diffusion
is substantial at coarse resolution, especially across contact discontinuities;
this is a bounded initial solver rather than a high-resolution aircraft database.

``CompressibleEquationOfState`` provides pressure, temperature and acoustic speed
from density/internal energy, and an inverse pressure initializer.
``CaloricallyPerfectGasEOS`` supplies a constant-R/gamma ideal gas with explicit
identifier, source and temperature coverage. Every reconstructed/updated state is
checked against the EOS. There is no energy/pressure/density clipping to hide a
failed state. A failed EOS, positivity check or substep budget leaves the original
grid unchanged. `workBudgetExceeded` is distinct from physical `outsideCoverage`.

```swift
let gas = try CaloricallyPerfectGasEOS(identifier: "synthetic-gas",
    specificGasConstantJoulesPerKilogramKelvin: 1, specificHeatRatio: 1.4,
    temperatureRangeKelvin: 1e-4...1e4, source: "Sod numerical benchmark")
let geometry = try CartesianFluidGrid(dimensions: SIMD3(160, 1, 1),
    minimumCornerMetres: .zero, cellSizeMetres: 1.0 / 160)
var solver = try CompressibleFlowGrid(configuration:
    CompressibleFlowConfiguration(grid: geometry, equationOfState: gas)) { p in
    try ConservedFluidState(densityKilogramsPerCubicMetre: p.x < 0.5 ? 1 : 0.125,
        velocityMetresPerSecond: .zero, pressurePascals: p.x < 0.5 ? 1 : 0.1,
        equationOfState: gas)
}
let report = try solver.advance(timeStepSeconds: 0.2, maximumSubsteps: 10000)
```

Use the labelled `initialState:` closure if also supplying `solidAtCellCentre:`.
Numerical benchmark units are synthetic SI; the 1 K state is not an Earth-air preset.
Dimensions may include one-cell directions for 1-D/2-D reference problems, with a
maximum of 262,144 total cells. Sampling is piecewise constant. Physical continuum
applicability, microscopic transport and material validity remain caller/provider
responsibilities; positive density is not proof of continuum coverage.

The acoustic CFL condition uses a conservative sum of maxima `|u_d| + c` over all
three axes, including prescribed boundary states. The Courant factor is at most
0.45 (default 0.4). Sound waves therefore impose work even at zero bulk speed, and
hypersonic flow can require many substeps for a rendering/vehicle timestep.
The caller sets a work cap and update cadence; there is no automatic real-time
fallback. Cost is O(cells × substeps), with O(cells) storage, and up to three face
sweeps per substep. This initial implementation is serial CPU Swift.

## Boundaries and coupled loads

Low/high boundaries are supplied in X/Y/Z order. Periodic boundaries must be
paired; fixed-state boundaries prescribe a ghost state; transmissive boundaries
copy the adjacent state and are not a nonreflecting far field. Reflecting boundaries
reverse normal momentum while retaining density and energy: stationary adiabatic
slip walls. Static solid voxels use the same reflected face flux, and may not touch
the outer cell layer. They do not represent moving geometry or no-slip walls.

``CompressibleStepReport`` returns net outward mass/momentum/energy integrals and
step-average solid reaction loads from the **same** face momentum flux used by the
fluid update. Lever arms are at voxel-face centres about the requested reference.
Use that contribution once; do not also add a pressure/whole-body force for the
same contribution. Tests balance fluid momentum change, solid impulse and boundary
flux, and verify periodic mass/energy conservation. General hybrid coupling and
angular-momentum accuracy still depend on discretization and model ownership.

The grid owns no vehicle state or clock. `state(at:)` returns density, velocity,
pressure, temperature and sound speed. `relativeFlow(bodyVelocityMetresPerSecond:)`
constructs shared flow inputs after subtracting the sampled medium velocity, with
optional separately supplied viscosity and mean free path. Positions/velocities
must share a nonrotating inertial or constant-translating frame; moving grids and
rotating/accelerating-frame source terms are not inferred.

## Entry use and limits

This solver provides compressible pressure/shock physics. The energy equation
includes shock-induced internal-energy changes, but the initial Euler formulation
has no viscous stress, thermal conduction, radiation transport or wall heat flux.
A reflected adiabatic boundary does not establish zero physical reentry heating.
Constant gamma is not a high-enthalpy chemistry model. More elaborate equilibrium
EOS providers can use the closure; finite-rate species/vibrational chemistry needs
additional evolved states and source terms, not just another EOS implementation.

Aerodynamics retains the separate bounded stagnation/leading-edge/heating-data
providers. Stagnation correlations require their declared freestream conditions;
do not feed an already compressed stagnation cell to a freestream correlation and
count the shock twice. Local enthalpy/Stanton data needs the appropriate boundary-
layer-edge state and its own geometry/transport coverage. No automatic coupling
claims to turn an unresolved voxel cell into a validated boundary-layer heat flux.

Tests compare the Sod solution with independently specified exact star/fan/shock
states, demonstrate grid convergence, preserve a uniform supersonic state, check
conservative solid reaction and exercise budget/positivity/EOS failures.
Run `swift run FluidPhysicsDemo --compressible` for the reference harness.

Background: [Clawpack Euler equations and Riemann solvers](https://www.clawpack.org/riemann_book/html/Euler_approximate.html).

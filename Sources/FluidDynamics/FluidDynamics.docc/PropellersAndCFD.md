# Optional propwash and local CFD

Choose a prescribed slipstream or evolve a local incompressible flow field, and sample either through the existing aerodynamic evaluators.

## Shared velocity convention

``FluidVelocityField`` returns physical **fluid** velocity in caller-selected common
axes and frame. ``InducedFluidVelocityField`` marks an increment; the actuator-disk
provider conforms to it. ``SuperposedFluidVelocityField`` adds increments to one
background. A ``IncompressibleFlowGrid`` returns total fluid velocity, including its boundary
background, and is deliberately not an induced-increment provider.

FluidFlow stores body motion relative to the medium, the opposite sign for a
stationary body in a moving-air grid. Subtract fluid velocity from body velocity in
the same frame before expressing the difference in section/body axes:

```swift
let localFlow = try ambientFlow.replacingRelativeVelocity(
    bodyVelocityInGridAxes - field.velocityMetresPerSecond(at: surfacePoint))
```

Do not add a total CFD velocity to the same background again. Atmosphere rotation,
vehicle translation and frame conversion remain explicit as in [fluid conventions](/documentation/fluiddynamics).

## Actuator-disk propwash

``ActuatorDiskSlipstream`` accepts supplied positive thrust, signed fluid torque,
radius, centre, downstream axis, density and nonnegative aligned upstream speed.
For disk area A it solves `T = 2 rho A vi (V + vi)`, with disk speed `V + vi`, far
wake speed `V + 2vi` and mass flow `rho A (V + vi)`. Downstream radius follows
constant mass flow as axial velocity increases. These are ideal axial momentum
relations; thrust and torque are inputs, not calculated from blade geometry or RPM.

Intermediate downstream development uses the explicitly chosen approximation
`vi(z) = vi_disk (1 + z / sqrt(z² + R²))`, with a uniform top-hat cross-section.
Solid-body swirl carries supplied torque: `omega = 2 Q / (massFlow Rwake²)` and
`tangentialVelocity = omega axis × radialOffset`. This predicts a spiral slipstream,
not individual tip-vortex filaments or blade-passage fluctuations. The returned
field is zero upstream and outside the wake, with an abrupt disk/wake edge. It is
not a divergence-free solution around obstacles. Swirl does not feed back into
the axial momentum calculation; no shaft-power/efficiency closure is inferred.

```swift
let disk = try ActuatorDiskSlipstream(centreMetres: .zero,
    downstreamAxis: SIMD3(-1, 0, 0), radiusMetres: 1,
    thrustNewtons: 500, fluidTorqueNewtonMetres: 40,
    densityKilogramsPerCubicMetre: 1.2,
    upstreamAxialSpeedMetresPerSecond: 30)
let field = try SuperposedFluidVelocityField(
    backgroundMetresPerSecond: SIMD3(-30, 0, 0), increments: [disk])
let airAtTail = try field.velocityMetresPerSecond(at: SIMD3(-3, 0.4, 0))
```

The axis points along force **on the fluid**; aircraft thrust and torque have opposite
signs. ``ActuatorDiskSlipstream/reactionLoads(momentReferencePointMetres:)`` returns
that reaction once, including its lever arm. Do not add it to the same propeller
loads already returned by a propulsion provider. With an existing propulsion model,
use its supplied thrust/torque only to construct the wake.

This idealization assumes aligned incompressible normal propeller operation.
Crossflow, reversed thrust, windmilling, vortex-ring states, upstream induction,
wake breakdown and body deflection need other providers. Sound-speed and gas/flow
validity must be established by the caller; density alone is not coverage.
Evaluation is O(1) per sample. Multiple disks can be superposed as a prescribed
approximation, without solving mutual induction.

Reference: [NASA ideal propeller momentum theory](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/propeller-thrust/).
Tests independently integrate axial/angular momentum, check hover/forward-flight
limits, axis/torque reversal, reaction moments and wing response to local wash.

## Local incompressible CFD

``IncompressibleFlowGrid`` is an optional experimental CPU solver, not a replacement for every
coefficient or heating model. ``IncompressibleFlowConfiguration`` accepts rectangular grids
of 4–64 cells per axis, with a default of 16³ (4,096 cells); powers of two are not
required. Spacing, origin, constant density, sound speed, kinematic
viscosity and outer boundary velocity are explicit. It conservatively restricts
its component-wise speed bound to Mach ≤ 0.3 throughout the domain.
This is a continuum, constant-density approximation without shock/energy equations.

Velocity lives on staggered faces and gauge pressure at cell centres. Each step
uses first-order semi-Lagrangian advection, explicit viscosity, prescribed fluid
acceleration and pressure projection with a Jacobi-preconditioned conjugate-gradient
solve. The pressure operator uses the same face gradient/divergence discretization.
Projection reports iterations, relative L2 residual and maximum divergence before
and after the solve. A relative residual tolerance and absolute-divergence floor
control convergence; the actual reported maximum divergence should be inspected.

Static solids occupy whole cells and impose zero normal velocity with a free-slip
viscous treatment. They may not touch the outermost cell layer. Outer faces have
prescribed normal velocity with balanced total box flux; tangential interpolation
uses edge values. These are finite-domain boundary conditions, **not** a nonreflecting
far-field or unlimited downstream wake. Increase domain extent independently of
resolution to assess boundary influence. Thin wings, boundary layers and detailed
separation cannot be resolved on a handful of voxels. No RANS/LES closure or artificial
vorticity-confinement force is added. Semi-Lagrangian advection dissipates vorticity.

The grid axes must not rotate and its origin must translate at constant velocity.
No planetary/vehicle rotation or fictitious acceleration is inferred. Fixed geometry
is supplied at cell centres; moving bodies and moving-grid/remeshing support remain
future extensions. The application owns frame transformations and motion coupling.

```swift
let config = try IncompressibleFlowConfiguration(minimumCornerMetres: SIMD3(repeating: -4),
    cellSizeMetres: 0.5, densityKilogramsPerCubicMetre: 1.2,
    speedOfSoundMetresPerSecond: 340,
    kinematicViscositySquareMetresPerSecond: 1.5e-5,
    boundaryVelocityMetresPerSecond: SIMD3(-10, 0, 0))
var grid = try IncompressibleFlowGrid(configuration: config,
    solidAtCellCentre: { p in (p*p).sum() < 1 })
let initialProjection = try grid.project(timeStepSeconds: 0.01)
let report = try grid.advance(timeStepSeconds: 0.01,
    maximumSubsteps: 8, maximumPressureIterations: 300,
    relativePressureTolerance: 1e-6)
let air = try grid.velocityMetresPerSecond(at: SIMD3(-2, 0, 0))
```

Arbitrary initial velocity may contain divergence: project it before treating it
as incompressible. Sampling outside the box or inside a solid throws. Projection
and advance operate transactionally: invalid forcing, coverage failure, pressure
nonconvergence or substep-budget exhaustion leave the original state unchanged.
There is no unconverged fallback. Snapshots are Swift values; thread ownership,
update cadence, asynchronous work and interpolation between updates are external.

Substeps enforce an advective distance limit of 0.75 cell (including a conservative
forcing bound) and the explicit diffusion limit `dt ≤ h²/(6 nu)`. The caller caps
substeps and pressure iterations. A call's supplied acceleration field is held
fixed over its substeps; it is acceleration of the **fluid**, not a requested velocity.
Cost is O(substeps × cells × pressureIterations), plus O(cells) interpolation/stencil
work, with O(cells) storage. The first implementation is serial CPU Swift; SIMD,
parallel/GPU backends and more efficient pressure solvers remain possible improvements.

Numerical formulation background:
[Bridson and Müller-Fischer, Fluid Simulation course notes](https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf).
Tests use a manufactured discrete Poisson/gradient solution, a divergence-free
vortex, an independent viscous eigenmode, uniform transport, hydrostatic buoyancy,
obstacle symmetry, rotor forcing and transactional failure checks. These establish
numerical contracts, not validated aircraft drag/stall or turbulence predictions.

## Rotor forcing and load coupling

``ActuatorDiskAccelerationField`` spreads prescribed disk thrust and fluid torque over
a finite cylinder. Its continuous force and torque integrals match the inputs;
face-sampled grid quadrature has discretization error. Use the same density as the
grid, resolve radius/thickness with multiple cells, keep the source away from outer
boundaries/solids and check resolution convergence. Clipped/blocked source faces
lose that source contribution. This initial source is not discretely normalized to
guarantee exact prescribed total thrust/torque on an arbitrary coarse grid.

```swift
let cfdDisk = try ActuatorDiskSlipstream(centreMetres: SIMD3(2, 0, 0),
    downstreamAxis: SIMD3(-1, 0, 0), radiusMetres: 0.75,
    thrustNewtons: 80, fluidTorqueNewtonMetres: 2,
    densityKilogramsPerCubicMetre: 1.2, upstreamAxialSpeedMetresPerSecond: 10)
let source = try ActuatorDiskAccelerationField(disk: cfdDisk, thicknessMetres: 1)
let report = try grid.advance(timeStepSeconds: 0.01,
    accelerationMetresPerSecondSquared: source.accelerationMetresPerSecondSquared)
```

Use this **instead of** adding analytic propwash to that CFD field. The pressure
solve and advection then determine how the forced flow meets voxel bodies. This is
one-way prescribed rotor forcing; it does not infer rotor loading from the altered
inflow or close an engine/propeller operating-point solve.

Sample the resulting velocity for existing wing/tail/body/drag-area evaluators.
``IncompressibleFlowGrid/pressureLoadsOnSolids(momentReferencePointMetres:)`` optionally
integrates gauge pressure on voxel solids using one-sided linear reconstruction
(falling back to adjacent-cell pressure where necessary). It includes lever arms,
but no viscous shear or unresolved thin-wing lift. Choose which model owns each
pressure/lift/wake contribution. Do not also add the same whole-body coefficient
loads, buoyancy from the same pressure gradient, or a second wake correction.
Sampled-flow coupling alone does not deposit wing forces back into the fluid;
mutual coupling must supply consistent fluid forcing and account for discretization
error. The library does not claim automatic conservation for arbitrary hybrid models.

## Budget measurement

Run the opt-in workload with an optimized build:

```sh
swift run --disable-sandbox -c release -debug-info-format none FluidPhysicsDemo --cfd
```

It reports 16³, 20³ and 32³ timings, pressure iterations, divergence and voxel
pressure force for a fixed 8 m box, free-slip sphere and rotor forcing. Each case
measures three 0.01 s steps after an initial projection, with relative pressure
tolerance 1e-8. This short serial workload excludes setup and rigid-body/render work;
it is not a sustained real-time guarantee or a validated sphere drag coefficient.
Choose resolution, domain, tolerance, cadence and maximum work from measurements
on the target hardware. Propwash and CFD are both optional; neither is scheduled
automatically by other aerodynamic evaluators.

The common solver does not require gas molecular properties for a liquid. Continuum
applicability is a material/model responsibility. Hydrodynamics adds liquid-property
and absolute-pressure reference helpers; see [Hydrodynamics](/documentation/hydrodynamics).
Compressible flow uses a distinct conservation solver described in <doc:CompressibleFlow>.

## RPM-driven coefficient rotor

`PropellerRotor` converts caller-supplied CT, profile CP, radius, angular speed,
collective and density to thrust/shaft torque. It uses `ActuatorDiskSlipstream`
for ideal induced power and accepts added coaxial inflow. CT uses rho A (omega R)^2;
CP uses rho A (omega R)^3. Spin inertia is metadata for a caller's coupled motor
integration. No blade-element lift/stall, dynamic inflow, vortex-ring descent,
windmilling, ground effect, blade flexibility or measured propeller map is supplied.
Do not add its reaction thrust a second time through a separate actuator disk.

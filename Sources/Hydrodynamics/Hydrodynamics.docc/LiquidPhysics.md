# Liquid buoyancy, hulls and submersibles

Evaluate dedicated liquid models first, and use the shared incompressible solver for fully filled submerged-flow domains.

## Material and pressure inputs

``LiquidPropertyModel`` evaluates material state at supplied temperature and absolute
pressure. ``LiquidProperties`` carries identifier/provenance, density, dynamic
viscosity, sound speed and vapor pressure. ``ConstantLiquidModel`` is an explicit
bounded approximation; its declared range must remain above its supplied vapor
pressure. There is no default ocean composition, stratification, salinity correction,
thermal evolution or phase-equilibrium solver. ``TabulatedLiquidModel`` accepts
matched scalar data with `temperature` (K) and absolute `pressure` (Pa) axes for
density (kg/m³), viscosity (Pa·s), sound speed (m/s) and vapor pressure (Pa). It
validates common grids/units, positive properties and liquid-phase coverage at
all nodes, and preserves source/version data without extrapolation. Composition
is an explicit dataset identity; no lake mixture is selected automatically.
User-supplied providers can add richer phase/thermodynamic models.

Titan's seas require composition/temperature/pressure-dependent hydrocarbon
properties rather than a water model with a different density. The demo uses
explicitly labelled constant-property illustrations for Earth and Titan, not a
calibrated cryogenic property catalogue. NASA's conceptual Titan submarine work is
one useful source of declared design assumptions:
[NASA/TM-2015-218831](https://ntrs.nasa.gov/api/citations/20150014581/downloads/20150014581.pdf).

``LiquidProperties/cavitationNumber(absolutePressurePascals:speedMetresPerSecond:)``
provides `(p_absolute - p_vapor) / (0.5 rho V²)` as a diagnostic. A positive value
alone does not establish safe propeller/hydrofoil operation; blade-local pressure
coefficients and appropriate reference speed are still needed. Cavitation/ventilation,
bubbles and mass transfer are not simulated by the initial models.

## Calm-surface displacement

``HydrostaticLiquid`` defines a flat surface, constant density, physical gravity,
surface absolute pressure and a declared maximum depth. Up must oppose gravity.
The hydrostatic pressure is `p_surface + rho |g| depth`; queries above the surface
return no liquid pressure. Gravity is physical gravity, not vehicle acceleration.
The depth bound is additional to material-property coverage, not a density EOS.

``DisplacementGeometry`` allows custom supplied geometry to return a checked
``DisplacedLiquidVolume``. ``ConvexDisplacementHull`` validates closed outward-wound
planar convex faces, clips them against the liquid plane and integrates exact
polyhedral volume/centroid. Its clipping cap closes the displaced volume but is
excluded from wetted hull area. ``SphereDisplacement`` evaluates spherical caps
analytically. Dry shapes have zero displacement and no buoyancy centre.

```swift
let material = try LiquidProperties(identifier: "supplied-liquid", source: "caller data",
    densityKilogramsPerCubicMetre: 1000, dynamicViscosityPascalSeconds: 0.001,
    speedOfSoundMetresPerSecond: 1500, vaporPressurePascals: 2000)
let liquid = try HydrostaticLiquid(properties: material, up: SIMD3(0, 0, -1),
    gravityMetresPerSecondSquared: SIMD3(0, 0, 9.80665),
    surfaceAbsolutePressurePascals: 101325, maximumDepthMetres: 100)
let hull = try ConvexDisplacementHull.box(halfExtentsMetres: SIMD3(2, 1, 0.5))
let displacement = try hull.displacement(in: liquid)
let loads = try displacement.buoyancy(in: liquid,
    momentReferencePointMetres: suppliedCentreOfMass)
```

Rotate/translate the liquid plane and gravity into body axes externally for a
heeled/pitched hull. The buoyancy lever arm then supplies restoring or overturning
moments without a separate invented stability force. Weight acts at the actual
centre of mass through the consumer's gravity model, once. Convex components may
approximate a concave hull only if their displaced volumes do not overlap; no
vehicle assembly or automatic constructive-solid-geometry union is owned here.
Mesh construction checks convex face ordering and paired edges with a relative
geometry tolerance. Construction is O(faces × vertices); clipping includes cap
sorting/uniqueness work and is intended for modest dedicated buoyancy geometry.

This is calm-surface hydrostatics, without wave evolution, slamming, planing,
capillarity or air displacement above the waterline. Add atmospheric buoyancy for
the dry volume separately where needed. Pressure-datum contributions must not be
double counted across gas/liquid models.

## Submerged loads

``HydrodynamicDamping`` accepts nonnegative body-axis linear/quadratic force and
moment coefficients. Relative velocity must be sampled at the damping reference
point, and relative rates must match the coefficient convention. The model is
passive: force/torque remove energy from relative motion. Coefficients are supplied
data, not universal submarine drag or a substitute for Reynolds-dependent tables.
Shared `DragAreaForces` and surface traction integration remain available.

``HydrodynamicAddedInertia`` provides diagonal translational/rotational added inertia
and a fixed-axis small-motion acceleration-reaction approximation. Include that
in effective inertia **or** apply the acceleration reaction; never both. General
rotating-frame marine Coriolis terms and translation/rotation coupling are not
inferred. ``SubmergedBodyHydrodynamics`` combines supplied displaced-volume
buoyancy with damping for a fully wetted body, without subtracting vehicle weight.

Ballast/gas inventory, flooding, pumps, seals, actuators, control logic, failure and
mass/inertia evolution belong to external systems/vehicle consumers. This module
evaluates the supplied physical state and does not run a vessel controller.

## Liquid CFD and shared propellers

The shared `IncompressibleFlowGrid` accepts a liquid configuration through
`IncompressibleFlowConfiguration(grid:liquid:boundaryVelocityMetresPerSecond:)`.
It uses the liquid density, sound speed and kinematic viscosity without requiring
gas mean-free-path inputs. Actuator-disk slipstream and volume-force models in
FluidDynamics accept the same supplied thrust/torque with liquid density. Source
quadrature, low-Mach coverage, fixed voxel/free-slip boundaries and duplicate-force
rules remain as documented in [FluidDynamics](/documentation/fluiddynamics/propellersandcfd).

``LiquidCFDPressureReference`` calibrates the arbitrary CFD gauge pressure to a known
absolute pressure at a known fluid point in the same connected fluid region. Include gravity/hydrostatic pressure consistently
before evaluating depth/phase/cavitation. An uncalibrated zero-mean pressure array
is not an absolute pressure measurement. If CFD pressure already contains hydrostatic
displacement or an unsteady added-inertia contribution, do not add that contribution
again from the dedicated models. A consistent decomposition can instead combine
dynamic-only CFD pressure with separately evaluated hydrostatic buoyancy.

The current grid is a **fully filled single-phase domain**. It does not evolve an
air/liquid interface. Sailing/floating can already use dedicated hull displacement
with atmospheric sail loads; CFD-resolved waves, waterline evolution and breaking
wakes need a later free-surface/volume-fraction extension. Do not voxelize the entire
box as liquid and interpret its top as a simulated lake surface.

Tests cover exact clipped tetrahedron/box/sphere volumes and centres, a metacentric
restoring slope, density/gravity scaling, passive damping, added inertia and a
calibrated hydrostatic CFD pressure field. `FluidPhysicsDemo` runs illustrative
Earth/Titan cases; its synthetic properties and loads are not validated vessel data.

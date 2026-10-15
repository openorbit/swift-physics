# ``Hydrodynamics``

Liquid properties, displacement, hydrostatic loads and supplied-state submerged-body physics.

Hydrodynamics depends only on FluidDynamics. Airfoils and gas-specific entry models
remain in Aerodynamics; orbital/body frames and rigid-body dynamics remain in CelMek.
Both domains return the same `FluidLoads`. All positions, velocities and gravity
vectors use common caller axes; water/air currents and rotating-body frames must
be handled before local force evaluation.

## Topics

### Floating and submerged physics
- <doc:LiquidPhysics>
- ``LiquidProperties``
- ``LiquidPropertyModel``
- ``ConstantLiquidModel``
- ``TabulatedLiquidModel``
- ``HydrostaticLiquid``
- ``DisplacementGeometry``
- ``DisplacedLiquidVolume``
- ``ConvexDisplacementHull``
- ``SphereDisplacement``
- ``HydrodynamicDamping``
- ``HydrodynamicAddedInertia``
- ``SubmergedBodyHydrodynamics``

### Shared CFD coupling
- ``LiquidCFDPressureReference``

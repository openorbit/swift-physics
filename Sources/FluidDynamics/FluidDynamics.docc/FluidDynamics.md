# ``FluidDynamics``

Shared gas/liquid flow contracts, physical loads, propeller slipstreams and numerical fluid solvers.

FluidDynamics has no dependency on Aerodynamics, Hydrodynamics, Atmospheres, CelMek
or Yams. Those consumers supply material, geometry and regime-specific physics.
All values use SI in caller-selected common axes. `FluidFlow` is body motion
relative to the medium; a `FluidVelocityField` is physical medium velocity, so
subtract it from body velocity before rotating the difference into body axes.
Body conventions are X forward, Y right, Z down. Forces/moments and application
points share axes and an explicit reference point.

The physics owns field state and numerical timesteps where necessary; applications
own scheduling, frames, rigid-body integration and system/configuration evolution.
No old-name aliases or module re-exports preserve the pre-extraction API.

## Topics

### Shared contracts
- ``FluidError``
- ``FluidFlow``
- ``FluidLoads``
- ``FluidCoefficients``
- ``FluidReferenceGeometry``
- ``SurfaceHeatFlux``
- ``FluidScalarTable``
- ``FluidTableAxis``
- ``CartesianFluidGrid``

### Surface integration and displacement
- ``FluidSurfacePanel``
- ``SurfacePressureCoefficientModel``
- ``SurfacePressureIntegration``
- ``SurfacePressureResult``
- ``SurfaceTraction``
- ``SurfaceTractionModel``
- ``SurfaceTractionIntegration``
- ``SurfaceTractionResult``
- ``PressureOnlyTraction``
- ``DragCoefficientQuery``
- ``DragCoefficientModel``
- ``DragCoefficientSample``
- ``TabulatedDragCoefficients``
- ``DragAreaForces``
- ``Buoyancy``
- ``BuoyancyField``

### Materials and thermodynamics
- ``ThermalTransportProperties``
- ``PerfectGasSpecies``
- ``PerfectGasMassFraction``
- ``PerfectGasMixture``
- ``CompressibleEquationOfState``
- ``CaloricallyPerfectGasEOS``
- ``FluidThermodynamicState``
- ``ConservedFluidState``
- ``CompressiblePrimitiveState``

### Propellers and incompressible flow
- <doc:PropellersAndCFD>
- ``FluidVelocityField``
- ``InducedFluidVelocityField``
- ``SuperposedFluidVelocityField``
- ``ActuatorDiskSlipstream``
- ``ActuatorDiskAccelerationField``
- ``IncompressibleFlowConfiguration``
- ``IncompressibleFlowGrid``
- ``PressureProjectionReport``
- ``IncompressibleStepReport``

### Compressible flow
- <doc:CompressibleFlow>
- ``CompressibleFlowConfiguration``
- ``CompressibleBoundary``
- ``CompressibleFlowGrid``
- ``ConservedFluidTotals``
- ``CompressibleStepReport``

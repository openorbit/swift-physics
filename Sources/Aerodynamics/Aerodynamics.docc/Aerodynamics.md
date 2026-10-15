# ``Aerodynamics``

Local flow, coefficient data, airfoil geometry and physical aerodynamic loads.

Aerodynamics provides parametric and thin-airfoil models, NACA geometry and YAML
coefficient input. Atmospheric profiles live in Atmospheres. Shared flow, load, propeller and CFD
contracts live in [FluidDynamics](/documentation/fluiddynamics); import those modules explicitly.

## Topics

### Physics foundations
- <doc:ModelCoverage>
- ``GasFlowCoverage``
- ``CoveredBodyCoefficients``
- ``MachRegimeBodyCoefficients``
- <doc:Foundations>

### Checked coefficient models
- ``AirfoilCoefficientModel``
- ``AirfoilPolar``
- ``TabulatedAirfoil``
- ``AirfoilTableLoader``
- ``CoefficientBoundaryPolicy``

### Wing and whole-body evaluation
- <doc:WingsAndBodies>
- <doc:FiniteWing>
- ``FiniteWingAirfoil``
- ``FiniteWingSolution``
- <doc:NonlinearWakes>
- ``NonlinearLiftingLine``
- ``NonlinearLiftingLineResult``
- <doc:LiftingLine>
- ``LiftingLineWing``
- ``LiftingLineSection``
- ``LiftingLineResult``
- ``LiftingLineStationResult``
- ``SegmentedWing``
- ``WingSegment``
- ``WingSectionAxes``
- ``WingSectionCoefficientModel``
- ``AirfoilWingSectionModel``
- ``WholeBodyAerodynamics``
- ``BodyAerodynamicCoefficientModel``
- ``TabulatedBodyCoefficients``
- ``WindForceCoefficients``
- <doc:BodyRates>
- ``BodyAerodynamicRateReference``
- ``BodyAerodynamicRateDerivatives``
- ``BodyAerodynamicRateDerivativeModel``
- ``ConstantBodyRateDerivatives``
- ``LinearBodyRateAerodynamics``
- <doc:ControlSurfaces>
- ``ControlSurfaceSetting``
- ``TabulatedControlSurface``
- ``ControlSurfaceTableLoader``

### Airfoil models
- ``Airfoil``
- ``AirfoilCoefficients``
- ``ParametricAirfoil``
- ``ThinSymetricAirfoil``
- ``ThinCamberedAirfoil``
- ``NACA4``

### Hypersonic surface pressure
- <doc:SurfacePressure>
- ``ModifiedNewtonianPressure``
- <doc:RarefiedFlow>
- ``FreeMolecularGas``
- ``GasSurfaceReflection``
- ``MaxwellianFreeMolecularTraction``
- ``TabulatedRarefiedTraction``

### Aerothermal physics
- <doc:SurfaceHeating>
- ``SurfaceHeatingQuery``
- ``SurfaceConvectiveHeatingModel``
- ``SurfaceConvectionResult``
- ``LaminarFlatPlateHeating``
- ``TurbulentFlatPlateHeating``
- <doc:EntryHeating>
- ``AerothermalBoundaryLayer``
- ``EnthalpySurfaceHeatingCoverage``
- ``EnthalpySurfaceHeatingQuery``
- ``EnthalpySurfaceHeatingResult``
- ``EnthalpySurfaceHeatingModel``
- ``TabulatedEnthalpySurfaceHeating``
- ``LeadingEdgeHeatingQuery``
- ``LeadingEdgeConvectiveHeatingModel``
- ``SweptCylinderStagnationHeating``
- ``TabulatedRadiativeHeating``
- <doc:Aerothermal>
- ``StagnationHeatingQuery``
- ``AerothermalWallState``
- ``ConvectiveHeatingModel``
- ``RadiativeHeatingModel``
- ``SuttonGravesParameters``
- ``StagnationHeatingCoverage``
- ``SuttonGravesHeating``
- <doc:PlanetaryHeating>
- ``StagnationGasState``
- ``SuttonGravesBaseGas``
- ``SuttonGravesMixtureParameters``
- ``SuttonGravesMixtureHeating``
- ``HydrogenHeliumHeatingMixture``
- ``HydrogenHeliumColdWallHeating``
- ``HydrogenHeliumHeatingResult``

### Optional airflow simulation
- [shared propwash/CFD guide](/documentation/fluiddynamics/propellersandcfd)

### Parachutes and balloons
- <doc:ParachutesAndBalloons>
- ``ParachuteAerodynamics``
- ``BalloonAerodynamics``

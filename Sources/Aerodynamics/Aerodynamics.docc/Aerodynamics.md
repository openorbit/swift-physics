# ``Aerodynamics``

Local flow, coefficient data, airfoil geometry and physical aerodynamic loads.

Aerodynamics provides parametric and thin-airfoil models, NACA geometry and YAML
coefficient input. Atmospheric profiles live in the Atmospheres module, which
Aerodynamics re-exports for compatibility.

## Topics

### Physics foundations
- <doc:Foundations>
- ``AerodynamicFlow``
- ``AerodynamicReferenceGeometry``
- ``AerodynamicCoefficients``
- ``AerodynamicLoads``
- ``AerodynamicHeatFlux``

### Checked coefficient models
- ``AirfoilCoefficientModel``
- ``AirfoilPolar``
- ``TabulatedAirfoil``
- ``AirfoilTableLoader``
- ``CoefficientBoundaryPolicy``

### Wing and whole-body evaluation
- <doc:WingsAndBodies>
- ``SegmentedWing``
- ``WingSegment``
- ``WingSectionAxes``
- ``WingSectionCoefficientModel``
- ``AirfoilWingSectionModel``
- ``WholeBodyAerodynamics``
- ``BodyAerodynamicCoefficientModel``
- ``TabulatedBodyCoefficients``
- ``WindForceCoefficients``

### Airfoil models
- ``Airfoil``
- ``AirfoilCoefficients``
- ``ParametricAirfoil``
- ``ThinSymetricAirfoil``
- ``ThinCamberedAirfoil``
- ``NACA4``

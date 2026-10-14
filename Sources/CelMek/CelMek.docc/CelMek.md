# ``CelMek``

Astronomical time, celestial bodies, ephemerides, reference frames and spacecraft dynamics.

CelMek is one module in the `swift-physics` package. All shared physical states use
explicit units, epochs, origins and frames. Select the module products your application
needs, or the `Physics` product to include all five modules.

## Topics

### Getting started
- <doc:PackageOverview>
- <doc:CONVENTIONS>
- <doc:STATES>
- <doc:BODIES>

### Time and frames
- <doc:TimeGuide>
- <doc:FRAMES>
- ``Epoch``
- ``CartesianState``
- ``BodyIdentifier``

### Ephemerides and data
- <doc:EPHEMERIDES>
- <doc:DE441>
- <doc:SPICE_KERNELS>
- <doc:ACCURACY>
- <doc:Algorithms>

### Spacecraft simulation
- <doc:ORBITAL_MECHANICS>
- <doc:SPACECRAFT_DYNAMICS>
- <doc:ATTITUDE_DYNAMICS>
- <doc:GAME_INTEGRATION>
- ``BodyMagneticProvider``
- ``SphericalBodyAtmosphere``
- <doc:AerodynamicIntegration>
- ``AerodynamicWrench``

### Accuracy, performance and reproducibility
- <doc:VSOP87_ACCURACY_CHARACTERIZATION>
- <doc:NUTATION_ACCURACY>
- <doc:GAME_PERFORMANCE>
- <doc:THREAD_SAFETY>
- <doc:REPRODUCIBILITY>
- <doc:SERIALIZATION>
- <doc:Attribution>

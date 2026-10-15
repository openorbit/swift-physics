# swift-physics

Swift physics libraries for celestial mechanics, gas and liquid dynamics, atmospheres and magnetic fields.

The celmek library implements several astronomical algorithms.

| Algorithm | Status
| --- | --- 
| VSOP87 | Done
| ELP2000 82b | Done
| Goffin2000 | Done

## Status

CelMek is very much in development, and the scope has grown to implement more astronomical algorithms.
It is currently NOT in a 1.0 state, so expect that API breakage happens.

A companion app for interactivelly working with this library is developed here:
[https://github.com/lorrden/CelMek-Viewer]

## Evolution

The package has grown beyond the original CelMek port. Shared fluid physics now
lives in FluidDynamics, with dedicated Aerodynamics and Hydrodynamics modules.

## Modules

Use the `CelMek`, `Math`, `FluidDynamics`, `Aerodynamics`, `Hydrodynamics`,
`Atmospheres`, or `MagneticFields` product for the modules your application needs.
The `Physics` product bundles all seven.
Import the corresponding module names. See [MODULES.md](MODULES.md) for responsibilities and dependency rules.

## Documentation

The public guides and API references use DocC. See [documentation build instructions](Documentation/README.md)
for the seven module catalogs and the combined package archive.

```sh
python3 Tools/Documentation/build.py
python3 Tools/Documentation/preview.py
```

Browse the result at [localhost:8080/documentation](http://127.0.0.1:8080/documentation).

## Time API

Precise epochs, durations, UTC leap seconds, and calendar readings with explicit
TAI/TT/TDB/UTC scales are documented in [time guide](Sources/CelMek/CelMek.docc/TimeGuide.md).

Coordinate handedness, quaternion and angle conventions, transform direction, and
error handling are documented in [CONVENTIONS.md](Sources/CelMek/CelMek.docc/CONVENTIONS.md).

## Testing

The tests use Swift Testing and require a Swift 6 or newer toolchain. Run all
test targets with:

```sh
swift test
```

## License

Apache 2.0 for the package's original code. ERFA-derived Earth-orientation and
Fairhead TDB numerical code is isolated in the internal `ERFASwift` module under
its [BSD terms and attribution](Sources/ERFASwift/LICENSE.txt). See the
[module README](Sources/ERFASwift/README.md) for its scope.

Shared Cartesian states use metres and metres per second, with explicit unit conversion and a VSOP87 provider adapter. See [STATES.md](Sources/CelMek/CelMek.docc/STATES.md).

Time-dependent frame trees, rotating transforms, and their conventions are documented in [FRAMES.md](Sources/CelMek/CelMek.docc/FRAMES.md).

The IAU coefficient review, updated data, and independent SPICE reference tests are documented in [IAU_ROTATION_AUDIT.md](IAU_ROTATION_AUDIT.md).

Body identities, sourced SI properties, and the injectable JPL GM catalogue are documented in [BODIES.md](Sources/CelMek/CelMek.docc/BODIES.md).

The common ephemeris contract, validating state-provider adapter, and body-query collection are documented in [EPHEMERIDES.md](Sources/CelMek/CelMek.docc/EPHEMERIDES.md).

Shared accuracy intent and each adapter's supported levels are documented in [ACCURACY.md](Sources/CelMek/CelMek.docc/ACCURACY.md).

Reduced VSOP87 error and query-cost measurements, with a reproducible release tool, are in [VSOP87_ACCURACY_CHARACTERIZATION.md](Sources/CelMek/CelMek.docc/VSOP87_ACCURACY_CHARACTERIZATION.md).

Atmospheric profiles are provided by the [Atmospheres target](Sources/Atmospheres/Atmospheres.docc/Guide.md),
shared by CelMek and Aerodynamics. It includes Earth ISA, nominal planetary/moon
profiles, and gas-giant pressure during descent below the 1-bar level.

The [MagneticFields target](Sources/MagneticFields/MagneticFields.docc/Guide.md) supplies
pluggable magnetic-field models, ideal dipoles and CelMek body/orientation adapters.
It also includes Swift WMM/IGRF evaluators, coefficient-file loaders and bundled
WMM2025/IGRF-14 datasets with explicit date coverage.

## Fluid physics examples

```sh
swift run FluidPhysicsDemo
swift run -c release -debug-info-format none FluidPhysicsDemo --cfd --compressible
```

The example includes dedicated Earth/Titan liquid displacement and pressure cases,
optional incompressible rotor/body CFD and a compressible Sod reference problem.
Liquid properties and vehicle coefficients in these examples are illustrative.
See the FluidDynamics and Hydrodynamics DocC guides for coverage, numerical cost,
pressure-reference rules and the remaining real-gas/free-surface fidelity work.

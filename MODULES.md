# Module boundaries

Keep the current three Swift modules and package name. Add focused library products
for selective consumption, while retaining the existing `Physics` product. Module
responsibilities are defined below; existing public types stay in their current modules.

## Package products and imports

| SwiftPM library product | Included module targets | Target dependencies |
| --- | --- | --- |
| `CelMek` | `CelMek` | `Math` |
| `Math` | `Math` | None |
| `Aerodynamics` | `Aerodynamics` | `Yams` |
| `Physics` (existing bundle) | `CelMek`, `Aerodynamics`, `Math` | Dependencies of those targets |

`Physics` is a product name, not a Swift module: there is no `import Physics`.
Imports remain `import CelMek`, `import Math`, and `import Aerodynamics`. Applications
that directly use Math types should declare the Math product explicitly as well.
Existing consumers selecting Physics can keep their manifest and imports.

For a local development dependency:

```swift
// In the application's Package.swift:
.package(path: "/path/to/swift-celmek")

// In the application's target dependencies:
.product(name: "CelMek", package: "swift-celmek"),
.product(name: "Math", package: "swift-celmek")
```

Selecting CelMek avoids linking the Aerodynamics/Yams targets. The package manifest
still declares Yams and the documentation plugin, so SwiftPM may resolve or fetch
those package dependencies even when only CelMek is selected. Product selection is
not a promise of offline resolution or separately downloadable packages. CelMek still
compiles its built-in coefficient tables; optional model packaging is future work.

## Responsibilities and dependency direction

**Math** owns general numerical and geometric operations: quaternions, matrices,
coordinate representations, angle utilities, and the existing frame identity/parent
primitive. It does not depend on CelMek or Aerodynamics. Astronomical frame names,
epochs, units of physical state, ephemerides, and orientation models belong above it.
Existing astronomical coordinate structs and the legacy `Transformation` remain
here for compatibility; that is not a requirement to place new astronomy models here.

**CelMek** owns astronomical time, shared physical states and conventions, the frame
tree and epoch-dependent transforms, celestial properties, orientations, ephemerides,
and orbital mechanics. New force, propagation and attitude interfaces initially live
here because they share epochs, states and frames. Keep the following conceptual
layers within the module rather than immediately introducing more targets:

| Layer | Responsibilities / allowed dependencies |
| --- | --- |
| Time | Epochs, durations, time-scale and calendar adapters, injected time data; independent of body catalogues and dynamics. |
| State and frames | SI Cartesian states, provider interfaces, frame composition and named physical frames; depends on time and Math. |
| Bodies and models | Body identities/properties, ephemeris and orientation implementations, coefficient/file evaluation; uses shared states and frames. |
| Orbits and dynamics | Propagators, force/torque interfaces, trajectories and maneuvers; consumes prescribed models through interfaces. |
| Data adapters | Parsers/loaders and storage; produce validated model/provider data with coverage and provenance. |

These are development boundaries, not separately importable submodules. Provider
interfaces may reference shared types without depending on concrete implementations.
For example, a third-body force should consume an ephemeris interface; it should not
require a particular DE reader or analytic series. Ephemerides must not depend on the
spacecraft integrator. Frames consume origin/orientation providers without owning a
simulation clock. Data loading remains separate from evaluating a provider.

**Aerodynamics** owns airfoil and atmosphere models, including its current YAML input.
It remains independently usable and does not depend on CelMek. When spacecraft drag
is implemented, an optional adapter can translate atmosphere results into CelMek force
inputs after Aerodynamics units and validity have been reviewed. CelMek should not
acquire a Yams dependency merely to expose a drag-force interface.

**Game/renderer adapters** belong in consumer applications or future optional targets.
They consume CelMek and optionally Aerodynamics. Engine clocks, rendering axes,
origin shifting, Float conversion and engine objects must not become core dependencies.
All shared physical coordinates remain right-handed; engine conversion is an explicit
adapter boundary.

## When to extract another target

Extract a target when there is a concrete consumer or build/data-cost benefit and
its interface and dependency direction are understood. Likely candidates are large
optional ephemeris assets/file readers, spacecraft dynamics, game-engine bridges,
and atmosphere-to-force adapters. Such a split must preserve one shared definition
of epoch, state and frame identity; it must not duplicate frame singletons or require
mutual imports.

A future split can move shared contracts to a lower target before moving implementations,
with an explicit compatibility/migration plan. Until then, keep numerical implementations
in Swift and avoid adding runtime dependencies to compensate for missing algorithms.
Reference tools may use external implementations for independent fixture generation;
they are developer tooling rather than linked library dependencies.

No package rename, source relocation, broad public-visibility expansion or namespace
migration is part of this decision. New public APIs should expose validated contracts;
legacy internal routines can be adapted as needed rather than exported wholesale.

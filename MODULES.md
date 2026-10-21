# Module boundaries

The package is `swift-physics`. Eight public modules share physical contracts
without importing vehicle systems or rendering. There are no compatibility
re-exports or old aerodynamic names for the shared fluid APIs.

## Package products and imports

| SwiftPM product | Module dependencies |
| --- | --- |
| `FluidDynamics` | None |
| `Hydrodynamics` | `FluidDynamics` |
| `Aerodynamics` | `FluidDynamics`, `Atmospheres`, `Yams` |
| `Mechanics` | `CelMek`, `Math` |
| `CelMek` | `Math`, `ERFASwift` (internal), `Atmospheres`, `MagneticFields`, `FluidDynamics` |
| `Math` | None |
| `Atmospheres` | None |
| `MagneticFields` | None |
| `Physics` | Bundles all eight module targets |

`Physics` is a product name, not a module. Import each module whose APIs are used;
importing Aerodynamics does not re-export FluidDynamics or Atmospheres. A consumer
using hull and sail models imports Hydrodynamics, Aerodynamics and FluidDynamics.
CelMek's `FluidWrench` accepts shared loads without linking Aerodynamics or Yams.
SwiftPM still resolves package-level dependencies independently of product selection.

For a local development dependency:

```swift
.package(path: "/path/to/swift-celmek")
// Target dependencies; the package identity here follows the local folder name.
.product(name: "CelMek", package: "swift-celmek"),
.product(name: "FluidDynamics", package: "swift-celmek"),
.product(name: "Hydrodynamics", package: "swift-celmek")
```

## Responsibilities and dependency direction

**Math** owns general numerical and geometric operations: quaternions, matrices,
coordinate representations, angle utilities, and the existing frame identity/parent
primitive. It does not depend on CelMek or Aerodynamics. Astronomical frame names,
epochs, units of physical state, ephemerides, and orientation models belong above it.
Existing astronomical coordinate structs and the legacy `Transformation` remain
here for compatibility; that is not a requirement to place new astronomy models here.

**Mechanics** owns independent mechanical bodies, compliant ground contact,
suspended driven wheels, tension-only tethers, gas cushions and harpoon anchors.
These are not celestial/orbital models. It currently consumes CelMek's public
epoch, inertial-frame and inertia primitives and Math's rotations; CelMek does not
depend on Mechanics. Import Mechanics explicitly; no CelMek compatibility re-export
is provided. Vehicle-specific controllers and rendering stay in the application.

**ERFASwift** is an internal target, with no standalone product or public API.
It isolates the ERFA-derived Earth-orientation and full Fairhead TDB routines and
tables under their BSD terms, retained in `Sources/ERFASwift/LICENSE.txt`.
It depends only on Math; CelMek adapts epochs, coverage, frames, and errors at the
boundary using package-scoped entry points.

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

**Atmospheres** owns atmospheric profiles, body lookup, units and coverage. It is a
dependency-free target exposed through the `Atmospheres` product.

**MagneticFields** owns magnetic-field queries, vectors, coverage, provenance,
dipole models and nominal body-centre lookup. It is a target without
dependencies. WMM/IGRF coefficient loaders and Swift harmonic evaluators implement the same provider interface.
CelMek adds body identity, epoch-to-calendar injection and orientation adapters;
magnetic evaluators do not depend on ephemerides or spacecraft integration.

**FluidDynamics** owns medium-neutral flow, force/moment/reference geometry,
surface traction integration, displacement buoyancy, scalar tables, propeller
slipstreams and shared Cartesian numerical infrastructure. Incompressible projection
and compressible conservative solvers are separate formulations. Gas EOS and
transport contracts are providers rather than atmosphere/body catalogues. It has
no dependency on other package modules, YAML, clocks or rigid-body state.

**Aerodynamics** owns gas/airfoil/wing physics, aerodynamic coefficient data,
finite-wing models, rarefied gas traction and entry heating correlations. It consumes
FluidDynamics and explicit Atmospheres samples, without owning generic fluid contracts.

**Hydrodynamics** owns liquid property providers, calm-surface partial-immersion
geometry, hydrostatic/restoring loads, submerged-body damping/added inertia and
liquid-specific CFD pressure calibration. It reuses FluidDynamics propeller/CFD
machinery without importing Aerodynamics. Fully filled liquid CFD is distinct from
future evolving free-surface/wave/two-phase support.

Vehicle assembly, deployment, ballast/flooding/gas management, actuators, failures
and coupled subsystem evolution belong to a separate systems package. Fluid models
may own fluid field state, geometric interfaces and their numerical timesteps;
consumers own scheduling, rigid-body integration and frame/state orchestration.

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

Extract shared contracts before implementations and preserve acyclic dependencies. Until then, keep numerical implementations
in Swift and avoid adding runtime dependencies to compensate for missing algorithms.
Reference tools may use external implementations for independent fixture generation;
they are developer tooling rather than linked library dependencies.

Shared vector helpers use package visibility; implementation helpers are not promoted
to public APIs solely to cross a target boundary. New public APIs expose validated
physical contracts rather than vehicle-system ownership.

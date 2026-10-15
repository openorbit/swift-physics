# Package overview

Choose the modules and library products that match your application.

| Product and module | Responsibility | Target dependencies |
| --- | --- | --- |
| `CelMek` | Time, states, frames, bodies, ephemerides, orbits and spacecraft dynamics | Math, Atmospheres, MagneticFields |
| `Math` | General numerical and geometric operations | None |
| `FluidDynamics` | Shared fluid loads, propellers and CFD | None |
| `Hydrodynamics` | Liquid properties, buoyancy and submerged loads | FluidDynamics |
| `Aerodynamics` | Airfoils, gas flow and heating | FluidDynamics, Atmospheres, Yams |
| `Atmospheres` | Atmospheric profiles, reference data and body lookup | None |
| `MagneticFields` | Dipoles, WMM/IGRF, coefficient loaders and body lookup | None |

The `Physics` product bundles all seven modules. It is a product name, so import
individual modules rather than `Physics`.

```swift
import CelMek
import Math
import Atmospheres
import MagneticFields
```

Applications that use Math types directly should declare the Math product explicitly.
Select `CelMek` to avoid linking Aerodynamics and Yams; SwiftPM still resolves the
package's declared dependencies. Body environment providers are explicit model
choices, with provenance and coverage, independent of ephemerides and orientation.

Atmospheres and MagneticFields have their own catalogs in the combined package
archive. Import Atmospheres and FluidDynamics explicitly where their types are used;
Aerodynamics has no compatibility re-export layer.
The developer architecture policy remains in the repository's `MODULES.md`.

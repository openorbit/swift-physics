# Flow, coefficient data and physical loads

Evaluate aerodynamic physics from supplied local flow and geometry, independently of vehicle systems.

Aerodynamics owns coefficients, section geometry and physical loads. Vehicle assembly,
actuator commands, deployment/reefing sequences, gas management and thermal/structural
state evolution belong to external consumers. Wing strips and whole-body coefficient evaluators build on this foundation; see
<doc:WingsAndBodies>. Finite-wing induced flow, vehicle-specific entry datasets and
heating correlations remain follow-up work.

## Coefficient datasets

``AirfoilPolar`` preserves explicit angle samples at one Reynolds/Mach condition.
``TabulatedAirfoil`` interpolates linearly in angle, Reynolds number and Mach number.
The Reynolds/Mach grid must be rectangular and contain exactly one polar at every
condition. Different polars may have different angle spacing. For a query, only
contributing polars must cover the requested angle. A singleton axis covers only its
one value under the default rejecting policy; it does not imply Reynolds/Mach independence.

```swift
import Aerodynamics
let samples = [
    AirfoilSample(angleOfAttackRadians: -0.1,
        coefficients: AirfoilCoefficients(lift: -0.5, drag: 0.02, moment: 0)),
    AirfoilSample(angleOfAttackRadians: 0.1,
        coefficients: AirfoilCoefficients(lift: 0.5, drag: 0.02, moment: 0))
]
let polar = try AirfoilPolar(reynoldsNumber: 1_000_000, machNumber: 0.2, samples: samples)
let model = try TabulatedAirfoil(polars: [polar], source: "Example synthetic polar", version: "1")
let coefficients = try model.coefficients(
    angleOfAttackRadians: 0.05, reynoldsNumber: 1_000_000, machNumber: 0.2)
```

Tables reject malformed/nonfinite coefficients, nonpositive Reynolds numbers,
negative Mach numbers, duplicate or unordered angles, duplicate conditions and missing
grid cells. Coefficient signs are preserved; validation does not certify scientific
accuracy or a nonnegative measured drag. Codable decoding revalidates the data.
Out-of-coverage queries throw by default. ``CoefficientBoundaryPolicy/clamp`` explicitly
uses the nearest boundary value instead; this is not a physical extrapolation model.

``AirfoilTableLoader`` reads the existing XFOIL-style YAML columns: angle in degrees,
CL, CD, CDp and CM. Every `data` entry is retained; optional `mach-number` defaults to
zero. File I/O is separate from evaluation and source/version metadata is explicit.
Malformed rows or conditions fail before constructing the evaluator.

The legacy ``ParametricAirfoil`` remains available. Its scalar queries use radians,
query Mach zero and explicitly clamp to dataset boundaries. Invalid queries return
NaN. It now interpolates between loaded Reynolds datasets and preserves irregular
angle samples. Prefer the throwing table API when coverage matters. Existing 361-sample
Codable degree grids are readable, as are older 362-sample grids with −181° padding;
new encoding retains the validated table, coordinate samples and provenance.

## Flow conventions and missing properties

``AerodynamicFlow`` uses right-handed body or section axes: X forward, Y right,
Z down. Its velocity is the body's motion relative to the atmosphere. Positive
angle of attack is `atan2(vz, vx)`; sideslip is `atan2(vy, hypot(vx, vz))`.
Angle of attack is absent when the X/Z velocity projection is zero; sideslip is absent
at zero speed. All vectors must be expressed in the same axes by the caller.

Dynamic pressure uses `q = ρ V² / 2`. Mach, Reynolds and Knudsen numbers require
speed of sound, dynamic viscosity and mean free path respectively. Missing properties
remain unavailable; they never become zero. Reference lengths must be positive.
Atmospheres profiles can directly supply ambient properties:

```swift
import Aerodynamics
import Atmospheres
let ambient = try ISAAtmosphere().properties(at: 0)
let flow = try AerodynamicFlow(
    relativeVelocityMetresPerSecond: SIMD3(100, 0, 0), properties: ambient)
let reynolds = try flow.reynoldsNumber(referenceLengthMetres: 1)
let local = try flow.atLocalPoint(
    offsetMetres: SIMD3(0, 2, 0), angularVelocityRadiansPerSecond: SIMD3(0, 0, 0.1))
```

Local velocity adds `ω×r` in a spatially uniform medium. Callers supply actual winds,
atmospheric gradients and coordinate transformations; no orientation or ephemeris is
inferred. The units and definitions follow [NASA dynamic pressure](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/dynamic-pressure-2/)
and [Reynolds number](https://www.grc.nasa.gov/WWW/K-12/airplane/reynolds.html).

## Forces and moments

``AerodynamicCoefficients`` carries three signed force and three moment components.
``AerodynamicReferenceGeometry`` states area, three moment reference lengths and a
moment reference point. Aircraft moment lengths typically use span/chord/span.
``AerodynamicLoads`` contains force in N and moment in N m in the same caller-selected
axes. Coefficient conversion uses `F = q S C` and `M = q S L Cmoment`.
The moment already includes the full coefficient moment about its stated reference
point; do not add a second centre-of-pressure moment for the same load.

Transferring the moment to another point uses `Mnew = Mold + (rold−rnew)×F`.
Load addition first transfers the other moment to the receiving load's reference
point. This operation aggregates physics without owning a vehicle assembly.
Reference points and force/moment vectors must use the same origin and axes;
the types do not perform frame transforms or certify coordinate compatibility.

``AerodynamicHeatFlux`` provides a surface position and positive inward heat flux in
W/m² for future correlations and external thermal integration. It does not compute
heating, evolve temperature or manage ablation.

## NACA geometry

``NACA4`` computes thickness offsets normal to the actual camber-line tangent,
including a zero-camber limit without division by zero. Coordinates are normalized
by chord; evaluation covers `0...1`. The conventional thickness polynomial retains
a finite trailing-edge gap. Use the throwing geometry initializer and
`surfaceCoordinates(atChordFraction:)` for recoverable input validation; legacy
methods retain preconditions. The basis is [NASA-TM-X-3284](https://ntrs.nasa.gov/citations/19760003945).

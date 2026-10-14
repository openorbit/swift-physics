# Celestial body metadata

`BodyIdentifier` is a namespaced identity, independent of a display name, frame,
or ephemeris model. Named Solar System identifiers use [NAIF body codes](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html).
Planet centres and system barycentres are distinct even if a model approximates
them with the same coordinates. Additional NAIF codes use `.naif(Int32)`; game
objects use `.custom(namespace:name:)`. Custom names and namespaces are case-sensitive
and preserved exactly; empty or whitespace-only components fail when constructing
a `CelestialBody`. No automatic body-name parsing or code inference occurs.

`CelestialBody` is immutable metadata: identifier, display name, kind, optional
physical properties, and provenance. Use its identifier as the cross-model lookup
key; value equality also compares the metadata. It does not own an ephemeris,
orientation, shape evaluator, frame, or gravity model. Legacy `PlanetaryBody` and
`IrregularPlanetaryBody` remain separate, with their original API and conventions.

## Properties and units

| Property | Units / meaning |
| --- | --- |
| `gravitationalParameter.cubicMetresPerSecondSquared` | μ = GM, m³/s² |
| `massKilograms.value` | Independently supplied mass, kg |
| `referenceRadii.metres` | Reference X/Y/Z radii, m, with an explicit axes description |

Values must be positive and finite. Missing properties are `nil`, not zero.
Each supplied property carries a source and version, plus an optional reference.
The body's own provenance identifies its definition separately from individual
property sources. Reference radii are geometric metadata, not a shape or orientation
model; a barycentre cannot have reference radii.

A gravitational parameter also requires an explicit coordinate-time convention:
TDB-compatible, TT-compatible, TCB-compatible, or a caller-defined Newtonian model.
Converting km³/s² to m³/s² does not change that convention. No relativistic rescaling
or mass-to-GM conversion is performed automatically. Later dynamics must use a
compatible parameter and coordinate-time formulation; SI state units alone do not
establish compatibility.

## Bundled catalogue and injection

```swift
let earth = try BodyCatalogue.jplGM2022.body(.earth)
let mu = earth.properties.gravitationalParameter!
// mu.cubicMetresPerSecondSquared is in m³/s²; convention is TDB-compatible.

let source = try PhysicalDataProvenance(source: "Game world", version: "1")
let parameter = try GravitationalParameter(
  cubicMetresPerSecondSquared: 1e12, convention: .newtonian, provenance: source)
let body = try CelestialBody(
  identifier: .custom(namespace: "my-game", name: "terra"),
  name: "Terra", kind: .custom,
  properties: BodyPhysicalProperties(gravitationalParameter: parameter),
  provenance: source)
let catalogue = try BodyCatalogue(bodies: [body])
```

`jplGM2022` is an explicit snapshot of selected [NAIF gm_de440.tpc values](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc),
dated 2022-12-14. It includes the Sun, nine planetary system barycentres, nine
planet/dwarf-planet centres, and the Moon, plus the Solar System barycentre without
a GM. The subset contains DE440 values and separately labelled satellite-solution
values. JPL's [ephemeris convention](https://ssd.jpl.nasa.gov/planets/eph_export.html)
uses TDB; these parameters retain the TDB-compatible convention.
No masses or radii are inferred and no downloads occur during evaluation.
The snapshot name is deliberate; it does not promise the latest constants.

A catalogue rejects duplicate identifiers, permits duplicate display names, preserves
input order in `bodies`, and throws `unknownBody` on absent identities. Inject a new
catalogue to select another solution or game configuration; there is no mutable
global registry. Metadata lookup does not establish ephemeris coverage or availability.

The kernel subset fixture is extracted from the official file independently of the
Swift constants. Regenerate it with `Tools/generate_body_gm_fixtures.py`. Swift
Testing checks every converted value, centre/system distinctions, the Earth–Moon
GM relationship, custom identities, injection and invalid inputs.

## Independent body models

`CelestialBodyModels` selects optional ephemeris, orientation, shape and gravity
providers for one body. Every supplied provider must declare the same body identity.
Absent models stay absent even when metadata contains GM or reference radii. Each
provider retains its own provenance, coverage and coordinate conventions. A matching
identity does not prove compatible model versions or matching body-fixed axes.

| Interface | Responsibility |
| --- | --- |
| `EphemerisProvider` | Position/velocity relative to an explicit origin, native frame and epoch coverage |
| `BodyOrientationProvider` | Reference-to-body rotation and angular derivatives, body identity and source |
| `BodyShapeProvider` | Surface geometry in described body-fixed axes |
| `BodyGravityProvider` | Gravity field at a supplied body-centred displacement, with explicit time convention |

`IdentifiedBodyOrientationProvider` labels an existing `OrientationProvider`, such
as an IAU or prescribed spin provider. It preserves rotations and derivatives and
checks returned epochs and reference frames. The caller must identify the actual
body represented by that orientation; assigning an identity does not change a model.

```swift
let body = try BodyCatalogue.jplGM2022.body(.earth)
let source = try PhysicalDataProvenance(source: "Game geometry", version: "1")
let shape = try EllipsoidBodyShape(
  bodyIdentifier: .earth,
  referenceRadii: BodyReferenceRadii(
    metres: SIMD3(6_378_000, 6_378_000, 6_357_000),
    axesDescription: "Selected Earth body-fixed XYZ; illustrative game dimensions",
    provenance: source))
let gravity = try PointMassBodyGravity(
  bodyIdentifier: .earth, parameter: body.properties.gravitationalParameter!)
let models = try CelestialBodyModels(body: body, shape: shape, gravity: gravity)
let surface = try shape.surface(inDirection: SIMD3(1, 1, 0))
let field = try gravity.gravity(at: .j2000, displacementMetres: SIMD3(7e6, 0, 0))
```

`EllipsoidBodyShape` intersects a ray from the centre with a positive triaxial
ellipsoid. Equal radii define a sphere. The result is a surface position in metres
and an outward unit gradient normal. The normal generally differs from the radial
direction. Input direction magnitude is irrelevant; zero and nonfinite directions
fail. This query does not compute geodetic coordinates, closest points, signed
distance, altitude, arbitrary mesh intersections, or collision response.

`PointMassBodyGravity` supplies acceleration −μr/|r|³ and potential −μ/|r|, with zero
potential at infinity. It is spherically symmetric, so any orthonormal body-centred
axes are valid. It does not inspect the shape or clamp at a surface: its mathematical
field continues inside the supplied ellipsoid, with a singularity at zero radius.
Nonfinite positions and results fail explicitly. Double underflow/overflow and
extreme ellipsoid aspect ratios limit representable queries.

Gravity seconds follow the selected gravitational parameter's convention. The
sample carries that convention; it does not silently reinterpret TDB-compatible
parameters as atomic-time acceleration. A later dynamics adapter must select a
compatible formulation. No rigid-body force application or numerical integration
is performed by this field interface.

Convert world states into the intended body-centred axes explicitly through the
frame system before calling a shape or nonspherical gravity provider. Shape and
gravity providers do not fetch ephemerides, create frames, or infer orientation.
Shape axes descriptions must match the caller's selected body-fixed frame; the
model bundle validates identity only. This separation permits multiple gravity or
shape approximations for the same body without changing its motion or attitude.

Tests cover sphere and ellipsoid reference intersections/normals, direction scaling,
point-mass inverse-square values and independently differentiated potential,
singular inputs, absent models, mismatched identities, and orientation wrappers.

# MagneticFields

Swift target with no dependencies, available via this repository's
`MagneticFields` and `Physics` products. All fields are **magnetic flux
density B in tesla**, positions in metres, in the model's stated body-fixed axes.

```swift
import MagneticFields
let earth = BodyMagneticField.model(for: .earth)
let query = try MagneticFieldQuery(
    positionMetres: SIMD3(earth.dipole.referenceRadiusMetres, 0, 0),
    decimalYear: 2026.5)
let field = try earth.checkedField(at: query)
// Earth magnetic equator: nominal 30.6 µT in the +Z direction.
```

`MagneticFieldModel` is the interchangeable interface. It accepts a body-centred
Cartesian position and calendar decimal year and returns Cartesian B in the same
axes. Models expose source, version, assumptions, axis convention, geocentric-radius
coverage and optional date coverage. `checkedField` enforces coverage for any provider.
A nil date range explicitly denotes a time-independent model. Spatial coverage is a
chosen evaluation bound, not a surface, collision boundary or magnetopause.

`DipoleMagneticField` supports arbitrary tilt, signed polarity and displacement.
Its moment-axis vector is `μ₀ m/(4π a³)` in tesla at reference radius a. An equatorial
field is opposite that vector; the axial field has twice its magnitude. Displaced
dipoles use displacement from the dipole for evaluation while radius coverage remains
relative to the body centre. Evaluation at the dipole centre throws.
`ConstantMagneticField` supplies prescribed fields for controlled scenarios.

Presets: Mercury, Earth, Jupiter, Saturn, Uranus and Neptune. These are deliberately
centred, aligned, static dipoles with nominal strengths from the
[NASA fact sheets](https://nssdc.gsfc.nasa.gov/planetary/factsheet/).
They omit measured tilts and offsets, especially important for Uranus and Neptune,
as well as multipoles, secular variation and external magnetospheric fields.
Signed moments are -Z for Earth/Mercury and prescribed +Z for giant presets;
these are simulation conventions, not full reconstructions of measured polarities.
Coverage is 1...10 reference radii. Unsupported bodies and system barycentres return
nil, meaning no supplied intrinsic model. Venus, Mars and moons can still experience
external or induced fields; absence of a preset never implies zero total field.

## WMM and IGRF

`WMMMagneticField` reads standard six-column WMM `.COF` files. `IGRFMagneticField`
reads traditional IGRF g/h epoch tables. Both use an immutable Swift internal-field
harmonic evaluator with local scratch arrays and no C/Fortran runtime dependency.
WMM applies secular variation from its reference epoch; IGRF interpolates between
epochs and uses the final SV column after the last epoch. Dates outside the file's
coverage throw. Harmonic degrees 1...64 are supported; WMMHR is not supported.

```swift
import Foundation
import MagneticFields
let wmm = try WMMMagneticField.wmm2025() // Bundled NOAA WMM2025, 2025...2030
let igrf = try IGRFMagneticField.igrf14() // Bundled IAGA IGRF-14, 1900...2030
let updated = try WMMMagneticField(coefficientFile: URL(fileURLWithPath: "/path/to/WMM.COF"))
let local = try wmm.elements(at: GeodeticMagneticQuery(
    latitudeRadians: 52.0 * .pi / 180, longitudeRadians: 5.0 * .pi / 180,
    ellipsoidHeightMetres: 0, decimalYear: 2026.5))
// local.northEastDownTesla, totalIntensityTesla, declinationRadians, inclinationRadians
```

File loaders are explicit I/O; evaluators never fetch data. New coefficient editions
in these formats can be loaded without modifying the evaluator. Parsers reject
missing, duplicate, malformed and nonfinite coefficients. `.COF` files have no end
date, so WMM loading accepts an explicit `validityYears` (default five).
Geodetic queries use WGS84 ellipsoid heights, radians and metres; supported heights
are −1 to 850 km. Geocentric queries retain the common XYZ interface. The 6000...7500 km
radius bound is numerical coverage, not a promise of physical accuracy throughout.
There is no geoid correction, external-field model, crustal anomaly model, grid
variation output or uncertainty estimator. Compass declination at geographic poles
has no unique meaning even though Cartesian fields remain well defined.

The bundled IGRF-14 table is attributed to IAGA:
[IGRF-14 coefficients](https://doi.org/10.5281/zenodo.14012302),
[NOAA coefficient distribution](https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt).
WMM2025 coefficients/test data come from NOAA's public-domain download.
Tests compare the full published WMM2025 vector table and 48 independent IGRF-14
Fortran reference vectors, including historical interpolation and polar cases.

[NOAA WMM documentation](https://www.ncei.noaa.gov/products/world-magnetic-model)
describes five-year coefficient validity and calendar decimal years. Calendar decimal
year means year + fraction of that calendar year, including its actual length; it is
not `(JulianDate-J2000)/365.25 + 2000`. The geographic helper explicitly uses WGS84 ellipsoid height; other
coordinate adapters must state their ellipsoid and altitude datum.

## CelMek integration

`BodyIdentifier.earth.magneticModel` provides lookup. `ModelBodyMagneticField` wraps
any provider with body identity and an explicit `(Epoch) throws -> Double` date mapping:

```swift
import CelMek
import MagneticFields
let model = try ModelBodyMagneticField(
    bodyIdentifier: .earth, model: BodyMagneticField.model(for: .earth),
    decimalYear: { _ in 2026.5 }) // Fixed date is suitable here only because the dipole is static.
let sample = try model.magneticField(at: .j2000,
    displacementMetres: SIMD3(7_000_000, 0, 0))
```

`CelestialBodyModels` accepts an optional `BodyMagneticProvider`, enforcing body
identity. `OrientedBodyMagneticField` combines that provider with a matching
`BodyOrientationProvider`: it rotates a body-centred reference position into the
model's axes, evaluates B, then rotates B back into reference axes. Epoch and frame
identities are checked. This adapter also implements CelMek's existing
`MagneticFieldProvider`, so it can be passed directly to `MagnetorquerTorque`.
Axes must agree with the supplied orientation; no ephemeris
translation, calendar conversion or implicit orientation selection occurs.

Run `swift test --filter 'MagneticFieldTests|GeomagneticTests'` from the repository root.

# Frame transforms

`FrameTree` connects `Math.Frame` identities with time-dependent `TransformProvider` edges. Each edge maps parent coordinates into child coordinates. A tree can adopt a frame such as `CelestialFrames.vsop87E` as its local root, excluding its external ancestors. `CelestialFrames.makeTree()` registers the connected astronomical frames under ICRF. Queries to frames outside the registered tree throw `disconnectedFrames`.

```swift
import CelMek
import Math

let tree = try FrameTree(root: CelestialFrames.vsop87E)
let rotating = try tree.makeFrame(
  parent: tree.root, name: "Rotating ecliptic axes",
  provider: UniformRotationProvider(
    axis: SIMD3<Double>(0, 0, 1), rate: 0.001))
let state = try earth.state(at: .j2000)
let rotated = try tree.convert(state, to: rotating)
```

This example rotates barycentric axes about a common origin. Body-fixed frames combine an origin state with an orientation provider as described below; validated planetary orientation models remain follow-up work.

## Convention and derivatives

A `FrameTransform` is evaluated at one `Epoch`, using right-handed Cartesian coordinates and normalized Hamilton quaternions `(x, y, z, w)` acting on column vectors:

```
p_to = R p_from + t
v_to = R v_from + ω × (R p_from) + ṫ
a_to = R a_from + 2ω × (R v_from)
       + α × (R p_from) + ω × (ω × (R p_from)) + ẗ
```

Translation and its derivatives use metres, metres per second, and metres per second squared. Angular velocity ω and angular acceleration α use radians per second and radians per second squared, expressed in destination axes. All derivatives use elapsed atomic seconds, matching `Epoch`. Translation is the mapped source origin in destination coordinates, rather than the child origin in parent coordinates.

`followed(by:)` composes mappings in application order, and `inverted()` reverses the full mapping including derivatives. Epoch mismatches throw. Invalid rotations and nonfinite transform components throw; a valid non-unit quaternion is normalized. State conversion preserves unknown acceleration as nil.

`FixedTransformProvider` supplies a constant rotation/translation with zero derivatives. Its defaults provide an identity edge. `UniformRotationProvider` turns child axes by the right-hand rule about a fixed axis expressed in parent coordinates. Its coordinate mapping rotates with the opposite sign. `initialRotation` maps the parent axes to the child at `referenceEpoch`; the resulting angular rate is expressed in child axes. It has a stationary common origin. Very long elapsed intervals still have Double phase precision limits.

Custom providers implement `transform(at:)`, returning a parent-to-child mapping at the requested epoch. They can supply translation, angular acceleration, or data-backed motion, and throw their own coverage/data errors. No fallback or cache is implicit. Queries traverse only the edges between the two frames and their common ancestor, so an unavailable edge above that ancestor does not block the query.

`FrameTransform.applying(to:destination:)` is a low-level operation: callers are responsible for matching the source coordinates to their mapping. `FrameTree.convert` selects the mapping from the state's frame identity and handles this association.

## Ownership

Frame parents are immutable strong references: a surviving child retains its ancestry. Parents do not retain children, so the parent chain does not form a reference cycle. `FrameTree` owns its registered frames and providers, which keeps frame identities stable. Frames do not retain the tree, and a child can outlive it, though provider queries require a retained tree.

A custom provider must not strongly capture its owning tree: that would create a tree → provider → tree cycle. Use a weak reference when such access is necessary. The built-in providers hold no tree references. These types do not promise concurrent mutation or Sendable conformance.

The legacy `Math.Transformation` remains separate from the epoch-aware `FrameTransform`; its Julian-date field and misspelled properties are not used by this API.

## Body-fixed frames

`OrientationProvider` returns a `BodyOrientation`: the mapping from a specified reference frame's axes into body axes, plus angular velocity and angular acceleration in body axes. It contains no translation. `FixedAxisOrientationProvider` supplies prescribed constant-axis spin; its initial rotation and rate follow `UniformRotationProvider` conventions. This is useful for custom game bodies, without claiming a planetary model.

`BodyFixedTransformProvider` combines a centre `StateProvider` with an orientation provider. Their reference frame identities must match. It constructs `p_body = R * (p_parent - centre)`, including origin velocity/acceleration and all rotational cross terms. `FrameTree.makeBodyFixedFrame` registers this mapping under that common reference frame. The returned frame is not marked pseudo-inertial.

```swift
let tree = try FrameTree(root: CelestialFrames.vsop87E)
let spin = try FixedAxisOrientationProvider(
  referenceFrame: tree.root, axis: SIMD3<Double>(0, 0, 1), rate: 0.001)
let body = try tree.makeBodyFixedFrame(
  name: "Earth-centred prescribed spin", centre: earth, orientation: spin,
  accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 60)))
let relative = try tree.convert(earth.state(at: .j2000), to: body)
// Position and velocity are zero, since the queried object is the origin.
```

The example's axes are prescribed in the VSOP87 ecliptic frame; they do not represent Earth's geographic orientation or ITRF. Real planetary orientations must declare their reference axes and be transformed into the origin provider's frame before combination. The IAU tables and evaluator have been updated and checked against independent CSPICE fixtures; see [IAU_ROTATION_AUDIT.md](IAU_ROTATION_AUDIT.md). `IAUOrientationProvider` now exposes those orientations in ICRF axes with explicit caller-selected coverage, as described below.

Origin acceleration defaults to `requireProvided`. When it is absent, the query throws `accelerationUnavailable`. Callers may explicitly choose `assumeZero` or `estimate(step:)`; the latter centrally differences velocities at neighbouring epochs using a positive elapsed-second interval. Supplied acceleration always takes precedence. Estimation introduces numerical error and requires coverage at both neighbouring epochs; its interval is caller-selected rather than an accuracy guarantee. Provider errors propagate, and returned epochs and frame identities are validated at every sample. Unknown acceleration in the input state remains unknown after conversion, even if origin acceleration was estimated.

Choose a local body-centred parent where practical when operating on small distances: subtracting large barycentric coordinates still incurs ordinary Double cancellation.

## Supported astronomical frames

`CelestialFrames.makeTree()` registers shared immutable frame identities and fixed edges. Models retain their native frame. Convert on query:

```swift
let tree = try CelestialFrames.makeTree()
let state = try earth.state(at: .j2000, in: CelestialFrames.icrf, using: tree)
```

| Frame | Origin | Axes / connection |
| --- | --- | --- |
| `icrf` | Solar-system barycentre | ICRF-aligned equatorial axes; the tree root |
| `meanEquatorJ2000` | Solar-system barycentre | Mean dynamical equator/equinox J2000 using the IAU 2000 fixed frame-bias angles |
| `eclipticJ2000` | Solar-system barycentre | Ecliptic/equinox J2000 defined from EME2000 with the IAU 1976 obliquity, 84381.448 arcsec |
| `fk5J2000` | Solar-system barycentre | FK5 axes frozen at J2000, connected through the fixed Hipparcos/ICRS orientation tie |
| `vsop87E` | Solar-system barycentre | Native VSOP87E dynamical ecliptic/equinox J2000, connected through its published FK5 matrix |
| `makeGCRF(earth:)` | Earth's centre | Newtonian geocentric frame with ICRF-aligned axes; caller supplies geocentre motion |
| `makeGeocentricFrames(earth:).eme2000` | Earth's centre | EME2000 mean equator/equinox J2000, sharing GCRF's moving origin |

The [original VSOP87 description](https://raw.githubusercontent.com/ctdk/vsop87/master/vsop87.txt) publishes its FK5 conversion matrix. A quaternion reconstructed from that rounded matrix preserves orthogonality; reference tests agree within 8e-13 per matrix entry. The native frame is distinct from `eclipticJ2000`.

The fixed frame bias follows the IAU 2000 convention described by [ERFA BP00](https://github.com/liberfa/erfa/blob/master/src/bp00.c). The separate FK5-to-ICRF connection uses the fixed J2000 orientation in [ERFA FK5HIP](https://github.com/liberfa/erfa/blob/master/src/fk5hip.c). This is an approximate frame tie through the Hipparcos realization, without its time-dependent stellar catalogue spin or zonal-error corrections; it does not imply a precision realization of ICRF3. Matrix tests use the independent ERFA reference values.

`meanEquatorJ2000` and `eclipticJ2000` have barycentric origins. `makeGeocentricFrames` creates GCRF and EME2000 with Earth-centred origins, sharing the same translation edge. `makeGCRF` supplies only Newtonian translation, not the full relativistic BCRS/GCRS transformation. Its acceleration policy has the same explicit choices as body-fixed frames.

```swift
let geocentric = try tree.makeGeocentricFrames(
  earth: earth, accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 60)))
let earthCentre = try earth.state(at: .j2000, in: geocentric.eme2000, using: tree)
// Position and velocity are zero: EME2000 has Earth's moving origin.
```

`FrameConvertedStateProvider` presents a native provider in a connected destination frame, for composition with origin/orientation providers. It holds the tree weakly to avoid ownership cycles; the caller must retain the tree. A query after tree release throws `frameTreeUnavailable`. `registerFrame(_:provider:)` registers an existing identity under its immutable parent, without reparenting. Re-registering an identity or connecting it before its parent throws.

## IAU orientations in the frame tree

`IAUOrientationProvider` declares ICRF axes directly. [NAIF's frame documentation](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/frames.html) explains that generic text-PCK IAU orientations use ICRF axes despite the SPICE label J2000; that label is not used here to alias EME2000.

```swift
let start = try Epoch(julianDate: 2451545, scale: .tdb)
let end = try start.advanced(by: CelMek.Duration(seconds: 86400))
let orientation = try IAUOrientationProvider(body: "mars", coverage: start...end)
let marsICRF = try FrameConvertedStateProvider(mars, in: CelestialFrames.icrf, using: tree)
let marsFixed = try tree.makeBodyFixedFrame(
  name: "IAU Mars", centre: marsICRF, orientation: orientation,
  accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 60)))
let bodyState = try mars.state(at: start, in: marsFixed, using: tree)
```

Coverage is required and caller-selected; it is not a certified physical-validity interval. Origin acceleration estimation also needs neighbouring ephemeris samples. `availableBodies` lists the existing 73 models. Earth and Moon require `allowLegacyModel: true` to expose the labelled older approximations; this does not provide ITRF or a precision lunar attitude.

The provider converts TDB angular derivatives using both first and second derivatives of the provisional `Epoch` time adapter. Public-provider regression tests use all 511 CSPICE samples, including rotation rate and acceleration; the physical accuracy limits of that time adapter remain documented in [TIME.md](TIME.md).

## Earth frames

`makeEarthFrames` adds both IAU 2006/2000A frame chains, using native Swift calculations. The algorithms and coefficient series are adapted
from [ERFA 2.0.1](https://github.com/liberfa/erfa/tree/v2.0.1); attribution and BSD terms
are retained in `Documentation/licenses/ERFA.txt`. There is no C target or ERFA
runtime dependency. ERFA is used as an external reference for regression fixtures.

| Frame | Parent | Transform |
| --- | --- | --- |
| CIRF | GCRF | Bias/precession/nutation and CIO locator |
| TIRF | CIRF | Earth rotation angle from UT1 |
| ITRF (caller realization label) | TIRF | Polar motion and terrestrial origin locator s′ |
| TOD | GCRF | Bias/precession/nutation, true equator/equinox of date |
| GTOD | TOD | Greenwich apparent sidereal time consistent with the CIO chain |

```swift
let eop = try ConstantEarthOrientationProvider(
  ut1MinusTT: CelMek.Duration(seconds: -64.184),
  polarMotionX: 0, polarMotionY: 0,
  coverage: start...end, provenance: "Explicit constant approximation")
let frames = try tree.makeEarthFrames(
  earth: earth, earthOrientation: eop, terrestrialRealization: "ITRF2020",
  accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 60)))
let terrestrial = try earth.state(at: .j2000, in: frames.itrf, using: tree)
```

The example offset is illustrative, not a current Earth-orientation measurement.
Implement `EarthOrientationProvider` for measured EOP with your selected interpolation,
coverage and provenance. Samples supply UT1−TT in seconds and polar motion in radians;
this avoids an implicit leap-second table. No UTC/DUT1 file ingestion is included yet.
The constant provider is an explicit approximation and rejects extrapolation.

Angular velocity and acceleration use a five-point matrix stencil in atomic SI
seconds, with a configurable 60-second default step. Provider coverage must include
±2 steps around every query, so nominal coverage endpoints cannot be queried with
that stencil. EOP interpolation discontinuities affect derivatives; choose a smooth
provider when acceleration matters. No derivative accuracy guarantee is implied
for arbitrary step sizes or interpolators.

The initial implementation uses the theoretical celestial intermediate pole;
observed dX/dY corrections remain future work. ITRF naming records the caller's
realization: these rotations do not implement station coordinates, tectonic motion,
or Helmert ties between realizations. GCRF retains the Newtonian translation limitation
above. Tests check published ERFA matrices, CIO/equinox route consistency, rotating
state derivatives and explicit coverage failures.

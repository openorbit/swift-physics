# Frame transforms

`FrameTree` connects `Math.Frame` identities with time-dependent `TransformProvider` edges. Each edge maps parent coordinates into child coordinates. A tree can adopt a celestial root such as `CelestialFrames.vsop87E`; queries between unrelated trees throw `disconnectedFrames`.

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

The example's axes are prescribed in the VSOP87 ecliptic frame; they do not represent Earth's geographic orientation or ITRF. Real planetary orientations must declare their reference axes and be transformed into the origin provider's frame before combination. The IAU tables and evaluator have been updated and checked against independent CSPICE fixtures; see [IAU_ROTATION_AUDIT.md](IAU_ROTATION_AUDIT.md). Public planetary orientation integration and coverage policies remain follow-up work.

Origin acceleration defaults to `requireProvided`. When it is absent, the query throws `accelerationUnavailable`. Callers may explicitly choose `assumeZero` or `estimate(step:)`; the latter centrally differences velocities at neighbouring epochs using a positive elapsed-second interval. Supplied acceleration always takes precedence. Estimation introduces numerical error and requires coverage at both neighbouring epochs; its interval is caller-selected rather than an accuracy guarantee. Provider errors propagate, and returned epochs and frame identities are validated at every sample. Unknown acceleration in the input state remains unknown after conversion, even if origin acceleration was estimated.

Choose a local body-centred parent where practical when operating on small distances: subtracting large barycentric coordinates still incurs ordinary Double cancellation.

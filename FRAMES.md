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

This example rotates barycentric axes about a common origin. A body-fixed frame also needs an origin translation and a physical orientation provider; those astronomical models remain follow-up work.

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

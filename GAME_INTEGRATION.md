# Rendering through a moving origin

Keep simulation states in SI `Double` coordinates. `RendererAdapter` captures a
`CartesianState` origin and a positive number of metres per render unit. It
subtracts positions and derivatives in `Double`, scales them, and only then
converts to `Float`. This preserves nearby detail at solar-system distances;
precision still depends on the resulting relative distance and render scale.
Detail already lost in the input `Double` values cannot be recovered.

```swift
let origin = try spacecraft.state(at: clock.epoch)
let planet = try planetProvider.state(at: clock.epoch)
// First use tree.convert(planet, to: origin.frame) if frames differ.
let renderer = try RendererAdapter(origin: origin, metresPerRenderUnit: 1e9)
let rendered = try renderer.adapt(planet)
// Pass rendered.position to the renderer; spacecraft is at zero.
```

Both states must have the same epoch and identical frame object. Use the frame
system for orientation and physical centre conversion before adapting. The
adapter preserves the frame's right-handed axes and translates a moving origin;
it does not rotate axes. Rebuild it at each new epoch or origin change, and use
that snapshot for all objects in a rendered scene. Historical orbit samples
require an explicit choice of origin at each sample epoch; stale origins are
rejected.

Velocity is relative to the origin, in render units per simulation second.
Acceleration is relative only when both states provide it; otherwise it is nil.
Neither derivative incorporates playback warp. These render values do not carry
frame or epoch metadata and must remain presentation data, rather than inputs to
physical propagation. Scale errors and nonfinite differences, scaling results,
or Float overflow throw; small values may round or underflow to zero.

The [Voyager example](Examples/Voyager2GameDemo/README.md) exports five planet
positions per playback snapshot using the spacecraft as a moving origin and
one billion metres per render unit. No graphics framework is required.

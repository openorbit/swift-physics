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

## Cache invalidation

Cache validity follows the entire configuration that produced a result. There is
no automatic detection of dependency edits. Treat configured providers and assets
as immutable snapshots; replace a snapshot when changing its configuration.
`CachedEphemerisProvider` supplies exact-epoch native-state caching, with private
scopes that separate wrappers even when metadata looks identical. See
[ephemeris cache contracts](EPHEMERIDES.md#batch-queries-and-explicit-caches).

| Change | Required action |
| --- | --- |
| Initial state, central mass, force model, burn schedule, attitude inputs or solver options | Construct the new propagator/planner and discard predictions, sampled paths and mission plans derived from the old inputs. |
| Provider selection, frame, physical origin, coverage, time model or requested accuracy | Build a new configured provider and a new cache wrapper. Rebuild dependent frame providers and derived buffers. A lower accuracy hint may still return higher fidelity; a changed hint is still a configuration change. |
| Mutable dependency updated while provider metadata remains fixed | Serialize the update with provider evaluation, call `cached.invalidate()`, and discard downstream results. Every affected wrapper needs invalidation. |
| SPICE load, unload, replacement or clear | Obtain a new kernel snapshot and rebuild affected providers. Namespace external snapshot caches by `cacheScope` plus provider configuration. Existing providers retain their old snapshot; a kernel-set mutation does not update them. |
| DE441 or SPICE binary asset replacement | Open a new reader/provider for the new immutable asset, wrap it anew and rebuild dependents. Do not modify an asset in place while a reader is alive. |
| Frame-tree bindings, orientation, centre providers or leap-second/time-conversion inputs | Rebuild affected configured components. Discard converted states and samples; invalidate native states only if their evaluation uses the changed dependency. |
| Playback epoch, pause, rate or seek | Query the new exact epoch and update presentation. Unchanged native providers can keep their cached states. Discard pending predictions selected under the old playback request. |
| Render origin, axes or metres per render unit | Rebuild the render snapshot and rendered buffers. Native SI states remain usable when their physical configuration is unchanged. |

`TrajectorySampler`, `SimulationClock`, `RendererAdapter` and frame-tree queries
introduce no internal result caches. Application-owned path arrays, search results,
converted states and GPU buffers still need explicit replacement. A native state
cache cannot invalidate these downstream products for the application.

For an atomic scene change, stage the replacement configuration first, then publish
it together with a new application-owned revision token. Tag asynchronous jobs with
that token and accept their output only if it matches the current revision. An
opaque `EphemerisCacheScope()` can serve as an in-process token. Do not reuse it for
a different configuration or serialize it as a durable version. Stop or synchronize
old evaluations before mutating a shared dependency; a token check alone does not
make a mutable provider safe to read concurrently.

`cached.invalidate()` rotates its generation before removing old entries. A query
already in flight can return an old result, but cannot repopulate the old generation
after invalidation. Discard that result at the scene revision boundary. If custom
storage cleanup throws, the new generation still applies; handle the error and
release the stale storage separately. Clearing `MemoryEphemerisCache.removeAll()`
only releases storage and does not rotate wrapper generations, so it is not a
substitute for invalidation during a dependency update.

DE441's `clearRecordCache()` evicts decoded records from its current reader and
preserves cumulative statistics. It does not reload headers or establish a new data
version. For replacement wrappers sharing bounded memory storage, old entries can
remain until eviction; retain and invalidate the old wrappers when prompt release
is required, or clear the scene's shared storage after old jobs finish. Independent
scenes should generally use separate storage when their lifetimes differ.

Representative release workloads and measured tradeoffs are documented in
[GAME_PERFORMANCE.md](GAME_PERFORMANCE.md).

For component synchronization, frame/asset lifetimes and actor isolation, see
[THREAD_SAFETY.md](THREAD_SAFETY.md).

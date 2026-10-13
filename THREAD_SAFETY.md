# Thread safety and data ownership

The library has no global simulation clock, task scheduler or mutable kernel pool.
Applications own their configured models and choose where to evaluate them. Swift
`Sendable` conformance and runtime synchronization are separate contracts: a lock
inside one dependency does not make the enclosing provider graph safe to share.
`StateProvider`, `EphemerisProvider`, `TransformProvider`, `TDBTimeModel`,
`ForceModel`, `TorqueModel` and attitude-provider protocols do not require Sendable.
A struct holding one of these providers can still share reference-backed state.

## Component contracts

| Component | Ownership and concurrent evaluation |
| --- | --- |
| `Epoch`, `Duration`, `SimulationClock`, `RenderState` | Sendable value types. Copies of a clock advance independently. Serialize mutation of one shared clock variable. Passing a value does not synchronize shared storage holding it. |
| `CartesianState` and spacecraft/attitude states | Numeric fields are snapshots, but frame identities are retained references. They are not blanket Sendable transfer objects. Copying a state does not clone its frame or provider graph. |
| `Frame` / `FrameTree` | Parent links are immutable and strong. The tree retains registered frame identities and edge providers; registration mutates arrays/dictionaries without locking. Serialize construction and registration with all queries. A frozen tree's traversal reads registration storage, but safe overlapping queries also require every edge provider and dependency to support them. The tree has no Sendable conformance. |
| Built-in analytical models and two-body propagation | Configured coefficients/options are retained and evaluations use local working storage. Overlapping evaluation requires immutable configuration and safe injected time/frame/state dependencies. Provider types are not automatically Sendable merely because their own arithmetic is local. |
| Numerical translation, attitude and coupled propagation | Integration state, steps and diagnostics are local to each query. Repeated queries begin from the configured initial state; they do not advance a shared integrator cursor. Force, torque, guidance, environment and event callbacks must be repeatable and separately synchronized if shared. Never mutate captured configuration while evaluating. |
| `TrajectorySampler` and batch queries | Synchronous ordered calls on the invoking executor; no worker tasks are created. Caller-owned returned arrays retain their state/frame references. Independent calls inherit the queried provider's concurrency contract. |
| `MemoryEphemerisCache` | Storage, count, LRU and clearing are lock-protected. Capacity bounds all scopes in that instance. This does not serialize model evaluations. |
| `CachedEphemerisProvider` | Generation updates are synchronized. Concurrent misses may evaluate the same epoch repeatedly; no single-flight or duplicate suppression is promised. The wrapped provider must support overlap or be confined to one executor. Custom-cache callbacks must not reenter their invoking wrapper. |
| `DE441Ephemeris` | Random-access evaluation, record-cache mutations and statistics use the reader's lock. The open asset must remain unchanged. The reader owns its handle until deinitialization. Provider-level time conversion still inherits the supplied time model's contract. No blanket Sendable conformance is supplied. |
| `SpiceDAFFile` | Read-only metadata and synchronized seek/read operations; explicitly unchecked Sendable at this boundary. The asset must remain immutable for the reader and all snapshots retaining it. Its handle closes on deinitialization. |
| `SpiceKernelSet` | Load, unload, clear and snapshot operations serialize under a lock; explicitly unchecked Sendable at this boundary. Loads can hold the lock during filesystem IO, so do not run large loads on the UI executor. A successful mutation publishes a new cache scope. |
| `SpiceKernelSnapshot` | Sendable immutable pool/load metadata and retained synchronized binary readers. Providers retain their snapshot; unloading from the set does not update existing providers or close readers still retained elsewhere. Concurrent high-level queries depend on their injected time models and frame providers. |
| `RendererAdapter` | Immutable origin/scale snapshot; adaptation uses local arithmetic. The origin retains a frame identity, so the adapter itself is not a blanket Sendable container. Its numeric `RenderState` output is Sendable. It creates no GPU resources. |

For protocol implementations, document whether overlapping calls are supported,
whether callbacks mutate state, and what external resources must stay alive. A
lock around a coefficient file cannot protect a mutable guidance closure, frame
tree registration or application configuration accessed by the same evaluation.
The same applies to value-type wrappers around classes or captured closures.

## Frame and asset lifetimes

A child frame retains its parent. Parents have no child list. `FrameTree` owns
its edge providers and registered identities, so releasing a caller's child handle
does not allow identifier reuse inside the tree. A state retains its frame and
ancestors, but that does not retain the tree's provider registrations.

`FrameConvertedStateProvider` and SPICE SPK, binary-PCK and CK frame adapters hold
weak tree references to avoid ownership cycles. Keep a strong scene-level tree
reference while using them. A missing tree fails explicitly rather than reconstructing
registrations from a frame handle. Custom edge providers must not strongly capture
their owning tree. Rotating dynamics can retain a tree in their configuration; do
not install a provider that captures that dynamics object back into the same tree.

Configured providers retain their dependencies according to their own fields;
releasing a variable does not necessarily close a file or free coefficient arrays.
Cache entries retain states and frame identities until eviction or clearing. SPICE
snapshots retain binary readers after unload. Release old jobs, providers, snapshots
and scene storage when retiring an asset version. Replace assets using new immutable
files/readers rather than editing a live reader's file in place. For invalidation
and old jobs, see [game integration](GAME_INTEGRATION.md#cache-invalidation).

## A single-owner actor

The simplest concurrency boundary is an actor that constructs and owns its provider
graph. Other tasks pass Sendable input values and receive Sendable numeric outputs.
This example evaluates Earth relative to the Sun in ICRF, with explicit render scale:

```swift
import CelMek

actor SolarScene {
    private let tree: FrameTree
    private let earthModel: VSOP87Ephemeris
    private let sunModel: VSOP87Ephemeris

    init() throws {
        tree = try CelestialFrames.makeTree()
        let timeModel = try FairheadTDBModel()
        let coverage = FairheadTDBModel.referenceCoverage
        earthModel = try VSOP87Ephemeris(body: earth, coverage: coverage,
                                      tdbModel: timeModel, accuracy: .balanced)
        sunModel = try VSOP87Ephemeris(body: sun, coverage: coverage,
                                    tdbModel: timeModel, accuracy: .balanced)
    }

    func earthForRendering(at epoch: Epoch) throws -> RenderState {
        let origin = try tree.convert(sunModel.state(at: epoch), to: CelestialFrames.icrf)
        let target = try tree.convert(earthModel.state(at: epoch), to: CelestialFrames.icrf)
        return try RendererAdapter(origin: origin, metresPerRenderUnit: 1e9).adapt(target)
    }
}

// From another task: let rendered = try await scene.earthForRendering(at: epoch)
```

The query contains no await, so another actor message cannot interleave its
synchronous evaluation. Keep expensive path searches and asset loading off the main
actor. If an owner method suspends, another message can change its configuration
before it resumes; capture a revision and check it before publishing results.
For parallel numerical jobs, create independently owned graphs in each worker, or
share only components whose complete dependency graph has an explicit overlapping
query contract. Avoid adding unchecked Sendable to a whole scene just to bypass
compiler diagnostics. Frame identity comparisons require deliberately shared frame
objects within each graph; independently constructed graphs are not interchangeable.

## Validation scope

Existing Swift Testing cases exercise synchronized DAF reads, cache generation
changes during evaluation, independent clock copies and deterministic query order.
They do not certify arbitrary custom providers, callbacks, tree mutation or every
possible race. No new broad Sendable promises or locking are introduced by this
policy. Runtime read-only behavior is not a substitute for Swift's static isolation
requirements; keep non-Sendable provider graphs inside their owner.

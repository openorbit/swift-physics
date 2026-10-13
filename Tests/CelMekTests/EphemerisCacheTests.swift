import Math
import Testing
@testable import CelMek

@Suite struct EphemerisCacheTests {
  private enum TestError: Error { case transient, cacheFailure }
  private final class Model: EphemerisProvider {
    let metadata: EphemerisMetadata
    var calls = 0, base: Double = 0
    var fail = false
    var hook: (() throws -> Void)?
    init(frame: Frame, body: BodyIdentifier = .earth) throws {
      metadata = try EphemerisMetadata(body: body,origin: .solarSystemBarycentre,frame: frame,
        coverage: Epoch.j2000.advanced(by: CelMek.Duration(seconds: -10))...Epoch.j2000.advanced(by: CelMek.Duration(seconds: 10)),
        coverageBasis: .callerSelected,modelTimeScale: .tt,velocityMethod: .analytic,
        accelerationMethod: .unavailable,provenance: PhysicalDataProvenance(source: "deterministic test trajectory",version: "1"))
    }
    func state(at epoch: Epoch) throws -> CartesianState {
      calls += 1
      if fail { throw TestError.transient }
      try hook?()
      let t = try epoch.duration(since: .j2000).seconds
      return try CartesianState(epoch: epoch,frame: frame,position: SIMD3(base+t,0,0),velocity: SIMD3(1,0,0))
    }
  }
  private func epoch(_ seconds: Double) throws -> Epoch {
    try Epoch.j2000.advanced(by: CelMek.Duration(seconds: seconds))
  }
  @Test func exactEpochsAndOrderedBatches() throws {
    let model = try Model(frame: Frame(name: "reference",pseudoInertial: true))
    let cache = try MemoryEphemerisCache(capacity: 10)
    let cached = try CachedEphemerisProvider(provider: model,cache: cache)
    let epochs = try [epoch(1),epoch(-1),epoch(1),epoch(0),epoch(1e-9)]
    let states = try cached.states(at: epochs)
    #expect(states.map(\.epoch) == epochs)
    #expect(model.calls == 4 && cache.count == 4)
    #expect(states[3].position != states[4].position)
    let reversed = try cached.states(at: Array(epochs.reversed()))
    #expect(reversed == Array(states.reversed()))
    #expect(model.calls == 4)
    #expect(try cached.states(at: []).isEmpty)
    let collection = try EphemerisCollection(providers: [cached])
    #expect(try collection.states(of: .earth,at: epochs) == states)
    let requests = epochs.map { EphemerisRequest(body: .earth,epoch: $0) }
    #expect(try collection.states(for: requests) == states)
    #expect(throws: EphemerisError.outsideCoverage) { try cached.state(at: epoch(11)) }
    #expect(model.calls == 4)
  }
  @Test func leastRecentlyUsedEvictionIsBounded() throws {
    let model = try Model(frame: Frame(name: "reference",pseudoInertial: true))
    let cache = try MemoryEphemerisCache(capacity: 2)
    let cached = try CachedEphemerisProvider(provider: model,cache: cache)
    for t in [0.0,1,0,2,1] { _ = try cached.state(at: epoch(t)) }
    #expect(model.calls == 4 && cache.count == 2)
    cache.removeAll()
    #expect(cache.count == 0)
    _ = try cached.state(at: .j2000)
    #expect(model.calls == 5)
    #expect(throws: EphemerisCacheError.invalidCapacity) { try MemoryEphemerisCache(capacity: 0) }
  }
  @Test func sharedStorageAndExplicitInvalidation() throws {
    let frame = Frame(name: "reference",pseudoInertial: true)
    let a = try Model(frame: frame), b = try Model(frame: frame)
    a.base = 10; b.base = 20
    let cache = try MemoryEphemerisCache(capacity: 10)
    let first = try CachedEphemerisProvider(provider: a,cache: cache)
    let second = try CachedEphemerisProvider(provider: b,cache: cache)
    #expect(try first.state(at: .j2000).position.x == 10)
    #expect(try second.state(at: .j2000).position.x == 20)
    a.base = 30
    try first.invalidate()
    #expect(cache.count == 1)
    #expect(try second.state(at: .j2000).position.x == 20)
    #expect(b.calls == 1)
    #expect(try first.state(at: .j2000).position.x == 30)
    #expect(a.calls == 2)
    // Invalidation during provider evaluation must prevent the old generation from being stored.
    try first.invalidate()
    a.hook = { try first.invalidate() }
    _ = try first.state(at: .j2000)
    a.hook = nil
    #expect(cache.count == 1)
    _ = try first.state(at: .j2000)
    #expect(a.calls == 4 && cache.count == 2)
  }
  private final class CleanupFailingCache: EphemerisStateCache {
    let storage = try! MemoryEphemerisCache(capacity: 10)
    func state(for key: EphemerisCacheKey) throws -> CartesianState? { storage.state(for: key) }
    func store(_ state: CartesianState, for key: EphemerisCacheKey) { storage.store(state, for: key) }
    func removeAll(in scope: EphemerisCacheScope) throws { throw TestError.cacheFailure }
  }
  @Test func failedCleanupStillSeparatesNewGeneration() throws {
    let model = try Model(frame: Frame(name: "reference", pseudoInertial: true))
    let cache = CleanupFailingCache()
    let cached = try CachedEphemerisProvider(provider: model, cache: cache)
    #expect(try cached.state(at: .j2000).position.x == 0)
    model.base = 42
    #expect(throws: TestError.cacheFailure) { try cached.invalidate() }
    #expect(cache.storage.count == 1) // Old entries remain, but cannot be reused.
    #expect(try cached.state(at: .j2000).position.x == 42)
    #expect(try cached.state(at: .j2000).position.x == 42)
    #expect(model.calls == 2 && cache.storage.count == 2)
    cache.storage.removeAll()
  }
  private final class BadCache: EphemerisStateCache {
    var hit: CartesianState?
    var fails = false
    func state(for key: EphemerisCacheKey) throws -> CartesianState? {
      if fails { throw TestError.cacheFailure }; return hit
    }
    func store(_ state: CartesianState,for key: EphemerisCacheKey) {}
    func removeAll(in scope: EphemerisCacheScope) {}
  }
  @Test func failedEvaluationsAndInvalidCacheHitsAreNeverAccepted() throws {
    let model = try Model(frame: Frame(name: "reference",pseudoInertial: true))
    let cache = try MemoryEphemerisCache(capacity: 10)
    let cached = try CachedEphemerisProvider(provider: model,cache: cache)
    model.fail = true
    #expect(throws: TestError.transient) { try cached.state(at: .j2000) }
    #expect(cache.count == 0)
    model.fail = false
    _ = try cached.state(at: .j2000); _ = try cached.state(at: .j2000)
    #expect(model.calls == 2)
    let bad = BadCache(), untrusted = try CachedEphemerisProvider(provider: model,cache: bad)
    bad.hit = try model.state(at: epoch(1))
    #expect(throws: FrameTransformError.epochMismatch) { try untrusted.state(at: .j2000) }
    bad.hit = try CartesianState(epoch: .j2000,frame: Frame(name: "wrong",pseudoInertial: true),position: .zero,velocity: .zero)
    #expect(throws: EphemerisError.frameMismatch) { try untrusted.state(at: .j2000) }
    bad.fails = true
    #expect(throws: TestError.cacheFailure) { try untrusted.state(at: .j2000) }
  }
  @Test func multiBodyBatchesPreserveOrderAndFailFast() throws {
    let frame = Frame(name: "reference",pseudoInertial: true)
    let earth = try Model(frame: frame), moon = try Model(frame: frame,body: .moon)
    moon.base = 10
    let collection = try EphemerisCollection(providers: [earth,moon])
    let requests = try [EphemerisRequest(body: .moon,epoch: epoch(2)),EphemerisRequest(body: .earth,epoch: epoch(-1))]
    let states = try collection.states(for: requests)
    #expect(states.map { $0.position.x } == [12,-1])
    let before = earth.calls
    #expect(throws: EphemerisError.unsupportedBody) {
      try collection.states(for: [EphemerisRequest(body: .sun,epoch: .j2000),EphemerisRequest(body: .earth,epoch: .j2000)])
    }
    #expect(earth.calls == before)
  }
}

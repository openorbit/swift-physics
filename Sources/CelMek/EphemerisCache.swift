// SPDX-License-Identifier: Apache-2.0
import Foundation

public enum EphemerisCacheError: Error, Equatable { case invalidCapacity }

/// Opaque identity for one configured provider and invalidation generation.
public struct EphemerisCacheScope: Hashable, Sendable {
  private let identifier: UUID
  public init() { identifier = UUID() }
}

/// Exact atomic epoch, without rounding or interpolation. Scope prevents model,
/// origin, frame, accuracy and data-version collisions between provider wrappers.
public struct EphemerisCacheKey: Hashable, Sendable {
  public let scope: EphemerisCacheScope
  public let epoch: Epoch
  public init(scope: EphemerisCacheScope, epoch: Epoch) { self.scope = scope; self.epoch = epoch }
}

/// Injectable native-state cache. Custom implementations own storage/concurrency;
/// callbacks must not reenter the wrapper that invoked them. Errors propagate.
public protocol EphemerisStateCache: AnyObject {
  func state(for key: EphemerisCacheKey) throws -> CartesianState?
  func store(_ state: CartesianState, for key: EphemerisCacheKey) throws
  func removeAll(in scope: EphemerisCacheScope) throws
}

/// Bounded O(1) lookup and LRU updates, with a single lock around storage only.
/// Cache storage is synchronized; this does not make wrapped providers thread-safe.
public final class MemoryEphemerisCache: EphemerisStateCache {
  public let capacity: Int
  private final class Entry {
    let key: EphemerisCacheKey
    var state: CartesianState
    weak var previous: Entry?
    var next: Entry?
    init(key: EphemerisCacheKey, state: CartesianState) { self.key = key; self.state = state }
  }
  private let lock = NSLock()
  private var entries: [EphemerisCacheKey: Entry] = [:]
  private var mostRecent: Entry?
  private var leastRecent: Entry?
  public init(capacity: Int) throws {
    guard capacity > 0 else { throw EphemerisCacheError.invalidCapacity }
    self.capacity = capacity
  }
  public var count: Int {
    lock.lock(); defer { lock.unlock() }
    return entries.count
  }
  private func detach(_ entry: Entry) {
    if let previous = entry.previous { previous.next = entry.next } else { mostRecent = entry.next }
    if let next = entry.next { next.previous = entry.previous } else { leastRecent = entry.previous }
    entry.previous = nil; entry.next = nil
  }
  private func prepend(_ entry: Entry) {
    entry.next = mostRecent; mostRecent?.previous = entry
    mostRecent = entry
    if leastRecent == nil { leastRecent = entry }
  }
  public func state(for key: EphemerisCacheKey) -> CartesianState? {
    lock.lock(); defer { lock.unlock() }
    guard let entry = entries[key] else { return nil }
    detach(entry); prepend(entry)
    return entry.state
  }
  public func store(_ state: CartesianState, for key: EphemerisCacheKey) {
    lock.lock(); defer { lock.unlock() }
    if let entry = entries[key] {
      entry.state = state; detach(entry); prepend(entry)
    } else {
      let entry = Entry(key: key,state: state)
      entries[key] = entry; prepend(entry)
      if entries.count > capacity, let oldest = leastRecent {
        detach(oldest); entries.removeValue(forKey: oldest.key)
      }
    }
  }
  public func removeAll(in scope: EphemerisCacheScope) {
    lock.lock(); defer { lock.unlock() }
    let keys = entries.keys.filter { $0.scope == scope }
    for key in keys {
      if let entry = entries.removeValue(forKey: key) { detach(entry) }
    }
  }
  /// Clear all scopes, e.g. when releasing a game scene's shared storage.
  public func removeAll() {
    lock.lock(); defer { lock.unlock() }
    // Break next links iteratively rather than recursively releasing a long chain.
    while let entry = mostRecent { detach(entry) }
    entries.removeAll(); leastRecent = nil
  }
  deinit { removeAll() }
}

/// Exact native-state memoization for a fixed provider configuration. Construct a
/// new wrapper when metadata/frame/settings change; invalidate after data updates
/// that retain metadata. Error results are never cached. No fallback or interpolation.
public final class CachedEphemerisProvider: EphemerisProvider {
  public let metadata: EphemerisMetadata
  private let provider: any EphemerisProvider
  private let cache: any EphemerisStateCache
  private let lock = NSLock()
  private var scope = EphemerisCacheScope()
  public init(provider: any EphemerisProvider, cache: any EphemerisStateCache) throws {
    guard provider.frame === provider.metadata.frame else { throw EphemerisError.frameMismatch }
    self.provider = provider; self.cache = cache; metadata = provider.metadata
  }
  private func validate(_ state: CartesianState, at epoch: Epoch) throws {
    guard state.epoch == epoch else { throw FrameTransformError.epochMismatch }
    guard state.frame === metadata.frame else { throw EphemerisError.frameMismatch }
    guard (state.acceleration == nil) == (metadata.accelerationMethod == .unavailable) else {
      throw EphemerisError.accelerationMismatch
    }
  }
  public func state(at epoch: Epoch) throws -> CartesianState {
    guard metadata.coverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
    lock.lock(); let generation = scope; lock.unlock()
    let key = EphemerisCacheKey(scope: generation,epoch: epoch)
    if let hit = try cache.state(for: key) {
      try validate(hit,at: epoch)
      return hit
    }
    let result = try provider.checkedState(at: epoch)
    try validate(result,at: epoch)
    lock.lock(); defer { lock.unlock() }
    // A query in flight during invalidation may finish, but cannot repopulate the old scope.
    if scope == generation { try cache.store(result,for: key) }
    return result
  }
  /// Rotate generation before clearing storage. In-flight queries can complete with
  /// their original result; subsequent queries always use the new generation.
  public func invalidate() throws {
    lock.lock(); let old = scope; scope = EphemerisCacheScope(); lock.unlock()
    try cache.removeAll(in: old)
  }
}

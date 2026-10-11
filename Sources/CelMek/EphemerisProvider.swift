// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum EphemerisError: Error, Equatable {
  case invalidMetadata, invalidDerivativeMethod, outsideCoverage, frameMismatch
  case accelerationMismatch, duplicateBody, unsupportedBody
}

/// How a derivative is obtained. All returned derivatives still use atomic SI seconds.
public enum EphemerisDerivativeMethod: Hashable, Sendable {
  case analytic
  /// Numerical derivative; the interval and neighbouring coverage belong to the model.
  case finiteDifference(step: Duration)
  /// Model-supplied derivative, e.g. an integrated or tabulated state. Describe its origin.
  case supplied(description: String)
  case unavailable

  fileprivate func validate() throws {
    switch self {
    case .finiteDifference(let step):
      guard step > .zero else { throw EphemerisError.invalidDerivativeMethod }
    case .supplied(let description):
      guard !description.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
        throw EphemerisError.invalidDerivativeMethod
      }
    default: break
    }
  }
}

public enum EphemerisCoverageBasis: Sendable {
  /// A caller's usable interval, without a claim of certified physical accuracy.
  case callerSelected
  /// Coverage published by the model or data source, recorded in provenance.
  case modelDocumented
}

/// Query contract for one target relative to an explicit origin. The frame identity
/// must actually represent that origin; metadata does not translate or rotate states.
public struct EphemerisMetadata {
  public let body: BodyIdentifier
  public let origin: BodyIdentifier
  public let frame: Frame
  public let coverage: ClosedRange<Epoch>
  public let coverageBasis: EphemerisCoverageBasis
  /// Time argument consumed by the underlying model. State derivatives are converted
  /// to elapsed atomic seconds regardless of this coordinate-time representation.
  public let modelTimeScale: TimeScale
  /// All coordinate-time arguments used by a composed model.
  public let modelTimeScales: Set<TimeScale>
  public let velocityMethod: EphemerisDerivativeMethod
  public let accelerationMethod: EphemerisDerivativeMethod
  public let provenance: PhysicalDataProvenance
  public init(body: BodyIdentifier, origin: BodyIdentifier, frame: Frame,
              coverage: ClosedRange<Epoch>, coverageBasis: EphemerisCoverageBasis,
              modelTimeScale: TimeScale, velocityMethod: EphemerisDerivativeMethod,
              accelerationMethod: EphemerisDerivativeMethod, provenance: PhysicalDataProvenance,
              modelTimeScales: Set<TimeScale>? = nil) throws {
    try body.validate(); try origin.validate()
    try velocityMethod.validate(); try accelerationMethod.validate()
    // CartesianState always supplies velocity. An unavailable velocity cannot be represented.
    guard velocityMethod != .unavailable else { throw EphemerisError.invalidMetadata }
    let scales = modelTimeScales ?? [modelTimeScale]
    guard scales.contains(modelTimeScale) else { throw EphemerisError.invalidMetadata }
    self.modelTimeScales = scales
    self.body = body; self.origin = origin; self.frame = frame; self.coverage = coverage
    self.coverageBasis = coverageBasis; self.modelTimeScale = modelTimeScale
    self.velocityMethod = velocityMethod; self.accelerationMethod = accelerationMethod
    self.provenance = provenance
  }
}

/// Prescribed ephemeris for one body, independent of force-driven propagators.
/// Implementations must enforce metadata and return SI states at the requested epoch.
public protocol EphemerisProvider: StateProvider {
  var metadata: EphemerisMetadata { get }
}

extension EphemerisProvider {
  public var frame: Frame { metadata.frame }

  /// Validates custom implementations too. Coverage is checked before evaluation;
  /// this method never extrapolates, applies fallback, or manufactures derivatives.
  public func checkedState(at epoch: Epoch) throws -> CartesianState {
    guard metadata.coverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
    guard frame === metadata.frame else { throw EphemerisError.frameMismatch }
    let result = try state(at: epoch)
    guard result.epoch == epoch else { throw FrameTransformError.epochMismatch }
    guard result.frame === metadata.frame else { throw EphemerisError.frameMismatch }
    guard (result.acceleration == nil) == (metadata.accelerationMethod == .unavailable) else {
      throw EphemerisError.accelerationMismatch
    }
    return result
  }
}

/// Adapter for an existing StateProvider and explicit model metadata. Coverage is
/// caller-owned: restrict it to the intersection of all model/time-data dependencies.
/// A narrower provider failure propagates; metadata never extends physical coverage.
public struct StateProviderEphemeris: EphemerisProvider {
  public let metadata: EphemerisMetadata
  private let provider: any StateProvider
  public init(provider: any StateProvider, metadata: EphemerisMetadata) throws {
    guard provider.frame === metadata.frame else { throw EphemerisError.frameMismatch }
    self.provider = provider; self.metadata = metadata
  }
  public func state(at epoch: Epoch) throws -> CartesianState {
    guard metadata.coverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
    let state = try provider.state(at: epoch)
    guard state.epoch == epoch else { throw FrameTransformError.epochMismatch }
    guard state.frame === metadata.frame else { throw EphemerisError.frameMismatch }
    guard (state.acceleration == nil) == (metadata.accelerationMethod == .unavailable) else {
      throw EphemerisError.accelerationMismatch
    }
    return state
  }
}

/// Immutable selection of one ephemeris per body. Origins, frames and coverage may
/// differ per body; no implicit origin conversion or provider fallback takes place.
/// No cache or hidden dependency on query order. Provider thread safety remains its own contract.
public struct EphemerisCollection {
  public let bodies: [BodyIdentifier]
  private let providers: [BodyIdentifier: any EphemerisProvider]
  public init(providers: [any EphemerisProvider]) throws {
    var lookup: [BodyIdentifier: any EphemerisProvider] = [:]
    for provider in providers {
      guard provider.frame === provider.metadata.frame else { throw EphemerisError.frameMismatch }
      guard lookup.updateValue(provider,forKey: provider.metadata.body) == nil else { throw EphemerisError.duplicateBody }
    }
    self.providers = lookup; bodies = providers.map { $0.metadata.body }
  }
  public func provider(for body: BodyIdentifier) throws -> any EphemerisProvider {
    guard let provider = providers[body] else { throw EphemerisError.unsupportedBody }
    return provider
  }
  public func metadata(for body: BodyIdentifier) throws -> EphemerisMetadata {
    try provider(for: body).metadata
  }
  public func state(of body: BodyIdentifier, at epoch: Epoch) throws -> CartesianState {
    try provider(for: body).checkedState(at: epoch)
  }
}

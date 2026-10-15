//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2026 Mattias Holm
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
import Math

public enum BodyFixedFrameError: Error, Equatable {
  case frameMismatch, accelerationUnavailable, invalidDifferenceStep
}

/// Orientation from reference coordinates into body-fixed axes, without an
/// origin translation. Angular derivatives are expressed in body-fixed axes.
public struct BodyOrientation {
  public let referenceFrame: Frame
  public let epoch: Epoch
  public let rotation: Quaternion
  public let angularVelocity: SIMD3<Double>
  public let angularAcceleration: SIMD3<Double>

  public init(referenceFrame: Frame, epoch: Epoch, rotation: Quaternion,
              angularVelocity: SIMD3<Double> = .zero,
              angularAcceleration: SIMD3<Double> = .zero) throws {
    let validated = try FrameTransform(epoch: epoch, rotation: rotation,
      angularVelocity: angularVelocity, angularAcceleration: angularAcceleration)
    self.referenceFrame = referenceFrame
    self.epoch = epoch
    self.rotation = validated.rotation
    self.angularVelocity = angularVelocity
    self.angularAcceleration = angularAcceleration
  }
}

public protocol OrientationProvider {
  var referenceFrame: Frame { get }
  func orientation(at epoch: Epoch) throws -> BodyOrientation
}

/// Prescribed fixed-axis spin; no planetary orientation model is implied.
/// Axis, rate, and initial rotation follow UniformRotationProvider conventions.
public struct FixedAxisOrientationProvider: OrientationProvider {
  public let referenceFrame: Frame
  private let spin: UniformRotationProvider
  public init(referenceFrame: Frame, referenceEpoch: Epoch = .j2000,
              axis: SIMD3<Double>, rate: Double, initialRotation: Quaternion = .identity) throws {
    self.referenceFrame = referenceFrame
    spin = try UniformRotationProvider(referenceEpoch: referenceEpoch, axis: axis,
      rate: rate, initialRotation: initialRotation)
  }
  public func orientation(at epoch: Epoch) throws -> BodyOrientation {
    let transform = try spin.transform(at: epoch)
    return try BodyOrientation(referenceFrame: referenceFrame, epoch: epoch,
      rotation: transform.rotation, angularVelocity: transform.angularVelocity,
      angularAcceleration: transform.angularAcceleration)
  }
}

/// Handling of missing acceleration at the body's centre. A supplied value
/// always takes precedence; approximation is never selected implicitly.
public enum OriginAccelerationPolicy {
  case requireProvided
  /// Explicitly treat origin acceleration as zero, suitable for prescribed
  /// constant-velocity motion but approximate for orbiting bodies.
  case assumeZero
  /// Central difference of provider velocity over elapsed atomic seconds.
  /// Both neighbouring samples must be available in the same frame.
  case estimate(step: CelMek.Duration)
}

/// Parent-to-body-fixed mapping: p_body = R * (p_parent - centre).
/// The centre and orientation must use the exact same reference frame identity.
public struct BodyFixedTransformProvider: TransformProvider {
  public let referenceFrame: Frame
  private let centre: any StateProvider
  private let orientationProvider: any OrientationProvider
  private let accelerationPolicy: OriginAccelerationPolicy

  public init(centre: any StateProvider, orientation: any OrientationProvider,
              accelerationPolicy: OriginAccelerationPolicy = .requireProvided) throws {
    guard centre.frame === orientation.referenceFrame else { throw BodyFixedFrameError.frameMismatch }
    if case .estimate(let step) = accelerationPolicy {
      guard step > .zero, step.seconds.isFinite else { throw BodyFixedFrameError.invalidDifferenceStep }
    }
    referenceFrame = centre.frame
    self.centre = centre
    orientationProvider = orientation
    self.accelerationPolicy = accelerationPolicy
  }

  private func centreState(at epoch: Epoch) throws -> CartesianState {
    let state = try centre.state(at: epoch)
    guard state.epoch == epoch else { throw FrameTransformError.epochMismatch }
    guard state.frame === referenceFrame else { throw BodyFixedFrameError.frameMismatch }
    return state
  }

  public func transform(at epoch: Epoch) throws -> FrameTransform {
    let state = try centreState(at: epoch)
    let orientation = try orientationProvider.orientation(at: epoch)
    guard orientation.epoch == epoch else { throw FrameTransformError.epochMismatch }
    guard orientation.referenceFrame === referenceFrame else { throw BodyFixedFrameError.frameMismatch }
    let acceleration: SIMD3<Double>
    if let supplied = state.acceleration { acceleration = supplied }
    else {
      switch accelerationPolicy {
      case .requireProvided: throw BodyFixedFrameError.accelerationUnavailable
      case .assumeZero: acceleration = .zero
      case .estimate(let step):
        let before = try epoch.advanced(by: step.negated())
        let after = try epoch.advanced(by: step)
        let beforeState = try centreState(at: before)
        let afterState = try centreState(at: after)
        acceleration = (afterState.velocity - beforeState.velocity) / (2 * step.seconds)
      }
    }
    let centring = try FrameTransform(epoch: epoch, translation: -state.position,
      velocity: -state.velocity, acceleration: -acceleration)
    let axes = try FrameTransform(epoch: epoch, rotation: orientation.rotation,
      angularVelocity: orientation.angularVelocity, angularAcceleration: orientation.angularAcceleration)
    return try centring.followed(by: axes)
  }
}

extension FrameTree {
  /// Adds a body-fixed child to the providers' common reference frame.
  public func makeBodyFixedFrame(name: String, centre: any StateProvider,
                                orientation: any OrientationProvider,
                                accelerationPolicy: OriginAccelerationPolicy = .requireProvided) throws -> Frame {
    let provider = try BodyFixedTransformProvider(centre: centre, orientation: orientation,
      accelerationPolicy: accelerationPolicy)
    return try makeFrame(parent: provider.referenceFrame, name: name, provider: provider)
  }
}

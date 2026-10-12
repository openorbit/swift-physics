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
import Foundation
import Math

public enum FrameTransformError: Error, Equatable {
  case invalidRotation, nonFinite, epochMismatch, disconnectedFrames
}

/// Rigid mapping p_destination = R * p_source + translation.
/// Derivatives use elapsed atomic seconds. Angular derivatives are expressed
/// in destination axes and satisfy d(R*p)/dt = angularVelocity × (R*p).
public struct FrameTransform {
  public let epoch: Epoch
  public let rotation: Quaternion
  public let translation: SIMD3<Double>
  public let velocity: SIMD3<Double>
  public let acceleration: SIMD3<Double>
  public let angularVelocity: SIMD3<Double>
  public let angularAcceleration: SIMD3<Double>

  public init(epoch: Epoch, rotation: Quaternion = .identity,
              translation: SIMD3<Double> = .zero, velocity: SIMD3<Double> = .zero,
              acceleration: SIMD3<Double> = .zero, angularVelocity: SIMD3<Double> = .zero,
              angularAcceleration: SIMD3<Double> = .zero) throws {
    let scale = max(abs(rotation.x), abs(rotation.y), abs(rotation.z), abs(rotation.w))
    guard scale.isFinite, scale > 0,
          [rotation.x, rotation.y, rotation.z, rotation.w].allSatisfy({ $0.isFinite }) else {
      throw FrameTransformError.invalidRotation
    }
    let q = rotation / scale
    let norm = sqrt(Quaternion.dot(lhs: q, rhs: q))
    guard [translation, velocity, acceleration, angularVelocity, angularAcceleration]
      .allSatisfy({ $0.x.isFinite && $0.y.isFinite && $0.z.isFinite }) else {
      throw FrameTransformError.nonFinite
    }
    self.epoch = epoch
    self.rotation = q / norm
    self.translation = translation
    self.velocity = velocity
    self.acceleration = acceleration
    self.angularVelocity = angularVelocity
    self.angularAcceleration = angularAcceleration
  }

  public static func identity(at epoch: Epoch) -> FrameTransform {
    try! FrameTransform(epoch: epoch)
  }

  public func applying(to state: CartesianState, destination: Frame) throws -> CartesianState {
    guard state.epoch == epoch else { throw FrameTransformError.epochMismatch }
    let r = rotation.matrix * state.position
    let v = rotation.matrix * state.velocity
    let a = state.acceleration.map {
      rotation.matrix * $0 + 2 * cross(angularVelocity, v)
        + cross(angularAcceleration, r) + cross(angularVelocity, cross(angularVelocity, r)) + acceleration
    }
    return try CartesianState(epoch: epoch, frame: destination, position: r + translation,
      velocity: v + cross(angularVelocity, r) + velocity, acceleration: a)
  }

  /// Applies this mapping followed by next; both must describe the same epoch.
  public func followed(by next: FrameTransform) throws -> FrameTransform {
    guard epoch == next.epoch else { throw FrameTransformError.epochMismatch }
    let r = next.rotation.matrix * translation
    let v = next.rotation.matrix * velocity
    let w = next.rotation.matrix * angularVelocity
    return try FrameTransform(epoch: epoch, rotation: next.rotation * rotation,
      translation: r + next.translation,
      velocity: v + cross(next.angularVelocity, r) + next.velocity,
      acceleration: next.rotation.matrix * acceleration + 2 * cross(next.angularVelocity, v)
        + cross(next.angularAcceleration, r) + cross(next.angularVelocity, cross(next.angularVelocity, r))
        + next.acceleration,
      angularVelocity: w + next.angularVelocity,
      angularAcceleration: next.rotation.matrix * angularAcceleration + next.angularAcceleration
        + cross(next.angularVelocity, w))
  }

  public func inverted() throws -> FrameTransform {
    let inverse = rotation.conj.matrix
    return try FrameTransform(epoch: epoch, rotation: rotation.conj,
      translation: inverse * -translation,
      velocity: inverse * (-velocity + cross(angularVelocity, translation)),
      acceleration: inverse * (-acceleration + 2 * cross(angularVelocity, velocity)
        + cross(angularAcceleration, translation) - cross(angularVelocity, cross(angularVelocity, translation))),
      angularVelocity: inverse * -angularVelocity, angularAcceleration: inverse * -angularAcceleration)
  }
}

private func cross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
  SIMD3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x)
}

public protocol TransformProvider {
  /// Returns the parent-to-child mapping at the requested epoch.
  func transform(at epoch: Epoch) throws -> FrameTransform
}

/// A constant rigid mapping with zero derivatives.
public struct FixedTransformProvider: TransformProvider {
  private let fixed: FrameTransform
  public init(rotation: Quaternion = .identity, translation: SIMD3<Double> = .zero) throws {
    fixed = try FrameTransform(epoch: .j2000, rotation: rotation, translation: translation)
  }
  public func transform(at epoch: Epoch) throws -> FrameTransform {
    try FrameTransform(epoch: epoch, rotation: fixed.rotation, translation: fixed.translation)
  }
}

/// Rotation of coordinate axes about a fixed axis, with a stationary common origin.
/// Positive rate turns the child axes by the right-hand rule in parent coordinates;
/// the parent-to-child coordinate mapping therefore has negative angular velocity.
public struct UniformRotationProvider: TransformProvider {
  public let referenceEpoch: Epoch
  public let rate: Double
  private let axis: SIMD3<Double>
  private let initial: Quaternion
  public init(referenceEpoch: Epoch = .j2000, axis: SIMD3<Double>, rate: Double,
              initialRotation: Quaternion = .identity) throws {
    let scale = max(abs(axis.x), abs(axis.y), abs(axis.z))
    guard scale.isFinite, scale > 0, axis.x.isFinite, axis.y.isFinite, axis.z.isFinite, rate.isFinite else {
      throw FrameTransformError.nonFinite
    }
    let scaled = axis / scale
    self.axis = scaled / sqrt((scaled * scaled).sum())
    self.referenceEpoch = referenceEpoch
    self.rate = rate
    self.initial = try FrameTransform(epoch: referenceEpoch, rotation: initialRotation).rotation
  }
  public func transform(at epoch: Epoch) throws -> FrameTransform {
    let angle = -rate * (try epoch.duration(since: referenceEpoch)).seconds
    guard angle.isFinite else { throw FrameTransformError.nonFinite }
    return try FrameTransform(epoch: epoch,
      rotation: initial * Quaternion(angle: angle.truncatingRemainder(dividingBy: 2 * .pi), axis: axis),
      angularVelocity: initial.matrix * (-rate * axis))
  }
}

/// Owns the providers for a frame tree. Parent links are immutable and strong.
/// Providers may throw for coverage or data errors; queries do not cache results.
/// Custom providers must not strongly capture their owning tree.
public final class FrameTree {
  public let root: Frame
  private var edges: [ObjectIdentifier: any TransformProvider] = [:]
  public init(name: String = "Inertial", pseudoInertial: Bool = true) {
    root = Frame(name: name, pseudoInertial: pseudoInertial)
  }
  /// Adopts an existing frame as a local tree root. Its external ancestors
  /// are outside this tree; queries remain within the registered subtree.
  public init(root: Frame) throws {
    self.root = root
  }
  public func makeFrame(parent: Frame, name: String, pseudoInertial: Bool = false,
                        provider: any TransformProvider) throws -> Frame {
    guard contains(parent) else { throw FrameTransformError.disconnectedFrames }
    let child = Frame(parent: parent, name: name, pseudoInertial: pseudoInertial)
    edges[ObjectIdentifier(child)] = provider
    // Retain identities to prevent ObjectIdentifier reuse after callers release a child.
    frames.append(child)
    return child
  }
  /// Registers an existing identity under its immutable parent. The provider
  /// must map that parent into this frame. No reparenting or relabelling occurs.
  public func registerFrame(_ frame: Frame, provider: any TransformProvider) throws {
    guard frame !== root, !contains(frame), let parent = frame.parent, contains(parent) else {
      throw FrameTransformError.disconnectedFrames
    }
    edges[ObjectIdentifier(frame)] = provider
    frames.append(frame)
  }
  private var frames: [Frame] = []
  internal func isRegistered(_ frame: Frame) -> Bool { contains(frame) }
  private func contains(_ frame: Frame) -> Bool {
    frame === root || edges[ObjectIdentifier(frame)] != nil
  }
  public func transform(from source: Frame, to destination: Frame, at epoch: Epoch) throws -> FrameTransform {
    guard contains(source), contains(destination) else { throw FrameTransformError.disconnectedFrames }
    // Find the common ancestor so unrelated edge providers need not be evaluated.
    var ancestors: [Frame] = []
    var node: Frame? = source
    while let current = node { ancestors.append(current); node = current.parent }
    node = destination
    var down: [Frame] = []
    while let current = node, !ancestors.contains(where: { $0 === current }) {
      down.append(current); node = current.parent
    }
    guard let common = node else { throw FrameTransformError.disconnectedFrames }
    var result = FrameTransform.identity(at: epoch)
    var up = source
    while up !== common {
      let edge = try edges[ObjectIdentifier(up)]!.transform(at: epoch)
      result = try result.followed(by: edge.inverted())
      up = up.parent!
    }
    for child in down.reversed() {
      result = try result.followed(by: edges[ObjectIdentifier(child)]!.transform(at: epoch))
    }
    return result
  }
  public func convert(_ state: CartesianState, to destination: Frame) throws -> CartesianState {
    try transform(from: state.frame, to: destination, at: state.epoch)
      .applying(to: state, destination: destination)
  }
}

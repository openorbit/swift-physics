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

/// Input conventions for normalizing Cartesian coordinates to SI units.
public enum CartesianUnits: Sendable {
  /// Position m, velocity m/s, acceleration m/s².
  case metresAndSeconds
  /// Position km, velocity km/s, acceleration km/s².
  case kilometresAndSeconds
  /// Position AU, velocity AU/day, acceleration AU/day².
  case astronomicalUnitsAndDays

  public static let metresPerAstronomicalUnit = 149_597_870_700.0
  public static let secondsPerDay = 86_400.0

  fileprivate var scales: (position: Double, velocity: Double, acceleration: Double) {
    switch self {
    case .metresAndSeconds: return (1, 1, 1)
    case .kilometresAndSeconds: return (1000, 1000, 1000)
    case .astronomicalUnitsAndDays:
      let au = Self.metresPerAstronomicalUnit
      let day = Self.secondsPerDay
      return (au, au / day, au / (day * day))
    }
  }
}

public enum CartesianStateError: Error, Equatable, Sendable {
  case nonFinitePosition
  case nonFiniteVelocity
  case nonFiniteAcceleration
}

/// Position and its time derivatives at an epoch, expressed in one reference frame.
/// Stored units are always metres, metres/second, and metres/second².
/// A nil acceleration means unavailable, not zero. Frame identity is object identity.
public struct CartesianState: Equatable {
  public let epoch: Epoch
  public let frame: Frame
  public let position: SIMD3<Double>
  public let velocity: SIMD3<Double>
  public let acceleration: SIMD3<Double>?

  /// Normalizes the supplied Cartesian input units and rejects nonfinite values,
  /// including overflow during conversion. Spherical coordinates must first be
  /// converted to Cartesian coordinates in the intended frame.
  public init(epoch: Epoch, frame: Frame, position: SIMD3<Double>, velocity: SIMD3<Double>,
              acceleration: SIMD3<Double>? = nil, units: CartesianUnits = .metresAndSeconds) throws {
    let scales = units.scales
    let position = position * scales.position
    let velocity = velocity * scales.velocity
    let acceleration = acceleration.map { $0 * scales.acceleration }
    guard position.isFiniteVector else { throw CartesianStateError.nonFinitePosition }
    guard velocity.isFiniteVector else { throw CartesianStateError.nonFiniteVelocity }
    guard acceleration?.isFiniteVector != false else { throw CartesianStateError.nonFiniteAcceleration }
    self.epoch = epoch
    self.frame = frame
    self.position = position
    self.velocity = velocity
    self.acceleration = acceleration
  }

  public static func == (lhs: CartesianState, rhs: CartesianState) -> Bool {
    lhs.epoch == rhs.epoch && lhs.frame === rhs.frame && lhs.position == rhs.position
      && lhs.velocity == rhs.velocity && lhs.acceleration == rhs.acceleration
  }
}

/// A common query surface for prescribed ephemerides and future propagators.
/// A result must carry the requested epoch and this provider's native frame.
public protocol StateProvider {
  var frame: Frame { get }
  func state(at epoch: Epoch) throws -> CartesianState
}

/// Named astronomical frame identities that can be adopted as FrameTree roots.

private extension SIMD3 where Scalar == Double {
  var isFiniteVector: Bool { x.isFinite && y.isFinite && z.isFinite }
}

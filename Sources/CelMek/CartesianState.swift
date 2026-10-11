// SPDX-License-Identifier: Apache-2.0
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
public enum CelestialFrames {
  /// Solar-system barycentre origin; dynamical ecliptic and equinox J2000 axes.
  /// This is the native VSOP87E frame, not an ICRF/equatorial frame alias.
  public static let vsop87E = Frame(name: "VSOP87E barycentric ecliptic J2000", pseudoInertial: true)
}

private extension SIMD3 where Scalar == Double {
  var isFiniteVector: Bool { x.isFinite && y.isFinite && z.isFinite }
}

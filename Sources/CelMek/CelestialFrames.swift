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

/// Fixed astronomical frame identities. Origins are barycentric unless named
/// otherwise. Use makeTree() to register their physical transform connections.
public enum CelestialFrames {
  public static let icrf = Frame(name: "ICRF barycentric", pseudoInertial: true)
  /// Mean dynamical equator and equinox J2000, with IAU 2000 frame bias.
  public static let meanEquatorJ2000 = Frame(parent: icrf, name: "Mean equator J2000 barycentric", pseudoInertial: true)
  /// Mean ecliptic/equinox J2000 using the IAU 1976 obliquity 84381.448 arcsec.
  public static let eclipticJ2000 = Frame(parent: meanEquatorJ2000, name: "Ecliptic J2000 barycentric", pseudoInertial: true)
  /// FK5 axes frozen at J2000. Separate from the mean dynamical EME2000 frame.
  public static let fk5J2000 = Frame(parent: icrf, name: "FK5 J2000 barycentric", pseudoInertial: true)
  /// Native VSOP87E dynamical ecliptic/equinox J2000; not an alias for the
  /// conventional ecliptic frame. Connected through the published FK5 rotation.
  public static let vsop87E = Frame(parent: fk5J2000, name: "VSOP87E barycentric ecliptic J2000", pseudoInertial: true)

  private static let arcsecond = Double.pi / (180 * 3600)
  private static let obliquity = 84381.448 * arcsecond

  // Independent implementation of the IAU 2000 fixed frame-bias rotations.
  // SOFA/ERFA BP00 defines GCRS -> mean J2000 as Rx(-deps) Ry(dpsi*sin(eps)) Rz(dra)
  // with passive rotations; our Hamilton coordinate mappings have opposite signs.
  internal static var icrfToEME2000: Quaternion {
    let z = Quaternion(angle: 0.0146 * arcsecond, axis: SIMD3(0, 0, 1))
    let y = Quaternion(angle: 0.041775 * arcsecond * sin(obliquity), axis: SIMD3(0, 1, 0))
    let x = Quaternion(angle: -0.0068192 * arcsecond, axis: SIMD3(1, 0, 0))
    return x * y * z
  }

  // FK5 orientation at J2000 relative to the Hipparcos/ICRS realization.
  // This is a fixed-epoch tie, not a time-dependent FK5 stellar catalogue model.
  internal static var fk5ToICRF: Quaternion {
    let vector = SIMD3<Double>(-19.9, -9.1, 22.9) * (0.001 * arcsecond)
    return Quaternion(angle: -sqrt((vector * vector).sum()), axis: vector)
  }

  internal static var vsopToFK5: Quaternion {
    // Reconstruct an orthogonal rotation from the rounded published matrix;
    // avoids introducing a small scale/shear when converting states.
    let tilt = atan2(0.397776982902, 0.917482137087)
    let node = atan2(-0.000000479966, 1.0)
    return Quaternion(angle: node, axis: SIMD3(0, 0, 1))
      * Quaternion(angle: tilt, axis: SIMD3(1, 0, 0))
  }

  public static func makeTree() throws -> FrameTree {
    let tree = try FrameTree(root: icrf)
    try tree.registerFrame(meanEquatorJ2000, provider: FixedTransformProvider(rotation: icrfToEME2000))
    try tree.registerFrame(eclipticJ2000,
      provider: FixedTransformProvider(rotation: Quaternion(angle: -obliquity, axis: SIMD3(1, 0, 0))))
    try tree.registerFrame(fk5J2000, provider: FixedTransformProvider(rotation: fk5ToICRF.conj))
    try tree.registerFrame(vsop87E, provider: FixedTransformProvider(rotation: vsopToFK5.conj))
    return tree
  }
}

public enum FrameQueryError: Error, Equatable {
  case providerFrameMismatch, frameTreeUnavailable
}

extension StateProvider {
  /// Query the native state, then convert through the explicitly supplied tree.
  public func state(at epoch: Epoch, in destination: Frame, using tree: FrameTree) throws -> CartesianState {
    let native = try state(at: epoch)
    guard native.epoch == epoch else { throw FrameTransformError.epochMismatch }
    guard native.frame === frame else { throw FrameQueryError.providerFrameMismatch }
    return try tree.convert(native, to: destination)
  }
}

/// A provider viewed in another connected frame. The tree is weakly referenced
/// to avoid tree -> body-frame provider -> converted origin provider -> tree cycles.
/// Callers must retain the tree while using this adapter.
public final class FrameConvertedStateProvider: StateProvider {
  public let frame: Frame
  private let native: any StateProvider
  private weak var tree: FrameTree?
  public init(_ native: any StateProvider, in frame: Frame, using tree: FrameTree) throws {
    // Check connectivity without querying epoch-dependent provider data.
    guard tree.isRegistered(native.frame), tree.isRegistered(frame) else {
      throw FrameTransformError.disconnectedFrames
    }
    self.native = native
    self.frame = frame
    self.tree = tree
  }
  public func state(at epoch: Epoch) throws -> CartesianState {
    guard let tree else { throw FrameQueryError.frameTreeUnavailable }
    return try native.state(at: epoch, in: frame, using: tree)
  }
}

extension FrameTree {
  /// Newtonian geocentric frame with axes aligned to ICRF. This supplies origin
  /// translation only, not the relativistic BCRS/GCRS coordinate transformation.
  /// The caller supplies a geocentre ephemeris, not an Earth-Moon barycentre.
  public func makeGCRF(earth: any StateProvider,
                       accelerationPolicy: OriginAccelerationPolicy = .requireProvided) throws -> Frame {
    let origin = try FrameConvertedStateProvider(earth, in: CelestialFrames.icrf, using: self)
    let axes = try FixedAxisOrientationProvider(referenceFrame: CelestialFrames.icrf,
      axis: SIMD3(0, 0, 1), rate: 0)
    let provider = try BodyFixedTransformProvider(centre: origin, orientation: axes,
      accelerationPolicy: accelerationPolicy)
    return try makeFrame(parent: CelestialFrames.icrf, name: "GCRF", pseudoInertial: true, provider: provider)
  }
}

/// Geocentric ICRF-aligned and mean-J2000 frames sharing one moving origin.
public struct GeocentricFrames {
  public let gcrf: Frame
  public let eme2000: Frame
}

extension FrameTree {
  public func makeGeocentricFrames(earth: any StateProvider,
                                  accelerationPolicy: OriginAccelerationPolicy = .requireProvided) throws -> GeocentricFrames {
    let gcrf = try makeGCRF(earth: earth, accelerationPolicy: accelerationPolicy)
    let eme2000 = try makeFrame(parent: gcrf, name: "EME2000", pseudoInertial: true,
      provider: FixedTransformProvider(rotation: CelestialFrames.icrfToEME2000))
    return GeocentricFrames(gcrf: gcrf, eme2000: eme2000)
  }
}

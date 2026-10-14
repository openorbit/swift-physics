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

public enum BodyModelError: Error, Equatable {
  case bodyMismatch, invalidDirection, invalidPosition, singularGravity, nonFiniteResult
}

/// Orientation with an explicit body identity, independent of geometric and gravity models.
public protocol BodyOrientationProvider: OrientationProvider {
  var bodyIdentifier: BodyIdentifier { get }
  var provenance: PhysicalDataProvenance { get }
}

/// Labels an existing attitude provider without changing its frames or rotations.
public struct IdentifiedBodyOrientationProvider: BodyOrientationProvider {
  public let bodyIdentifier: BodyIdentifier
  public let provenance: PhysicalDataProvenance
  public var referenceFrame: Frame { provider.referenceFrame }
  private let provider: any OrientationProvider
  public init(bodyIdentifier: BodyIdentifier, provider: any OrientationProvider,
              provenance: PhysicalDataProvenance) throws {
    try bodyIdentifier.validate()
    self.bodyIdentifier = bodyIdentifier; self.provider = provider; self.provenance = provenance
  }
  public func orientation(at epoch: Epoch) throws -> BodyOrientation {
    let result = try provider.orientation(at: epoch)
    guard result.epoch == epoch else { throw FrameTransformError.epochMismatch }
    guard result.referenceFrame === referenceFrame else { throw BodyFixedFrameError.frameMismatch }
    return result
  }
}

/// Surface point and outward unit normal in the shape's stated body-fixed axes.
public struct BodySurfacePoint: Sendable {
  public let positionMetres: SIMD3<Double>
  public let outwardNormal: SIMD3<Double>
  public init(positionMetres: SIMD3<Double>, outwardNormal: SIMD3<Double>) throws {
    guard [positionMetres.x,positionMetres.y,positionMetres.z].allSatisfy(\.isFinite) else {
      throw BodyModelError.nonFiniteResult
    }
    self.positionMetres = positionMetres
    self.outwardNormal = try bodyModelUnit(outwardNormal)
  }
}

/// Geometry relative to the body's centre. This initial query describes a surface
/// reachable along a radial direction; it does not promise arbitrary mesh intersections.
public protocol BodyShapeProvider {
  var bodyIdentifier: BodyIdentifier { get }
  var axesDescription: String { get }
  var provenance: PhysicalDataProvenance { get }
  func surface(inDirection direction: SIMD3<Double>) throws -> BodySurfacePoint
}

/// Positive reference ellipsoid, including a sphere when all three radii match.
/// Radii and axes come from supplied shape data, independently of body metadata.
public struct EllipsoidBodyShape: BodyShapeProvider {
  public let bodyIdentifier: BodyIdentifier
  public let referenceRadii: BodyReferenceRadii
  public var axesDescription: String { referenceRadii.axesDescription }
  public var provenance: PhysicalDataProvenance { referenceRadii.provenance }
  public init(bodyIdentifier: BodyIdentifier, referenceRadii: BodyReferenceRadii) throws {
    try bodyIdentifier.validate()
    self.bodyIdentifier = bodyIdentifier; self.referenceRadii = referenceRadii
  }
  public func surface(inDirection direction: SIMD3<Double>) throws -> BodySurfacePoint {
    let u = try bodyModelUnit(direction), radii = referenceRadii.metres
    // Scale by the smallest radius, so intermediate ratios cannot overflow.
    let scale = min(radii.x,radii.y,radii.z), inverse = SIMD3(scale/radii.x,scale/radii.y,scale/radii.z)
    let weighted = u*inverse
    let denominator = sqrt((weighted*weighted).sum())
    guard denominator > 0 else { throw BodyModelError.nonFiniteResult }
    let point = u*(scale/denominator)
    // Gradient of sum((x_i/a_i)^2): direction is proportional to u_i/a_i^2.
    return try BodySurfacePoint(positionMetres: point,outwardNormal: u*inverse*inverse)
  }
}

/// Local gravity field. Units use the model's explicit coordinate-time convention;
/// conversion to atomic-time dynamics is a responsibility of the consuming dynamics model.
public struct BodyGravitySample: Sendable {
  public let epoch: Epoch
  public let acceleration: SIMD3<Double>
  /// Gravitational potential in m²/s² with zero at infinity for point-mass gravity.
  public let potential: Double
  public let convention: GravitationalParameterConvention
  public init(epoch: Epoch, acceleration: SIMD3<Double>, potential: Double,
              convention: GravitationalParameterConvention) throws {
    guard [acceleration.x,acceleration.y,acceleration.z,potential].allSatisfy(\.isFinite) else {
      throw BodyModelError.nonFiniteResult
    }
    self.epoch = epoch; self.acceleration = acceleration; self.potential = potential; self.convention = convention
  }
}

/// Field evaluated at a body-centred displacement in the model's described axes.
/// Providers neither fetch ephemerides nor transform world coordinates implicitly.
public protocol BodyGravityProvider {
  var bodyIdentifier: BodyIdentifier { get }
  var axesDescription: String { get }
  var provenance: PhysicalDataProvenance { get }
  var convention: GravitationalParameterConvention { get }
  func gravity(at epoch: Epoch, displacementMetres: SIMD3<Double>) throws -> BodyGravitySample
}

/// Prescribed Newtonian point-mass field: acceleration = −μ r/|r|³, potential = −μ/|r|.
/// No surface, collision, nonspherical field or relativistic correction is implied.
public struct PointMassBodyGravity: BodyGravityProvider {
  public let bodyIdentifier: BodyIdentifier
  public let parameter: GravitationalParameter
  public let axesDescription = "Any orthonormal body-centred Cartesian axes (spherically symmetric field)"
  public var provenance: PhysicalDataProvenance { parameter.provenance }
  public var convention: GravitationalParameterConvention { parameter.convention }
  public init(bodyIdentifier: BodyIdentifier, parameter: GravitationalParameter) throws {
    try bodyIdentifier.validate()
    self.bodyIdentifier = bodyIdentifier; self.parameter = parameter
  }
  public func gravity(at epoch: Epoch, displacementMetres r: SIMD3<Double>) throws -> BodyGravitySample {
    guard [r.x,r.y,r.z].allSatisfy(\.isFinite) else { throw BodyModelError.invalidPosition }
    let scale = max(abs(r.x),abs(r.y),abs(r.z))
    guard scale > 0 else { throw BodyModelError.singularGravity }
    let scaled = r/scale, length = sqrt((scaled*scaled).sum())
    let inverseRadius = (1/scale)/length
    let potential = -parameter.cubicMetresPerSecondSquared*inverseRadius
    let magnitude = potential*inverseRadius
    return try BodyGravitySample(epoch: epoch,acceleration: (scaled/length)*magnitude,
      potential: potential,convention: convention)
  }
}

/// Optional independent models selected for one body, including magnetic fields. Missing models stay absent;
/// GM/radii in metadata never silently create gravity or shape providers.
/// Identity validation does not certify consistency of model versions or axis descriptions.
public struct CelestialBodyModels {
  public let body: CelestialBody
  public let ephemeris: (any EphemerisProvider)?
  public let orientation: (any BodyOrientationProvider)?
  public let shape: (any BodyShapeProvider)?
  public let gravity: (any BodyGravityProvider)?
  public let magneticField: (any BodyMagneticProvider)?
  public init(body: CelestialBody, ephemeris: (any EphemerisProvider)? = nil,
              orientation: (any BodyOrientationProvider)? = nil,
              shape: (any BodyShapeProvider)? = nil, gravity: (any BodyGravityProvider)? = nil,
              magneticField: (any BodyMagneticProvider)? = nil) throws {
    let ids = [ephemeris?.metadata.body,orientation?.bodyIdentifier,shape?.bodyIdentifier,gravity?.bodyIdentifier,magneticField?.bodyIdentifier]
    guard ids.compactMap({ $0 }).allSatisfy({ $0 == body.identifier }) else { throw BodyModelError.bodyMismatch }
    self.body = body; self.ephemeris = ephemeris; self.orientation = orientation; self.shape = shape; self.gravity = gravity
    self.magneticField = magneticField
  }
}

private func bodyModelUnit(_ vector: SIMD3<Double>) throws -> SIMD3<Double> {
  let scale = max(abs(vector.x),abs(vector.y),abs(vector.z))
  guard [vector.x,vector.y,vector.z].allSatisfy(\.isFinite), scale > 0 else { throw BodyModelError.invalidDirection }
  let scaled = vector/scale
  return scaled/sqrt((scaled*scaled).sum())
}

// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum LocalOrbitalFrameError: Error, Equatable {
  case nonInertialReference, frameMismatch, zeroPosition, zeroVelocity, undefinedOrbitalPlane
  case invalidConfiguration
}

/// Right-handed axes expressed in the orbit's inertial reference frame.
public enum LocalOrbitalAxes: Sendable {
  /// X = radial outward, Z = r×v, Y = Z×X (also RSW/QSW).
  case rtn
  /// X = velocity, Z = r×v, Y = Z×X.
  case tnw
  /// Z = nadir (−radial), Y = −orbital normal, X = Y×Z.
  case lvlhCCSDS
}

/// Orbital attitude from motion relative to a supplied central body. Both providers
/// must use the same pseudo-inertial frame. Without a central provider, the reference
/// frame's origin defines the orbit centre. No force model or propagation is implied.
/// Five-point angular derivatives use atomic seconds and require ±2*derivativeStep
/// state coverage. The step must resolve motion; no accuracy guarantee is implied.
public struct LocalOrbitalOrientationProvider: OrientationProvider {
  public let referenceFrame: Frame
  public let axes: LocalOrbitalAxes
  public let derivativeStep: Duration
  public let minimumPlaneSine: Double
  private let spacecraft: any StateProvider
  private let centralBody: (any StateProvider)?
  private let mounting: FrameTransform

  /// mountingRotation maps orbital coordinates into spacecraft coordinates.
  /// minimumPlaneSine rejects nearly parallel r/v using |r̂×v̂|, independently of units.
  public init(spacecraft: any StateProvider, centralBody: (any StateProvider)? = nil,
              axes: LocalOrbitalAxes, derivativeStep: Duration,
              minimumPlaneSine: Double = 1e-12, mountingRotation: Quaternion = .identity) throws {
    guard spacecraft.frame.pseudoInertial else { throw LocalOrbitalFrameError.nonInertialReference }
    guard centralBody == nil || centralBody!.frame === spacecraft.frame else {
      throw LocalOrbitalFrameError.frameMismatch
    }
    guard derivativeStep > .zero, minimumPlaneSine.isFinite,
          minimumPlaneSine > 0, minimumPlaneSine < 1 else { throw LocalOrbitalFrameError.invalidConfiguration }
    self.spacecraft = spacecraft; self.centralBody = centralBody
    referenceFrame = spacecraft.frame; self.axes = axes; self.derivativeStep = derivativeStep
    self.minimumPlaneSine = minimumPlaneSine
    mounting = try FrameTransform(epoch: .j2000, rotation: mountingRotation)
  }
  private func state(_ provider: any StateProvider, at epoch: Epoch) throws -> CartesianState {
    let state = try provider.state(at: epoch)
    guard state.epoch == epoch else { throw FrameTransformError.epochMismatch }
    guard state.frame === referenceFrame else { throw LocalOrbitalFrameError.frameMismatch }
    return state
  }
  private func unit(_ vector: SIMD3<Double>, zero: LocalOrbitalFrameError) throws -> SIMD3<Double> {
    let scale = max(abs(vector.x),abs(vector.y),abs(vector.z))
    guard scale.isFinite else { throw FrameTransformError.nonFinite }
    guard scale > 0 else { throw zero }
    let scaled = vector/scale
    return scaled/sqrt((scaled*scaled).sum())
  }
  private func cross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
    SIMD3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
  }
  internal func matrix(at epoch: Epoch) throws -> [Double] {
    let s = try state(spacecraft,at: epoch)
    let centre = try centralBody.map { try state($0,at: epoch) }
    let radial = try unit(s.position-(centre?.position ?? .zero),zero: .zeroPosition)
    let velocity = try unit(s.velocity-(centre?.velocity ?? .zero),zero: .zeroVelocity)
    let normal = cross(radial,velocity)
    guard sqrt((normal*normal).sum()) >= minimumPlaneSine else { throw LocalOrbitalFrameError.undefinedOrbitalPlane }
    let zNormal = try unit(normal,zero: .undefinedOrbitalPlane)
    let x: SIMD3<Double>, y: SIMD3<Double>, z: SIMD3<Double>
    switch axes {
    case .rtn: x = radial; z = zNormal; y = cross(z,x)
    case .tnw: x = velocity; z = zNormal; y = cross(z,x)
    case .lvlhCCSDS: z = -radial; y = -zNormal; x = cross(y,z)
    }
    // Rows are child axes in parent coordinates: p_child = R p_parent.
    return [x.x,x.y,x.z,y.x,y.y,y.z,z.x,z.y,z.z]
  }
  public func orientation(at epoch: Epoch) throws -> BodyOrientation {
    let h = derivativeStep.seconds, r = try matrix(at: epoch)
    let samples = try [-2.0,-1,1,2].map {
      try matrix(at: epoch.advanced(by: Duration(seconds: $0*h)))
    }
    let d = (0..<9).map { (samples[0][$0]-8*samples[1][$0]+8*samples[2][$0]-samples[3][$0])/(12*h) }
    let dd = (0..<9).map { (-samples[3][$0]+16*samples[2][$0]-30*r[$0]+16*samples[1][$0]-samples[0][$0])/(12*h*h) }
    func axial(_ derivative: [Double]) -> SIMD3<Double> {
      func entry(_ i: Int,_ j: Int) -> Double { (0..<3).reduce(0) { $0+derivative[3*i+$1]*r[3*j+$1] } }
      return SIMD3((entry(2,1)-entry(1,2))/2,(entry(0,2)-entry(2,0))/2,(entry(1,0)-entry(0,1))/2)
    }
    let orbital = try FrameTransform(epoch: epoch, rotation: quaternion(r),
      angularVelocity: axial(d), angularAcceleration: axial(dd))
    let mounted = try orbital.followed(by: FrameTransform(epoch: epoch, rotation: mounting.rotation))
    return try BodyOrientation(referenceFrame: referenceFrame, epoch: epoch, rotation: mounted.rotation,
      angularVelocity: mounted.angularVelocity, angularAcceleration: mounted.angularAcceleration)
  }
  private func quaternion(_ m: [Double]) -> Quaternion {
    let candidates = [1+m[0]+m[4]+m[8],1+m[0]-m[4]-m[8],1-m[0]+m[4]-m[8],1-m[0]-m[4]+m[8]]
    let k = candidates.indices.max { candidates[$0] < candidates[$1] }!
    let v = sqrt(candidates[k])/2, f = 1/(4*v)
    switch k {
    case 0: return Quaternion(x:(m[7]-m[5])*f,y:(m[2]-m[6])*f,z:(m[3]-m[1])*f,w:v)
    case 1: return Quaternion(x:v,y:(m[1]+m[3])*f,z:(m[2]+m[6])*f,w:(m[7]-m[5])*f)
    case 2: return Quaternion(x:(m[1]+m[3])*f,y:v,z:(m[5]+m[7])*f,w:(m[2]-m[6])*f)
    default: return Quaternion(x:(m[2]+m[6])*f,y:(m[5]+m[7])*f,z:v,w:(m[3]-m[1])*f)
    }
  }
}

extension FrameTree {
  /// Spacecraft-centred orbital frame, with orbit defined relative to centralBody
  /// or the parent origin. Translation uses the spacecraft's absolute parent state.
  public func makeLocalOrbitalFrame(name: String, spacecraft: any StateProvider,
    centralBody: (any StateProvider)? = nil, axes: LocalOrbitalAxes, derivativeStep: Duration,
    minimumPlaneSine: Double = 1e-12,
    accelerationPolicy: OriginAccelerationPolicy = .requireProvided) throws -> Frame {
    let orientation = try LocalOrbitalOrientationProvider(spacecraft: spacecraft,centralBody: centralBody,
      axes: axes,derivativeStep: derivativeStep,minimumPlaneSine: minimumPlaneSine)
    return try makeSpacecraftFrame(name: name,spacecraft: spacecraft,attitude: orientation,
      accelerationPolicy: accelerationPolicy)
  }
  /// Spacecraft-centred axes prescribed by any attitude provider, including a
  /// LocalOrbitalOrientationProvider with a fixed mounting rotation. Attitude rates
  /// are coordinate-map rates in spacecraft axes, as in BodyOrientation.
  public func makeSpacecraftFrame(name: String, spacecraft: any StateProvider,
    attitude: any OrientationProvider,
    accelerationPolicy: OriginAccelerationPolicy = .requireProvided) throws -> Frame {
    try makeBodyFixedFrame(name: name,centre: spacecraft,orientation: attitude,
      accelerationPolicy: accelerationPolicy)
  }
}

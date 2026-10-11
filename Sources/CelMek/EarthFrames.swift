// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum EarthOrientationError: Error, Equatable {
  case outsideCoverage, invalidParameters, invalidDerivativeStep, missingRealization
  case invalidTable, excessiveGap, invalidFile(line: Int), correctionDidNotConverge
}

/// Model to which measured celestial pole offsets refer.
public enum CelestialPoleReference: Sendable { case iau2000A, iau2006A }

/// Earth orientation at an atomic epoch. Angles are radians. UT1−TT is supplied
/// directly, avoiding implicit UTC/leap-second assumptions.
public struct EarthOrientationSample {
  public let epoch: Epoch
  public let ut1MinusTT: Duration
  public let polarMotionX: Double
  public let polarMotionY: Double
  public let celestialPoleX: Double
  public let celestialPoleY: Double
  public let poleReference: CelestialPoleReference
  public init(epoch: Epoch, ut1MinusTT: Duration, polarMotionX: Double, polarMotionY: Double,
              celestialPoleX: Double = 0, celestialPoleY: Double = 0,
              poleReference: CelestialPoleReference = .iau2006A) throws {
    guard [polarMotionX, polarMotionY, celestialPoleX, celestialPoleY].allSatisfy(\.isFinite) else { throw EarthOrientationError.invalidParameters }
    self.epoch = epoch; self.ut1MinusTT = ut1MinusTT
    self.polarMotionX = polarMotionX; self.polarMotionY = polarMotionY
    self.celestialPoleX = celestialPoleX; self.celestialPoleY = celestialPoleY
    self.poleReference = poleReference
  }
}

/// Implementations own coverage, interpolation and data provenance. Samples must
/// match the requested epoch; neighbouring samples also determine frame rates.
public protocol EarthOrientationProvider {
  var provenance: String { get }
  func sample(at epoch: Epoch) throws -> EarthOrientationSample
}

/// Explicit approximation for games or reference cases. No extrapolation beyond
/// caller-selected coverage; constant UT1−TT implies one UT1 second per TT second.
public struct ConstantEarthOrientationProvider: EarthOrientationProvider {
  public let provenance: String
  public let coverage: ClosedRange<Epoch>
  private let offset: Duration
  private let x: Double
  private let y: Double
  private let dx: Double
  private let dy: Double
  private let reference: CelestialPoleReference
  public init(ut1MinusTT: Duration, polarMotionX: Double, polarMotionY: Double,
              coverage: ClosedRange<Epoch>, provenance: String,
              celestialPoleX: Double = 0, celestialPoleY: Double = 0,
              poleReference: CelestialPoleReference = .iau2006A) throws {
    _ = try EarthOrientationSample(epoch: coverage.lowerBound, ut1MinusTT: ut1MinusTT,
                                  polarMotionX: polarMotionX, polarMotionY: polarMotionY,
                                  celestialPoleX: celestialPoleX, celestialPoleY: celestialPoleY, poleReference: poleReference)
    dx = celestialPoleX; dy = celestialPoleY; reference = poleReference
    self.offset = ut1MinusTT; x = polarMotionX; y = polarMotionY
    self.coverage = coverage; self.provenance = provenance
  }
  public func sample(at epoch: Epoch) throws -> EarthOrientationSample {
    guard coverage.contains(epoch) else { throw EarthOrientationError.outsideCoverage }
    return try EarthOrientationSample(epoch: epoch, ut1MinusTT: offset, polarMotionX: x, polarMotionY: y,
      celestialPoleX: dx, celestialPoleY: dy, poleReference: reference)
  }
}

/// IAU 2006 precession, IAU 2000A nutation with 2006 adjustments, and IAU 2000
/// Earth rotation/polar motion conventions. Numerical core: native Swift, adapted from ERFA 2.0.1.
/// Angular derivatives use a five-point stencil in atomic SI seconds. Data must
/// cover ±2*derivativeStep, including at the requested coverage boundaries.
public struct EarthFrameTransformProvider: TransformProvider {
  public enum Edge: Int32 { case celestialIntermediate, earthRotation, polarMotion, trueOfDate, apparentSidereal }
  public let edge: Edge
  public let earthOrientation: any EarthOrientationProvider
  public let derivativeStep: Duration
  public init(edge: Edge, earthOrientation: any EarthOrientationProvider,
              derivativeStep: Duration = try! Duration(seconds: 60)) throws {
    guard derivativeStep.seconds > 0, derivativeStep.seconds.isFinite else {
      throw EarthOrientationError.invalidDerivativeStep
    }
    self.edge = edge; self.earthOrientation = earthOrientation; self.derivativeStep = derivativeStep
  }
  internal func matrix(at epoch: Epoch) throws -> [Double] {
    let tt = try epoch.julianDate(in: .tt)
    let eop = try earthOrientation.sample(at: epoch)
    guard eop.epoch == epoch else { throw FrameTransformError.epochMismatch }
    // Preserve the split date; ERFA accepts arbitrary two-part Julian dates.
    return try EarthOrientationModel.matrix(edge: edge, tt1: Double(tt.wholeDay), tt2: tt.fractionalDay,
      ut11: Double(tt.wholeDay), ut12: tt.fractionalDay + eop.ut1MinusTT.seconds / 86400,
      xp: eop.polarMotionX, yp: eop.polarMotionY,
      dx: eop.celestialPoleX, dy: eop.celestialPoleY, poleReference: eop.poleReference)
  }
  public func transform(at epoch: Epoch) throws -> FrameTransform {
    let h = derivativeStep.seconds
    let r = try matrix(at: epoch)
    let samples = try [-2.0, -1, 1, 2].map { try matrix(at: epoch.advanced(by: Duration(seconds: $0*h))) }
    let d = (0..<9).map { (samples[0][$0] - 8*samples[1][$0] + 8*samples[2][$0] - samples[3][$0]) / (12*h) }
    let dd = (0..<9).map { (-samples[3][$0] + 16*samples[2][$0] - 30*r[$0] + 16*samples[1][$0] - samples[0][$0]) / (12*h*h) }
    func axial(_ derivative: [Double]) -> SIMD3<Double> {
      func entry(_ i: Int, _ j: Int) -> Double { (0..<3).reduce(0) { $0 + derivative[3*i+$1]*r[3*j+$1] } }
      return SIMD3((entry(2,1)-entry(1,2))/2, (entry(0,2)-entry(2,0))/2, (entry(1,0)-entry(0,1))/2)
    }
    return try FrameTransform(epoch: epoch, rotation: quaternion(r), angularVelocity: axial(d), angularAcceleration: axial(dd))
  }
  private func quaternion(_ m: [Double]) -> Quaternion {
    // Largest-component construction is stable through half-turn rotations.
    let candidates = [1+m[0]+m[4]+m[8], 1+m[0]-m[4]-m[8], 1-m[0]+m[4]-m[8], 1-m[0]-m[4]+m[8]]
    let k = candidates.indices.max { candidates[$0] < candidates[$1] }!
    let v = sqrt(candidates[k]) / 2, f = 1/(4*v)
    switch k {
    case 0: return Quaternion(x: (m[7]-m[5])*f, y: (m[2]-m[6])*f, z: (m[3]-m[1])*f, w: v)
    case 1: return Quaternion(x: v, y: (m[1]+m[3])*f, z: (m[2]+m[6])*f, w: (m[7]-m[5])*f)
    case 2: return Quaternion(x: (m[1]+m[3])*f, y: v, z: (m[5]+m[7])*f, w: (m[2]-m[6])*f)
    default: return Quaternion(x: (m[2]+m[6])*f, y: (m[5]+m[7])*f, z: v, w: (m[3]-m[1])*f)
    }
  }
}

public struct EarthFrames {
  public let gcrf: Frame
  public let eme2000: Frame
  public let cirf: Frame
  public let tirf: Frame
  public let itrf: Frame
  public let tod: Frame
  public let gtod: Frame
  public let terrestrialRealization: String
  public let earthOrientationProvenance: String
}

extension FrameTree {
  /// Earth-centred frame families. The realization label describes caller EOP;
  /// no station coordinates, tectonic motion or realization Helmert ties are applied.
  public func makeEarthFrames(earth: any StateProvider, earthOrientation: any EarthOrientationProvider,
    terrestrialRealization: String, derivativeStep: Duration = try! Duration(seconds: 60),
    accelerationPolicy: OriginAccelerationPolicy = .requireProvided) throws -> EarthFrames {
    guard !terrestrialRealization.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
      throw EarthOrientationError.missingRealization
    }
    // Validate configuration before registering any frames.
    func provider(_ edge: EarthFrameTransformProvider.Edge) throws -> EarthFrameTransformProvider {
      try EarthFrameTransformProvider(edge: edge, earthOrientation: earthOrientation, derivativeStep: derivativeStep)
    }
    let cirfProvider = try provider(.celestialIntermediate)
    let geo = try makeGeocentricFrames(earth: earth, accelerationPolicy: accelerationPolicy)
    let cirf = try makeFrame(parent: geo.gcrf, name: "CIRF", provider: cirfProvider)
    let tirf = try makeFrame(parent: cirf, name: "TIRF", provider: provider(.earthRotation))
    let itrf = try makeFrame(parent: tirf, name: terrestrialRealization, provider: provider(.polarMotion))
    let tod = try makeFrame(parent: geo.gcrf, name: "TOD", provider: provider(.trueOfDate))
    let gtod = try makeFrame(parent: tod, name: "GTOD", provider: provider(.apparentSidereal))
    return EarthFrames(gcrf: geo.gcrf, eme2000: geo.eme2000, cirf: cirf, tirf: tirf, itrf: itrf,
      tod: tod, gtod: gtod, terrestrialRealization: terrestrialRealization,
      earthOrientationProvenance: earthOrientation.provenance)
  }
}

import Foundation
import Testing
@testable import CelMek
import Math

@Suite struct EarthFrameTests {
  private struct Reference: Decodable {
    let day: Int64
    let fraction: Double
    let edge: Int32
    let matrix: [Double]
  }
  @Test func nativeSwiftMatchesExternalReferenceAcrossCenturies() throws {
    let url = try #require(Bundle.module.url(forResource: "earth-orientation-erfa", withExtension: "json", subdirectory: "Fixtures"))
    let references = try JSONDecoder().decode([Reference].self, from: Data(contentsOf: url))
    #expect(references.count == 150)
    for reference in references {
      let epoch = try Epoch(julianDate: JulianDateComponents(wholeDay: reference.day,
        fractionalDay: reference.fraction, scale: .tt))
      let edge = try #require(EarthFrameTransformProvider.Edge(rawValue: reference.edge))
      let p = try EarthFrameTransformProvider(edge: edge, earthOrientation: eop(at: epoch, offset: -64.184))
      let actual = try p.transform(at: epoch).rotation.matrix
      // Long-baseline ERA arithmetic differs slightly with compiler contraction.
      let tolerance = edge == .earthRotation || edge == .apparentSidereal ? 2e-12 : 2e-14
      for i in 0..<3 { for j in 0..<3 {
        #expect(abs(actual[i,j]-reference.matrix[3*i+j]) < tolerance, "day=\(reference.day) f=\(reference.fraction) edge=\(reference.edge) error=\(actual[i,j]-reference.matrix[3*i+j])")
      }}
    }
  }
  func eop(at epoch: Epoch, offset: Double = 0) throws -> ConstantEarthOrientationProvider {
    try ConstantEarthOrientationProvider(ut1MinusTT: CelMek.Duration(seconds: offset),
      polarMotionX: 2.55060238e-7, polarMotionY: 1.860359247e-6,
      coverage: epoch.advanced(by: CelMek.Duration(seconds: -1000))...epoch.advanced(by: CelMek.Duration(seconds: 1000)),
      provenance: "ERFA validation case; constant approximation")
  }
  @Test func celestialToTerrestrialReferenceAndRouteAgreement() throws {
    let epoch = try Epoch(julianDate: 2453736.5, scale: .tt)
    let data = try eop(at: epoch)
    func transform(_ edge: EarthFrameTransformProvider.Edge) throws -> FrameTransform {
      try EarthFrameTransformProvider(edge: edge, earthOrientation: data).transform(at: epoch)
    }
    let cirf = try transform(.celestialIntermediate)
    let tirf = try cirf.followed(by: transform(.earthRotation))
    let itrf = try tirf.followed(by: transform(.polarMotion))
    // Published ERFA t_c2t06a reference, independent of the Swift composition.
    let expected = [-0.1810332128305897282, 0.9834769806938592296, 0.6555550962998436505e-4,
      -0.9834768134136214897, -0.1810332203649130832, 0.5749800844905594110e-3,
      0.5773474024748545878e-3, 0.3961816829632690581e-4, 0.9999998325501747785]
    for i in 0..<3 { for j in 0..<3 { #expect(abs(itrf.rotation.matrix[i,j]-expected[3*i+j]) < 1e-12) } }
    let gtod = try transform(.trueOfDate).followed(by: transform(.apparentSidereal))
    for i in 0..<3 { for j in 0..<3 { #expect(abs(gtod.rotation.matrix[i,j]-tirf.rotation.matrix[i,j]) < 1e-14) } }
    for i in 0..<3 {
      #expect(abs(gtod.angularVelocity[i]-tirf.angularVelocity[i]) < 1e-16)
      #expect(abs(gtod.angularAcceleration[i]-tirf.angularAcceleration[i]) < 1e-18)
    }
  }
  @Test func trueOfDateReference() throws {
    let epoch = try Epoch(julianDate: JulianDateComponents(wholeDay: 2450124, fractionalDay: 0.4999, scale: .tt))
    let t = try EarthFrameTransformProvider(edge: .trueOfDate, earthOrientation: eop(at: epoch)).transform(at: epoch)
    let expected = [0.9999995832794205484, 0.8372382772630962111e-3, 0.3639684771140623099e-3,
      -0.8372533744743683605e-3, 0.9999996486492861646, 0.4132905944611019498e-4,
      -0.3639337469629464969e-3, -0.4163377605910663999e-4, 0.9999999329094260057]
    for i in 0..<3 { for j in 0..<3 { #expect(abs(t.rotation.matrix[i,j]-expected[3*i+j]) < 1e-14) } }
  }
  @Test func rotationRatesAndStateDerivatives() throws {
    let epoch = Epoch.j2000
    let p = try EarthFrameTransformProvider(edge: .earthRotation, earthOrientation: eop(at: epoch, offset: -64.184))
    let t = try p.transform(at: epoch)
    #expect(abs(t.angularVelocity.z + 7.29211514670698e-5) < 2e-14)
    #expect(abs(t.angularAcceleration.z) < 1e-18)
    let frame = Frame(name: "reference", pseudoInertial: true)
    let state = try CartesianState(epoch: epoch, frame: frame, position: SIMD3(7e6, 2e6, 1e6), velocity: .zero, acceleration: .zero)
    let result = try t.applying(to: state, destination: frame)
    let delta = 0.5
    let plus = try p.transform(at: epoch.advanced(by: CelMek.Duration(seconds: delta))).rotation.matrix * state.position
    let minus = try p.transform(at: epoch.advanced(by: CelMek.Duration(seconds: -delta))).rotation.matrix * state.position
    for i in 0..<3 {
      #expect(abs((plus[i]-minus[i])/(2*delta)-result.velocity[i]) < 1e-5)
      let second = (plus[i] - 2.0 * result.position[i] + minus[i]) / (delta * delta)
      #expect(abs(second - result.acceleration![i]) < 1e-7)
    }
  }
  @Test func frameTreeGeocentreAndCoverage() throws {
    let tree = try CelestialFrames.makeTree()
    let data = try eop(at: .j2000, offset: -64.184)
    let frames = try tree.makeEarthFrames(earth: earth, earthOrientation: data, terrestrialRealization: "ITRF2020",
      accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 60)))
    let centre = try earth.state(at: .j2000, in: frames.itrf, using: tree)
    for i in 0..<3 { #expect(abs(centre.position[i]) < 0.0002); #expect(abs(centre.velocity[i]) < 1e-8) }
    #expect(frames.tirf.parent === frames.cirf && frames.itrf.parent === frames.tirf)
    #expect(frames.gtod.parent === frames.tod && frames.tod.parent === frames.gcrf)
    #expect(frames.earthOrientationProvenance == data.provenance)
    #expect(throws: EarthOrientationError.outsideCoverage) {
      try tree.transform(from: frames.gcrf, to: frames.itrf, at: data.coverage.upperBound)
    }
    #expect(throws: EarthOrientationError.invalidDerivativeStep) {
      try EarthFrameTransformProvider(edge: .earthRotation, earthOrientation: data, derivativeStep: .zero)
    }
    #expect(throws: EarthOrientationError.invalidParameters) {
      try EarthOrientationSample(epoch: .j2000, ut1MinusTT: .zero, polarMotionX: .nan, polarMotionY: 0)
    }
  }
  private struct DriftingUT1: EarthOrientationProvider {
    let provenance = "Synthetic UT1 rate validation"
    func sample(at epoch: Epoch) throws -> EarthOrientationSample {
      let t = epoch.offsetFromJ2000.seconds
      return try EarthOrientationSample(epoch: epoch,
        ut1MinusTT: CelMek.Duration(seconds: -64.184 + 0.001*t + 1e-7*t*t),
        polarMotionX: 0, polarMotionY: 0)
    }
  }
  @Test func eopTimeVariationDrivesBothAngularDerivatives() throws {
    let provider = try EarthFrameTransformProvider(edge: .earthRotation, earthOrientation: DriftingUT1())
    let t = try provider.transform(at: .j2000)
    let rate = 7.29211514670698e-5
    #expect(abs(t.angularVelocity.z + rate*1.001) < 2e-14)
    #expect(abs(t.angularAcceleration.z + rate*2e-7) < 1e-18)
  }
}

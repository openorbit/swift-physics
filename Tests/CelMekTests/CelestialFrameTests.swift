// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import CelMek
import Math

@Suite
struct CelestialFrameTests {
  func matrixClose(_ actual: Matrix3x3, _ expected: [Double], tolerance: Double) {
    for row in 0..<3 { for col in 0..<3 {
      #expect(abs(actual[row,col] - expected[3*row+col]) < tolerance)
    }}
  }
  @Test
  func biasAndFK5TieMatchERFAReferenceMatrices() throws {
    let tree = try CelestialFrames.makeTree()
    let bias = try tree.transform(from: CelestialFrames.icrf, to: CelestialFrames.meanEquatorJ2000, at: .j2000)
    matrixClose(bias.rotation.matrix, [
      0.9999999999999942498, -0.7078279744199196626e-7, 0.8056217146976134152e-7,
      0.7078279477857337206e-7, 0.9999999999999969484, 0.3306041454222136517e-7,
      -0.8056217380986972157e-7, -0.3306040883980552500e-7, 0.9999999999999962084], tolerance: 3e-16)
    let tie = try tree.transform(from: CelestialFrames.fk5J2000, to: CelestialFrames.icrf, at: .j2000)
    matrixClose(tie.rotation.matrix, [
      0.9999999999999928638, 0.1110223351022919694e-6, 0.4411803962536558154e-7,
      -0.1110223308458746430e-6, 0.9999999999999891830, -0.9647792498984142358e-7,
      -0.4411805033656962252e-7, 0.9647792009175314354e-7, 0.9999999999999943728], tolerance: 3e-16)
    #expect(CelestialFrames.fk5J2000 !== CelestialFrames.meanEquatorJ2000)
  }

  @Test
  func nativeVSOPRotationMatchesPublishedMatrixAndRoundTrips() throws {
    let tree = try CelestialFrames.makeTree()
    let rotation = try tree.transform(from: CelestialFrames.vsop87E, to: CelestialFrames.fk5J2000, at: .j2000)
    matrixClose(rotation.rotation.matrix, [1, 0.000000440360, -0.000000190919,
      -0.000000479966, 0.917482137087, -0.397776982902,
      0, 0.397776982902, 0.917482137087], tolerance: 8e-13)
    let native = try earth.state(at: .j2000)
    let icrf = try earth.state(at: .j2000, in: CelestialFrames.icrf, using: tree)
    #expect(icrf.frame === CelestialFrames.icrf && icrf.epoch == native.epoch)
    let restored = try tree.convert(icrf, to: native.frame)
    for axis in 0..<3 {
      #expect(abs(restored.position[axis]-native.position[axis]) < 0.0002)
      #expect(abs(restored.velocity[axis]-native.velocity[axis]) < 1e-10)
    }
    #expect(icrf.acceleration == nil)
    let ecliptic = try tree.transform(from: CelestialFrames.eclipticJ2000, to: CelestialFrames.meanEquatorJ2000, at: .j2000)
    let vector = ecliptic.rotation.matrix * SIMD3<Double>(0, 1, 0)
    #expect(abs(vector.y - cos(84381.448 * .pi/(180*3600))) < 1e-15)
    #expect(vector.z > 0)
  }

  @Test
  func geocentreAndIAUBodyFramesComposeThroughNativeModel() throws {
    let tree = try CelestialFrames.makeTree()
    let geocentric = try tree.makeGeocentricFrames(earth: earth, accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 60)))
    let gcrf = geocentric.gcrf
    let centre = try earth.state(at: .j2000, in: gcrf, using: tree)
    for axis in 0..<3 {
      #expect(abs(centre.position[axis]) < 0.0002)
      #expect(abs(centre.velocity[axis]) < 1e-10)
    }
    #expect(gcrf.pseudoInertial && gcrf.parent === CelestialFrames.icrf)
    #expect(geocentric.eme2000.parent === gcrf && geocentric.eme2000.pseudoInertial)
    let emeCentre = try earth.state(at: .j2000, in: geocentric.eme2000, using: tree)
    for axis in 0..<3 { #expect(abs(emeCentre.position[axis]) < 0.0002) }
    let bias = try tree.transform(from: gcrf, to: geocentric.eme2000, at: .j2000)
    let barycentricBias = try tree.transform(from: CelestialFrames.icrf, to: CelestialFrames.meanEquatorJ2000, at: .j2000)
    for row in 0..<3 { for col in 0..<3 {
      #expect(bias.rotation.matrix[row,col] == barycentricBias.rotation.matrix[row,col])
    }}
    let start = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: -100))
    let end = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: 100))
    let orientation = try IAUOrientationProvider(body: "mars", coverage: start...end)
    let origin = try FrameConvertedStateProvider(mars, in: CelestialFrames.icrf, using: tree)
    let bodyFrame = try tree.makeBodyFixedFrame(name: "IAU Mars", centre: origin, orientation: orientation,
      accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 60)))
    let result = try mars.state(at: .j2000, in: bodyFrame, using: tree)
    for axis in 0..<3 {
      #expect(abs(result.position[axis]) < 0.0002)
      #expect(abs(result.velocity[axis]) < 1e-7)
    }
  }

  @Test
  func convertedProvidersDoNotCreateOwnershipCycles() throws {
    var tree: FrameTree? = try CelestialFrames.makeTree()
    weak let observed = tree
    let adapter = try FrameConvertedStateProvider(earth, in: CelestialFrames.icrf, using: tree!)
    _ = try tree!.makeGCRF(earth: earth, accelerationPolicy: .assumeZero)
    tree = nil
    #expect(observed == nil)
    #expect(throws: FrameQueryError.frameTreeUnavailable) { try adapter.state(at: .j2000) }
    let local = try FrameTree(root: CelestialFrames.vsop87E)
    #expect(throws: FrameTransformError.disconnectedFrames) {
      try FrameConvertedStateProvider(earth, in: CelestialFrames.icrf, using: local)
    }
  }

  @Test
  func orientationCoverageAndLegacyOptInAreEnforced() throws {
    let epoch = Epoch.j2000
    #expect(IAUOrientationProvider.availableBodies.count == 73)
    #expect(throws: IAUOrientationError.legacyModelRequiresOptIn) {
      try IAUOrientationProvider(body: "earth", coverage: epoch...epoch)
    }
    let legacy = try IAUOrientationProvider(body: "moon", coverage: epoch...epoch, allowLegacyModel: true)
    #expect(legacy.usesLegacyModel)
    let mars = try IAUOrientationProvider(body: "mars", coverage: epoch...epoch)
    #expect(mars.referenceFrame === CelestialFrames.icrf)
    #expect(throws: IAUOrientationError.outsideCoverage) {
      try mars.orientation(at: epoch.advanced(by: CelMek.Duration(seconds: 1)))
    }
    #expect(throws: IAUOrientationError.unknownBody) {
      try IAUOrientationProvider(body: "missing", coverage: epoch...epoch)
    }
  }
}

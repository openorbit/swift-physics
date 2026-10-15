import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct BodyModelTests {
  private func source() throws -> PhysicalDataProvenance {
    try PhysicalDataProvenance(source: "test field and geometry",version: "1")
  }
  private func shape(_ radii: SIMD3<Double>, body: BodyIdentifier = .earth) throws -> EllipsoidBodyShape {
    try EllipsoidBodyShape(bodyIdentifier: body,referenceRadii: BodyReferenceRadii(metres: radii,
      axesDescription: "body-fixed XYZ",provenance: source()))
  }
  @Test func sphereAndEllipsoidSurfaceNormals() throws {
    let sphere = try shape(SIMD3(repeating: 10))
    let p = try sphere.surface(inDirection: SIMD3(3,4,0))
    #expect(p.positionMetres == SIMD3(6,8,0))
    #expect(abs(p.outwardNormal.x-0.6) < 1e-15)
    #expect(abs(p.outwardNormal.y-0.8) < 1e-15)
    let ellipsoid = try shape(SIMD3(3,2,1))
    let diagonal = try ellipsoid.surface(inDirection: SIMD3(1,1,1))
    for i in 0..<3 { #expect(abs(diagonal.positionMetres[i]-6.0/7) < 1e-15) }
    #expect(abs(diagonal.outwardNormal.y/diagonal.outwardNormal.x-9.0/4) < 1e-14)
    #expect(abs(diagonal.outwardNormal.z/diagonal.outwardNormal.x-9) < 1e-14)
    #expect(try ellipsoid.surface(inDirection: SIMD3(0,-1,0)).positionMetres == SIMD3(0,-2,0))
    let scaled = try sphere.surface(inDirection: SIMD3(3e200,4e200,0))
    for i in 0..<3 { #expect(abs(scaled.positionMetres[i]-p.positionMetres[i]) < 1e-14) }
    #expect(throws: BodyModelError.invalidDirection) { try sphere.surface(inDirection: .zero) }
    #expect(throws: BodyModelError.invalidDirection) { try sphere.surface(inDirection: SIMD3(.nan,1,0)) }
  }
  @Test func pointMassFieldAndPotentialGradient() throws {
    let gm = try GravitationalParameter(cubicMetresPerSecondSquared: 4e14,convention: .newtonian,provenance: source())
    let gravity = try PointMassBodyGravity(bodyIdentifier: .earth,parameter: gm)
    let s = try gravity.gravity(at: .j2000,displacementMetres: SIMD3(1e7,0,0))
    #expect(s.acceleration == SIMD3(-4,0,0))
    #expect(s.potential == -4e7)
    #expect(s.convention == .newtonian && s.epoch == .j2000)
    let far = try gravity.gravity(at: .j2000,displacementMetres: SIMD3(2e7,0,0))
    #expect(far.acceleration == SIMD3(-1,0,0))
    let r = SIMD3<Double>(3e6,-4e6,5e6), h = 10.0
    let field = try gravity.gravity(at: .j2000,displacementMetres: r)
    for i in 0..<3 {
      var plus = r, minus = r; plus[i] += h; minus[i] -= h
      let p = try gravity.gravity(at: .j2000,displacementMetres: plus).potential
      let m = try gravity.gravity(at: .j2000,displacementMetres: minus).potential
      #expect(abs(field.acceleration[i]+(p-m)/(2*h)) < 1e-8)
    }
    #expect(throws: BodyModelError.singularGravity) { try gravity.gravity(at: .j2000,displacementMetres: .zero) }
    #expect(throws: BodyModelError.invalidPosition) { try gravity.gravity(at: .j2000,displacementMetres: SIMD3(.infinity,0,0)) }
  }
  @Test func modelsStayIndependentAndRejectIdentityMismatch() throws {
    let body = try BodyCatalogue.jplGM2022.body(.earth)
    let empty = try CelestialBodyModels(body: body)
    #expect(empty.gravity == nil && empty.shape == nil && empty.ephemeris == nil && empty.orientation == nil)
    let gravity = try PointMassBodyGravity(bodyIdentifier: .earth,parameter: #require(body.properties.gravitationalParameter))
    let geometry = try shape(SIMD3(repeating: 10))
    let models = try CelestialBodyModels(body: body,shape: geometry,gravity: gravity)
    #expect(models.gravity?.convention == .tdbCompatible)
    #expect(models.body.properties.referenceRadii == nil)
    #expect(throws: BodyModelError.bodyMismatch) {
      try CelestialBodyModels(body: body,shape: shape(SIMD3(repeating: 10),body: .moon))
    }
    let lunarGravity = try PointMassBodyGravity(bodyIdentifier: .moon,parameter: gravity.parameter)
    #expect(throws: BodyModelError.bodyMismatch) { try CelestialBodyModels(body: body,gravity: lunarGravity) }
  }
  private struct WrongOrientation: OrientationProvider {
    let referenceFrame: Frame
    let wrongFrame: Frame
    var wrongEpoch: Bool
    func orientation(at epoch: Epoch) throws -> BodyOrientation {
      try BodyOrientation(referenceFrame: wrongEpoch ? referenceFrame : wrongFrame,
        epoch: wrongEpoch ? epoch.advanced(by: CelMek.Duration(seconds: 1)) : epoch,rotation: .identity)
    }
  }
  @Test func identifiedAttitudePreservesKinematicsAndValidatesResults() throws {
    let frame = Frame(name: "reference",pseudoInertial: true)
    let spin = try FixedAxisOrientationProvider(referenceFrame: frame,axis: SIMD3(0,0,1),rate: 0.1)
    let wrapper = try IdentifiedBodyOrientationProvider(bodyIdentifier: .earth,provider: spin,provenance: source())
    let actual = try wrapper.orientation(at: .j2000), expected = try spin.orientation(at: .j2000)
    #expect(actual.angularVelocity == expected.angularVelocity)
    #expect(actual.referenceFrame === frame)
    let body = try BodyCatalogue.jplGM2022.body(.earth)
    #expect(try CelestialBodyModels(body: body,orientation: wrapper).orientation?.bodyIdentifier == .earth)
    let wrongBody = try IdentifiedBodyOrientationProvider(bodyIdentifier: .moon,provider: spin,provenance: source())
    #expect(throws: BodyModelError.bodyMismatch) { try CelestialBodyModels(body: body,orientation: wrongBody) }
    for wrongEpoch in [true,false] {
      let p = try IdentifiedBodyOrientationProvider(bodyIdentifier: .earth,
        provider: WrongOrientation(referenceFrame: frame,wrongFrame: Frame(name: "wrong",pseudoInertial: true),wrongEpoch: wrongEpoch),
        provenance: source())
      if wrongEpoch { #expect(throws: FrameTransformError.epochMismatch) { try p.orientation(at: .j2000) } }
      else { #expect(throws: BodyFixedFrameError.frameMismatch) { try p.orientation(at: .j2000) } }
    }
  }
}

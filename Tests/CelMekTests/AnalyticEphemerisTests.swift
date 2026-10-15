import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct AnalyticEphemerisTests {
  private struct Centre: StateProvider {
    let frame = CelestialFrames.icrf
    let position: SIMD3<Double>
    func state(at epoch: Epoch) throws -> CartesianState {
      try CartesianState(epoch: epoch,frame: frame,position: position,velocity: .zero,acceleration: .zero)
    }
  }
  private func setup() throws -> (FrameTree,AnalyticEphemerisFrames) {
    let tree = try CelestialFrames.makeTree()
    let frames = try tree.makeAnalyticEphemerisFrames(earth: Centre(position: SIMD3(1e11,2e11,3e10)),
      sun: Centre(position: SIMD3(1e8,2e8,3e8)))
    return (tree,frames)
  }
  private func adapter(_ model: SphericalAnalyticEphemeris.Model, epoch: Epoch,
    frames: AnalyticEphemerisFrames, step: Double = 30) throws -> SphericalAnalyticEphemeris {
    try SphericalAnalyticEphemeris(model: model,frames: frames,
      positionCoverage: epoch.advanced(by: CelMek.Duration(seconds: -1000))...epoch.advanced(by: CelMek.Duration(seconds: 1000)),
      velocityStep: CelMek.Duration(seconds: step))
  }
  @Test func sphericalPublishedCasesAndUnits() throws {
    let (_,frames) = try setup()
    for (model,jd,longitude,latitude,distance) in [
      (SphericalAnalyticEphemeris.Model.meeusMoon,2448724.5,133.162659,-3.229127,368409700.0),
      (.goffinPluto,2448908.5,232.74009,14.58769,29.711383*CartesianUnits.metresPerAstronomicalUnit)] {
      let epoch = try Epoch(julianDate: jd,scale: .tt)
      let p = try adapter(model,epoch: epoch,frames: frames)
      let s = try p.checkedState(at: epoch)
      let r = sqrt((s.position*s.position).sum())
      let lon = atan2(s.position.y,s.position.x)*180/Double.pi
      let lat = asin(s.position.z/r)*180/Double.pi
      let expected = longitude > 180 ? longitude-360 : longitude
      #expect(abs(lon-expected) < 3e-5)
      #expect(abs(lat-latitude) < 3e-5)
      #expect(abs(r-distance) < (model == .meeusMoon ? 100 : 150000))
      #expect(s.acceleration == nil)
      #expect(p.metadata.origin == (model == .meeusMoon ? .earth : .sun))
      #expect(p.metadata.body == (model == .meeusMoon ? .moon : .plutoBarycentre))
    }
  }
  @Test func frameConversionsAndDifferentiatedICRFStates() throws {
    let (tree,frames) = try setup(), epoch = Epoch.j2000
    for model in [SphericalAnalyticEphemeris.Model.meeusMoon,.goffinPluto] {
      let p = try adapter(model,epoch: epoch,frames: frames)
      let other = try adapter(model,epoch: epoch,frames: frames,step: 60)
      let s = try p.state(at: epoch), second = try other.state(at: epoch)
      for i in 0..<3 { #expect(abs(s.velocity[i]-second.velocity[i]) < 0.001) }
      let converted = try p.state(at: epoch,in: CelestialFrames.icrf,using: tree)
      let before = try p.state(at: epoch.advanced(by: CelMek.Duration(seconds: -120)),in: CelestialFrames.icrf,using: tree)
      let after = try p.state(at: epoch.advanced(by: CelMek.Duration(seconds: 120)),in: CelestialFrames.icrf,using: tree)
      let numerical = (after.position-before.position)/240
      for i in 0..<3 { #expect(abs(converted.velocity[i]-numerical[i]) < 0.01) }
      let back = try tree.transform(from: converted.frame,to: p.frame,at: epoch).applying(to: converted,destination: p.frame)
      for i in 0..<3 {
        #expect(abs(back.position[i]-s.position[i]) < 0.01)
        #expect(abs(back.velocity[i]-s.velocity[i]) < 1e-8)
      }
    }
  }
  @Test func coverageStencilAndVSOPIntegration() throws {
    let (_,frames) = try setup()
    let p = try adapter(.meeusMoon,epoch: .j2000,frames: frames)
    #expect(try p.metadata.coverage.lowerBound.duration(since: p.positionCoverage.lowerBound).seconds == 60)
    _ = try p.state(at: p.metadata.coverage.lowerBound)
    _ = try p.state(at: p.metadata.coverage.upperBound)
    #expect(throws: EphemerisError.outsideCoverage) { try p.state(at: p.positionCoverage.lowerBound) }
    #expect(throws: EphemerisError.invalidDerivativeMethod) { try adapter(.meeusMoon,epoch: .j2000,frames: frames,step: 0) }
    let tdb = try FairheadTDBModel()
    let vsop = try VSOP87Ephemeris(body: earth,coverage: p.metadata.coverage,tdbModel: tdb)
    #expect(vsop.accuracy == .full && vsop.metadata.accuracy == .full)
    #expect(vsop.availableAccuracyLevels == Set(AccuracyLevel.allCases))
    #expect(p.accuracy == .full && p.metadata.accuracy == .full)
    let direct = try earth.state(at: .j2000,using: tdb)
    #expect(try vsop.state(at: .j2000) == direct)
    let collection = try EphemerisCollection(providers: [vsop,p])
    #expect(collection.bodies == [.earth,.moon])
    #expect(try collection.state(of: .moon,at: .j2000).frame === frames.moonEclipticOfDate)
  }

  @Test func meanDatePrecessionMatchesPublishedERFA() throws {
    // ERFA 2.0.1 t_pmat76; BSD attribution in Documentation/licenses/ERFA.txt.
    let epoch = try Epoch(julianDate: JulianDateComponents(wholeDay: 2450124,fractionalDay: 0.4999,scale: .tt))
    let t = try centuriesTT(epoch)
    let eps = (84381.448-46.815*t-0.00059*t*t+0.001813*t*t*t)*Double.pi/(180*3600)
    let q = try MeanEclipticOfDateProvider(step: CelMek.Duration(seconds: 60)).rotation(at: epoch)
    let precession = (Quaternion(angle: eps,axis: SIMD3(1,0,0))*q).matrix
    let expected = [0.9999995504328350733,0.8696632209480960785e-3,0.3779153474959888345e-3,
      -0.8696632209485112192e-3,0.9999996218428560614,-0.1643284776111886407e-6,
      -0.3779153474950335077e-3,-0.1643306746147366896e-6,0.9999999285899790119]
    for i in 0..<3 { for j in 0..<3 { #expect(abs(precession[i,j]-expected[3*i+j]) < 2e-14) } }
  }


  @Test(arguments: [AccuracyLevel.low,.balanced,.high])
  func accuracyHintsMayUseFullEvaluation(_ level: AccuracyLevel) throws {
    let (_,frames) = try setup()
    let coverage = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: -1000))...Epoch.j2000.advanced(by: CelMek.Duration(seconds: 1000))
    for model in [SphericalAnalyticEphemeris.Model.meeusMoon,.goffinPluto] {
      let p = try SphericalAnalyticEphemeris(model: model,frames: frames,positionCoverage: coverage,
        velocityStep: CelMek.Duration(seconds: 30),accuracy: level)
      let full = try adapter(model,epoch: .j2000,frames: frames)
      let actual = try p.checkedState(at: .j2000), expected = try full.state(at: .j2000)
      #expect(actual == expected)
      #expect(p.accuracy == level && p.effectiveAccuracy == .full)
      #expect(p.metadata.accuracy == level && p.metadata.effectiveAccuracy == .full)
    }
  }
  @Test func preparedVSOPSubsetsBoundPositionErrorAndRetainAnalyticVelocity() throws {
    let tdb = try FairheadTDBModel()
    let coverage = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: -86400*365))...Epoch.j2000.advanced(by: CelMek.Duration(seconds: 86400*365))
    let full = try VSOP87Ephemeris(body: earth,coverage: coverage,tdbModel: tdb)
    var lastCount = 0
    for level in AccuracyLevel.allCases {
      let p = try VSOP87Ephemeris(body: earth,coverage: coverage,tdbModel: tdb,accuracy: level)
      let selection = p.coefficientSelection
      #expect(selection.retainedTermCount >= lastCount)
      lastCount = selection.retainedTermCount
      #expect(p.metadata.accuracy == level)
      if level != .full { #expect(selection.retainedTermCount < selection.originalTermCount) }
      else { #expect(selection.positionErrorBoundMetres == .zero) }
      for day in stride(from: -365,through: 365,by: 73) {
        let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: Double(day)*86400))
        let reduced = try p.state(at: epoch), reference = try full.state(at: epoch)
        for i in 0..<3 {
          #expect(abs(reduced.position[i]-reference.position[i]) <= selection.positionErrorBoundMetres[i]+0.05)
        }
      }
      let centre = try p.state(at: .j2000)
      let before = try p.state(at: Epoch.j2000.advanced(by: CelMek.Duration(seconds: -60)))
      let after = try p.state(at: Epoch.j2000.advanced(by: CelMek.Duration(seconds: 60)))
      let differentiated = (after.position-before.position)/120
      for i in 0..<3 { #expect(abs(centre.velocity[i]-differentiated[i]) < 0.001) }
    }
  }
}

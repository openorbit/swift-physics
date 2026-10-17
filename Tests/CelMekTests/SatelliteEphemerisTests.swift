import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct SatelliteEphemerisTests {
  @Test func gustMatchesCorrectedOriginalCAtTheoryEpoch() throws {
    // Original C evaluator compiled independently with z_amp phase lookup fixed.
    let expected: [[Double]] = [
      [4.4434305872331761,-0.21457821013446771,0.0010282835208842648,0.0007921473134023762,0.031705921704991716,-0.020747425290356083],
      [2.4925197958844794,3.0961598451027394,-0.0012138340729567966,0.0015778323839911394,0.00043640965732884968,0.00043270271719131326],
      [1.5159500367588301,2.2861519168688176,-0.0010298527412352335,0.003306501979786579,0.00072396580836163373,0.00084300027895201231],
      [0.72175895283551439,0.85618990354019331,0.0018067067743402058,0.0016838216635544118,-0.00042154997875563162,0.00026667997829451416],
      [0.46663157355102919,-0.91534670559029885,0.00073564126336227588,-0.00080516225419183609,-8.4627810102214058e-05,-0.00096175050800333681]
    ]
    let epoch = try Epoch(julianDate: 2444239.5,scale: .tt)
    for moon in GUST86.Moon.allCases {
      let e = try GUST86.elements(of: moon,at: epoch)
      let actual = [e.meanMotion,e.meanLongitude,e.eccentricity.x,e.eccentricity.y,e.inclination.x,e.inclination.y]
      for (a,b) in zip(actual,expected[moon.rawValue]) { #expect(abs(a-b) < 1e-14) }
    }
  }
  @Test func gustFrequenciesMatchCorrectedOriginalCAwayFromEpoch() throws {
    let expected: [[Double]] = [
      [4.4433252076136478,5487.3428064320633,-0.00048016903537524718,0.0014483297719270387,-0.0077727440457283385,-0.037033873821789284],
      [2.4924208140321422,3080.6483352711339,-0.0015618761678382199,0.0012493936011243251,0.00062449684990241992,0.00031495709711976973],
      [1.5159211957705372,1873.9699926481237,-0.0022686309217629717,0.0030027745082965157,0.00085508877232730101,0.00072150526402348948],
      [0.72172968324651177,891.81807466407588,0.0018147184140608363,0.00015384089948809539,-0.00036256556318418941,0.00029999771696105051],
      [0.46667665098529676,575.21569717109014,0.00070448720477993287,0.0011119586260887966,-0.00015153185676850569,-0.0009712805695904534],
    ]
    let epoch = try Epoch(julianDate: 2444239.5+1234.5,scale: .tt)
    for moon in GUST86.Moon.allCases {
      let e = try GUST86.elements(of: moon,at: epoch)
      let actual = [e.meanMotion,e.meanLongitude,e.eccentricity.x,e.eccentricity.y,e.inclination.x,e.inclination.y]
      for (a,b) in zip(actual,expected[moon.rawValue]) { #expect(abs(a-b) < 1e-11) }
    }
  }
  @Test func allMoonsIntegrateWithCollectionAndRespectCoverage() throws {
    let frames = SatelliteEphemerisFrames()
    let range = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: -1000))...Epoch.j2000.advanced(by: CelMek.Duration(seconds: 1000))
    let models = LieskeE2x3.Moon.allCases.map { SatelliteEphemeris.Model.lieskeE2x3($0) }
      + GUST86.Moon.allCases.map { SatelliteEphemeris.Model.gust86($0) }
    let providers = try models.map { try SatelliteEphemeris(model: $0,frames: frames,positionCoverage: range) }
    let collection = try EphemerisCollection(providers: providers)
    #expect(collection.bodies.count == 9)
    for p in providers {
      let state = try collection.state(of: p.metadata.body,at: .j2000)
      let r = sqrt((state.position*state.position).sum())
      #expect(r > 1e8 && r < 2e9)
      let speed = sqrt((state.velocity*state.velocity).sum())
      #expect(speed > 2000 && speed < 20000)
      #expect(state.acceleration == nil)
      let other = try SatelliteEphemeris(model: p.model,frames: frames,positionCoverage: range,velocityStep: CelMek.Duration(seconds: 60))
      let check = try other.state(at: .j2000)
      for k in 0..<3 { #expect(abs(state.velocity[k]-check.velocity[k]) < 0.001) }
      _ = try p.checkedState(at: p.metadata.coverage.lowerBound)
      _ = try p.checkedState(at: p.metadata.coverage.upperBound)
      #expect(throws: EphemerisError.outsideCoverage) { try p.state(at: range.lowerBound) }
      #expect(try p.state(at: .j2000) == state)
    }
    #expect(throws: EphemerisError.invalidDerivativeMethod) {
      try SatelliteEphemeris(model: .lieskeE2x3(.io),frames: frames,positionCoverage: range,velocityStep: .zero)
    }
    #expect(throws: EphemerisError.outsideCoverage) {
      try SatelliteEphemeris(model: .gust86(.miranda),frames: frames,positionCoverage: Epoch.j2000...Epoch.j2000)
    }
  }
  @Test func lieskeMatchesCorrectedOriginalCAtJ2000() throws {
    // Independent original-C evaluation with documented port corrections.
    let expected: [SIMD3<Double>] = [
      SIMD3(206001650.07642466,366911152.16646528,-60819.928830089586),
      SIMD3(-170110037.92347211,-643621535.06931674,-864082.29922565143),
      SIMD3(-137176168.12921852,-1062835871.0664721,3459344.4942878205),
      SIMD3(-1021262555.564065,1583354181.5214057,-3129269.90261812)
    ]
    for moon in LieskeE2x3.Moon.allCases {
      let actual = try LieskeE2x3.position(of: moon,at: .j2000)
      for k in 0..<3 { #expect(abs(actual[k]-expected[moon.rawValue][k]) < 0.001) }
    }
  }
  @Test func ioMovesAndQueriesHaveNoSharedMutableState() throws {
    let start = try LieskeE2x3.position(of: .io,at: .j2000)
    let later = try LieskeE2x3.position(of: .io,at: .j2000.advanced(by: CelMek.Duration(seconds: 1000)))
    #expect(sqrt(((later-start)*(later-start)).sum()) > 1e7)
    _ = try LieskeE2x3.position(of: .callisto,at: .j2000)
    #expect(try LieskeE2x3.position(of: .io,at: .j2000) == start)
  }
}

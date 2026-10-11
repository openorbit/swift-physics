import Foundation
import Testing
import CelMek
import Math

@Suite struct TDBTimeModelTests {
  private struct Reference: Decodable {
    let daysSinceJ2000: Double
    let tdbMinusTT: Double
    let rate: Double
    let acceleration: Double
  }
  @Test func fullSeriesMatchesIndependentReferenceAndDerivatives() throws {
    let model = try FairheadTDBModel()
    let url = try #require(Bundle.module.url(forResource: "tdb-erfa",withExtension: "json",subdirectory: "Fixtures"))
    let rows = try JSONDecoder().decode([Reference].self,from: Data(contentsOf: url))
    #expect(rows.count == 4874)
    for row in rows {
      let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: row.daysSinceJ2000*86400))
      let sample = try epoch.tdbCorrection(using: model)
      #expect(abs(sample.tdbMinusTT.seconds-row.tdbMinusTT) < 2e-15)
      // Independent centred differences have a finite 300-second stencil error.
      #expect(abs((sample.rate-1)-row.rate) < 3e-16)
      #expect(abs(sample.acceleration-row.acceleration) < 2e-21)
      do {
        let date = try epoch.julianDate(in: .tdb,using: model)
        let recovered = try Epoch(julianDate: date,using: model)
        #expect(abs(try recovered.duration(since: epoch).seconds) < 3e-11)
      }
    }
  }
  @Test func coverageAndModelContract() throws {
    let model = try FairheadTDBModel()
    for boundary in [model.coverage.lowerBound,model.coverage.upperBound] {
      #expect(try model.correction(at: boundary).epoch == boundary)
    }
    let outside = try model.coverage.lowerBound.advanced(by: CelMek.Duration(seconds: -1))
    #expect(throws: TDBModelError.outsideCoverage) { try outside.julianDate(in: .tdb,using: model) }
    #expect(throws: TDBModelError.outsideCoverage) { try FairheadTDBModel(coverage: outside...model.coverage.upperBound) }
    #expect(throws: TDBModelError.invalidCorrection) {
      try TDBCorrection(epoch: .j2000,tdbMinusTT: .zero,rate: .nan,acceleration: 0)
    }
    let cheap = TwoTermTDBModel(coverage: model.coverage)
    let defaultDate = try Epoch.j2000.julianDate(in: .tdb)
    let selectedDate = try Epoch.j2000.julianDate(in: .tdb,using: cheap)
    #expect(defaultDate == selectedDate)
  }
  private struct LinearModel: TDBTimeModel {
    let coverage: ClosedRange<Epoch>
    let provenance = "Synthetic exact linear model"
    func correction(at epoch: Epoch) throws -> TDBCorrection {
      let seconds = epoch.offsetFromJ2000.seconds
      return try TDBCorrection(epoch: epoch,tdbMinusTT: CelMek.Duration(seconds: 0.125+0.001*seconds),rate: 1.001,acceleration: 0)
    }
  }
  @Test func inversionUsesSelectedDerivativeAndHandlesCoverageBoundaries() throws {
    let limit = try CelMek.Duration(seconds: 100)
    let lower = try Epoch.j2000.advanced(by: limit.negated())
    let upper = try Epoch.j2000.advanced(by: limit)
    let model = LinearModel(coverage: lower...upper)
    for epoch in [model.coverage.lowerBound,Epoch.j2000,model.coverage.upperBound] {
      let date = try epoch.julianDate(in: .tdb,using: model)
      let recovered = try Epoch(julianDate: date,using: model)
      #expect(abs(try recovered.duration(since: epoch).seconds) < 3e-11)
    }
    let reading = try Epoch.j2000.calendar(in: .tdb,tdbModel: model)
    #expect(abs(try reading.epoch(tdbModel: model).duration(since: .j2000).seconds) < 3e-11)
  }
  private struct FaultyModel: TDBTimeModel {
    let coverage = FairheadTDBModel.referenceCoverage
    let provenance = "Synthetic invalid contract"
    let wrongEpoch: Bool
    func correction(at epoch: Epoch) throws -> TDBCorrection {
      try TDBCorrection(epoch: wrongEpoch ? .j2000 : epoch,
        tdbMinusTT: CelMek.Duration(seconds: 0.125+epoch.offsetFromJ2000.seconds),
        rate: 1,acceleration: 0) // Incorrect derivative intentionally prevents Newton convergence.
    }
  }
  @Test func rejectsWrongEpochAndNonconvergentCustomModels() throws {
    let later = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: 1))
    #expect(throws: TDBModelError.epochMismatch) { try later.tdbCorrection(using: FaultyModel(wrongEpoch: true)) }
    #expect(throws: TDBModelError.noConvergence) {
      try Epoch(julianDate: 2451545,scale: .tdb,using: FaultyModel(wrongEpoch: false))
    }
  }
  @Test func selectedModelFlowsThroughVSOPAndIAUProviders() throws {
    let model = try FairheadTDBModel()
    let provider = VSOP87StateProvider(body: earth,tdbModel: model)
    let tree = try CelestialFrames.makeTree()
    let centre = try provider.state(at: .j2000,in: CelestialFrames.icrf,using: tree)
    #expect(centre.epoch == .j2000 && centre.frame === CelestialFrames.icrf)
    let direct = try earth.state(at: .j2000,using: model)
    #expect(try provider.state(at: .j2000) == direct)
    let delta = 300.0
    let plus = try provider.state(at: Epoch.j2000.advanced(by: CelMek.Duration(seconds: delta)))
    let minus = try provider.state(at: Epoch.j2000.advanced(by: CelMek.Duration(seconds: -delta)))
    for i in 0..<3 { #expect(abs((plus.position[i]-minus.position[i])/(2*delta)-direct.velocity[i]) < 0.0001) }
    let orientation = try IAUOrientationProvider(body: "mars",coverage: model.coverage,tdbModel: model)
    #expect(try orientation.orientation(at: .j2000).epoch == .j2000)
    let before = try model.coverage.lowerBound.advanced(by: CelMek.Duration(seconds: -1))
    #expect(throws: TDBModelError.outsideCoverage) { try provider.state(at: before) }
  }
}

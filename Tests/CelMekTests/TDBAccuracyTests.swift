import Foundation
import Testing
@testable import CelMek

@Suite struct TDBAccuracyTests {
  private struct Reference: Decodable {
    let daysSinceJ2000: Double
    let tdbMinusTT: Double
    let rate: Double
    let acceleration: Double
  }
  @Test func independentGeocentricReferenceFrom1950To2050() throws {
    let url = try #require(Bundle.module.url(forResource: "tdb-erfa", withExtension: "json", subdirectory: "Fixtures"))
    let rows = try JSONDecoder().decode([Reference].self, from: Data(contentsOf: url))
    #expect(rows.count == 4874)
    var maximum = 0.0, squareSum = 0.0
    for row in rows {
      let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: row.daysSinceJ2000*86400))
      let tdb = try epoch.julianDate(in: .tdb)
      // Read the representation's clock offset through split dates, avoiding
      // subtraction of large scalar Julian dates that would mask this error.
      let clock = try Epoch(julianDate: JulianDateComponents(wholeDay: tdb.wholeDay,
        fractionalDay: tdb.fractionalDay, scale: .tt))
      let offset = try clock.duration(since: epoch).seconds
      let error = abs(offset-row.tdbMinusTT)
      maximum = max(maximum,error); squareSum += error*error
      #expect(error < 37e-6)
      #expect(abs((epoch.coordinateTimeRate(in: .tdb)-1)-row.rate) < 1.1e-11)
      #expect(abs(epoch.coordinateTimeAcceleration(in: .tdb)-row.acceleration) < 1.3e-17)
      let restored = try Epoch(julianDate: tdb)
      #expect(abs(try restored.duration(since: epoch).seconds) < 3e-11)
    }
    #expect(maximum > 30e-6)
    #expect(sqrt(squareSum/Double(rows.count)) < 18e-6)
  }
}

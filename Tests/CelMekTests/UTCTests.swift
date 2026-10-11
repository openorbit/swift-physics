import Testing
import CelMek

@Suite struct UTCTests {
  let table = LeapSecondTable.iersBulletinC72
  func utc(_ day: Int64, _ second: Double) throws -> UTCDate {
    try UTCDate(modifiedJulianDay: day, secondOfDay: CelMek.Duration(seconds: second))
  }
  @Test func j2000AndKnownTAIReference() throws {
    let date = try Epoch.j2000.utc(using: table)
    #expect(date.modifiedJulianDay == 51544)
    #expect(abs(date.secondOfDay.seconds-43135.816) < 1e-11)
    let midnight = try Epoch(utc: utc(57754,0), using: table)
    let tai = try Epoch(julianDate: JulianDateComponents(wholeDay: 2457754,
      fractionalDay: 0.5+37.0/86400, scale: .tai))
    #expect(abs(try midnight.duration(since: tai).seconds) < 1e-11)
  }
  @Test func allHistoricalPositiveLeapBoundaries() throws {
    for entry in table.entries.dropFirst() {
      let day = entry.effectiveModifiedJulianDay
      #expect(try table.dayLength(on: day-1) == 86401)
      let before = try Epoch(utc: utc(day-1,86399), using: table)
      let leap = try Epoch(utc: utc(day-1,86400), using: table)
      let after = try Epoch(utc: utc(day,0), using: table)
      #expect(abs(try leap.duration(since: before).seconds-1) < 1e-12)
      #expect(abs(try after.duration(since: leap).seconds-1) < 1e-12)
      for second in [86399.0,86400,86400.25,86400.999999] {
        let original = try utc(day-1,second)
        let restored = try Epoch(utc: original, using: table).utc(using: table)
        #expect(restored.modifiedJulianDay == original.modifiedJulianDay)
        #expect(abs(restored.secondOfDay.seconds-second) < 1e-11)
      }
      #expect(try after.utc(using: table).modifiedJulianDay == day)
    }
  }
  @Test func coverageAndInvalidLeapLabels() throws {
    #expect(throws: UTCError.outsideCoverage) { try Epoch(utc: utc(41316,0), using: table) }
    #expect(throws: UTCError.outsideCoverage) { try Epoch(utc: utc(61584,0), using: table) }
    #expect(throws: UTCError.invalidSecondOfDay) { try Epoch(utc: utc(57754,86400), using: table) }
    #expect(throws: UTCError.invalidSecondOfDay) { try utc(57753,-1) }
    #expect(throws: UTCError.invalidSecondOfDay) { try utc(57753,86401) }
    let fineLabel = try UTCDate(modifiedJulianDay: 57753, secondOfDay: CelMek.Duration(wholeSeconds: 86400, fractionalSecond: 0.999999999999999))
    #expect(fineLabel.secondOfDay.wholeSeconds == 86400)
    let first = try Epoch(utc: utc(41317,0), using: table)
    #expect(throws: UTCError.outsideCoverage) { try first.advanced(by: CelMek.Duration(seconds: -0.001)).utc(using: table) }
    let last = try Epoch(utc: utc(61583,86399.5), using: table)
    #expect(try last.utc(using: table).modifiedJulianDay == 61583)
    let end = try Epoch(utc: utc(61583,0), using: table).advanced(by: CelMek.Duration(seconds: 86400))
    #expect(throws: UTCError.outsideCoverage) { try end.utc(using: table) }
  }
  @Test func syntheticNegativeLeapSecond() throws {
    let negative = try LeapSecondTable(entries: [
      .init(effectiveModifiedJulianDay: 60000,taiMinusUTC: 37),
      .init(effectiveModifiedJulianDay: 60002,taiMinusUTC: 36)],
      validUntilModifiedJulianDay: 60004, version: "synthetic",source: "test")
    #expect(try negative.dayLength(on: 60001) == 86399)
    #expect(throws: UTCError.invalidSecondOfDay) { try Epoch(utc: utc(60001,86399), using: negative) }
    let before = try Epoch(utc: utc(60001,86398.75), using: negative)
    let after = try Epoch(utc: utc(60002,0), using: negative)
    #expect(abs(try after.duration(since: before).seconds-0.25) < 1e-12)
    #expect(try before.utc(using: negative).modifiedJulianDay == 60001)
    #expect(try after.utc(using: negative).modifiedJulianDay == 60002)
  }
  @Test func injectableIERSDataAndValidation() throws {
    let parsed = try LeapSecondTable(iersData: "# version supplied by caller\n41317.0 1 1 1972 10\n41499.0 1 7 1972 11 # comment",
      validUntilModifiedJulianDay: 41683, version: "fixture",source: "test")
    #expect(try parsed.taiMinusUTC(on: 41498) == 10)
    #expect(try parsed.taiMinusUTC(on: 41499) == 11)
    #expect(throws: UTCError.invalidTable) {
      try LeapSecondTable(entries: [.init(effectiveModifiedJulianDay: 41317,taiMinusUTC: 10),
        .init(effectiveModifiedJulianDay: 41318,taiMinusUTC: 12)],validUntilModifiedJulianDay: 41400,version: "bad",source: "test")
    }
    #expect(throws: UTCError.malformedData) {
      try LeapSecondTable(iersData: "41317.5 1 1 1972 10",validUntilModifiedJulianDay: 41400,version: "bad",source: "test")
    }
  }
}

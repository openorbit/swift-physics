// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
@testable import CelMek

@Suite
struct TimeAndCalendarRegressionTests {
  @Test(arguments: ["", "2024", "2024-02", "2024-02-30", "2023-02-29",
                    "2024-01-00", "2024-13-01", "2024-01-01T", "2024-01-01T12",
                    "2024-01-01T12:30", "2024-01-01T12:30:00:01", "2024-01-01T24:00:00",
                    "2024-01-01T12:60:00", "2024-01-01T12:00:-1", "2024-01-01T12:00:nan",
                    "2024-01-01T12:00:inf", "2024-01-01T12:00:1e1", "2024-01-01T12:00:01.",
                    "2024-01-01T12:00:.5", "2024-01-01T12:00:01Z", "2024-01-01T12:00:01+01:00",
                    "2016-12-31T23:59:60", "2024-01-01T12:00:00T13:00:00"])
  func rejectsMalformedOrUnsupportedCalendarReading(input: String) {
    #expect(DateTime(isoDate: input) == nil)
  }

  @Test(arguments: ["2024-02-29T23:59:59.125", "2024-02-29 23:59:59.125"])
  func parsesFractionalCalendarReading(input: String) throws {
    let value = try #require(DateTime(isoDate: input))
    #expect(value.isValid())
    #expect(value.date.year == 2024)
    #expect(value.date.month == .february)
    #expect(value.date.day == 29)
    #expect(value.time.hh == 23)
    #expect(value.time.mm == 59)
    #expect(value.time.sec == 59.125)
    let dateOnly = try #require(DateTime(isoDate: "2024-02-29"))
    #expect(dateOnly.time.hh == 0 && dateOnly.time.mm == 0 && dateOnly.time.sec == 0)
  }

  @Test
  func validatesRawCalendarFields() {
    #expect(!Date(year: 2024, month: .january, day: 0).isValid())
    for second in [-1.0, 60, .nan, .infinity] {
      #expect(!Time(hh: 23, mm: 59, sec: second).isValid())
    }
  }

  @Test(arguments: [0.0, -1, 367, .nan, .infinity, -.infinity])
  func rejectsInvalidOrdinals(day: Double) {
    #expect(GregorianDate(year: 2024, dayInYear: day) == nil)
    #expect(JulianDate(year: 2024, dayInYear: day) == nil)
  }

  @Test(arguments: [1900, 2000, 2023, 2024])
  func ordinalAndMonthBoundaries(year: Int) throws {
    let commonLengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    var gregorianOrdinal = 1
    var julianOrdinal = 1
    for (index, commonLength) in commonLengths.enumerated() {
      let month = try #require(Month(rawValue: index + 1))
      let gregorianLength = commonLength + (index == 1 && isGregorianLeapYear(year) ? 1 : 0)
      let julianLength = commonLength + (index == 1 && isJulianLeapYear(year) ? 1 : 0)
      for day in 1...gregorianLength {
        let ordinal = Double(gregorianOrdinal) + 0.75
        let date = try #require(GregorianDate(year: year, dayInYear: ordinal))
        #expect(date.month == month)
        #expect(date.day == Double(day) + 0.75)
        #expect(date.dayOfYear() == ordinal)
        gregorianOrdinal += 1
      }
      for day in 1...julianLength {
        let ordinal = Double(julianOrdinal) + 0.75
        let date = try #require(JulianDate(year: year, dayInYear: ordinal))
        #expect(date.month == month)
        #expect(date.day == Double(day) + 0.75)
        #expect(date.dayOfYear() == ordinal)
        julianOrdinal += 1
      }
    }
    #expect(GregorianDate(year: year, dayInYear: Double(gregorianOrdinal)) == nil)
    #expect(JulianDate(year: year, dayInYear: Double(julianOrdinal)) == nil)
  }

  @Test
  func weekdayChangesAtMidnightRatherThanNoon() {
    #expect(2451544.499.weekday == .friday)
    #expect(2451544.5.weekday == .saturday)
    #expect(2451545.0.weekday == .saturday)
    #expect(2451545.499.weekday == .saturday)
    #expect(2451545.5.weekday == .sunday)
    #expect((-0.5).weekday == .monday)
    #expect((-1.5).weekday == .sunday)
    #expect((-2.5).weekday == .saturday)
  }

  @Test(arguments: [2415020.0, 2451545.0, 2488070.0])
  func timeConversionsRoundTripWithinJulianDateResolution(jd: Double) {
    #expect(abs(cm_tdbToTt(tdb: cm_ttToTdb(tt: jd)) - jd) <= jd.ulp)
    #expect(abs(cm_ttToTdb(tt: cm_tdbToTt(tdb: jd)) - jd) <= jd.ulp)
    #expect(abs(cm_ttToTdb(tt: jd) - jd) * 86400 < 0.002)
    #expect(cm_tdtToTdb(tdt: jd) == cm_ttToTdb(tt: jd))
    #expect(cm_tdbToTdt(tdb: jd) == cm_tdbToTt(tdb: jd))
    #expect(abs(cm_ttToTai(tt: cm_taiToTt(tai: jd)) - jd) <= jd.ulp)
    #expect(abs(cm_taiToGps(tai: cm_gpsToTai(gps: jd)) - jd) <= jd.ulp)
    #expect(abs(cm_tcgToTt(tcg: cm_ttToTcg(tt: jd)) - jd) <= 2 * jd.ulp)
    #expect(abs(cm_tcbToTdb(tcb: cm_tdbToTcb(tdb: jd)) - jd) <= 2 * jd.ulp)
  }

  @Test
  func tdbApproximationHasMillisecondRatherThanDayScale() {
    let tt = 2451545 + (90 - 357.53) / 0.9856003
    let offsetSeconds = (cm_ttToTdb(tt: tt) - tt) * 86400
    #expect(abs(offsetSeconds - 0.001658) <= tt.ulp * 86400)
  }

  @Test
  func deltaTTableEndpointsAndInterpolation() {
    #expect(deltaT(year: 1620) == 121)
    #expect(deltaT(year: 1621) == 116.5)
    #expect(abs(deltaT(year: 1997.5) - 62.65) < 1e-12)
    #expect(deltaT(year: 1998) == 63)
    #expect(deltaT(year: 1998.001).isFinite)
    for year in [948.0, 1619.999, 2100, 2100.001] {
      #expect(deltaT(year: year).isFinite)
    }
    for year in [Double.nan, .infinity, -.infinity] {
      #expect(deltaT(year: year).isNaN)
    }
  }
}

import Testing
import CelMek

@Suite struct CalendarDateTimeTests {
  @Test func knownDayNumbersAndNegativeYears() throws {
    for (year,month,day,mjd) in [(1858,11,17,0),(1970,1,1,40587),(2000,1,1,51544),
      (2017,1,1,57754),(0,1,1,-678941),(-1,1,1,-679306)] {
      let date = try GregorianDay(year: year,month: month,day: day)
      #expect(date.modifiedJulianDay == Int64(mjd))
      #expect(try GregorianDay(modifiedJulianDay: Int64(mjd)) == date)
    }
    for year in [-1_000_000,-400,-100,0,100,400,1900,2000,1_000_000] {
      for month in 1...12 {
        let date = try GregorianDay(year: year,month: month,day: 1)
        #expect(try GregorianDay(modifiedJulianDay: date.modifiedJulianDay) == date)
      }
    }
  }
  @Test func j2000ClockScalesAndNoonBoundary() throws {
    let tt = try Epoch.j2000.calendar(in: .tt)
    #expect(try tt.date == GregorianDay(year: 2000,month: 1,day: 1))
    #expect(tt.hour == 12 && tt.minute == 0 && tt.second == .zero)
    let tai = try Epoch.j2000.calendar(in: .tai)
    #expect(tai.hour == 11 && tai.minute == 59)
    #expect(abs(tai.second.seconds-27.816) < 2e-11)
    for scale in [CalendarTimeScale.tai,.tt,.tdb] {
      for seconds in [0.0,0.125,43199.125,43200,86399.999999] {
        let reading = try CalendarDateTime(date: tt.date,hour: Int(seconds)/3600,
          minute: Int(seconds)%3600/60,second: CelMek.Duration(seconds: seconds.truncatingRemainder(dividingBy: 60)),scale: scale)
        let epoch = try reading.epoch()
        let restored = try epoch.calendar(in: scale)
        #expect(restored.date == reading.date && restored.hour == reading.hour && restored.minute == reading.minute)
        #expect(abs(restored.second.seconds-reading.second.seconds) < 3e-11)
      }
    }
  }
  @Test func utcLeapLabelAndTableRequirements() throws {
    let reading = try CalendarDateTime(date: GregorianDay(year: 2016,month: 12,day: 31),hour: 23,minute: 59,
      second: CelMek.Duration(seconds: 60.25),scale: .utc)
    #expect(throws: CalendarConversionError.leapSecondTableRequired) { try reading.epoch() }
    let epoch = try Epoch(calendar: reading,using: .iersBulletinC72)
    #expect(try epoch.calendar(in: .utc,using: .iersBulletinC72) == reading)
    let next = try epoch.advanced(by: CelMek.Duration(seconds: 1)).calendar(in: .utc,using: .iersBulletinC72)
    #expect(try next.date == GregorianDay(year: 2017,month: 1,day: 1))
    #expect(next.hour == 0 && next.minute == 0 && abs(next.second.seconds-0.25) < 1e-12)
    let invalid = try CalendarDateTime(date: next.date,hour: 23,minute: 59,second: CelMek.Duration(seconds: 60),scale: .utc)
    #expect(throws: UTCError.invalidSecondOfDay) { try invalid.epoch(using: .iersBulletinC72) }
  }
  @Test func fullGregorianEraAndDistantEpochs() throws {
    let start = try GregorianDay(year: 2000,month: 3,day: 1).modifiedJulianDay
    for mjd in start..<(start+146097) {
      let date = try GregorianDay(modifiedJulianDay: mjd)
      #expect(date.modifiedJulianDay == mjd)
    }
    for year in [-1_000_000, -16000, 0, 16000, 1_000_000] {
      let reading = try CalendarDateTime(date: GregorianDay(year: year,month: 12,day: 31),
        hour: 23,minute: 59,second: CelMek.Duration(seconds: 59.125),scale: .tt)
      let restored = try reading.epoch().calendar(in: .tt)
      #expect(restored.date == reading.date && restored.hour == 23 && restored.minute == 59)
      #expect(abs(restored.second.seconds-59.125) < 3e-11)
    }
  }
  @Test func invalidAndRangeInputs() throws {
    #expect(throws: CalendarConversionError.invalidDate) { try GregorianDay(year: 1900,month: 2,day: 29) }
    #expect(throws: CalendarConversionError.invalidDate) { try GregorianDay(year: 2000,month: 13,day: 1) }
    #expect(throws: CalendarConversionError.outsideCalendarRange) { try GregorianDay(year: Int.max,month: 1,day: 1) }
    #expect(throws: CalendarConversionError.outsideCalendarRange) { try GregorianDay(modifiedJulianDay: Int64.min) }
    let date = try GregorianDay(year: 2000,month: 2,day: 29)
    #expect(throws: CalendarConversionError.invalidClock) { try CalendarDateTime(date: date,hour: 24,minute: 0,second: .zero,scale: .tt) }
    #expect(throws: CalendarConversionError.invalidClock) { try CalendarDateTime(date: date,hour: 23,minute: 59,second: CelMek.Duration(seconds: 60),scale: .tai) }
  }
}

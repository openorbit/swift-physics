// SPDX-License-Identifier: Apache-2.0
import Foundation

public enum CalendarConversionError: Error, Equatable, Sendable {
  case invalidDate, invalidClock, outsideCalendarRange, leapSecondTableRequired
}

/// Proleptic Gregorian day with astronomical year numbering (year 0 = 1 BCE).
/// Supported years: −1,000,000...1,000,000. No historical Julian/Gregorian cutover.
public struct GregorianDay: Hashable, Sendable {
  public let year: Int
  public let month: Int
  public let day: Int
  public init(year: Int, month: Int, day: Int) throws {
    guard (-1_000_000...1_000_000).contains(year) else { throw CalendarConversionError.outsideCalendarRange }
    guard (1...12).contains(month) else { throw CalendarConversionError.invalidDate }
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0)
    let lengths = [31,leap ? 29 : 28,31,30,31,30,31,31,30,31,30,31]
    guard (1...lengths[month-1]).contains(day) else { throw CalendarConversionError.invalidDate }
    self.year = year; self.month = month; self.day = day
  }
  /// Integer MJD at midnight, independent of the time scale of a clock reading.
  public var modifiedJulianDay: Int64 {
    // March-based 400-year era decomposition described by Howard Hinnant:
    // https://howardhinnant.github.io/date_algorithms.html
    let y = Int64(year) - (month <= 2 ? 1 : 0)
    let era = (y >= 0 ? y : y-399)/400
    let within = y-era*400
    let doy = (153*Int64(month > 2 ? month-3 : month+9)+2)/5+Int64(day)-1
    return era*146097 + within*365 + within/4 - within/100 + doy - 719468 + 40587
  }
  public init(modifiedJulianDay: Int64) throws {
    let first = try GregorianDay(year: -1_000_000,month: 1,day: 1).modifiedJulianDay
    let last = try GregorianDay(year: 1_000_000,month: 12,day: 31).modifiedJulianDay
    guard (first...last).contains(modifiedJulianDay) else { throw CalendarConversionError.outsideCalendarRange }
    let z = modifiedJulianDay-40587+719468
    let era = (z >= 0 ? z : z-146096)/146097
    let doe = z-era*146097
    let yoe = (doe-doe/1460+doe/36524-doe/146096)/365
    let doy = doe-(365*yoe+yoe/4-yoe/100)
    let mp = (5*doy+2)/153
    let d = doy-(153*mp+2)/5+1
    let m = mp+(mp < 10 ? 3 : -9)
    let y = yoe+era*400+(m <= 2 ? 1 : 0)
    try self.init(year: Int(y),month: Int(m),day: Int(d))
  }
}

/// Scale of a Gregorian clock reading. UTC requires a leap-second table;
/// continuous Julian dates continue to use TimeScale (TAI/TT/TDB).
public enum CalendarTimeScale: Hashable, Sendable { case tai, tt, tdb, utc }

public struct CalendarDateTime: Hashable, Sendable {
  public let date: GregorianDay
  public let hour: Int
  public let minute: Int
  public let second: Duration
  public let scale: CalendarTimeScale
  public init(date: GregorianDay, hour: Int, minute: Int, second: Duration, scale: CalendarTimeScale) throws {
    guard (0..<24).contains(hour), (0..<60).contains(minute), second >= .zero,
          second.wholeSeconds < 61 else { throw CalendarConversionError.invalidClock }
    if second.wholeSeconds == 60 {
      guard scale == .utc, hour == 23, minute == 59 else { throw CalendarConversionError.invalidClock }
    }
    self.date = date; self.hour = hour; self.minute = minute; self.second = second; self.scale = scale
  }
  public func epoch(using table: LeapSecondTable? = nil) throws -> Epoch {
    let seconds = try Duration(wholeSeconds: Int64(hour*3600+minute*60)).adding(second)
    if scale == .utc {
      guard let table else { throw CalendarConversionError.leapSecondTableRequired }
      return try Epoch(utc: UTCDate(modifiedJulianDay: date.modifiedJulianDay, secondOfDay: seconds), using: table)
    }
    let continuous: TimeScale = scale == .tai ? .tai : (scale == .tt ? .tt : .tdb)
    let jd = try JulianDateComponents(wholeDay: date.modifiedJulianDay+2400000,
      fractionalDay: 0.5+seconds.seconds/86400,scale: continuous)
    return try Epoch(julianDate: jd)
  }
}

extension Epoch {
  public init(calendar: CalendarDateTime, using table: LeapSecondTable? = nil) throws {
    self = try calendar.epoch(using: table)
  }
  public func calendar(in scale: CalendarTimeScale, using table: LeapSecondTable? = nil) throws -> CalendarDateTime {
    let mjd: Int64
    let seconds: Duration
    if scale == .utc {
      guard let table else { throw CalendarConversionError.leapSecondTableRequired }
      let date = try utc(using: table)
      mjd = date.modifiedJulianDay; seconds = date.secondOfDay
    } else {
      let continuous: TimeScale = scale == .tai ? .tai : (scale == .tt ? .tt : .tdb)
      let jd = try julianDate(in: continuous)
      // JD turns at noon; Gregorian clock days turn at midnight.
      let shifted = try JulianDateComponents(wholeDay: jd.wholeDay-2400001,
        fractionalDay: jd.fractionalDay+0.5,scale: continuous)
      mjd = shifted.wholeDay
      let raw = try Duration(seconds: shifted.fractionalDay*86400)
      // Avoid a rounded clock label of 24:00 for a date already chosen below midnight.
      seconds = raw.wholeSeconds < 86400 ? raw : try Duration(wholeSeconds: 86399,fractionalSecond: Double(1).nextDown)
    }
    let day = try GregorianDay(modifiedJulianDay: mjd)
    if seconds.wholeSeconds == 86400 {
      return try CalendarDateTime(date: day,hour: 23,minute: 59,
        second: Duration(wholeSeconds: 60,fractionalSecond: seconds.fractionalSecond),scale: scale)
    }
    return try CalendarDateTime(date: day,hour: Int(seconds.wholeSeconds/3600),
      minute: Int(seconds.wholeSeconds%3600/60),
      second: Duration(wholeSeconds: seconds.wholeSeconds%60,fractionalSecond: seconds.fractionalSecond),scale: scale)
  }
}

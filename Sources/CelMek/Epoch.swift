// SPDX-License-Identifier: Apache-2.0
import Foundation

public enum TimeArithmeticError: Error, Equatable, Sendable {
  case nonFinite
  case outOfRange
}

/// Continuous time scales supported by the initial epoch adapters.
/// UTC uses UTCDate and an explicit LeapSecondTable rather than a continuous JD.
public enum TimeScale: Hashable, Sendable {
  case tai
  case tt
  /// The library's low-order periodic approximation to TDB.
  case tdb
}

/// A finite SI-second interval, stored without combining its large and small parts.
/// Fractions are canonicalized to 0..<1, including for negative intervals.
/// The supported total range is ±(2^53 - 1) seconds, about 285 million years.
public struct Duration: Hashable, Comparable, Sendable {
  public static let maximumWholeSeconds: Int64 = 9_007_199_254_740_991
  public let wholeSeconds: Int64
  public let fractionalSecond: Double

  public init(wholeSeconds: Int64, fractionalSecond: Double = 0) throws {
    let parts = try normalizeTimeParts(whole: wholeSeconds, fraction: fractionalSecond)
    guard parts.whole >= -Self.maximumWholeSeconds,
          parts.whole <= Self.maximumWholeSeconds,
          parts.whole != Self.maximumWholeSeconds || parts.fraction == 0 else {
      throw TimeArithmeticError.outOfRange
    }
    self.wholeSeconds = parts.whole
    self.fractionalSecond = parts.fraction
  }

  public init(seconds: Double) throws {
    try self.init(wholeSeconds: 0, fractionalSecond: seconds)
  }

  /// A potentially lossy scalar conversion. Use the split fields for precision.
  public var seconds: Double { Double(wholeSeconds) + fractionalSecond }
  public static let zero = Duration(uncheckedWhole: 0, fraction: 0)

  public func adding(_ other: Duration) throws -> Duration {
    try Duration(wholeSeconds: wholeSeconds + other.wholeSeconds,
                 fractionalSecond: fractionalSecond + other.fractionalSecond)
  }

  public func subtracting(_ other: Duration) throws -> Duration {
    try Duration(wholeSeconds: wholeSeconds - other.wholeSeconds,
                 fractionalSecond: fractionalSecond - other.fractionalSecond)
  }

  public func negated() throws -> Duration {
    try Duration(wholeSeconds: -wholeSeconds, fractionalSecond: -fractionalSecond)
  }

  public static func < (lhs: Duration, rhs: Duration) -> Bool {
    lhs.wholeSeconds == rhs.wholeSeconds
      ? lhs.fractionalSecond < rhs.fractionalSecond
      : lhs.wholeSeconds < rhs.wholeSeconds
  }

  private init(uncheckedWhole: Int64, fraction: Double) {
    wholeSeconds = uncheckedWhole
    fractionalSecond = fraction
  }
}

/// A split Julian date in a specified scale, with a canonical fraction in 0..<1.
/// Julian date integer boundaries occur at noon, not civil midnight.
public struct JulianDateComponents: Hashable, Sendable {
  public let wholeDay: Int64
  public let fractionalDay: Double
  public let scale: TimeScale

  public init(wholeDay: Int64, fractionalDay: Double = 0, scale: TimeScale) throws {
    let parts = try normalizeTimeParts(whole: wholeDay, fraction: fractionalDay)
    self.wholeDay = parts.whole
    self.fractionalDay = parts.fraction
    self.scale = scale
  }

  public init(value: Double, scale: TimeScale) throws {
    try self.init(wholeDay: 0, fractionalDay: value, scale: scale)
  }

  /// A potentially lossy conversion suitable for legacy astronomical APIs.
  public var value: Double { Double(wholeDay) + fractionalDay }
}

/// An instant stored as elapsed continuous atomic seconds from J2000:
/// Julian date 2451545.0 TT (2000-01-01 12:00 TT).
/// Time scale belongs to a representation, not to the identity of the instant.
public struct Epoch: Hashable, Comparable, Sendable {
  public let offsetFromJ2000: Duration
  public static let j2000 = Epoch(offsetFromJ2000: .zero)

  public init(offsetFromJ2000: Duration) {
    self.offsetFromJ2000 = offsetFromJ2000
  }

  public init(julianDate: Double, scale: TimeScale) throws {
    try self.init(julianDate: JulianDateComponents(value: julianDate, scale: scale))
  }

  public init(julianDate: JulianDateComponents) throws {
    let (days, overflow) = julianDate.wholeDay.subtractingReportingOverflow(2_451_545)
    guard !overflow else { throw TimeArithmeticError.outOfRange }
    let (seconds, multiplyOverflow) = days.multipliedReportingOverflow(by: 86_400)
    guard !multiplyOverflow else { throw TimeArithmeticError.outOfRange }
    let represented = try Duration(wholeSeconds: seconds,
                                   fractionalSecond: julianDate.fractionalDay * 86_400)
    switch julianDate.scale {
    case .tt:
      offsetFromJ2000 = represented
    case .tai:
      offsetFromJ2000 = try represented.adding(Self.ttMinusTAI)
    case .tdb:
      var tt = represented
      for _ in 0..<4 {
        let correction = try Duration(seconds: Self.tdbMinusTT(at: tt))
        let next = try represented.subtracting(correction)
        if next == tt { break }
        tt = next
      }
      offsetFromJ2000 = tt
    }
  }

  public func advanced(by duration: Duration) throws -> Epoch {
    Epoch(offsetFromJ2000: try offsetFromJ2000.adding(duration))
  }

  public func duration(since other: Epoch) throws -> Duration {
    try offsetFromJ2000.subtracting(other.offsetFromJ2000)
  }

  public func julianDate(in scale: TimeScale) throws -> JulianDateComponents {
    let represented: Duration
    switch scale {
    case .tt: represented = offsetFromJ2000
    case .tai: represented = try offsetFromJ2000.subtracting(Self.ttMinusTAI)
    case .tdb:
      represented = try offsetFromJ2000.adding(Duration(seconds: Self.tdbMinusTT(at: offsetFromJ2000)))
    }
    // Floor division keeps the remainder positive for dates before J2000.
    var days = represented.wholeSeconds / 86_400
    var seconds = represented.wholeSeconds % 86_400
    if seconds < 0 { days -= 1; seconds += 86_400 }
    let fraction = (Double(seconds) + represented.fractionalSecond) / 86_400
    return try JulianDateComponents(wholeDay: 2_451_545 + days,
                                    fractionalDay: fraction, scale: scale)
  }

  public static func < (lhs: Epoch, rhs: Epoch) -> Bool {
    lhs.offsetFromJ2000 < rhs.offsetFromJ2000
  }

  // Derivative of a scale's coordinate seconds with respect to the elapsed
  // atomic seconds stored by Epoch. Used when adapting provider velocities.
  internal func coordinateTimeRate(in scale: TimeScale) -> Double {
    switch scale {
    case .tai, .tt: return 1
    case .tdb:
      let g = Self.tdbMeanAnomaly(at: offsetFromJ2000)
      let angularRate = 0.9856003 * Double.pi / (180 * 86_400)
      return 1 + (0.001658 * cos(g) + 0.000028 * cos(2 * g)) * angularRate
    }
  }

  // Second derivative of coordinate seconds with respect to atomic seconds.
  internal func coordinateTimeAcceleration(in scale: TimeScale) -> Double {
    switch scale {
    case .tai, .tt: return 0
    case .tdb:
      let g = Self.tdbMeanAnomaly(at: offsetFromJ2000)
      let angularRate = 0.9856003 * Double.pi / (180 * 86_400)
      return -(0.001658 * sin(g) + 0.000056 * sin(2*g)) * angularRate * angularRate
    }
  }

  // Both constants are exactly the intended decimal offsets at Double precision.
  private static let ttMinusTAI = try! Duration(wholeSeconds: 32, fractionalSecond: 0.184)

  /// Same approximation as the legacy TT/TDB conversion, evaluated from a split
  /// offset to avoid subtracting large Julian dates. Its physical error is not
  /// characterized; it is not a precision relativistic time-scale model.
  private static func tdbMinusTT(at tt: Duration) -> Double {
    let g = tdbMeanAnomaly(at: tt)
    return 0.001658 * sin(g) + 0.000014 * sin(2 * g)
  }

  private static func tdbMeanAnomaly(at tt: Duration) -> Double {
    let days = Double(tt.wholeSeconds) / 86_400 + tt.fractionalSecond / 86_400
    let degrees = (357.53 + 0.9856003 * days).truncatingRemainder(dividingBy: 360)
    return degrees * .pi / 180
  }
}

// Bound both input parts before any integer conversion or addition. At this
// bound their sum fits Int64, and the Double integer part converts exactly.
private func normalizeTimeParts(whole: Int64, fraction: Double) throws -> (whole: Int64, fraction: Double) {
  guard fraction.isFinite else { throw TimeArithmeticError.nonFinite }
  let limit = Duration.maximumWholeSeconds
  let integral = floor(fraction)
  guard whole >= -2 * limit, whole <= 2 * limit,
        integral >= -Double(limit), integral <= Double(limit) else {
    throw TimeArithmeticError.outOfRange
  }
  var normalizedWhole = whole + Int64(integral)
  var normalizedFraction = fraction - integral
  // Rounding a tiny negative fraction can produce exactly 1 rather than <1.
  if normalizedFraction == 1 { normalizedWhole += 1; normalizedFraction = 0 }
  if normalizedFraction == 0 { normalizedFraction = 0 } // Canonical positive zero.
  guard normalizedWhole >= -limit, normalizedWhole <= limit else {
    throw TimeArithmeticError.outOfRange
  }
  return (normalizedWhole, normalizedFraction)
}

// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import CelMek

// Normal imports exercise the public API. Qualify Duration to avoid Swift.Duration.
@Suite
struct EpochTests {
  @Test
  func normalizesSignedIntervals() throws {
    let negative = try CelMek.Duration(seconds: -0.25)
    #expect(negative.wholeSeconds == -1)
    #expect(negative.fractionalSecond == 0.75)
    let carry = try CelMek.Duration(wholeSeconds: 4, fractionalSecond: 1.25)
    #expect(carry.wholeSeconds == 5 && carry.fractionalSecond == 0.25)
    let borrow = try CelMek.Duration(wholeSeconds: 4, fractionalSecond: -0.25)
    #expect(borrow.wholeSeconds == 3 && borrow.fractionalSecond == 0.75)
    #expect(try negative.negated() == CelMek.Duration(seconds: 0.25))
    #expect(try negative.adding(carry) == CelMek.Duration(seconds: 5))
    #expect(negative < .zero)
    #expect(Set([negative, try CelMek.Duration(wholeSeconds: -1, fractionalSecond: 0.75)]).count == 1)
  }

  @Test(arguments: [Double.nan, .infinity, -.infinity])
  func rejectsNonfiniteInput(value: Double) {
    #expect(throws: TimeArithmeticError.nonFinite) { try CelMek.Duration(seconds: value) }
    #expect(throws: TimeArithmeticError.nonFinite) { try JulianDateComponents(value: value, scale: .tt) }
    #expect(throws: TimeArithmeticError.nonFinite) { try Epoch(julianDate: value, scale: .tt) }
  }

  @Test
  func reportsRangeAndArithmeticOverflow() throws {
    let limit = CelMek.Duration.maximumWholeSeconds
    let maximum = try CelMek.Duration(wholeSeconds: limit)
    let minimum = try CelMek.Duration(wholeSeconds: -limit)
    #expect(try maximum.negated() == minimum)
    #expect(try maximum.adding(minimum) == .zero)
    #expect(throws: TimeArithmeticError.outOfRange) { try maximum.adding(CelMek.Duration(seconds: 1)) }
    #expect(throws: TimeArithmeticError.outOfRange) { try minimum.subtracting(CelMek.Duration(seconds: 1)) }
    #expect(throws: TimeArithmeticError.outOfRange) { try maximum.subtracting(minimum) }
    #expect(throws: TimeArithmeticError.outOfRange) { try CelMek.Duration(wholeSeconds: .max) }
    #expect(throws: TimeArithmeticError.outOfRange) { try CelMek.Duration(wholeSeconds: .min) }
    #expect(throws: TimeArithmeticError.outOfRange) { try CelMek.Duration(seconds: .greatestFiniteMagnitude) }
    #expect(throws: TimeArithmeticError.outOfRange) {
      try Epoch(julianDate: JulianDateComponents(wholeDay: limit, scale: .tt))
    }
  }

  @Test
  func j2000IsAnInstantWithDifferentScaleRepresentations() throws {
    let tt = try Epoch.j2000.julianDate(in: .tt)
    #expect(tt.wholeDay == 2_451_545 && tt.fractionalDay == 0)
    #expect(tt.scale == .tt)
    let tai = try Epoch.j2000.julianDate(in: .tai)
    #expect(tai.wholeDay == 2_451_544)
    #expect(abs((1 - tai.fractionalDay) * 86400 - 32.184) < 1e-11)
    #expect(abs(try Epoch(julianDate: tai).duration(since: .j2000).seconds) < 2e-11)
    #expect(try Epoch(julianDate: 2451545, scale: .tt) == .j2000)
    let taiNoon = try Epoch(julianDate: 2451545, scale: .tai)
    #expect(abs(try taiNoon.duration(since: .j2000).seconds - 32.184) < 1e-14)
  }

  @Test(arguments: [-500_000_000_000 as Int64, 500_000_000_000])
  func nanosecondStepsSurviveFarFromReference(whole: Int64) throws {
    // Approximately 16,000 years from J2000, representative of DE441-scale spans.
    let initial = Epoch(offsetFromJ2000: try CelMek.Duration(wholeSeconds: whole, fractionalSecond: 0.25))
    let step = try CelMek.Duration(seconds: 1e-9)
    var advanced = initial
    for _ in 0..<1000 { advanced = try advanced.advanced(by: step) }
    #expect(advanced > initial)
    #expect(abs(try advanced.duration(since: initial).seconds - 1e-6) < 3e-14)
    let earlier = try initial.advanced(by: step.negated())
    #expect(earlier < initial)
    #expect(abs(try initial.duration(since: earlier).seconds - 1e-9) < 1e-16)
  }

  @Test
  func splitJulianDatePreservesFineStepsThatScalarDatesLose() throws {
    let initial = Epoch.j2000
    let next = try initial.advanced(by: CelMek.Duration(seconds: 1e-9))
    let firstDate = try initial.julianDate(in: .tt)
    let nextDate = try next.julianDate(in: .tt)
    #expect(firstDate.value == nextDate.value)
    #expect(nextDate.fractionalDay > firstDate.fractionalDay)
    #expect(abs(try Epoch(julianDate: nextDate).duration(since: initial).seconds - 1e-9) < 1e-20)
  }

  @Test(arguments: [TimeScale.tai, .tt, .tdb])
  func splitDateConversionsRoundTrip(scale: TimeScale) throws {
    for whole in [-500_000_000_000 as Int64, -1, 0, 1, 500_000_000_000] {
      let initial = Epoch(offsetFromJ2000: try CelMek.Duration(wholeSeconds: whole, fractionalSecond: 0.123456789))
      let date = try initial.julianDate(in: scale)
      #expect(date.fractionalDay >= 0 && date.fractionalDay < 1)
      let recovered = try Epoch(julianDate: date)
      #expect(abs(try recovered.duration(since: initial).seconds) < 3e-11)
    }
  }

  @Test
  func approximateTDBHasMillisecondCorrection() throws {
    let tt = try Epoch.j2000.julianDate(in: .tt)
    let tdb = try Epoch.j2000.julianDate(in: .tdb)
    let differenceDays = Double(tdb.wholeDay - tt.wholeDay) + tdb.fractionalDay - tt.fractionalDay
    #expect(abs(differenceDays * 86400) < 0.002)
    #expect(differenceDays != 0)
    // Same numeric JD in TDB is a different instant from that JD in TT.
    let tdbNoon = try Epoch(julianDate: 2451545, scale: .tdb)
    #expect(abs(try tdbNoon.duration(since: .j2000).seconds) > 0.00001)
  }

  @Test
  func negativeJulianFractionsNormalizeAcrossNoon() throws {
    let date = try JulianDateComponents(wholeDay: 2451545, fractionalDay: -0.25, scale: .tt)
    #expect(date.wholeDay == 2451544 && date.fractionalDay == 0.75)
    #expect(try Epoch(julianDate: date).offsetFromJ2000.seconds == -21600)
  }
}

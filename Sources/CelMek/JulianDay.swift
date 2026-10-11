//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2022 Mattias Holm
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

import Foundation

@frozen
public struct GregorianDate {
  public var year: Int
  public var month: Month
  public var day: Double
  
  public init(year: Int, month: Month, day: Double) {
    self.year = year
    self.month = month
    self.day = day
  }
  public init?(year: Int, dayInYear: Double) {
    let leap = isGregorianLeapYear(year)
    guard dayInYear.isFinite, dayInYear >= 1,
          dayInYear < (leap ? 367 : 366) else { return nil }
    let (month, day) = calendarMonthAndDay(dayInYear: dayInYear, leap: leap)
    self.year = year
    self.month = month
    self.day = day
  }
}


extension GregorianDate: CustomStringConvertible {
  public var description: String {
    return "\(year)-\(month)-\(day)"
  }
}


@frozen
public struct JulianDate {
  public var year: Int
  public var month: Month
  public var day: Double

  public init(year: Int, month: Month, day: Double) {
    self.year = year
    self.month = month
    self.day = day
  }

  public init?(year: Int, dayInYear: Double) {
    let leap = isJulianLeapYear(year)
    guard dayInYear.isFinite, dayInYear >= 1,
          dayInYear < (leap ? 367 : 366) else { return nil }
    let (month, day) = calendarMonthAndDay(dayInYear: dayInYear, leap: leap)
    self.year = year
    self.month = month
    self.day = day
  }
}
extension JulianDate: CustomStringConvertible {
  public var description: String {
    return "\(year)-\(month)-\(day)"
  }
}


extension GregorianDate {
  public func toJD() -> Double {
    // Meeus, eq 7.1
    let y = month.rawValue > 2 ? year : year - 1
    
    let m = month.rawValue > 2 ? month.rawValue : month.rawValue + 12
    
    let A = y / 100
    let B = 2 - A + A / 4
    
    let JD = Double(Int32(365.25 * Double(y + 4716))) + Double(Int32(30.6001 * Double(m + 1))) + self.day + Double(B) - 1524.5
    return JD
  }
}

extension JulianDate {
  public func toJD() -> Double {
    // Meeus, eq 7.1
    let y = month.rawValue > 2 ? year : year - 1
    
    let m = month.rawValue > 2 ? month.rawValue : month.rawValue + 12
    
    let JD = Double(Int32(365.25 * Double(y + 4716))) + Double(Int32(30.6001 * Double(m + 1))) + self.day - 1524.5
    return JD
  }
}

public func
isGregorianLeapYear(_ y: Int) -> Bool
{
  if (y % 400) == 0 {
    return true
  }
  if (y % 100) == 0 {
    return false
  }
  if (y % 4) == 0 {
    return true
  }
  return false
}
public func
isJulianLeapYear(_ y: Int) -> Bool
{
  if (y % 4) == 0 {
    return true
  }
  return false
}


extension GregorianDate {
  public func dayOfYear() -> Double {
    // Meeus, p65
    let K = isGregorianLeapYear(year) ? 1 : 2
    let N1 = 275 * month.rawValue / 9
    let N2 = K * ((month.rawValue + 9) / 12)
    let N = Double(N1 - N2 - 30) + day
    return N
  }
}

extension JulianDate {
  public func dayOfYear() -> Double {
    // Meeus, p65
    let K = isJulianLeapYear(year) ? 1 : 2
    let N1 = 275 * month.rawValue / 9
    let N2 = K * ((month.rawValue + 9) / 12)
    let N = Double(N1 - N2 - 30) + day
    return N
  }
}

extension Double {
  public var asMJD : Double { return self - MJD_JD_DIFFERENCE }
  public var asJD : Double { return self + MJD_JD_DIFFERENCE }
}

extension Double {
  public func toGregorian() -> GregorianDate {
    assert(self >= 0.0)
 
    let tmp = self + 0.5
    let Z = modf(tmp).0
    let F = modf(tmp).1
    
    let alpha = Int((Z - 1867216.25) / 36524.25)
    let A = Int(Z) + 1 + alpha - alpha / 4
    let B = Double(A + 1524)
    let C = floor((Double(B) - 122.1) / 365.25)
    let D = floor(365.25 * C)
    let E = floor((B - D) / 30.6001)
    
    let dayInMonth = B - D - floor(30.6001 * E) + F
    let month = Int(E) < 14 ? Int(E) - 1 : Int(E) - 13
    let year = Int(C) - (month > 2 ? 4716 : 4715)
    let date = GregorianDate(year: year, month: Month(rawValue: month)!, day: dayInMonth)
    return date
  }

  public func toJulian() -> JulianDate {
    let tmp = self + 0.5
    let Z = modf(tmp).0
    let F = modf(tmp).1
    let A = Int(Z)
    let B = Double(A + 1524)
    let C = floor((Double(B) - 122.1) / 365.25)
    let D = floor(365.25 * C)
    let E = floor((B - D) / 30.6001)
    
    let dayInMonth = B - D - floor(30.6001 * E) + F
    let month = Int(E) < 14 ? Int(E) - 1 : Int(E) - 13
    let year = Int(C) - (month > 2 ? 4716 : 4715)
    let date = JulianDate(year: year, month: Month(rawValue: month)!, day: dayInMonth)
    return date

  }
}

extension Double {
  public var weekday: Weekday {
    precondition(isFinite, "Weekday requires a finite Julian date")
    // Julian dates change their integer part at noon; civil days start at midnight.
    let remainder = fmod(floor(self + 0.5), 7)
    let weekdayNumber = Int(fmod(remainder + 8, 7))
    return Weekday(rawValue: weekdayNumber)!
  }
}

public func jd0(year: Int) -> Double {
  let Y = Double(year - 1)
  let A = floor(Y/100)
  let JD0 = floor(365.25 * Y) - A + floor(A/4) + 1721424.5
  return JD0
}

public func julianCenturiesFromJ2000(jd: Double) -> Double {
  // Meeus, Eq. 12.1
  let T =  (jd - 2451545) / 36525
  return T
}

// The ordinal range is checked before this helper is called. Keeping the
// fractional day intact avoids shifting late hours into the following month.
private func calendarMonthAndDay(dayInYear: Double, leap: Bool) -> (Month, Double) {
  let lengths = [31, leap ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
  var remaining = dayInYear
  for (index, length) in lengths.enumerated() {
    if remaining < Double(length + 1) {
      return (Month(rawValue: index + 1)!, remaining)
    }
    remaining -= Double(length)
  }
  preconditionFailure("Ordinal day outside the validated calendar year")
}

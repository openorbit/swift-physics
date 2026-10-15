//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2020 Mattias Holm
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

fileprivate let month_length : [Int] = [31, 28, 31, 30, 31, 30,
                                        31, 31, 30, 31, 30, 31];
fileprivate let month_length_leap : [Int]  = [31, 29, 31, 30, 31, 30,
                                              31, 31, 30, 31, 30, 31];

fileprivate let month_offset : [Int]  = [  0,  31,  59,  90, 120, 151,
                                           181, 212, 243, 273, 304, 334];
fileprivate let month_offset_leap : [Int]  = [  0,  31,  60,  91, 121, 152,
                                                182, 213, 244, 274, 305, 335];

// Ordered in work week order, but numbered to simplify calculations
public enum Weekday: Int, Sendable {
  case monday = 1
  case tuesday = 2
  case wednesday = 3
  case thursday = 4
  case friday = 5
  case saturday = 6
  case sunday = 0

  @available(*, deprecated, renamed: "monday")
  public static var Monday: Self { .monday }
  @available(*, deprecated, renamed: "tuesday")
  public static var Tuesday: Self { .tuesday }
  @available(*, deprecated, renamed: "wednesday")
  public static var Wednesday: Self { .wednesday }
  @available(*, deprecated, renamed: "thursday")
  public static var Thursday: Self { .thursday }
  @available(*, deprecated, renamed: "friday")
  public static var Friday: Self { .friday }
  @available(*, deprecated, renamed: "saturday")
  public static var Saturday: Self { .saturday }
  @available(*, deprecated, renamed: "sunday")
  public static var Sunday: Self { .sunday }
}
public enum Month: Int, Sendable {
  case january = 1
  case february = 2
  case march = 3
  case april = 4
  case may = 5
  case june = 6
  case july = 7
  case august = 8
  case september = 9
  case october = 10
  case november = 11
  case december = 12

  @available(*, deprecated, renamed: "january")
  public static var January: Self { .january }
  @available(*, deprecated, renamed: "february")
  public static var February: Self { .february }
  @available(*, deprecated, renamed: "march")
  public static var March: Self { .march }
  @available(*, deprecated, renamed: "april")
  public static var April: Self { .april }
  @available(*, deprecated, renamed: "may")
  public static var May: Self { .may }
  @available(*, deprecated, renamed: "june")
  public static var June: Self { .june }
  @available(*, deprecated, renamed: "july")
  public static var July: Self { .july }
  @available(*, deprecated, renamed: "august")
  public static var August: Self { .august }
  @available(*, deprecated, renamed: "september")
  public static var September: Self { .september }
  @available(*, deprecated, renamed: "october")
  public static var October: Self { .october }
  @available(*, deprecated, renamed: "november")
  public static var November: Self { .november }
  @available(*, deprecated, renamed: "december")
  public static var December: Self { .december }
}
extension Month: CustomStringConvertible {
  public var description: String {
    switch self {
    case .january:
      return "January"
    case .february:
      return "February"
    case .march:
      return "March"
    case .april:
      return "April"
    case .may:
      return "May"
    case .june:
      return "June"
    case .july:
      return "July"
    case .august:
      return "August"
    case .september:
      return "September"
    case .october:
      return "October"
    case .november:
      return "November"
    case .december:
      return "December"
    }
  }
}

public enum JewishMonth: Int, Sendable {
  case tishri = 1
  case heshvan = 2
  case kislev = 3
  case tevet = 4
  case shevat = 5
  case adar = 6
  case veadar = 7
  case nisan = 8
  case iyar = 9
  case sivan = 10
  case tammuz = 11
  case av = 12
  case elul = 13

  @available(*, deprecated, renamed: "tishri")
  public static var Tishri: Self { .tishri }
  @available(*, deprecated, renamed: "heshvan")
  public static var Heshvan: Self { .heshvan }
  @available(*, deprecated, renamed: "kislev")
  public static var Kislev: Self { .kislev }
  @available(*, deprecated, renamed: "tevet")
  public static var Tevet: Self { .tevet }
  @available(*, deprecated, renamed: "shevat")
  public static var Shevat: Self { .shevat }
  @available(*, deprecated, renamed: "adar")
  public static var Adar: Self { .adar }
  @available(*, deprecated, renamed: "veadar")
  public static var Veadar: Self { .veadar }
  @available(*, deprecated, renamed: "nisan")
  public static var Nisan: Self { .nisan }
  @available(*, deprecated, renamed: "iyar")
  public static var Iyar: Self { .iyar }
  @available(*, deprecated, renamed: "sivan")
  public static var Sivan: Self { .sivan }
  @available(*, deprecated, renamed: "tammuz")
  public static var Tammuz: Self { .tammuz }
  @available(*, deprecated, renamed: "av")
  public static var Av: Self { .av }
  @available(*, deprecated, renamed: "elul")
  public static var Elul: Self { .elul }
}

public enum MoslemMonth: Int, Sendable {
  case muharram = 1
  case safar = 2
  case rabiAlAwwal = 3
  case rabiAthThani = 4
  case jumadaLUla = 5
  case jumadaTTania = 6
  case rajab = 7
  case shaBan = 8
  case ramadan = 9
  case shawwal = 10
  case dhuLQaDa = 11
  case dhuLHijja = 12

  @available(*, deprecated, renamed: "muharram")
  public static var Muharram: Self { .muharram }
  @available(*, deprecated, renamed: "safar")
  public static var Safar: Self { .safar }
  @available(*, deprecated, renamed: "rabiAlAwwal")
  public static var RabiAlAwwal: Self { .rabiAlAwwal }
  @available(*, deprecated, renamed: "rabiAthThani")
  public static var RabiAthThani: Self { .rabiAthThani }
  @available(*, deprecated, renamed: "jumadaLUla")
  public static var JumadaLUla: Self { .jumadaLUla }
  @available(*, deprecated, renamed: "jumadaTTania")
  public static var JumadaTTania: Self { .jumadaTTania }
  @available(*, deprecated, renamed: "rajab")
  public static var Rajab: Self { .rajab }
  @available(*, deprecated, renamed: "shaBan")
  public static var ShaBan: Self { .shaBan }
  @available(*, deprecated, renamed: "ramadan")
  public static var Ramadan: Self { .ramadan }
  @available(*, deprecated, renamed: "shawwal")
  public static var Shawwal: Self { .shawwal }
  @available(*, deprecated, renamed: "dhuLQaDa")
  public static var DhuLQaDa: Self { .dhuLQaDa }
  @available(*, deprecated, renamed: "dhuLHijja")
  public static var DhuLHijja: Self { .dhuLHijja }
}

/// A Gregorian calendar date; this legacy type does not identify a time scale.
public struct Date {
  var year: Int32
  var month: Month
  var day: UInt8

  func isValid() -> Bool {
    let lengths = isLeapYear(year) ? month_length_leap : month_length
    return day >= 1 && Int(day) <= lengths[month.rawValue - 1]
  }
}

/// A clock reading without timezone or leap-second information.
public struct Time {
  var hh: UInt8
  var mm: UInt8
  var sec: Double

  func isValid() -> Bool {
    return hh < 24 && mm < 60 && sec.isFinite && sec >= 0 && sec < 60
  }
}

public struct DateTime {
  var date: Date
  var time: Time

  func isValid() -> Bool {
    return date.isValid() && time.isValid()
  }

  /// Parses YYYY-MM-DD with optional HH:mm:ss[.fraction], separated by T or space.
  /// This is a calendar reading, not a UTC instant. Timezone suffixes and leap
  /// seconds are rejected until a time-scale-aware API can validate them.
  init?(isoDate str: String) {
    let components = str.components(separatedBy: str.contains("T") ? "T" : " ")
    guard components.count == 1 || components.count == 2 else { return nil }
    let dateComponents = components[0].components(separatedBy: "-")
    guard dateComponents.count == 3,
          dateComponents.allSatisfy({ !$0.isEmpty && $0.allSatisfy({ $0.isASCII && $0.isNumber }) }),
          let year = Int32(dateComponents[0]),
          let monthNumber = Int(dateComponents[1]),
          let month = Month(rawValue: monthNumber),
          let day = UInt8(dateComponents[2]) else { return nil }
    date = Date(year: year, month: month, day: day)
    guard date.isValid() else { return nil }

    if components.count == 2 {
      let clock = components[1].components(separatedBy: ":")
      guard clock.count == 3,
            clock[0].allSatisfy({ $0.isASCII && $0.isNumber }),
            clock[1].allSatisfy({ $0.isASCII && $0.isNumber }),
            let hour = UInt8(clock[0]), let minute = UInt8(clock[1]) else { return nil }
      let seconds = clock[2].components(separatedBy: ".")
      guard seconds.count == 1 || seconds.count == 2,
            seconds.allSatisfy({ !$0.isEmpty && $0.allSatisfy({ $0.isASCII && $0.isNumber }) }),
            let second = Double(clock[2]) else { return nil }
      time = Time(hh: hour, mm: minute, sec: second)
    } else {
      time = Time(hh: 0, mm: 0, sec: 0)
    }
    guard time.isValid() else { return nil }
  }
}

func
isLeapYear(_ y: Int32) -> Bool
{
  if (y % 400) == 0 {return true;}
  if (y % 100) == 0 {return false;}
  if (y % 4) == 0 {return true;}
  return false;
}


func
cm_julianCenturiesFromEpoch(jd:  Double, epoch: Double) -> Double
{
  return (jd-epoch)/36525.0;
}




func
cm_tcbToTdb(tcb: Double) -> Double
{
  return tcb - CM_LB * (tcb - CM_TCB_0) + CM_TDB0 / 86400.0;
}

func
cm_tdbToTcb(tdb: Double) -> Double
{
  return (tdb - CM_TDB0 / 86400.0 - CM_LB * CM_TCB_0) / (1.0-CM_LB);
}

/// Low-order periodic approximation to TDB - TT, in Julian days.
/// Both arguments/results of the conversion routines are Julian dates.
/// This is not a replacement for a precision relativistic time model.
private func tdbMinusTT(atTT tt: Double) -> Double {
  let g = 357.53.deg + 0.9856003.deg * (tt - 2451545.0)
  return (0.001658 * sin(g) + 0.000014 * sin(2 * g)) / SEC_PER_JD
}

func cm_ttToTdb(tt: Double) -> Double {
  return tt + tdbMinusTT(atTT: tt)
}

func
cm_tdbToTt(tdb: Double) -> Double
{
  // Solve the same forward approximation rather than using a different model.
  var tt = tdb
  for _ in 0..<3 {
    let next = tdb - tdbMinusTT(atTT: tt)
    if next == tt { break }
    tt = next
  }
  return tt
}

func
cm_gpsToTai(gps: Double)->Double
{
  return gps + 19.0/SEC_PER_JD;
}

func
cm_taiToGps(tai: Double)->Double
{
  return tai - 19.0/SEC_PER_JD;
}

func
cm_taiToTt(tai: Double)->Double
{
  return tai + 32.184/SEC_PER_JD;
}
func
cm_ttToTai(tt: Double)->Double
{
  return tt - 32.184/SEC_PER_JD;
}

func
cm_tcgToTt(tcg: Double)->Double
{
  let EJD = 2443144.5003725;
  return EJD + (tcg - EJD) * (1.0 - CM_LG);
}

func
cm_ttToTcg(tt: Double)->Double
{
  let EJD = 2443144.5003725;
  return (tt - EJD)/(1.0 - CM_LG) + EJD;
}

/// Legacy TDT spelling for TT; delegates to the shared approximation.
func cm_tdtToTdb(tdt: Double) -> Double {
  return cm_ttToTdb(tt: tdt)
}

/// Legacy TDT spelling for TT; delegates to the shared approximation inverse.
func cm_tdbToTdt(tdb: Double) -> Double {
  return cm_tdbToTt(tdb: tdb)
}

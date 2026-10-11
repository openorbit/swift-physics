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

#if os(Linux)
import Glibc
#else
import Darwin
#endif

fileprivate let month_length : [Int] = [31, 28, 31, 30, 31, 30,
                                        31, 31, 30, 31, 30, 31];
fileprivate let month_length_leap : [Int]  = [31, 29, 31, 30, 31, 30,
                                              31, 31, 30, 31, 30, 31];

fileprivate let month_offset : [Int]  = [  0,  31,  59,  90, 120, 151,
                                           181, 212, 243, 273, 304, 334];
fileprivate let month_offset_leap : [Int]  = [  0,  31,  60,  91, 121, 152,
                                                182, 213, 244, 274, 305, 335];

// Ordered in work week order, but numbered to simplify calculations
public enum Weekday : Int {
  case Monday = 1
  case Tuesday = 2
  case Wednesday = 3
  case Thursday = 4
  case Friday = 5
  case Saturday = 6
  case Sunday = 0

}
public enum Month : Int {
  case January = 1;
  case February = 2;
  case March = 3;
  case April = 4;
  case May = 5;
  case June = 6;
  case July = 7;
  case August = 8;
  case September = 9;
  case October = 10;
  case November = 11;
  case December = 12;
}
extension Month: CustomStringConvertible {
  public var description: String {
    switch self {
    case .January:
      return "January"
    case .February:
      return "February"
    case .March:
      return "March"
    case .April:
      return "April"
    case .May:
      return "May"
    case .June:
      return "June"
    case .July:
      return "July"
    case .August:
      return "August"
    case .September:
      return "September"
    case .October:
      return "October"
    case .November:
      return "November"
    case .December:
      return "December"
    }
  }
}

public enum JewishMonth : Int {
  case Tishri = 1;
  case Heshvan = 2;
  case Kislev = 3;
  case Tevet = 4;
  case Shevat = 5;
  case Adar = 6;
  case Veadar = 7;
  case Nisan = 8;
  case Iyar = 9;
  case Sivan = 10;
  case Tammuz = 11;
  case Av = 12;
  case Elul = 13;
}

public enum MoslemMonth : Int {
  case Muharram = 1;
  case Safar = 2;
  case RabiAlAwwal = 3;
  case RabiAthThani = 4;
  case JumadaLUla = 5;
  case JumadaTTania = 6;
  case Rajab = 7;
  case ShaBan = 8;
  case Ramadan = 9;
  case Shawwal = 10;
  case DhuLQaDa = 11;
  case DhuLHijja = 12;
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

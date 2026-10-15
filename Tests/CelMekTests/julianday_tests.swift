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

import Testing
@testable import CelMek;

@Suite
struct jd_tests {
  @Test
  func testGregorianToJD() {
    // Meeus, example 7.a
    let date = GregorianDate(year: 1957, month: .october, day: 4.81)
    #expect(date.toJD() == 2436116.31)
  }
  
  @Test
  func testJulianToJD() {
    // Meeus, example 7.b
    let date = JulianDate(year: 333, month: .january, day: 27.5)
    #expect(date.toJD() == 1842713.0)
  }
  
  @Test
  func testGregorianConversions() {
    // Meeus, p 62
    let dates = [
      (GregorianDate(year: 2000, month: .january, day: 1.5), 2451545.0),
      (GregorianDate(year: 1999, month: .january, day: 1.0), 2451179.5),
      (GregorianDate(year: 1987, month: .january, day: 27.0), 2446822.5),
      (GregorianDate(year: 1987, month: .june, day: 19.5), 2446966.0),
      (GregorianDate(year: 1988, month: .january, day: 27.0), 2447187.5),
      (GregorianDate(year: 1988, month: .june, day: 19.5), 2447332.0),
      (GregorianDate(year: 1900, month: .january, day: 1.0), 2415020.5),
      (GregorianDate(year: 1600, month: .january, day: 1.0), 2305447.5),
      (GregorianDate(year: 1600, month: .december, day: 31.0), 2305812.5),
    ]
    
    for (date, jd) in dates {
      #expect(date.toJD() == jd)
      
      let gdate = jd.toGregorian()
      #expect(gdate.year == date.year)
      #expect(gdate.month == date.month)
      #expect(gdate.day == date.day)
    }
  }
  
  @Test
  func testJulianConversion() {
    // Meeus, p 62
    let dates = [
      (JulianDate(year: 837, month: .april, day: 10.3), 2026871.8),
      (JulianDate(year: -123, month: .december, day: 31.0), 1676496.5),
      (JulianDate(year: -122, month: .january, day: 1.0), 1676497.5),
      (JulianDate(year: -1000, month: .july, day: 12.5), 1356001.0),
      (JulianDate(year: -1000, month: .february, day: 29.0), 1355866.5),
      (JulianDate(year: -1001, month: .august, day: 17.9), 1355671.4),
      (JulianDate(year: -4712, month: .january, day: 1.5), 0.0),
    ]
    
    for (date, jd) in dates {
      #expect(date.toJD() == jd)
    }
  }
  
  @Test
  func testMJDConversion() {
    // Meeus, p 63
    let date = GregorianDate(year: 1858, month: .november, day: 17.0)
    #expect(date.toJD().asMJD == 0.0)
  }
  
  @Test
  func testWeekday() {
    // Meeus, p65
    let date = GregorianDate(year: 1954, month: .june, day: 30.1)
    #expect(date.toJD().weekday == .wednesday)
  }
  
  @Test
  func testDayOfYear() {
    // Meeus, p65, examples 7.f and 7.g
    let date_7f = GregorianDate(year: 1978, month: .november, day: 14)
    #expect(date_7f.dayOfYear() == 318)
    
    let date_7g = GregorianDate(year: 1988, month: .april, day: 22)
    #expect(date_7g.dayOfYear() == 113)
  }

  @Test
  func testGregorianDateFromDayInYear() throws {
    // Meeus, p65, examples 7.f and 7.g in reverse
    let date_7f = try #require(GregorianDate(year: 1978, dayInYear: 318))
    #expect(date_7f.year == 1978)
    #expect(date_7f.month == .november)
    #expect(date_7f.day == 14)
    
    let date_7g = try #require(GregorianDate(year: 1988, dayInYear: 113))
    #expect(date_7g.year == 1988)
    #expect(date_7g.month == .april)
    #expect(date_7g.day == 22)
  }
  
  @Test
  func testJulianYears() {
    let date = JulianDate(year: 333, month: .february, day: 6.33333333333)
    #expect(abs((date.years()) - (333.1)) <= 0.01)
  }
  
  @Test
  func testGregorianYears() {
    let date = GregorianDate(year: 1999, month: .february, day: 6.33333333333)
    #expect(abs((date.years()) - (1999.1)) <= 0.01)
  }

  @Test
  func testJD0() {
    let jd0_1999 = jd0(year: 1999)
    let greg0_1999 = GregorianDate(year: 1999, month: .january, day: 0.0)
    let jd0_expect = greg0_1999.toJD()
    #expect(jd0_1999 == jd0_expect)
  }
  
  @Test
  func testGregorianToJulianConversion() {
    let gregDate = GregorianDate(year: 2023, month: .january, day: 4)
    let jd = gregDate.toJD()
    let julianDate = jd.toJulian()
    
    #expect(julianDate.year == 2022)
    #expect(julianDate.month == .december)
    #expect(julianDate.day == 22)

  }
}


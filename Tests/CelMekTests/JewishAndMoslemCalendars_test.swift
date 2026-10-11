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
struct JewishCalendar_tests {
  @Test
  func testJewishEaster() {
    // Meeus, example 9.a
    let date = jewishEaster(gregorianYear: 1990)
    
    #expect(date.year == 1990)
    #expect(date.month == .April)
    #expect(date.day == 10.0)
  }
  
  @Test
  func testJewishCommonYear() {
    // Meeu,s example 9.a
    let isCommon = jewishYearIsCommon(jewishYear: 5751)
    #expect(isCommon == true)
  }
  
  @Test
  func testJewishYearStart() {
    // Meeu,s example 9.a
    let date = jewishYearStart(jewishYear: 5751)
    #expect(date.year == 1990)
    #expect(date.month == .September)
    #expect(date.day == 20.0)
  }
  
  @Test
  func testDaysOfJewishYear() {
    // Meeu,s example 9.a
    let days = daysOfJewishYear(jewishYear: 5751)
    #expect(days == 354)
  }
  
  @Test
  func testMoslemToJulian() {
    // Meeu,s example 9.b
    let date = JulianDate(
      moslemDate: MoslemDate(year: 1421,
                             month: .Muharram,
                             day: 1))
    
    #expect(date.year == 2000)
    #expect(date.dayOfYear() == 84.0)

    let jd = date.toJD()
    let gregDate = jd.toGregorian()
    #expect(gregDate.year == 2000)
    #expect(gregDate.month == .April)
    #expect(gregDate.day == 6.0)
  }
  
  @Test
  func testJulianToMoslem() throws {
    // Meeu,s example 9.c
    let date = MoslemDate(
      julianDate: JulianDate(year: 1991,
                             month: .July,
                             day: 31))

    #expect(date.year == 1412)
    #expect(date.month == .Safar)
    #expect(date.day == 2)

    let julianDate = try #require(JulianDate(year: 2000, dayInYear: 84.0))
    #expect(julianDate.year == 2000)
    #expect(julianDate.month == .March)
    #expect(julianDate.day == 24)
    
    let date2 = MoslemDate(julianDate: julianDate)
    #expect(date2.year == 1421)
    #expect(date2.month == .Muharram)
    #expect(date2.day == 1)
    // Currently incorrect
    // Gregorian date 2023-01-04
    // Julian date 2022-12-22
    // Islamic date: 1444-06-11 According to Wiki
    let date3 = MoslemDate(julianDate: JulianDate(year: 2022, month: .December, day: 22))
    #expect(date3.year == 1444)
    #expect(date3.month == .JumadaTTania)
    #expect(date3.day == 11)

  }
}

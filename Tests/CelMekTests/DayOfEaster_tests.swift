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
struct dayOfEaster_tests {
  @Test
  func testGregorianDateOfEaster() {
    // Meeus, p68
    let date_0 = gregorianDateOfEaster(year: 1991)
    #expect(date_0.year == 1991)
    #expect(date_0.month == .march)
    #expect(date_0.day == 31)

    let date_1 = gregorianDateOfEaster(year: 1992)
    #expect(date_1.year == 1992)
    #expect(date_1.month == .april)
    #expect(date_1.day == 19)

    let date_2 = gregorianDateOfEaster(year: 1993)
    #expect(date_2.year == 1993)
    #expect(date_2.month == .april)
    #expect(date_2.day == 11)

    let date_3 = gregorianDateOfEaster(year: 1954)
    #expect(date_3.year == 1954)
    #expect(date_3.month == .april)
    #expect(date_3.day == 18)

    let date_4 = gregorianDateOfEaster(year: 2000)
    #expect(date_4.year == 2000)
    #expect(date_4.month == .april)
    #expect(date_4.day == 23)

    let date_5 = gregorianDateOfEaster(year: 1818)
    #expect(date_5.year == 1818)
    #expect(date_5.month == .march)
    #expect(date_5.day == 22)
  }
  
  @Test
  func testJulianDateOfEaster() {
    let date_0 = julianDateOfEaster(year: 179)
    #expect(date_0.year == 179)
    #expect(date_0.month == .april)
    #expect(date_0.day == 12)

    let date_1 = julianDateOfEaster(year: 711)
    #expect(date_1.year == 711)
    #expect(date_1.month == .april)
    #expect(date_1.day == 12)

    let date_2 = julianDateOfEaster(year: 1243)
    #expect(date_2.year == 1243)
    #expect(date_2.month == .april)
    #expect(date_2.day == 12)
  }
}


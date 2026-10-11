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
struct time_tests {
  @Test
  func testDateTime() throws {
    let dt0 = try #require(DateTime(isoDate: "2020-12-12"))
    #expect(dt0.date.year == 2020)
    #expect(dt0.date.month == Month.December)
    #expect(dt0.date.day == 12)
    
    let dt1 = DateTime(isoDate: "2020-0-12")
    #expect(dt1 == nil)
    
  }

  @Test
  func testIsLeap() {
    #expect(isLeapYear(2000))
    #expect(isLeapYear(1996))
    #expect(isLeapYear(2004))
    #expect(!(isLeapYear(1900)))
    #expect(!(isLeapYear(2001)))
    #expect(!(isLeapYear(1999)))
  }
}



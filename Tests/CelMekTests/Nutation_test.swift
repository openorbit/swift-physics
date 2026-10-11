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
import Math

@Suite
struct Nutation_tests {
  @Test
  func testFastNutation() {
    // Meeus, example 22.a
    let nut = fastNutation(jd: 2446895.5)
    let e0 = fastMeanObliquityOfTheEcliptic(jd: 2446895.5)
    #expect(abs(nut.nutationInLongitude - (-3.788.arcsec)) <= 0.35.arcsec)
    #expect(abs(nut.nutationInObliquity - 9.443.arcsec) <= 0.1.arcsec)
    #expect(abs((e0) - (AngleOfArc(degrees: 23, minutes: 26, seconds: 27.407).toRad())) <= 0.00001)

    let e = fastTrueObliquityOfTheEcliptic(jd: 2446895.5)
    #expect(abs(e - AngleOfArc(degrees: 23, minutes: 26, seconds: 36.850).toRad()) <= 0.1.arcsec)
   }

  @Test
  func testNutation() {
    // Meeus, example 22.a
    let nut = nutation(jd: 2446895.5)
    let e0 = meanObliquityOfTheEcliptic(jd: 2446895.5)
    #expect(abs((nut.nutationInLongitude) - (-3.788.arcsec)) <= 0.000000001)
    #expect(abs((nut.nutationInObliquity) - (9.443.arcsec)) <= 0.00000001)
    #expect(abs((e0) - (AngleOfArc(degrees: 23, minutes: 26, seconds: 27.407).toRad())) <= 0.00000001)
    let e = trueObliquityOfTheEcliptic(jd: 2446895.5)
    #expect(abs((e) - (AngleOfArc(degrees: 23, minutes: 26, seconds: 36.850).toRad())) <= 0.00000001)
  }
}

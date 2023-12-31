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

import XCTest
@testable import CelMek;
import Math

final class Nutation_tests : XCTestCase {
  func testFastNutation() {
    // Meeus, example 22.a
    let nut = fastNutation(jd: 2446895.5)
    let e0 = fastMeanObliquityOfTheEcliptic(jd: 2446895.5)
    XCTAssertEqual(nut.nutationInLongitude, -3.788.arcsec, accuracy: 0.00001)
    XCTAssertEqual(nut.nutationInObliquity, 9.443.arcsec, accuracy: 0.00001)
    XCTAssertEqual(e0, AngleOfArc(degrees: 23, minutes: 26, seconds: 27.407).toRad(),
                   accuracy: 0.00001)

    let e = fastTrueObliquityOfTheEcliptic(jd: 2446895.5)
    XCTAssertEqual(e, AngleOfArc(degrees: 23, minutes: 26, seconds: 36.850).toRad(),
                   accuracy: 0.00001)
   }

  func testNutation() {
    // Meeus, example 22.a
    let nut = nutation(jd: 2446895.5)
    let e0 = meanObliquityOfTheEcliptic(jd: 2446895.5)
    XCTAssertEqual(nut.nutationInLongitude, -3.788.arcsec, accuracy: 0.000000001)
    XCTAssertEqual(nut.nutationInObliquity, 9.443.arcsec, accuracy: 0.00000001)
    XCTAssertEqual(e0, AngleOfArc(degrees: 23, minutes: 26, seconds: 27.407).toRad(),
                   accuracy: 0.00000001)
    let e = trueObliquityOfTheEcliptic(jd: 2446895.5)
    XCTAssertEqual(e, AngleOfArc(degrees: 23, minutes: 26, seconds: 36.850).toRad(),
                   accuracy: 0.00000001)
  }
}

//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2023 Mattias Holm
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
struct angularSeparation_tests {
  @Test
  func testAngularSeparation() {
    // Meeus, Example 17.a
    let p0 = EquatorialCoordinate(rightAscension: 213.9154.deg,
                                  declination: 19.1825.deg)
    let p1 = EquatorialCoordinate(rightAscension: 201.2983.deg,
                                  declination: -11.1614.deg)
    let d = angularSeparation(p0, p1)

    #expect(abs((d) - (32.7930.deg)) <= 0.000005)
  }
}


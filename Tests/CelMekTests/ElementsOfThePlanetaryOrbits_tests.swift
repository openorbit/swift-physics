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
struct PlanetaryOrbits_tests {
  @Test
  func testMercuryMean() {
    // Meeus, Example 31.a
    let elements = elementsOfMercuryForMeanEquinox(jd: 2475460.5)
    #expect(abs((elements.meanLongitude) - (203.494701)) <= 0.0000005)
    #expect(elements.semimajorAxis == 0.387098310)
    #expect(abs((elements.eccentricity) - (0.20564510)) <= 0.0000005)
    #expect(abs((elements.inclinationOnEclipticPlane) - (7.006171)) <= 0.0000005)
    #expect(abs((elements.longitudeOfAscendingNode) - (49.107650)) <= 0.0000005)
    #expect(abs((elements.longitudeOfPerihelion) - (78.475382)) <= 0.0000005)
  }
}


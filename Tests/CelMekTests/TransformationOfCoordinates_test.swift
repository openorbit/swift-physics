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
struct TransformationOfCoordinates_tests {
  @Test
  func testEclipticalFromEquatorial() {
    // Meeus, p 95, example 13.a
    let eq = EquatorialCoordinate(rightAscension: 116.328942.deg, declination: 28.026183.deg)
    let ecl = EclipticCoordinate(equatorial: eq, obliquityOfEcliptic: J2000_OBLIQUITY_OF_ECLIPTIC)
    
    let 𝜆 = ecl.longitude
    let 𝛽 = ecl.latitude

    #expect(abs((𝜆) - (113.215630.deg)) <= 0.00000005)
    #expect(abs((𝛽) - (6.684170.deg)) <= 0.00000005)
  }
  
  @Test
  func testEquatorialFromEcliptical() {
    // Meeus, p 95, example 13.a
    let ecl = EclipticCoordinate(latitude: 6.684170.deg, longitude: 113.215630.deg)
    let eq = EquatorialCoordinate(ecliptic: ecl, obliquityOfEcliptic: J2000_OBLIQUITY_OF_ECLIPTIC)

    #expect(abs((eq.rightAscension) - (116.328942.deg)) <= 0.00000005)
    #expect(abs((eq.declination) - (28.026183.deg)) <= 0.00000005)
  }

  // TODO: Add example 13.b when nutations have been implemented
  @Test
  func testGalacticFromEquatorial() {
    // Meeus, p 96 exercise
    let eq = EquatorialCoordinate(
      rightAscension: HourAngle(hours: 17, minutes: 48, seconds: 59.74).radians,
      declination: AngleOfArc(degrees: -14, minutes: 43, seconds: 8.2).radians)
    let gal = GalacticCoordinate(equatorial: eq)
    
    
    #expect(abs((gal.latitude) - (6.0463.deg)) <= 0.00005)
    #expect(abs((gal.longitude) - (12.9593.deg)) <= 0.00005)
  }


}


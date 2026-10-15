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
struct EarthGlobe_tests {
  @Test
  func testPComp() {
    // Meeus, example 11.a
    let (psin, pcos) = pcomp(H: 1706, geographicLat: 33.356111.deg)
    #expect(abs((psin) - (0.546861)) <= 0.0000005)
    #expect(abs((pcos) - (0.836339)) <= 0.0000005)
  }
  
  @Test
  func testRadiusOfCurvatureOfMeridian() {
    // Meeus, p84
    let eqRad = radiusOfCurvatureOfMeridian(geographicalLat: 0.0)
    #expect(abs((eqRad) - (6335.44)) <= 0.005)
    let polRad = radiusOfCurvatureOfMeridian(geographicalLat: 90.0.deg)
    #expect(abs((polRad) - (6399.60)) <= 0.005)
    
    // Meeus, example 11.b
    let chicagoRad = radiusOfCurvatureOfMeridian(geographicalLat: 42.0.deg)
    #expect(abs((chicagoRad) - (6364.033)) <= 0.0005)
  }
  
  @Test
  func testParallelOfLatitudeRadius() {
    // Meeus, p84
    let eqRad = parallelOfLatitudeRadius(geographicalLat: 0.0)
    #expect(abs((eqRad) - (EARTH_EQUATORIAL_RADIUS_KM)) <= 0.000005)
    let polRad = parallelOfLatitudeRadius(geographicalLat: 90.0.deg)
    #expect(abs((polRad) - (0.0)) <= 0.00001)
    
    // Meeus, example 11.b
    let chicagoRad = parallelOfLatitudeRadius(geographicalLat: 42.0.deg)
    #expect(abs((chicagoRad) - (4747.001)) <= 0.0005)
  }
  
  @Test
  func testGeodesicDistanceFast() {
    // Meeus, Example 11.c
    let paris = GeographicCoordinate(
      latitude: AngleOfArc(degrees: 48, minutes: 50, seconds: 11).radians,
      longitude: AngleOfArc(degrees: -2, minutes: 20, seconds: 14).radians)
    
    let washington = GeographicCoordinate(
      latitude: AngleOfArc(degrees: 38, minutes: 55, seconds: 17).radians,
      longitude: AngleOfArc(degrees: 77, minutes: 3, seconds: 56).radians)
    let s = geodesicDistanceFast(p0: paris, p1: washington)
    #expect(abs((s) - (6166)) <= 0.5)
  }

  @Test
  func testGeodesicDistance() {
    // Meeus, Example 11.c
    let paris = GeographicCoordinate(
      latitude: AngleOfArc(degrees: 48, minutes: 50, seconds: 11).radians,
      longitude: AngleOfArc(degrees: -2, minutes: 20, seconds: 14).radians)
    
    let washington = GeographicCoordinate(
      latitude: AngleOfArc(degrees: 38, minutes: 55, seconds: 17).radians,
      longitude: AngleOfArc(degrees: 77, minutes: 3, seconds: 56).radians)
    let s = geodesicDistance(p0: paris, p1: washington)
    #expect(abs((s) - (6181.63)) <= 0.005)
  }
}

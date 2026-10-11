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

import Foundation
import Math
// Meeus, chapter 15
// L: geographic longitude
// 𝜑: geographic latitude
// dT: difference between TD - UT in seconds
// h0: standard altitude, geometric altitude of center of body
let MEAN_ATMOSPHERIC_REFRACTION_ARC = AngleOfArc(degrees: 0, minutes: 34, seconds: 0)
let MEAN_ATMOSPHERIC_REFRACTION = MEAN_ATMOSPHERIC_REFRACTION_ARC.toRad()
let MEAN_ATMOSPHERIC_REFRACTION_SUN_ARC = AngleOfArc(degrees: 0, minutes: 50, seconds: 0)
let MEAN_ATMOSPHERIC_REFRACTION_SUN = MEAN_ATMOSPHERIC_REFRACTION_SUN_ARC.toRad()
let MEAN_ATMOSPHERIC_REFRACTION_MOON = 0.125.deg

func angleForRiseOrSet(declination: Double, observerLatitude: Double) -> Double {
  let 𝛿 = declination
  let 𝜑 = observerLatitude

  let H0 = acos(-tan(𝜑) * tan(𝛿))
  return H0
}


func hourAngleForRiseOrSet(declination: Double, observerLatitude: Double) -> HourAngle {
  let H0 = angleForRiseOrSet(declination: declination, observerLatitude: observerLatitude)
  return HourAngle(degrees: H0.asDeg)
}

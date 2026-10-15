//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2026 Mattias Holm
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

// Internal coefficient evaluator used by the public IAUOrientationProvider.
// Coefficients come from the pinned NAIF PCK; all stored angles are in degrees.
enum IAURotationProvenance {
  case iau2015Corrected2019
  case legacy2009 // Earth/Moon approximations retained by NAIF, absent from 2015 report.
}

enum IAURotationError: Error, Equatable {
  case unknownBody, nonFinite
}

struct IAUAngle {
  let radians: Double
  let rate: Double // radians per TDB second
  let acceleration: Double // radians per TDB second squared
}

struct IAURotationalElements {
  let rightAscension: IAUAngle
  let declination: IAUAngle
  let primeMeridian: IAUAngle
}

struct IAURotationModelData {
  let naifID: Int
  let provenance: IAURotationProvenance
  let referenceJulianDate: Double
  let rightAscension: IAUAngleSeries
  let declination: IAUAngleSeries
  let primeMeridian: IAUAngleSeries
}

struct IAUHarmonic {
  enum Function { case sine, cosine }
  let amplitudeDegrees: Double
  // Polynomial in TDB centuries; includes any integer harmonic multipliers.
  let phaseDegrees: [Double]
  let function: Function
}

struct IAUAngleSeries {
  enum Variable { case days, centuries }
  let polynomialDegrees: [Double]
  let variable: Variable
  let harmonics: [IAUHarmonic]

  func evaluate(days: Double) throws -> IAUAngle {
    let century = days / 36_525
    let secondsPerUnit: Double = variable == .days ? 86_400 : 36_525 * 86_400
    let x = variable == .days ? days : century
    let base = polynomial(polynomialDegrees, at: x)
    var value = base.value
    var first = base.first / secondsPerUnit
    var second = base.second / (secondsPerUnit * secondsPerUnit)
    let radiansPerDegree = Double.pi / 180
    let centurySeconds: Double = 36_525 * 86_400
    for term in harmonics {
      let phase = polynomial(term.phaseDegrees, at: century)
      let angle = phase.value.truncatingRemainder(dividingBy: 360) * radiansPerDegree
      let phaseRate = phase.first * radiansPerDegree / centurySeconds
      let phaseAcceleration = phase.second * radiansPerDegree / (centurySeconds * centurySeconds)
      switch term.function {
      case .sine:
        value += term.amplitudeDegrees * sin(angle)
        first += term.amplitudeDegrees * cos(angle) * phaseRate
        second += term.amplitudeDegrees * (cos(angle) * phaseAcceleration - sin(angle) * phaseRate * phaseRate)
      case .cosine:
        value += term.amplitudeDegrees * cos(angle)
        first -= term.amplitudeDegrees * sin(angle) * phaseRate
        second -= term.amplitudeDegrees * (sin(angle) * phaseAcceleration + cos(angle) * phaseRate * phaseRate)
      }
    }
    guard value.isFinite, first.isFinite, second.isFinite else { throw IAURotationError.nonFinite }
    return IAUAngle(radians: value * radiansPerDegree,
      rate: first * radiansPerDegree, acceleration: second * radiansPerDegree)
  }
}

// Simultaneous Horner evaluation of a polynomial and both derivatives.
private func polynomial(_ coefficients: [Double], at x: Double) -> (value: Double, first: Double, second: Double) {
  var value = 0.0, first = 0.0, second = 0.0
  for coefficient in coefficients.reversed() {
    second = second * x + 2 * first
    first = first * x + value
    value = value * x + coefficient
  }
  return (value, first, second)
}

struct IAURotModel {
  /// Explicit TDB days since JD 2451545.0; no mutable shared argument cache.
  func elements(body: String, tdbDaysSinceJ2000 days: Double) throws -> IAURotationalElements {
    guard days.isFinite else { throw IAURotationError.nonFinite }
    guard let model = iauRotationModels[body] else { throw IAURotationError.unknownBody }
    let modelDays = days - (model.referenceJulianDate - 2_451_545)
    return try IAURotationalElements(rightAscension: model.rightAscension.evaluate(days: modelDays),
      declination: model.declination.evaluate(days: modelDays), primeMeridian: model.primeMeridian.evaluate(days: modelDays))
  }
}

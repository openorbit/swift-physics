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

// NACA: 4 digits describing airfoil

/// Four-digit geometry with the conventional finite trailing-edge thickness.
/// Source: NASA-TM-X-3284, https://ntrs.nasa.gov/citations/19760003945
public struct NACA4: Sendable {
  public let camberPercentage: Int
  public let maximumCamber: Int
  public let maximumThicknessPercentage: Int

  /// Four-digit NACA parameters; coordinates are normalized by chord.
  public init(camberPercentage: Int, maximumCamber: Int, maximumThicknessPercentage: Int) {
    precondition(Self.valid(camberPercentage, maximumCamber, maximumThicknessPercentage), "Invalid NACA four-digit geometry")
    self.camberPercentage = camberPercentage
    self.maximumCamber = maximumCamber
    self.maximumThicknessPercentage = maximumThicknessPercentage
  }
  public init(validatingCamberPercentage camberPercentage: Int, maximumCamber: Int, maximumThicknessPercentage: Int) throws {
    guard Self.valid(camberPercentage, maximumCamber, maximumThicknessPercentage) else { throw AerodynamicError.invalidGeometry }
    self.camberPercentage = camberPercentage; self.maximumCamber = maximumCamber
    self.maximumThicknessPercentage = maximumThicknessPercentage
  }
  private static func valid(_ m: Int, _ p: Int, _ t: Int) -> Bool {
    (0...9).contains(m) && (0...9).contains(p) && (0...99).contains(t) && (m == 0 || p > 0)
  }
  private func checkChord(_ x: Double) {
    precondition(x.isFinite && (0...1).contains(x), "Normalized chord coordinate must be in 0...1")
  }
  public func thickness(_ x: Double) -> Double {
    checkChord(x)
    let t = Double(maximumThicknessPercentage)/100.0
    let yt = 5 * t * (0.2969 * sqrt(x) -
                      0.1260 * x -
                      0.3516 * x * x +
                      0.2843 * x * x * x -
                      0.1015 * x * x * x * x)
    return yt
  }
  public func meanCamberLine(_ x: Double) -> Double {
    checkChord(x)
    if camberPercentage == 0 { return 0 }
    let m = Double(camberPercentage) / 100.0
    let p = Double(maximumCamber) / 10.0

    if 0 <= x && x <= p {
      let yc = m / (p * p) * (2 * p * x - x * x)
      return yc
    }

    let yc = m / ((1 - p) * (1 - p)) * ((1 - 2 * p) + 2 * p * x - x * x)
    return yc
  }

  /// Derivative of the normalized camber line with respect to normalized chord.
  public func meanCamberSlope(_ x: Double) -> Double {
    checkChord(x)
    guard camberPercentage != 0 else { return 0 }
    let m = Double(camberPercentage)/100, p = Double(maximumCamber)/10
    return x <= p ? 2*m*(p-x)/(p*p) : 2*m*(p-x)/((1-p)*(1-p))
  }

  /// Checked surface evaluation for caller-owned geometry input.
  public func surfaceCoordinates(atChordFraction x: Double) throws -> (SIMD2<Double>, SIMD2<Double>) {
    guard x.isFinite, (0...1).contains(x) else { throw AerodynamicError.invalidGeometry }
    return surfaceCoord(x: x)
  }

  public func surfaceCoord(x: Double) -> (SIMD2<Double>, SIMD2<Double>) {
    let yt = thickness(x)
    if camberPercentage == 0 && maximumCamber == 0 {
      let u = SIMD2(x, yt)
      let l = SIMD2(x, -yt)
      return (u, l)
    }

    let yc = meanCamberLine(x)
    let dyx = meanCamberSlope(x)
    let theta = atan(dyx)
    let xu = x - yt * sin(theta)
    let xl = x + yt * sin(theta)
    let yu = yc + yt * cos(theta)
    let yl = yc - yt * cos(theta)
    let u = SIMD2(xu, yu)
    let l = SIMD2(xl, yl)
    return (u, l)
  }
}

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

public struct NACA4 {
  public let camberPercentage: Int
  public let maximumCamber: Int
  public let maximumThicknessPercentage: Int

  public init(camberPercentage: Int, maximumCamber: Int, maximumThicknessPercentage: Int) {
    self.camberPercentage = camberPercentage
    self.maximumCamber = maximumCamber
    self.maximumThicknessPercentage = maximumThicknessPercentage
  }
  public func thickness(_ x: Double) -> Double {
    let t = Double(maximumThicknessPercentage)/100.0
    let yt = 5 * t * (0.2969 * sqrt(x) -
                      0.1260 * x -
                      0.3516 * x * x +
                      0.2843 * x * x * x -
                      0.1015 * x * x * x * x)
    return yt
  }
  public func meanCamberLine(_ x: Double) -> Double {
    let m = Double(camberPercentage) / 100.0
    let p = Double(maximumCamber) / 10.0

    if 0 <= x && x <= p {
      let yc = m / (p * p) * (2 * p * x - x * x)
      return yc
    }

    let yc = m / ((1 - p) * (1 - p)) * ((1 - 2 * p) + 2 * p * x - x * x)
    return yc
  }

  public func surfaceCoord(x: Double) -> (SIMD2<Double>, SIMD2<Double>) {
    let yt = thickness(x)
    if camberPercentage == 0 && maximumCamber == 0 {
      let u = SIMD2(x, yt)
      let l = SIMD2(x, -yt)
      return (u, l)
    }

    let yc = meanCamberLine(x)
    let dyx = 0.0
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

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

/// Supplied E2x3 series, with degree-aware trigonometry and local evaluation state.
/// Positions are metres in Jupiter's equatorial plane, x along the ascending
/// node on the ecliptic of date, y prograde. No observer/light-time corrections.
/// Legacy coefficients are retained; this is not the later E5 theory.
public enum LieskeE2x3 {
  public enum Moon: Int, CaseIterable, Sendable { case io, europa, ganymede, callisto }
  public static func position(of moon: Moon, at epoch: Epoch) throws -> SIMD3<Double> {
    let date = try epoch.julianDate(in: .tt)
    let t = Double(date.wholeDay - 2443000) + date.fractionalDay - 0.5
    let d = Double.pi / 180
    var p = [Double](repeating: 0, count: 26)
    for i in 0..<16 { p[i] = rules[i].0 + rules[i].1*t }
    p[14] += 0.33033*sin((163.679+0.0010512*t)*d)
      + 0.03439*sin((34.486-0.0161731*t)*d)
    p[24] = 1; p[25] = 13.469942
    let corrections = (0..<4).map { series(longitude[$0], p, sine: true) }
    for i in 0..<4 { p[20+i] = corrections[i]; p[16+i] = p[i]+corrections[i] }
    let i = moon.rawValue
    let latitude = atan(series(latitude[i], p, sine: true))
    let radius = distances[i] * (1 + series(radius[i], p, sine: false)) * 71_492_000
    let longitude = (p[16+i]-p[13])*d
    return SIMD3(radius*cos(latitude)*cos(longitude), radius*cos(latitude)*sin(longitude), radius*sin(latitude))
  }
  private enum Parameter: Int, Sendable { case ell1, ell2, ell3, ell4, pi1, pi2, pi3, pi4, w1, w2, w3, w4, phi, psi, g, gp, l1, l2, l3, l4, z1, z2, z3, z4, constant, piC }
  private struct Term: Sendable { let amplitude: Double; let argument: [(Double, Parameter)] }
  private static func series(_ terms: [Term], _ p: [Double], sine: Bool) -> Double {
    terms.reduce(0) { result, term in
      let angle = term.argument.reduce(0) { $0 + $1.0*p[$1.1.rawValue] } * .pi/180
      return result + term.amplitude * (sine ? sin(angle) : cos(angle))
    }
  }
  private static let distances = [5.90730, 9.39912, 14.99240, 26.36990]
  private static let rules: [(Double,Double)] = [(106.07947,203.488955432),(175.72938,101.374724550),(120.55434,50.317609110),(84.44868,21.571071314),(58.3329,0.16103936),(132.8959,0.04647985),(187.2887,0.00712740),(335.3418,0.00183998),(331.0793,0.13279430),(100.5099,0.03263047),(119.1688,0.00717704),(322.6729,0.00175934),(191.8132,0.17390023),(316.5182,0.00000208),(30.23756,0.0830925701),(31.97853,0.0334597339)]
  private static let longitude: [[Term]] = [
    [
      Term(amplitude: +0.47259, argument: [(2.0, .ell1), (-2.0, .ell2)]),
      Term(amplitude: -0.03480, argument: [(1.0, .pi1), (-1.0, .pi4)]),
      Term(amplitude: -0.01756, argument: [(1.0, .pi1), (1.0, .pi3), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: +0.01080, argument: [(1.0, .ell2), (-2.0, .ell3), (1.0, .pi3)]),
      Term(amplitude: +0.00757, argument: [(1.0, .phi)]),
      Term(amplitude: +0.00663, argument: [(1.0, .ell2), (-2.0, .ell3), (1.0, .pi4)]),
      Term(amplitude: +0.00453, argument: [(1.0, .ell1), (-1.0, .pi3)]),
      Term(amplitude: +0.00453, argument: [(1.0, .ell2), (-2.0, .ell3), (1.0, .pi2)]),
      Term(amplitude: -0.00354, argument: [(1.0, .ell1), (-1.0, .ell2)]),
      Term(amplitude: -0.00317, argument: [(2.0, .psi), (-2.0, .piC)]),
      Term(amplitude: -0.00269, argument: [(1.0, .ell2), (-2.0, .ell3), (1.0, .pi1)]),
      Term(amplitude: +0.00263, argument: [(1.0, .ell1), (-1.0, .pi4)]),
      Term(amplitude: +0.00186, argument: [(1.0, .ell1), (-1.0, .pi1)]),
      Term(amplitude: -0.00186, argument: [(1.0, .g)]),
      Term(amplitude: +0.00167, argument: [(1.0, .pi2), (-1.0, .pi3)]),
      Term(amplitude: +0.00158, argument: [(4.0, .ell1), (-4.0, .ell2)]),
      Term(amplitude: -0.00155, argument: [(1.0, .ell1), (-1.0, .ell3)]),
      Term(amplitude: -0.00142, argument: [(1.0, .psi), (1.0, .w3), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: -0.00115, argument: [(2.0, .ell1), (-4.0, .ell2), (2.0, .w2)]),
      Term(amplitude: +0.00089, argument: [(1.0, .pi2), (-1.0, .pi4)]),
      Term(amplitude: +0.00084, argument: [(1.0, .w2), (-1.0, .w3)]),
      Term(amplitude: +0.00084, argument: [(1.0, .ell1), (1.0, .pi3), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: +0.00053, argument: [(1.0, .psi), (-1.0, .w2)]),
    ],
    [
      Term(amplitude: +1.06476, argument: [(2.0, .ell2), (-2.0, .ell3)]),
      Term(amplitude: 0.04253, argument: [(1.0, .ell1), (-2.0, .ell2), (1.0, .pi3)]),
      Term(amplitude: 0.03579, argument: [(1.0, .ell2), (-1.0, .pi3)]),
      Term(amplitude: 0.02383, argument: [(1.0, .ell1), (-2.0, .ell2), (1.0, .pi4)]),
      Term(amplitude: 0.01977, argument: [(1.0, .ell2), (-1.0, .pi4)]),
      Term(amplitude: -0.01843, argument: [(1.0, .phi)]),
      Term(amplitude: 0.01299, argument: [(1.0, .pi3), (-1.0, .pi4)]),
      Term(amplitude: -0.01142, argument: [(1.0, .ell2), (-1.0, .ell3)]),
      Term(amplitude: 0.01078, argument: [(1.0, .ell2), (-1.0, .pi2)]),
      Term(amplitude: -0.01058, argument: [(1.0, .g)]),
      Term(amplitude: 0.00870, argument: [(1.0, .ell2), (-2.0, .ell3), (1.0, .pi2)]),
      Term(amplitude: -0.00775, argument: [(2.0, .psi), (-2.0, .piC)]),
      Term(amplitude: 0.00524, argument: [(2.0, .ell1), (-2.0, .ell2)]),
      Term(amplitude: -0.00460, argument: [(1.0, .ell1), (-1.0, .ell3)]),
      Term(amplitude: 0.00450, argument: [(1.0, .ell2), (-2.0, .ell3), (1.0, .pi1)]),
      Term(amplitude: 0.00327, argument: [(1.0, .psi), (-2.0, .g), (1.0, .w3), (-2.0, .piC)]),
      Term(amplitude: -0.00296, argument: [(1.0, .pi1), (1.0, .pi3), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: -0.00151, argument: [(2.0, .g)]),
      Term(amplitude: 0.00146, argument: [(1.0, .psi), (-1.0, .w3)]),
      Term(amplitude: 0.00125, argument: [(1.0, .psi), (-1.0, .w4)]),
      Term(amplitude: -0.00117, argument: [(1.0, .ell1), (-2.0, .ell3), (1.0, .pi3)]),
      Term(amplitude: -0.00095, argument: [(2.0, .ell2), (-2.0, .w2)]),
      Term(amplitude: 0.00086, argument: [(2.0, .ell1), (-4.0, .ell2), (2.0, .w2)]),
      Term(amplitude: -0.00086, argument: [(5.0, .gp), (-2.0, .g), (52.225, .constant)]),
      Term(amplitude: -0.00078, argument: [(1.0, .ell2), (-1.0, .ell4)]),
      Term(amplitude: -0.00064, argument: [(1.0, .ell1), (-2.0, .ell3), (1.0, .pi4)]),
      Term(amplitude: -0.00063, argument: [(3.0, .ell3), (-7.0, .ell4), (4.0, .pi4)]),
      Term(amplitude: 0.00061, argument: [(1.0, .pi1), (-1.0, .pi4)]),
      Term(amplitude: 0.00058, argument: [(2.0, .psi), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: 0.00058, argument: [(1.0, .w3), (-1.0, .w4)]),
      Term(amplitude: 0.00056, argument: [(2.0, .ell2), (-2.0, .ell4)]),
      Term(amplitude: 0.00055, argument: [(2.0, .ell1), (-2.0, .ell3)]),
      Term(amplitude: 0.00052, argument: [(3.0, .ell3), (-7.0, .ell4), (1.0, .pi3), (3.0, .pi4)]),
      Term(amplitude: -0.00043, argument: [(1.0, .ell1), (-1.0, .pi3)]),
      Term(amplitude: 0.00042, argument: [(1.0, .pi3), (-1.0, .pi2)]),
      Term(amplitude: 0.00041, argument: [(5.0, .ell2), (-5.0, .ell3)]),
      Term(amplitude: 0.00041, argument: [(1.0, .pi4), (-1.0, .piC)]),
      Term(amplitude: 0.00038, argument: [(1.0, .ell2), (-1.0, .pi1)]),
      Term(amplitude: 0.00032, argument: [(1.0, .w2), (-1.0, .w3)]),
      Term(amplitude: 0.00032, argument: [(2.0, .ell3), (-2.0, .g), (-2.0, .piC)]),
      Term(amplitude: 0.00029, argument: [(1.0, .pi1), (-1.0, .pi3)]),
    ],
    [
      Term(amplitude: 0.16477, argument: [(1.0, .ell3), (-1.0, .pi3)]),
      Term(amplitude: 0.09062, argument: [(1.0, .ell3), (-1.0, .pi4)]),
      Term(amplitude: -0.06907, argument: [(1.0, .ell2), (-1.0, .ell3)]),
      Term(amplitude: 0.03786, argument: [(1.0, .pi3), (-1.0, .pi4)]),
      Term(amplitude: 0.01844, argument: [(2.0, .ell3), (-2.0, .ell4)]),
      Term(amplitude: -0.01340, argument: [(1.0, .g)]),
      Term(amplitude: 0.00703, argument: [(1.0, .ell2), (-2.0, .ell3), (1.0, .pi3)]),
      Term(amplitude: -0.00670, argument: [(2.0, .psi), (-2.0, .piC)]),
      Term(amplitude: -0.00540, argument: [(1.0, .ell3), (-1.0, .ell4)]),
      Term(amplitude: 0.00481, argument: [(1.0, .pi1), (1.0, .pi2), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: -0.00409, argument: [(1.0, .ell2), (-2.0, .ell3), (1.0, .pi2)]),
      Term(amplitude: 0.00379, argument: [(1.0, .ell2), (-2.0, .ell3), (1.0, .pi4)]),
      Term(amplitude: 0.00235, argument: [(1.0, .psi), (-1.0, .w3)]),
      Term(amplitude: 0.00198, argument: [(1.0, .psi), (-1.0, .w4)]),
      Term(amplitude: 0.00180, argument: [(1.0, .phi)]),
      Term(amplitude: 0.00129, argument: [(3.0, .ell3), (-3.0, .ell4)]),
      Term(amplitude: 0.00124, argument: [(1.0, .ell1), (-1.0, .ell3)]),
      Term(amplitude: -0.00119, argument: [(5.0, .gp), (-2.0, .g), (52.225, .constant)]),
      Term(amplitude: 0.00109, argument: [(1.0, .ell1), (-1.0, .ell2)]),
      Term(amplitude: -0.00099, argument: [(3.0, .ell3), (-7.0, .ell4), (4.0, .pi4)]),
      Term(amplitude: 0.00091, argument: [(1.0, .w3), (-1.0, .w4)]),
      Term(amplitude: 0.00081, argument: [(3.0, .ell3), (-7.0, .ell4), (1.0, .pi3), (3.0, .pi4)]),
      Term(amplitude: -0.00076, argument: [(2.0, .ell2), (-3.0, .ell3), (1.0, .pi3)]),
      Term(amplitude: 0.00069, argument: [(1.0, .pi4), (-1.0, .piC)]),
      Term(amplitude: -0.00058, argument: [(2.0, .ell3), (-3.0, .ell4), (1.0, .pi4)]),
      Term(amplitude: 0.00057, argument: [(1.0, .ell3), (1.0, .pi3), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: -0.00057, argument: [(1.0, .ell3), (-2.0, .ell4), (1.0, .pi4)]),
      Term(amplitude: -0.00052, argument: [(1.0, .pi2), (-1.0, .pi3)]),
      Term(amplitude: -0.00052, argument: [(1.0, .ell2), (-2.0, .ell3), (1.0, .pi1)]),
      Term(amplitude: 0.00048, argument: [(1.0, .ell3), (-2.0, .ell4), (1.0, .pi3)]),
      Term(amplitude: -0.00045, argument: [(2.0, .ell2), (-3.0, .ell3), (1.0, .pi4)]),
      Term(amplitude: -0.00041, argument: [(1.0, .pi2), (-1.0, .pi4)]),
      Term(amplitude: -0.00038, argument: [(2.0, .g)]),
      Term(amplitude: -0.00033, argument: [(1.0, .pi3), (-1.0, .pi4), (1.0, .w3), (-1.0, .w4)]),
      Term(amplitude: -0.00032, argument: [(3.0, .ell3), (-7.0, .ell4), (2.0, .pi3), (2.0, .pi4)]),
      Term(amplitude: 0.00030, argument: [(4.0, .ell3), (-4.0, .ell4)]),
      Term(amplitude: -0.00029, argument: [(1.0, .w3), (1.0, .psi), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: 0.00029, argument: [(1.0, .ell3), (1.0, .pi4), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: 0.00026, argument: [(1.0, .ell3), (-1.0, .piC), (-2.0, .g)]),
      Term(amplitude: 0.00024, argument: [(1.0, .ell2), (-3.0, .ell3), (2.0, .ell4)]),
      Term(amplitude: 0.00021, argument: [(2.0, .ell3), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: -0.00021, argument: [(1.0, .ell3), (-1.0, .pi2)]),
      Term(amplitude: 0.00017, argument: [(2.0, .ell3), (-2.0, .pi3)]),
    ],
    [
      Term(amplitude: 0.84109, argument: [(1.0, .ell4), (-1.0, .pi4)]),
      Term(amplitude: 0.03429, argument: [(1.0, .pi4), (-1.0, .pi3)]),
      Term(amplitude: -0.03305, argument: [(2.0, .psi), (-2.0, .piC)]),
      Term(amplitude: -0.03211, argument: [(1.0, .g)]),
      Term(amplitude: -0.01860, argument: [(1.0, .ell4), (-1.0, .pi3)]),
      Term(amplitude: 0.01182, argument: [(1.0, .psi), (-1.0, .w4)]),
      Term(amplitude: 0.00622, argument: [(1.0, .ell4), (1.0, .pi4), (-3.0, .g), (-2.0, .piC)]),
      Term(amplitude: 0.00385, argument: [(2.0, .ell4), (-2.0, .pi4)]),
      Term(amplitude: -0.00284, argument: [(5.0, .gp), (-2.0, .g), (52.225, .constant)]),
      Term(amplitude: -0.00233, argument: [(2.0, .psi), (-2.0, .pi4)]),
      Term(amplitude: -0.00223, argument: [(1.0, .ell3), (-1.0, .ell4)]),
      Term(amplitude: -0.00208, argument: [(1.0, .ell4), (-1.0, .piC)]),
      Term(amplitude: 0.00177, argument: [(1.0, .psi), (1.0, .w4), (-2.0, .pi4)]),
      Term(amplitude: 0.00134, argument: [(1.0, .pi4), (-1.0, .piC)]),
      Term(amplitude: 0.00125, argument: [(2.0, .ell4), (-2.0, .g), (-2.0, .piC)]),
      Term(amplitude: -0.00117, argument: [(2.0, .g)]),
      Term(amplitude: -0.00112, argument: [(2.0, .ell3), (-2.0, .ell4)]),
      Term(amplitude: 0.00106, argument: [(3.0, .ell3), (-7.0, .ell4), (4.0, .pi4)]),
      Term(amplitude: 0.00102, argument: [(1.0, .ell4), (-1.0, .g), (-1.0, .piC)]),
      Term(amplitude: 0.00096, argument: [(2.0, .ell4), (-1.0, .psi), (-1.0, .w4)]),
      Term(amplitude: 0.00087, argument: [(2.0, .psi), (-2.0, .w4)]),
      Term(amplitude: -0.00087, argument: [(3.0, .ell3), (-7.0, .ell4), (1.0, .pi3), (3.0, .pi4)]),
      Term(amplitude: 0.00085, argument: [(1.0, .ell3), (-2.0, .ell4), (1.0, .pi4)]),
      Term(amplitude: -0.00081, argument: [(2.0, .ell4), (-2.0, .psi)]),
      Term(amplitude: 0.00071, argument: [(1.0, .ell4), (1.0, .pi4), (-2.0, .piC), (-3.0, .g)]),
      Term(amplitude: 0.00060, argument: [(1.0, .ell1), (-1.0, .ell4)]),
      Term(amplitude: -0.00056, argument: [(1.0, .psi), (-1.0, .w3)]),
      Term(amplitude: -0.00055, argument: [(1.0, .ell3), (-2.0, .ell4), (1.0, .pi3)]),
      Term(amplitude: 0.00051, argument: [(1.0, .ell2), (-1.0, .ell4)]),
      Term(amplitude: 0.00042, argument: [(2.0, .psi), (-2.0, .g), (-2.0, .piC)]),
      Term(amplitude: 0.00039, argument: [(2.0, .pi4), (-2.0, .w4)]),
      Term(amplitude: 0.00036, argument: [(1.0, .psi), (1.0, .piC), (-1.0, .pi4), (-1.0, .w4)]),
      Term(amplitude: 0.00035, argument: [(2.0, .gp), (-1.0, .g), (188.37, .constant)]),
      Term(amplitude: -0.00035, argument: [(1.0, .ell4), (-1.0, .pi4), (2.0, .piC), (-2.0, .psi)]),
      Term(amplitude: -0.00032, argument: [(1.0, .ell4), (1.0, .pi4), (-2.0, .piC), (-1.0, .g)]),
      Term(amplitude: 0.00030, argument: [(3.0, .ell3), (-7.0, .ell4), (2.0, .pi3), (2.0, .pi4)]),
      Term(amplitude: 0.00030, argument: [(2.0, .gp), (-2.0, .g), (149.15, .constant)]),
      Term(amplitude: 0.00028, argument: [(1.0, .ell4), (-1.0, .pi4), (2.0, .psi), (-2.0, .piC)]),
      Term(amplitude: -0.00028, argument: [(2.0, .ell4), (-2.0, .w4)]),
      Term(amplitude: -0.00027, argument: [(1.0, .pi3), (-1.0, .pi4), (1.0, .w3), (-1.0, .w4)]),
      Term(amplitude: -0.00026, argument: [(5.0, .gp), (-3.0, .g), (188.37, .constant)]),
      Term(amplitude: 0.00025, argument: [(1.0, .w4), (-1.0, .w3)]),
      Term(amplitude: -0.00025, argument: [(1.0, .ell2), (-3.0, .ell3), (2.0, .ell4)]),
      Term(amplitude: -0.00023, argument: [(3.0, .ell3), (-3.0, .ell4)]),
      Term(amplitude: 0.00021, argument: [(2.0, .ell4), (-2.0, .piC), (-3.0, .g)]),
      Term(amplitude: -0.00021, argument: [(2.0, .ell3), (-3.0, .ell4), (1.0, .pi4)]),
      Term(amplitude: 0.00019, argument: [(1.0, .ell4), (-1.0, .pi4), (-1.0, .g)]),
      Term(amplitude: -0.00019, argument: [(2.0, .ell4), (-1.0, .pi3), (-1.0, .pi4)]),
      Term(amplitude: -0.00018, argument: [(1.0, .ell4), (-1.0, .pi4), (1.0, .g)]),
      Term(amplitude: -0.00016, argument: [(1.0, .ell4), (1.0, .pi3), (-2.0, .piC), (-2.0, .g)]),
    ],
  ]
  private static let latitude: [[Term]] = [
    [
      Term(amplitude: 0.0006502, argument: [(1.0, .l1), (-1.0, .w1)]),
      Term(amplitude: 0.0001835, argument: [(1.0, .l1), (-1.0, .w2)]),
      Term(amplitude: 0.0000329, argument: [(1.0, .l1), (-1.0, .w3)]),
      Term(amplitude: -0.0000311, argument: [(1.0, .l1), (-1.0, .psi)]),
      Term(amplitude: 0.0000093, argument: [(1.0, .l1), (-1.0, .w4)]),
      Term(amplitude: 0.0000075, argument: [(3.0, .l1), (-4.0, .ell2), (-1.9927, .z1), (1.0, .w2)]),
      Term(amplitude: 0.0000046, argument: [(1.0, .l1), (1.0, .psi), (-2.0, .piC), (-2.0, .g)]),
    ],
    [
      Term(amplitude: 0.0081275, argument: [(1.0, .l2), (-1.0, .w2)]),
      Term(amplitude: 0.0004512, argument: [(1.0, .l2), (-1.0, .w3)]),
      Term(amplitude: -0.0003286, argument: [(1.0, .l2), (-1.0, .psi)]),
      Term(amplitude: 0.0001164, argument: [(1.0, .l2), (-1.0, .w4)]),
      Term(amplitude: 0.0000273, argument: [(1.0, .ell1), (-2.0, .ell3), (1.0146, .z2), (1.0, .w2)]),
      Term(amplitude: 0.0000143, argument: [(1.0, .l2), (1.0, .psi), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: -0.0000143, argument: [(1.0, .l2), (-1.0, .w1)]),
      Term(amplitude: 0.0000035, argument: [(1.0, .l2), (-1.0, .psi), (1.0, .g)]),
      Term(amplitude: -0.0000028, argument: [(1.0, .ell1), (-2.0, .ell3), (1.0146, .z2), (1.0, .w3)]),
    ],
    [
      Term(amplitude: 0.0032364, argument: [(1.0, .l3), (-1.0, .w3)]),
      Term(amplitude: -0.0016911, argument: [(1.0, .l3), (-1.0, .psi)]),
      Term(amplitude: 0.0006849, argument: [(1.0, .l3), (-1.0, .w4)]),
      Term(amplitude: -0.0002806, argument: [(1.0, .l3), (-1.0, .w2)]),
      Term(amplitude: 0.0000321, argument: [(1.0, .l3), (1.0, .psi), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: 0.0000051, argument: [(1.0, .l3), (-1.0, .psi), (1.0, .g)]),
      Term(amplitude: -0.0000045, argument: [(1.0, .l3), (-1.0, .psi), (-1.0, .g)]),
      Term(amplitude: -0.0000045, argument: [(1.0, .l3), (1.0, .psi), (-2.0, .piC)]),
      Term(amplitude: 0.0000037, argument: [(1.0, .l3), (1.0, .psi), (-2.0, .piC), (-3.0, .g)]),
      Term(amplitude: 0.0000030, argument: [(2.0, .ell2), (-3.0, .l3), (4.03, .z3), (1.0, .w2)]),
      Term(amplitude: -0.0000021, argument: [(2.0, .ell2), (-3.0, .l3), (4.03, .z3), (1.0, .w3)]),
    ],
    [
      Term(amplitude: -0.0076579, argument: [(1.0, .l4), (-1.0, .psi)]),
      Term(amplitude: 0.0044148, argument: [(1.0, .l4), (-1.0, .w4)]),
      Term(amplitude: -0.0005106, argument: [(1.0, .l4), (-1.0, .w3)]),
      Term(amplitude: 0.0000773, argument: [(1.0, .l4), (1.0, .psi), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: 0.0000104, argument: [(1.0, .l4), (-1.0, .psi), (1.0, .g)]),
      Term(amplitude: -0.0000102, argument: [(1.0, .l4), (-1.0, .psi), (-1.0, .g)]),
      Term(amplitude: 0.0000088, argument: [(1.0, .l4), (1.0, .psi), (-2.0, .piC), (-3.0, .g)]),
      Term(amplitude: -0.0000038, argument: [(1.0, .l4), (1.0, .psi), (-2.0, .piC), (-1.0, .g)]),
    ],
  ]
  private static let radius: [[Term]] = [
    [
      Term(amplitude: -0.0041339, argument: [(2.0, .ell1), (-2.0, .ell2)]),
      Term(amplitude: -0.0000395, argument: [(1.0, .ell1), (-1.0, .pi3)]),
      Term(amplitude: -0.0000214, argument: [(1.0, .ell1), (-1.0, .pi4)]),
      Term(amplitude: 0.0000170, argument: [(1.0, .ell1), (-1.0, .ell2)]),
      Term(amplitude: -0.0000162, argument: [(1.0, .ell1), (-1.0, .pi1)]),
      Term(amplitude: -0.0000130, argument: [(4.0, .ell1), (-4.0, .ell2)]),
      Term(amplitude: 0.0000106, argument: [(1.0, .ell1), (-1.0, .ell3)]),
      Term(amplitude: -0.0000063, argument: [(1.0, .ell1), (1.0, .pi3), (-2.0, .piC), (-2.0, .g)]),
    ],
    [
      Term(amplitude: 0.0093847, argument: [(1.0, .ell1), (-1.0, .ell2)]),
      Term(amplitude: -0.0003114, argument: [(1.0, .ell2), (-1.0, .pi3)]),
      Term(amplitude: -0.0001738, argument: [(1.0, .ell2), (-1.0, .pi4)]),
      Term(amplitude: -0.0000941, argument: [(1.0, .ell2), (-1.0, .pi2)]),
      Term(amplitude: 0.0000553, argument: [(1.0, .ell2), (-1.0, .ell3)]),
      Term(amplitude: 0.0000523, argument: [(1.0, .ell1), (-1.0, .ell3)]),
      Term(amplitude: -0.0000290, argument: [(2.0, .ell1), (-2.0, .ell2)]),
      Term(amplitude: 0.0000166, argument: [(2.0, .ell2), (-2.0, .w2)]),
      Term(amplitude: 0.0000107, argument: [(1.0, .ell1), (-2.0, .ell3), (1.0, .pi3)]),
      Term(amplitude: -0.0000102, argument: [(1.0, .ell2), (-1.0, .pi1)]),
      Term(amplitude: -0.0000091, argument: [(2.0, .ell1), (-2.0, .ell3)]),
    ],
    [
      Term(amplitude: -0.0014377, argument: [(1.0, .ell3), (-1.0, .pi3)]),
      Term(amplitude: -0.0007904, argument: [(1.0, .ell3), (-1.0, .pi4)]),
      Term(amplitude: 0.0006342, argument: [(1.0, .ell2), (-1.0, .ell3)]),
      Term(amplitude: -0.0001758, argument: [(2.0, .ell3), (-2.0, .ell4)]),
      Term(amplitude: 0.0000294, argument: [(1.0, .ell3), (-1.0, .ell4)]),
      Term(amplitude: -0.0000156, argument: [(3.0, .ell3), (-3.0, .ell4)]),
      Term(amplitude: 0.0000155, argument: [(1.0, .ell1), (-1.0, .ell3)]),
      Term(amplitude: -0.0000153, argument: [(1.0, .ell1), (-1.0, .ell2)]),
      Term(amplitude: 0.0000070, argument: [(2.0, .ell2), (-3.0, .ell3), (1.0, .pi3)]),
      Term(amplitude: -0.0000051, argument: [(1.0, .ell3), (1.0, .pi3), (-2.0, .piC), (-2.0, .g)]),
    ],
    [
      Term(amplitude: -0.0073391, argument: [(1.0, .ell4), (-1.0, .pi4)]),
      Term(amplitude: 0.0001620, argument: [(1.0, .ell4), (-1.0, .pi3)]),
      Term(amplitude: 0.0000974, argument: [(1.0, .ell3), (-1.0, .ell4)]),
      Term(amplitude: -0.0000541, argument: [(1.0, .ell4), (1.0, .pi4), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: -0.0000269, argument: [(2.0, .ell4), (-2.0, .pi4)]),
      Term(amplitude: 0.0000182, argument: [(1.0, .ell4), (-1.0, .piC)]),
      Term(amplitude: 0.0000177, argument: [(2.0, .ell3), (-2.0, .ell4)]),
      Term(amplitude: -0.0000167, argument: [(2.0, .ell4), (-1.0, .psi), (-1.0, .w4)]),
      Term(amplitude: 0.0000167, argument: [(1.0, .psi), (-1.0, .w4)]),
      Term(amplitude: -0.0000155, argument: [(2.0, .ell4), (-2.0, .piC), (-2.0, .g)]),
      Term(amplitude: 0.0000142, argument: [(2.0, .ell4), (-2.0, .psi)]),
      Term(amplitude: 0.0000104, argument: [(1.0, .ell1), (-1.0, .ell4)]),
      Term(amplitude: 0.0000092, argument: [(1.0, .ell2), (-1.0, .ell4)]),
      Term(amplitude: -0.0000089, argument: [(1.0, .ell4), (-1.0, .piC), (-1.0, .g)]),
      Term(amplitude: -0.0000062, argument: [(1.0, .ell4), (1.0, .pi4), (-2.0, .piC), (-3.0, .g)]),
      Term(amplitude: 0.0000048, argument: [(2.0, .ell4), (-2.0, .w4)]),
    ],
  ]
}

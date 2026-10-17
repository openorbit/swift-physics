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

/// Port of the supplied, unfinished GUST86 coefficient set. All original
/// amplitudes and argument combinations are retained. Not precision validated.
public enum GUST86 {
  public enum Moon: Int, CaseIterable, Sendable { case miranda, ariel, umbriel, titania, oberon }
  public struct Elements: Sendable {
    /// Radians/day and radians, respectively.
    public let meanMotion: Double
    public let meanLongitude: Double
    /// (e cos(varpi), e sin(varpi)).
    public let eccentricity: SIMD2<Double>
    /// (sin(i/2) cos(Omega), sin(i/2) sin(Omega)).
    public let inclination: SIMD2<Double>
  }
  public static func elements(of moon: Moon, at epoch: Epoch) throws -> Elements {
    let date = try epoch.julianDate(in: .tt)
    let t = Double(date.wholeDay-2444239)+date.fractionalDay-0.5
    let i = moon.rawValue
    func real(_ terms: [Term], sine: Bool) -> Double {
      terms.reduce(0) { $0 + $1.amplitude * (sine ? sin($1.frequency*t+$1.phase) : cos($1.frequency*t+$1.phase)) }
    }
    func complex(_ terms: [Term]) -> SIMD2<Double> {
      terms.reduce(.zero) { $0 + $1.amplitude * SIMD2(cos($1.frequency*t+$1.phase),sin($1.frequency*t+$1.phase)) }
    }
    return Elements(meanMotion: baseN[i]+real(n[i],sine: false),
      meanLongitude: baseLongitude[i]+rates[i]*t+real(longitude[i],sine: true),
      eccentricity: complex(z[i]), inclination: complex(inclination[i]))
  }
  /// Metres in the native GUST86 Uranus equatorial frame (B1950 pole).
  /// The missing scale conversion uses a0*(n0/n)^(2/3) from supplied table 5.
  public static func position(of moon: Moon, at epoch: Epoch) throws -> SIMD3<Double> {
    let e = try elements(of: moon, at: epoch)
    let ecc = hypot(e.eccentricity.x,e.eccentricity.y)
    let peri = atan2(e.eccentricity.y,e.eccentricity.x)
    let node = atan2(e.inclination.y,e.inclination.x)
    let inc = 2*asin(hypot(e.inclination.x,e.inclination.y))
    let m = (e.meanLongitude-peri).truncatingRemainder(dividingBy: 2 * .pi)
    var anomaly = m
    for _ in 0..<12 {
      let delta = (anomaly-ecc*sin(anomaly)-m)/(1-ecc*cos(anomaly))
      anomaly -= delta
      if abs(delta) < 1e-14 { break }
    }
    let i = moon.rawValue
    let a = radii[i]*1000*pow(baseN[i]/e.meanMotion,2.0/3)
    let x = a*(cos(anomaly)-ecc), y = a*sqrt(1-ecc*ecc)*sin(anomaly)
    let w = peri-node
    let u = x*cos(w)-y*sin(w), v = x*sin(w)+y*cos(w)
    return SIMD3(u*cos(node)-v*cos(inc)*sin(node),u*sin(node)+v*cos(inc)*cos(node),v*sin(inc))
  }
  private struct Term: Sendable { let amplitude: Double; let frequency: Double; let phase: Double }
  private static let baseN = [4.44352267,2.49254257,1.51595490,0.72166316,0.46658054]
  private static let baseLongitude = [-0.23805158,3.09804641,2.28540169,0.85635879,-0.91559180]
  private static let rates = [4.44519055,2.49295252,1.51614811,0.72171851,0.46669212]
  private static let radii = [129872.0,190945,265998,436298,583519]
  private static let n: [[Term]] = [
    [
      Term(amplitude: -34.92e-6, frequency: -0.0013707850000002075, phase: 1.3218003071795854),
      Term(amplitude: 8.47e-6, frequency: -0.002741570000000415, phase: 2.643600614359171),
      Term(amplitude: 1.31e-6, frequency: -0.004112354999998402, phase: 3.965400921538759),
      Term(amplitude: -52.28e-6, frequency: 1.9522380310000003, phase: 2.947088307179586),
      Term(amplitude: -136.65e-6, frequency: 3.9044760620000005, phase: 5.894176614359172),
    ],
    [
      Term(amplitude: 2.55e-6, frequency: -0.0013707850000002075, phase: 1.3218003071795854),
      Term(amplitude: -42.16e-6, frequency: 0.9768044080000002, phase: 0.8126440000000001),
      Term(amplitude: -102.56e-6, frequency: 1.9536088160000005, phase: 1.6252880000000003),
    ],
    [
      Term(amplitude: 9.74e-6, frequency: 0.07284798697405517, phase: 2.6404579999999997),
      Term(amplitude: -106.00e-6, frequency: 0.9768044080000002, phase: 0.8126440000000001),
      Term(amplitude: 54.16e-6, frequency: 1.9536088160000005, phase: 1.6252880000000003),
      Term(amplitude: -23.59e-6, frequency: 0.7944296019999999, phase: 1.429043),
      Term(amplitude: -70.70e-6, frequency: 1.5888592039999998, phase: 2.858086),
      Term(amplitude: -36.28e-6, frequency: 2.383288806, phase: 4.287129),
    ],
    [
      Term(amplitude: -2.64e-6, frequency: 0.07284798697405517, phase: 2.6404579999999997),
      Term(amplitude: -2.16e-6, frequency: 0.04337912418079259, phase: 4.886260999999999),
      Term(amplitude: 6.45e-6, frequency: 0.04345997953895008, phase: 5.194624999999999),
      Term(amplitude: -1.11e-6, frequency: 0.04349755197405539, phase: 0.24408269282041317),
      Term(amplitude: -62.23e-6, frequency: 1.7712340100000001, phase: 2.241687),
      Term(amplitude: -56.13e-6, frequency: 0.7944296019999999, phase: 1.429043),
      Term(amplitude: -39.94e-6, frequency: 0.255026389, phase: 1.771951),
      Term(amplitude: -91.85e-6, frequency: 0.510052778, phase: 3.543902),
      Term(amplitude: -58.31e-6, frequency: 0.7650791670000001, phase: 5.315853),
      Term(amplitude: -38.60e-6, frequency: 1.020105556, phase: 0.804618692820414),
      Term(amplitude: -26.18e-6, frequency: 1.275131945, phase: 2.5765696928204136),
      Term(amplitude: -18.06e-6, frequency: 1.5301583340000002, phase: 4.348520692820413),
    ],
    [
      Term(amplitude: 2.08e-6, frequency: 0.04337912418079259, phase: 4.886260999999999),
      Term(amplitude: -6.22e-6, frequency: 0.04345997953895008, phase: 5.194624999999999),
      Term(amplitude: 1.07e-6, frequency: 0.04349755197405539, phase: 0.24408269282041317),
      Term(amplitude: -43.10e-6, frequency: 2.0262603990000003, phase: 4.013638),
      Term(amplitude: -38.94e-6, frequency: 1.0494559909999999, phase: 3.2009939999999997),
      Term(amplitude: -80.11e-6, frequency: 0.255026389, phase: 1.771951),
      Term(amplitude: 59.06e-6, frequency: 0.510052778, phase: 3.543902),
      Term(amplitude: 37.49e-6, frequency: 0.7650791670000001, phase: 5.315853),
      Term(amplitude: 24.82e-6, frequency: 1.020105556, phase: 0.804618692820414),
      Term(amplitude: 16.84e-6, frequency: 1.275131945, phase: 2.5765696928204136),
    ],
  ]
  private static let longitude: [[Term]] = [
    [
      Term(amplitude: 25472.17e-6, frequency: -0.0013707850000002075, phase: 1.3218003071795854),
      Term(amplitude: -3088.31e-6, frequency: -0.002741570000000415, phase: 2.643600614359171),
      Term(amplitude: -318.10e-6, frequency: -0.004112354999998402, phase: 3.965400921538759),
      Term(amplitude: -37.49e-6, frequency: -0.00548314000000083, phase: 5.287201228718342),
      Term(amplitude: -57.85e-6, frequency: 1.9522380310000003, phase: 2.947088307179586),
      Term(amplitude: -62.32e-6, frequency: 3.9044760620000005, phase: 5.894176614359172),
      Term(amplitude: -27.95e-6, frequency: 5.856714093000002, phase: 2.5580796143591726),
    ],
    [
      Term(amplitude: -1860.50e-6, frequency: -0.0013707850000002075, phase: 1.3218003071795854),
      Term(amplitude: 219.99e-6, frequency: -0.002741570000000415, phase: 2.643600614359171),
      Term(amplitude: 23.10e-6, frequency: -0.004112354999998402, phase: 3.965400921538759),
      Term(amplitude: 4.30e-6, frequency: -0.00548314000000083, phase: 5.287201228718342),
      Term(amplitude: -90.11e-6, frequency: 0.9768044080000002, phase: 0.8126440000000001),
      Term(amplitude: -91.07e-6, frequency: 1.9536088160000005, phase: 1.6252880000000003),
      Term(amplitude: -42.75e-6, frequency: 2.9304132240000005, phase: 2.437932),
      Term(amplitude: -16.49e-6, frequency: 3.5424680200000003, phase: 4.483374),
    ],
    [
      Term(amplitude: 660.57e-6, frequency: -0.0013707850000002075, phase: 1.3218003071795854),
      Term(amplitude: -76.51e-6, frequency: -0.002741570000000415, phase: 2.643600614359171),
      Term(amplitude: -8.96e-6, frequency: -0.004112354999998402, phase: 3.965400921538759),
      Term(amplitude: -2.53e-6, frequency: -0.00548314000000083, phase: 5.287201228718342),
      Term(amplitude: -52.91e-6, frequency: 0.029350434999999786, phase: 2.3963753071795866),
      Term(amplitude: -7.34e-6, frequency: 0.07272955918079238, phase: 0.999451),
      Term(amplitude: -1.83e-6, frequency: 0.07281041453894986, phase: 1.307815),
      Term(amplitude: 147.91e-6, frequency: 0.07284798697405517, phase: 2.6404579999999997),
      Term(amplitude: -7.77e-6, frequency: 0.0730081895634291, phase: 2.980758),
      Term(amplitude: 97.76e-6, frequency: 0.9768044080000002, phase: 0.8126440000000001),
      Term(amplitude: 73.13e-6, frequency: 1.9536088160000005, phase: 1.6252880000000003),
      Term(amplitude: 34.71e-6, frequency: 2.9304132240000005, phase: 2.437932),
      Term(amplitude: 18.89e-6, frequency: 3.907217632000001, phase: 3.2505760000000006),
      Term(amplitude: -67.89e-6, frequency: 0.7944296019999999, phase: 1.429043),
      Term(amplitude: -82.86e-6, frequency: 1.5888592039999998, phase: 2.858086),
      Term(amplitude: -33.81e-6, frequency: 2.383288806, phase: 4.287129),
      Term(amplitude: -15.79e-6, frequency: 3.1777184079999996, phase: 5.716172),
      Term(amplitude: -10.21e-6, frequency: 1.0494559909999999, phase: 3.2009939999999997),
      Term(amplitude: -17.08e-6, frequency: 2.0989119819999997, phase: 0.11880269282041311),
    ],
    [
      Term(amplitude: 20.61e-6, frequency: 0.029350434999999786, phase: 2.3963753071795866),
      Term(amplitude: -2.07e-6, frequency: 0.07272955918079238, phase: 0.999451),
      Term(amplitude: -2.88e-6, frequency: 0.07281041453894986, phase: 1.307815),
      Term(amplitude: -40.79e-6, frequency: 0.07284798697405517, phase: 2.6404579999999997),
      Term(amplitude: 2.11e-6, frequency: 0.0730081895634291, phase: 2.980758),
      Term(amplitude: -51.83e-6, frequency: 0.04337912418079259, phase: 4.886260999999999),
      Term(amplitude: 159.87e-6, frequency: 0.04345997953895008, phase: 5.194624999999999),
      Term(amplitude: -35.05e-6, frequency: 0.04349755197405539, phase: 0.24408269282041317),
      Term(amplitude: -1.56e-6, frequency: 0.29840551318079256, phase: 0.37502669282041357),
      Term(amplitude: 40.54e-6, frequency: 1.7712340100000001, phase: 2.241687),
      Term(amplitude: 46.17e-6, frequency: 0.7944296019999999, phase: 1.429043),
      Term(amplitude: -317.76e-6, frequency: 0.255026389, phase: 1.771951),
      Term(amplitude: -305.59e-6, frequency: 0.510052778, phase: 3.543902),
      Term(amplitude: -148.36e-6, frequency: 0.7650791670000001, phase: 5.315853),
      Term(amplitude: -82.92e-6, frequency: 1.020105556, phase: 0.804618692820414),
      Term(amplitude: -49.98e-6, frequency: 1.275131945, phase: 2.5765696928204136),
      Term(amplitude: -31.56e-6, frequency: 1.5301583340000002, phase: 4.348520692820413),
      Term(amplitude: -20.56e-6, frequency: 1.7851847229999995, phase: 6.120471692820413),
      Term(amplitude: -13.69e-6, frequency: 2.040211112, phase: 1.609237385640828),
    ],
    [
      Term(amplitude: -7.82e-6, frequency: 0.029350434999999786, phase: 2.3963753071795866),
      Term(amplitude: 51.29e-6, frequency: 0.04337912418079259, phase: 4.886260999999999),
      Term(amplitude: -158.24e-6, frequency: 0.04345997953895008, phase: 5.194624999999999),
      Term(amplitude: 34.51e-6, frequency: 0.04349755197405539, phase: 0.24408269282041317),
      Term(amplitude: 47.51e-6, frequency: 2.0262603990000003, phase: 4.013638),
      Term(amplitude: 38.96e-6, frequency: 1.0494559909999999, phase: 3.2009939999999997),
      Term(amplitude: 359.73e-6, frequency: 0.255026389, phase: 1.771951),
      Term(amplitude: 282.78e-6, frequency: 0.510052778, phase: 3.543902),
      Term(amplitude: 138.60e-6, frequency: 0.7650791670000001, phase: 5.315853),
      Term(amplitude: 78.03e-6, frequency: 1.020105556, phase: 0.804618692820414),
      Term(amplitude: 47.29e-6, frequency: 1.275131945, phase: 2.5765696928204136),
      Term(amplitude: 30.00e-6, frequency: 1.5301583340000002, phase: 4.348520692820413),
      Term(amplitude: 19.62e-6, frequency: 1.7851847229999995, phase: 6.120471692820413),
      Term(amplitude: 13.11e-6, frequency: 2.040211112, phase: 1.609237385640828),
    ],
  ]
  private static let z: [[Term]] = [
    [
      Term(amplitude: 1312.38e-6, frequency: 0.0009595856065847822, phase: 0.611392),
      Term(amplitude: 71.81e-6, frequency: 0.0002970965634292017, phase: 2.408074),
      Term(amplitude: 69.77e-6, frequency: 0.0001368939740552811, phase: 2.067774),
      Term(amplitude: 6.75e-6, frequency: 9.932153894997327e-05, phase: 0.735131),
      Term(amplitude: 6.27e-6, frequency: 1.8466180792483302e-05, phase: 0.426767),
      Term(amplitude: -123.31e-6, frequency: 0.5407144879999999, phase: 0.15095769282041438),
      Term(amplitude: 39.52e-6, frequency: -1.4115235430000004, phase: 3.487054692820415),
      Term(amplitude: 194.10e-6, frequency: 4.44519055, phase: 6.045134307179586),
    ],
    [
      Term(amplitude: -3.35e-6, frequency: 0.0009595856065847822, phase: 0.611392),
      Term(amplitude: 1187.63e-6, frequency: 0.0002970965634292017, phase: 2.408074),
      Term(amplitude: 861.59e-6, frequency: 0.0001368939740552811, phase: 2.067774),
      Term(amplitude: 71.50e-6, frequency: 9.932153894997327e-05, phase: 0.735131),
      Term(amplitude: 55.59e-6, frequency: 1.8466180792483302e-05, phase: 0.426767),
      Term(amplitude: -84.60e-6, frequency: 0.5393437029999997, phase: 1.4727579999999998),
      Term(amplitude: 91.81e-6, frequency: -0.43746070500000034, phase: 0.6601140000000001),
      Term(amplitude: 20.03e-6, frequency: -1.0495155010000001, phase: 4.897857307179586),
      Term(amplitude: 89.77e-6, frequency: 2.492952519, phase: 3.098046),
    ],
    [
      Term(amplitude: -0.21e-6, frequency: 0.0009595856065847822, phase: 0.611392),
      Term(amplitude: -227.95e-6, frequency: 0.0002970965634292017, phase: 2.408074),
      Term(amplitude: 3904.69e-6, frequency: 0.0001368939740552811, phase: 2.067774),
      Term(amplitude: 309.17e-6, frequency: 9.932153894997327e-05, phase: 0.735131),
      Term(amplitude: 221.92e-6, frequency: 1.8466180792483302e-05, phase: 0.426767),
      Term(amplitude: 29.34e-6, frequency: 2.492952519, phase: 3.098046),
      Term(amplitude: 26.20e-6, frequency: 1.516148111, phase: 2.285402),
      Term(amplitude: 51.19e-6, frequency: 0.5393437029999997, phase: 1.4727579999999998),
      Term(amplitude: -103.86e-6, frequency: -0.43746070500000034, phase: 0.6601140000000001),
      Term(amplitude: -27.16e-6, frequency: -1.4142651130000008, phase: 6.130655307179586),
      Term(amplitude: -16.22e-6, frequency: 0.721718509, phase: 0.856359),
      Term(amplitude: 549.23e-6, frequency: -0.0727110929999999, phase: 5.7105013071795865),
      Term(amplitude: 34.70e-6, frequency: -0.8671406949999998, phase: 4.281458307179586),
      Term(amplitude: 12.81e-6, frequency: -1.661570297, phase: 2.852415307179586),
      Term(amplitude: 21.81e-6, frequency: -0.5827638709999999, phase: 2.1665993071795864),
      Term(amplitude: 46.25e-6, frequency: 1.516148111, phase: 2.285402),
    ],
    [
      Term(amplitude: -0.02e-6, frequency: 0.0009595856065847822, phase: 0.611392),
      Term(amplitude: -1.29e-6, frequency: 0.0002970965634292017, phase: 2.408074),
      Term(amplitude: -324.51e-6, frequency: 0.0001368939740552811, phase: 2.067774),
      Term(amplitude: 932.81e-6, frequency: 9.932153894997327e-05, phase: 0.735131),
      Term(amplitude: 1120.89e-6, frequency: 1.8466180792483302e-05, phase: 0.426767),
      Term(amplitude: 33.86e-6, frequency: 2.492952519, phase: 3.098046),
      Term(amplitude: 17.46e-6, frequency: 0.721718509, phase: 0.856359),
      Term(amplitude: 16.58e-6, frequency: -1.0495155010000001, phase: 4.897857307179586),
      Term(amplitude: 28.89e-6, frequency: 1.516148111, phase: 2.285402),
      Term(amplitude: -35.86e-6, frequency: -0.0727110929999999, phase: 5.7105013071795865),
      Term(amplitude: -17.86e-6, frequency: 0.721718509, phase: 0.856359),
      Term(amplitude: -32.10e-6, frequency: 0.46669212, phase: 5.367593307179586),
      Term(amplitude: -177.83e-6, frequency: 0.21166573099999997, phase: 3.5956423071795864),
      Term(amplitude: 793.43e-6, frequency: -0.04336065800000011, phase: 1.8236913071795868),
      Term(amplitude: 99.48e-6, frequency: -0.29838704700000007, phase: 0.05174030717958633),
      Term(amplitude: 44.83e-6, frequency: -0.553413436, phase: 4.562974614359172),
      Term(amplitude: 25.13e-6, frequency: -0.8084398250000002, phase: 2.7910236143591725),
      Term(amplitude: 15.43e-6, frequency: -1.063466214, phase: 1.0190726143591728),
    ],
    [
      Term(amplitude: 0.00e-6, frequency: 0.0009595856065847822, phase: 0.611392),
      Term(amplitude: -0.35e-6, frequency: 0.0002970965634292017, phase: 2.408074),
      Term(amplitude: 74.53e-6, frequency: 0.0001368939740552811, phase: 2.067774),
      Term(amplitude: -758.68e-6, frequency: 9.932153894997327e-05, phase: 0.735131),
      Term(amplitude: 1397.34e-6, frequency: 1.8466180792483302e-05, phase: 0.426767),
      Term(amplitude: 39.00e-6, frequency: 2.492952519, phase: 3.098046),
      Term(amplitude: 17.66e-6, frequency: -1.559568279, phase: 1.3539553071795858),
      Term(amplitude: 32.42e-6, frequency: 1.516148111, phase: 2.285402),
      Term(amplitude: 79.75e-6, frequency: 0.721718509, phase: 0.856359),
      Term(amplitude: 75.66e-6, frequency: 0.46669212, phase: 5.367593307179586),
      Term(amplitude: 134.04e-6, frequency: 0.21166573099999997, phase: 3.5956423071795864),
      Term(amplitude: -987.26e-6, frequency: -0.04336065800000011, phase: 1.8236913071795868),
      Term(amplitude: -126.09e-6, frequency: -0.29838704700000007, phase: 0.05174030717958633),
      Term(amplitude: -57.42e-6, frequency: -0.553413436, phase: 4.562974614359172),
      Term(amplitude: -32.41e-6, frequency: -0.8084398250000002, phase: 2.7910236143591725),
      Term(amplitude: -19.99e-6, frequency: -1.063466214, phase: 1.0190726143591728),
      Term(amplitude: -12.94e-6, frequency: -1.3184926029999997, phase: 5.5303069215387595),
    ],
  ]
  private static let inclination: [[Term]] = [
    [
      Term(amplitude: 37871.71e-6, frequency: -0.0009704523777238522, phase: 5.702313),
      Term(amplitude: 27.01e-6, frequency: -0.00030041892488663255, phase: 0.395757),
      Term(amplitude: 30.76e-6, frequency: -0.0001355115203772675, phase: 0.589326),
      Term(amplitude: 12.18e-6, frequency: -8.80488861773251e-05, phase: 1.746237),
      Term(amplitude: 5.37e-6, frequency: -1.2376188261917353e-05, phase: 4.206896),
    ],
    [
      Term(amplitude: -121.75e-6, frequency: -0.0009704523777238522, phase: 5.702313),
      Term(amplitude: 358.25e-6, frequency: -0.00030041892488663255, phase: 0.395757),
      Term(amplitude: 290.08e-6, frequency: -0.0001355115203772675, phase: 0.589326),
      Term(amplitude: 97.78e-6, frequency: -8.80488861773251e-05, phase: 1.746237),
      Term(amplitude: 33.97e-6, frequency: -1.2376188261917353e-05, phase: 4.206896),
    ],
    [
      Term(amplitude: -10.86e-6, frequency: -0.0009704523777238522, phase: 5.702313),
      Term(amplitude: -86.51e-6, frequency: -0.00030041892488663255, phase: 0.395757),
      Term(amplitude: 1113.36e-6, frequency: -0.0001355115203772675, phase: 0.589326),
      Term(amplitude: 350.14e-6, frequency: -8.80488861773251e-05, phase: 1.746237),
      Term(amplitude: 106.50e-6, frequency: -1.2376188261917353e-05, phase: 4.206896),
    ],
    [
      Term(amplitude: -1.43e-6, frequency: -0.0009704523777238522, phase: 5.702313),
      Term(amplitude: -1.06e-6, frequency: -0.00030041892488663255, phase: 0.395757),
      Term(amplitude: -140.13e-6, frequency: -0.0001355115203772675, phase: 0.589326),
      Term(amplitude: 685.72e-6, frequency: -8.80488861773251e-05, phase: 1.746237),
      Term(amplitude: 378.32e-6, frequency: -1.2376188261917353e-05, phase: 4.206896),
    ],
    [
      Term(amplitude: -0.44e-6, frequency: -0.0009704523777238522, phase: 5.702313),
      Term(amplitude: -0.31e-6, frequency: -0.00030041892488663255, phase: 0.395757),
      Term(amplitude: 36.89e-6, frequency: -0.0001355115203772675, phase: 0.589326),
      Term(amplitude: -596.33e-6, frequency: -8.80488861773251e-05, phase: 1.746237),
      Term(amplitude: 451.69e-6, frequency: -1.2376188261917353e-05, phase: 4.206896),
    ],
  ]
}

// Supporting source tables retained for future completion and coefficient audits.
// Values retain the original table units (including micro-radians where printed).
internal enum GUST86SourceTables {
  static let table_1_e: [[Double]] = [
    [0.999997, 0.059274, 0.017385, 0.005428, 0.003471],
    [-0.005223, 0.980346, 0.214697, 0.057494, 0.030783],
    [-0.000160,-0.188166, 0.972995, 0.248610, 0.122884],
    [-0.000013,-0.001061,-0.080862, 0.750105, 0.620678],
    [-0.000003,-0.000286, 0.018571,-0.610082, 0.773756],
  ]
  static let table_2a: [[Double]] = [
    [19.440, 1.033, 20.082,  1075.0e-6,  753.0e-6],
    [5.944, 1.046,  6.217,  -900.0e-6,  810.0e-6],
    [2.731, 1.049,  2.865, -1913.0e-6, 3527.0e-6],
    [1.735, 1.198,  2.078,   922.0e-6,  834.0e-6],
    [0.383, 1.009,  0.386,  1639.0e-6,  747.0e-6],
    [-19.453, 1.044, -20.309, 31660.0e-6, -20782.0e-6],
    [-5.999, 1.048,  -6.288,   340.0e-6,    142.0e-6],
    [-2.706, 1.048,  -2.836,   964.0e-6,    645.0e-6],
    [-1.828, 1.008,  -1.843,  -170.0e-6,    964.0e-6],
    [-0.259, 1.000,  -0.259,  -290.0e-6,   -524.0e-6],
  ]
  static let table_2b: [[Double]] = [
    [4.4,  -238051e-6, 4445190.550e-6],
    [86.1,  3098046e-6, 2492952.519e-6],
    [84.0,  2285402e-6, 1516148.111e-6],
    [230.0,   856359e-6,  721718.509e-6],
    [200.0,  -915592e-6,  466692.120e-6],
  ]
  static let table_3b: [[Double]] = [
    [4.8,  -237924, 4445190.857],
    [94.4,  3098578, 2492952.726],
    [71.2,  2285286, 1516148.050],
    [201.3,   856441,  721718.519],
    [233.6,  -915620,  466692.130],
  ]
  static let table_4: [[Double]] = [
    [1312.38e-6, 0.611392],
    [1211.43e-6, 2.408074],
    [4013.06e-6, 2.067774],
    [1243.57e-6, 0.735131],
    [1805.91e-6, 0.426767],
    [37871.91e-6, 5.702313],
    [368.39e-6, 0.395757],
    [1160.02e-6, 0.589326],
    [978.84e-6, 1.746237],
    [599.73e-6, 4.206896],
  ]
  static let table_5: [[Double]] = [
    [0.075, 4445190.550e-6, 4443522.67, 129872],
    [1.49 , 2492952.519e-6, 2492542.57, 190945],
    [1.45 , 1516148.111e-6, 1515954.90, 265998],
    [3.97 ,  721718.509e-6,  721663.16, 436298],
    [3.45 ,  466692.120e-6,  466580.54, 583519],
  ]
}

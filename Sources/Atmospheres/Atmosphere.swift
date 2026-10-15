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

/// Legacy scalar queries: altitude in metres, results in SI units.
public protocol Atmosphere {
  func pressure(at altitude: Double) -> Double
  func temperature(at altitude: Double) -> Double
  func density(at altitude: Double) -> Double
  func speedOfSound(at altitude: Double) -> Double
  func dynamicViscosity(at altitude: Double) -> Double
}

public extension Atmosphere {
  func kinematicViscosity(at altitude: Double) -> Double {
    return dynamicViscosity(at: altitude)/density(at: altitude)
  }
}

/// Compatibility facade for ISA. All values use SI (K, Pa, kg/m³, m/s, Pa s).
/// Out-of-coverage and nonfinite altitudes return NaN; use ISAAtmosphere for checked queries.
public struct EarthAtmosphere: Atmosphere {
  public init() {}
  private func sample(_ altitude: Double) -> AtmosphericProperties? {
    try? ISAAtmosphere().properties(at: altitude)
  }
  public func pressure(at altitude: Double) -> Double { sample(altitude)?.pressurePascals ?? .nan }
  public func temperature(at altitude: Double) -> Double { sample(altitude)?.temperatureKelvin ?? .nan }
  public func density(at altitude: Double) -> Double { sample(altitude)?.densityKilogramsPerCubicMetre ?? .nan }
  public func speedOfSound(at altitude: Double) -> Double { sample(altitude)?.speedOfSoundMetresPerSecond ?? .nan }
  public func dynamicViscosity(at altitude: Double) -> Double { sample(altitude)?.dynamicViscosityPascalSeconds ?? .nan }
}

extension Atmosphere {
  func reynoldsNumber(at altitude: Double, forLength length: Double, atVelocity velocity: Double) -> Double {
    // Re = density * velocity * length / viscosity
    return self.density(at: altitude) * velocity * length / self.dynamicViscosity(at: altitude)
  }
}

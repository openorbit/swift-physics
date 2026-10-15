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

import FluidDynamics
import Foundation
import Yams

public protocol Airfoil {
  func liftCoefficient(at alpha: Double, for reynoldsNumber: Double) -> Double
  func dragCoefficient(at alpha: Double, for reynoldsNumber: Double) -> Double
  func momentCoefficient(at alpha: Double, for reynoldsNumber: Double) -> Double
  func coefficients(at alpha: Double, for reynoldsNumber: Double) -> AirfoilCoefficients
}

public struct AirfoilCoefficients: Codable, Sendable {
  public let lift: Double
  public let drag: Double
  public let moment: Double
  public init(lift: Double, drag: Double, moment: Double) {
    self.lift = lift; self.drag = drag; self.moment = moment
  }
  internal var isFinite: Bool { [lift,drag,moment].allSatisfy(\.isFinite) }
}

/// Legacy scalar facade. Queries use radians and explicitly clamp to the loaded
/// angle/Reynolds/Mach coverage with a Mach-zero query. Invalid queries return NaN; prefer the
/// throwing TabulatedAirfoil API for coverage-sensitive physics.
public struct ParametricAirfoil: Airfoil, Codable, Sendable {
  private let model: TabulatedAirfoil
  private init(model: TabulatedAirfoil) { self.model = model }

  /// Coefficients at one-degree intervals from −180° through +180°, inclusive.
  /// Use validatingReynoldsNumber for recoverable validation errors.
  public init(reynoldsNumber: Double, coefficients: [AirfoilCoefficients]) {
    do { try self.init(validatingReynoldsNumber: reynoldsNumber, coefficients: coefficients) }
    catch { preconditionFailure("Invalid coefficient table: \(error)") }
  }

  public init(validatingReynoldsNumber reynoldsNumber: Double, coefficients: [AirfoilCoefficients]) throws {
    guard coefficients.count == 361 else { throw FluidError.invalidCoefficients }
    let samples = coefficients.enumerated().map {
      AirfoilSample(angleOfAttackRadians: Double($0.offset-180) * .pi/180, coefficients: $0.element)
    }
    model = try TabulatedAirfoil(polars: [AirfoilPolar(reynoldsNumber: reynoldsNumber,samples: samples)],
                                boundaryPolicy: .clamp,source: "Caller-supplied degree grid",version: "1")
  }

  public init?(withString string: String) {
    guard let parsed = try? AirfoilTableLoader.yaml(string,boundaryPolicy: .clamp) else { return nil }
    model = parsed
  }
  public init?(withFile file: String) {
    guard let parsed = try? AirfoilTableLoader.yamlFile(URL(fileURLWithPath: file),boundaryPolicy: .clamp) else { return nil }
    model = parsed
  }

  public func liftCoefficient(at alpha: Double, for reynoldsNumber: Double) -> Double {
    return coefficients(at: alpha,for: reynoldsNumber).lift
  }
  public func dragCoefficient(at alpha: Double, for reynoldsNumber: Double) -> Double {
    return coefficients(at: alpha,for: reynoldsNumber).drag
  }
  public func momentCoefficient(at alpha: Double, for reynoldsNumber: Double) -> Double {
    return coefficients(at: alpha,for: reynoldsNumber).moment
  }
  public func coefficients(at alpha: Double, for reynoldsNumber: Double) -> AirfoilCoefficients {
    guard alpha.isFinite, abs(alpha) <= .pi else {
      return AirfoilCoefficients(lift: .nan,drag: .nan,moment: .nan)
    }
    return (try? model.coefficients(angleOfAttackRadians: alpha,reynoldsNumber: reynoldsNumber,machNumber: 0))
      ?? AirfoilCoefficients(lift: .nan,drag: .nan,moment: .nan)
  }

  private enum CodingKeys: String, CodingKey { case model, reynoldsNumber, coefficients }
  public init(from decoder: Decoder) throws {
    let c = try decoder.container(keyedBy: CodingKeys.self)
    if c.contains(.model) { try self.init(model: c.decode(TabulatedAirfoil.self,forKey: .model)) }
    else {
      var values = try c.decode([AirfoilCoefficients].self,forKey: .coefficients)
      // The original YAML loader added an unused −181° padding sample.
      if values.count == 362 { values.removeFirst() }
      try self.init(validatingReynoldsNumber: c.decode(Double.self,forKey: .reynoldsNumber),coefficients: values)
    }
  }
  public func encode(to encoder: Encoder) throws {
    var c = encoder.container(keyedBy: CodingKeys.self)
    try c.encode(model,forKey: .model)
  }
}

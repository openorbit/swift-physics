// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Math
@testable import CelMek

@Suite
struct IAURotationTests {
  struct Reference: Decodable {
    let toolkit: String
    let kernelSHA256: String
    let samples: [Sample]
  }
  struct Sample: Decodable {
    let body: String
    let days: Double
    let rotation: [Double]
    let derivative: [Double]
    let secondDerivative: [Double]
  }

  // Coordinate mapping from corrected node convention (90° + alpha).
  // Each elemental mapping supplies the analytic derivatives of its Euler angle.
  func mapping(_ elements: IAURotationalElements) throws -> FrameTransform {
    let ra = elements.rightAscension, dec = elements.declination, w = elements.primeMeridian
    let node = try FrameTransform(epoch: .j2000,
      rotation: Quaternion(angle: -(ra.radians + .pi/2), axis: SIMD3(0, 0, 1)),
      angularVelocity: SIMD3(0, 0, -ra.rate), angularAcceleration: SIMD3(0, 0, -ra.acceleration))
    let tilt = try FrameTransform(epoch: .j2000,
      rotation: Quaternion(angle: dec.radians - .pi/2, axis: SIMD3(1, 0, 0)),
      angularVelocity: SIMD3(dec.rate, 0, 0), angularAcceleration: SIMD3(dec.acceleration, 0, 0))
    let spin = try FrameTransform(epoch: .j2000,
      rotation: Quaternion(angle: -w.radians, axis: SIMD3(0, 0, 1)),
      angularVelocity: SIMD3(0, 0, -w.rate), angularAcceleration: SIMD3(0, 0, -w.acceleration))
    return try node.followed(by: tilt).followed(by: spin)
  }

  @Test
  func allModelsMatchIndependentCSPICEReference() throws {
    let url = try #require(Bundle.module.url(forResource: "iau-pck00011", withExtension: "json", subdirectory: "Fixtures"))
    let reference = try JSONDecoder().decode(Reference.self, from: Data(contentsOf: url))
    #expect(reference.toolkit == "CSPICE_N0067")
    #expect(reference.kernelSHA256 == "3dff7b1dbeceaa01f25467767d3fa25816051c85d162d1edf04acb310ee28bb1")
    #expect(Set(reference.samples.map(\.body)) == Set(iauRotationModels.keys))
    #expect(reference.samples.count == 73 * 7)
    let frame = Frame(name: "Test", pseudoInertial: true)
    for sample in reference.samples {
      let transform = try mapping(IAURotModel().elements(body: sample.body, tdbDaysSinceJ2000: sample.days))
      for column in 0..<3 {
        var basis = SIMD3<Double>.zero
        basis[column] = 1
        let state = try CartesianState(epoch: .j2000, frame: frame, position: basis, velocity: .zero, acceleration: .zero)
        let result = try transform.applying(to: state, destination: frame)
        for row in 0..<3 {
          let index = 3*row + column
          #expect(abs(result.position[row] - sample.rotation[index]) < 2e-9,
                  "\(sample.body), day \(sample.days), matrix \(row),\(column)")
          #expect(abs(result.velocity[row] - sample.derivative[index]) < 2e-12,
                  "\(sample.body), day \(sample.days), derivative \(row),\(column)")
          #expect(abs(result.acceleration![row] - sample.secondDerivative[index]) < 2e-14,
                  "\(sample.body), day \(sample.days), acceleration \(row),\(column)")
        }
      }
    }
  }

  @Test
  func correctedPhobosAndUpdatedMercuryCoefficients() throws {
    let phobos = try #require(iauRotationModels["phobos"])
    #expect(phobos.primeMeridian.polynomialDegrees[0] == 35.18774440)
    #expect(phobos.primeMeridian.harmonics.last?.amplitudeDegrees == -1.143)
    #expect(phobos.primeMeridian.harmonics.last?.phaseDegrees.count == 3)
    #expect(iauRotationModels["mercury"]?.primeMeridian.polynomialDegrees[0] == 329.5988)
    #expect(iauRotationModels["ida"]?.declination.polynomialDegrees[0] == -87.12)
    #expect(iauRotationModels["moon"]?.primeMeridian.polynomialDegrees[2] == -1.4e-12)
    let tempel = try #require(iauRotationModels["9p_tempel_1"])
    #expect(tempel.referenceJulianDate == 2455607.694660)
    let atEncounter = try IAURotModel().elements(body: "9p_tempel_1",
      tdbDaysSinceJ2000: tempel.referenceJulianDate - 2451545)
    #expect(abs(atEncounter.primeMeridian.radians - 69.2 * .pi/180) < 1e-12)
  }

  @Test
  func provenanceAndInvalidInputAreExplicit() throws {
    #expect(iauRotationModels["earth"]?.provenance == .legacy2009)
    #expect(iauRotationModels["moon"]?.provenance == .legacy2009)
    #expect(iauRotationModels.filter { $0.value.provenance == .iau2015Corrected2019 }.count == 71)
    #expect(throws: IAURotationError.unknownBody) {
      try IAURotModel().elements(body: "unknown", tdbDaysSinceJ2000: 0)
    }
    #expect(throws: IAURotationError.nonFinite) {
      try IAURotModel().elements(body: "mars", tdbDaysSinceJ2000: .nan)
    }
    #expect(throws: IAURotationError.nonFinite) {
      try IAURotModel().elements(body: "phobos", tdbDaysSinceJ2000: Double.greatestFiniteMagnitude)
    }
  }

  @Test
  func quadraticPhasesAndCosineDerivativesFollowIndependentFormula() throws {
    let series = IAUAngleSeries(polynomialDegrees: [1, 2, 3], variable: .days,
      harmonics: [IAUHarmonic(amplitudeDegrees: 4, phaseDegrees: [20, 30, 40], function: .cosine)])
    let d = 36525.0, t = 1.0, rad = Double.pi/180, seconds = 36525.0*86400
    let result = try series.evaluate(days: d)
    let g = (20 + 30*t + 40*t*t)*rad
    let rate = (30 + 80*t)*rad/seconds
    let accel = 80*rad/(seconds*seconds)
    #expect(abs(result.rate - ((2+6*d)/86400 - 4*sin(g)*rate)*rad) < 1e-15)
    #expect(abs(result.acceleration - (6/(86400.0*86400) - 4*(cos(g)*rate*rate+sin(g)*accel))*rad) < 1e-22)
  }
}

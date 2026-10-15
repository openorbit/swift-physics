// SPDX-License-Identifier: Apache-2.0
import Testing
import Math

@Suite
struct PublicAPIContractTests {
  private func requireSendable<T: Sendable>(_ value: T) {}

  @Test
  func mathValuesCanCrossIsolationBoundaries() {
    requireSendable(Quaternion.identity)
    requireSendable(Matrix3x3.identity)
    requireSendable(Frame(name: "root", pseudoInertial: true))
    requireSendable(HourAngle(hours: 1, minutes: 0, seconds: 0))
    requireSendable(AngleOfArc(degrees: 1, minutes: 0, seconds: 0))
  }

  @Test
  func transformationCanBeConstructedFromAnotherModule() {
    let value = Transformation(date: 1, position: SIMD3(1, 2, 3), velocity: .zero,
      acceleration: SIMD3(4, 5, 6), quaternion: .identity,
      rotationRate: .zero, rotationAcceleration: SIMD3(7, 8, 9))
    requireSendable(value)
    #expect(value.acceleration == SIMD3(4, 5, 6))
    #expect(value.rotationAcceleration == SIMD3(7, 8, 9))
    #expect(value.rotationMatrix * SIMD3<Double>(1, 2, 3) == SIMD3(1, 2, 3))
  }

  @Test
  func anglePropertiesIncludeFractionalComponents() {
    let hourAngle = HourAngle(hours: 1, minutes: 30, seconds: 0)
    #expect(hourAngle.degrees == 22.5)
    #expect(abs(hourAngle.radians - .pi / 8) < 1e-14)
    let arc = AngleOfArc(degrees: -1, minutes: 30, seconds: 0)
    #expect(arc.totalDegrees == -1.5)
    #expect(abs(arc.radians + .pi / 120) < 1e-14)
  }
}

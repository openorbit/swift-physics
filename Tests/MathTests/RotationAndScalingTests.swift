// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Math

@Suite
struct RotationAndScalingTests {
  @Test(arguments: [1.0, 7.0, 1e200, 1e-200])
  func axisAngleNormalizesAxis(scale: Double) {
    let rotation = Quaternion(angle: .pi / 2, axis: SIMD3(0, 0, scale))
    #expect(abs(rotation.w - cos(Double.pi / 4)) < 1e-14)
    #expect(abs(Quaternion.dot(lhs: rotation, rhs: rotation) - 1) < 1e-14)
    let rotated = rotation.matrix * SIMD3<Double>(1, 0, 0)
    #expect(abs(rotated.x) < 1e-14)
    #expect(abs(rotated.y - 1) < 1e-14)
    #expect(abs(rotated.z) < 1e-14)
    #expect(abs(rotation.angle - .pi / 2) < 1e-14)
    #expect(rotation.axis == SIMD3<Double>(0, 0, 1))
  }

  @Test
  func extractsAngleAndAxisFromScaledQuaternion() {
    let axis = SIMD3<Double>(1, 2, 3) / sqrt(14)
    let rotation = Quaternion(angle: 1.2, axis: axis) * 5
    #expect(abs(rotation.angle - 1.2) < 1e-14)
    for component in 0..<3 {
      #expect(abs(rotation.axis[component] - axis[component]) < 1e-14)
    }
    #expect(Quaternion.identity.angle == 0)
    #expect(Quaternion.identity.axis == SIMD3<Double>(1, 0, 0))
    let zeroAngle = Quaternion(angle: 0, axis: .zero)
    #expect(zeroAngle.scalar == 1)
    #expect(zeroAngle.vector == .zero)
  }

  @Test
  func conjugationPreservesScalarAndReversesRotation() {
    let value = Quaternion(x: 1, y: -2, z: 3, w: 4)
    #expect(value.conj.vector == -value.vector)
    #expect(value.conj.scalar == value.scalar)
    let product = value * value.conj
    #expect(product.vector == .zero)
    #expect(product.scalar == 30)
    let rotation = Quaternion(angle: 0.8, axis: SIMD3(1, 2, -3))
    let original = SIMD3<Double>(2, -4, 7)
    let recovered = rotation.conj.matrix * (rotation.matrix * original)
    for component in 0..<3 {
      #expect(abs(recovered[component] - original[component]) < 1e-13)
    }
  }

  @Test(arguments: [0.0, 2.0, -3.0, 0.5])
  func scalarMultiplicationScalesEveryEntry(scale: Double) {
    let matrix = Matrix3x3(elements: [1, 2, 3, -4, 5, 6, 7, 8, -9])
    let scaled = matrix * scale
    for row in 0..<3 {
      for column in 0..<3 {
        #expect(scaled[row, column] == matrix[row, column] * scale)
      }
    }
  }
}

// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Extended modified divided differences, in km and TDB seconds.
/// Public format: NAIF SPK Required Reading and SPKW21/SPKE21 API documentation.
/// Implements Newton acceleration polynomials and their two integrals independently
/// of toolkit implementation source. Binary CSPICE is used only as a test oracle.
struct SpiceSPKType21Directory: Sendable {
    let size: Int
    let count: Int
    let dimension: Int
    let epochs: [Double]

    init(file: SpiceDAFFile, array: Int) throws {
        let metadata = file.metadata.arrays[array]
        let words = metadata.finalAddress - metadata.initialAddress + 1
        func malformed(_ reason: String) -> SpiceSPKError {
            .malformedSegment(path: file.url.path, array: array, reason: "type 21: " + reason)
        }
        guard words >= 2 else { throw malformed("missing directory") }
        let tail = try file.readWords(array: array, offset: words - 2, count: 2)
        guard tail[0].isFinite, tail[0] >= 15, tail[0] <= 25,
              tail[0].rounded(.towardZero) == tail[0],
              tail[1].isFinite, tail[1] >= 1, tail[1] <= Double(words),
              tail[1].rounded(.towardZero) == tail[1] else { throw malformed("invalid size/count") }
        // The binary writer stores MAXDIM in the penultimate word, despite the
        // Required Reading diagram labelling it as difference-line size.
        let size = 4 * Int(tail[0]) + 11, count = Int(tail[1])
        guard (size - 11) % 4 == 0, count <= (words - 2) / (size + 1),
              count * (size + 1) + count / 100 + 2 == words else { throw malformed("inconsistent layout") }
        guard size <= file.limits.maximumReadWords, count <= file.limits.maximumReadWords else {
            throw SpiceDAFError.resourceLimit("SPK type 21 record/epoch words")
        }
        let epochs = try file.readWords(array: array, offset: size * count, count: count)
        guard epochs.allSatisfy(\.isFinite), zip(epochs, epochs.dropFirst()).allSatisfy({ $0 < $1 }),
              let descriptor = metadata.spk, epochs.last! >= descriptor.endSecondsTDB else {
            throw malformed("unordered epochs or incomplete coverage")
        }
        if count >= 100 {
            let directory = try file.readWords(array: array, offset: (size + 1) * count, count: count / 100)
            for i in directory.indices {
                guard directory[i] == epochs[(i + 1) * 100 - 1] else { throw malformed("epoch directory mismatch") }
            }
        }
        self.size = size; self.count = count; self.dimension = (size - 11) / 4; self.epochs = epochs
    }

    /// First record ending at or after the query; the earlier record wins a shared endpoint.
    func evaluate(file: SpiceDAFFile, array: Int, seconds: Double) throws -> (Int, SIMD3<Double>, SIMD3<Double>) {
        var lower = 0, upper = count
        while lower < upper {
            let mid = lower + (upper - lower) / 2
            if epochs[mid] < seconds { lower = mid + 1 } else { upper = mid }
        }
        guard lower < count else { throw SpiceSPKError.outOfCoverage }
        let record = try file.readWords(array: array, offset: lower * size, count: size)
        func malformed(_ reason: String) -> SpiceSPKError {
            .malformedSegment(path: file.url.path, array: array, reason: "type 21: " + reason)
        }
        guard record.allSatisfy(\.isFinite) else { throw malformed("nonfinite difference line") }
        func integer(_ x: Double, max: Int) throws -> Int {
            guard x >= 1, x <= Double(max), x.rounded(.towardZero) == x else { throw malformed("invalid integration order") }
            return Int(x)
        }
        let orderOffset = 7 + 4 * dimension
        let orders = try (1...3).map { try integer(record[orderOffset + $0], max: dimension) }
        let maximumOrder = orders.max()!
        let maximumOrderPlusOne = try integer(record[orderOffset], max: dimension + 1)
        guard maximumOrderPlusOne >= max(3, maximumOrder + 1),
              record[1..<maximumOrder].allSatisfy({ $0 != 0 }) else { throw malformed("inconsistent order or zero step") }
        let dt = seconds - record[0]
        guard dt.isFinite else { throw SpiceSPKError.numericalRangeExceeded }
        var position = SIMD3<Double>.zero, velocity = SIMD3<Double>.zero
        for axis in 0..<3 {
            position[axis] = record[1 + dimension + 2 * axis]
            velocity[axis] = record[2 + dimension + 2 * axis]
        }
        if dt == 0 { return (lower, position, velocity) }
        // Acceleration has degree maximumOrder-1. Gauss-Legendre integrates both
        // acceleration and (1-u)*acceleration exactly in real arithmetic with this node count.
        // Evaluate in Newton form to avoid converting to ill-conditioned monomial coefficients.
        let quadrature = spk21Quadratures[(maximumOrder + 2) / 2 - 1]
        var velocityIntegral = SIMD3<Double>.zero, positionIntegral = SIMD3<Double>.zero
        for (u, weight) in quadrature {
            var acceleration = SIMD3<Double>.zero, basis = 1.0
            for j in 0..<maximumOrder {
                for axis in 0..<3 where j < orders[axis] {
                    acceleration[axis] += record[7 + dimension + axis * dimension + j] * basis
                }
                if j + 1 < maximumOrder {
                    let shift = j == 0 ? 0 : record[j]
                    basis *= (dt * u + shift) / record[j + 1]
                }
            }
            velocityIntegral += acceleration * weight
            positionIntegral += acceleration * (weight * (1 - u))
        }
        position += dt * (velocity + dt * positionIntegral)
        velocity += dt * velocityIntegral
        guard (0..<3).allSatisfy({ position[$0].isFinite && velocity[$0].isFinite }) else {
            throw SpiceSPKError.numericalRangeExceeded
        }
        return (lower, position, velocity)
    }
}

/// Orders 1...13, enough for the NAIF N0067 maximum difference table dimension of 25.
private let spk21Quadratures: [[(Double, Double)]] = (1...13).map { count in
    var result = [(Double, Double)]()
    for i in 0..<count {
        var x = cos(.pi * (Double(i) + 0.75) / (Double(count) + 0.5))
        var derivative = 0.0
        for _ in 0..<32 {
            var p = 1.0, previous = 0.0
            for j in 1...count {
                let next = ((2 * Double(j) - 1) * x * p - (Double(j) - 1) * previous) / Double(j)
                previous = p; p = next
            }
            derivative = Double(count) * (x * p - previous) / (x * x - 1)
            let step = p / derivative
            x -= step
            if abs(step) < 2e-16 { break }
        }
        result.append(((x + 1) / 2, 1 / ((1 - x * x) * derivative * derivative)))
    }
    return result
}

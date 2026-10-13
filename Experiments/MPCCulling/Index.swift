// SPDX-License-Identifier: Apache-2.0
import Foundation
import CelMek

struct Cell: Hashable { let x, y, z: Int }
struct Bounds {
    var low = SIMD3<Double>(repeating: .infinity)
    var high = SIMD3<Double>(repeating: -.infinity)
    mutating func include(_ p: SIMD3<Double>, padding: Double) {
        for axis in 0..<3 { low[axis] = min(low[axis], p[axis] - padding); high[axis] = max(high[axis], p[axis] + padding) }
    }
    func intersects(_ p: SIMD3<Double>, radius: Double) -> Bool {
        var squared = 0.0
        for axis in 0..<3 {
            let d = max(low[axis] - p[axis], max(0, p[axis] - high[axis]))
            squared += d * d
        }
        return squared <= radius * radius
    }
}
struct Bucket {
    var indices = [Int]()
    var bounds = Bounds()
    var radialLow = Double.infinity, radialHigh = -Double.infinity
}
struct CullingIndex {
    let centerEpoch: Epoch
    let halfWindowSeconds: Double
    let buckets: [Bucket]
    let radialBounds: [(Double, Double)]
    init(records: [OrbitRecord], epoch: Epoch, halfWindowSeconds: Double, cellWidth: Double) throws {
        guard halfWindowSeconds.isFinite, halfWindowSeconds >= 0, cellWidth.isFinite, cellWidth > 0 else { throw ExperimentError.invalidQuery }
        centerEpoch = epoch; self.halfWindowSeconds = halfWindowSeconds
        var cells = [Cell: Bucket](), radial = [(Double, Double)]()
        radial.reserveCapacity(records.count)
        for (index, record) in records.enumerated() {
            let position = try record.orbit.state(at: epoch).position
            // Perihelion speed is a global speed bound for elliptic, parabolic and hyperbolic two-body motion.
            // Numerical guard is experimental, not a certified floating-point error bound.
            let padding = record.maximumSpeed * halfWindowSeconds + 1000 + norm(position) * 1e-12
            let r = norm(position)
            let lower = max(0, max(min(record.perihelion, r) - 1000, r - padding))
            let upper = min(max(record.aphelion, r) + 1000, r + padding)
            radial.append((lower, upper))
            let coordinates = (0..<3).map { floor(position[$0] / cellWidth) }
            guard coordinates.allSatisfy({ $0.isFinite && $0 > Double(Int.min) && $0 < Double(Int.max) }) else { throw ExperimentError.invalidQuery }
            let cell = Cell(x: Int(coordinates[0]), y: Int(coordinates[1]), z: Int(coordinates[2]))
            // Mutate through Dictionary's default subscript to avoid copying the bucket's array on every insertion.
            cells[cell, default: Bucket()].indices.append(index)
            cells[cell, default: Bucket()].bounds.include(position, padding: padding)
            cells[cell, default: Bucket()].radialLow = min(cells[cell]!.radialLow, lower)
            cells[cell, default: Bucket()].radialHigh = max(cells[cell]!.radialHigh, upper)
        }
        buckets = Array(cells.values); radialBounds = radial
    }
    func candidates(observer: SIMD3<Double>, radius: Double, epoch: Epoch, spatial: Bool) throws -> [Int] {
        guard radius.isFinite, radius >= 0, (0..<3).allSatisfy({ observer[$0].isFinite }) else { throw ExperimentError.invalidQuery }
        guard abs(try epoch.duration(since: centerEpoch).seconds) <= halfWindowSeconds else { throw ExperimentError.outsideIndexWindow }
        let r = norm(observer), low = max(0, r - radius), high = r + radius
        var selected = [Int]()
        for bucket in buckets {
            guard bucket.radialLow <= high, bucket.radialHigh >= low else { continue }
            if spatial && !bucket.bounds.intersects(observer, radius: radius) { continue }
            for i in bucket.indices {
                let bounds = radialBounds[i]
                if bounds.0 <= high && bounds.1 >= low { selected.append(i) }
            }
        }
        return selected
    }
}
func norm(_ p: SIMD3<Double>) -> Double { hypot(hypot(p.x, p.y), p.z) }
func matches(_ records: [OrbitRecord], indices: [Int], observer: SIMD3<Double>, radius: Double, epoch: Epoch) throws -> Set<Int> {
    var result = Set<Int>()
    for i in indices {
        if norm(try records[i].orbit.state(at: epoch).position - observer) <= radius { result.insert(i) }
    }
    return result
}

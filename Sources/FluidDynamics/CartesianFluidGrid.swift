// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Shared uniform Cartesian geometry. X varies fastest; cells are cubic, dimensions
/// may be rectangular. Solver-specific resolution/coverage restrictions are additional.
public struct CartesianFluidGrid: Sendable {
    public let dimensions: SIMD3<Int>
    public let minimumCornerMetres: SIMD3<Double>
    public let cellSizeMetres: Double
    public var cellCount: Int { dimensions.x*dimensions.y*dimensions.z }
    public var cellVolumeCubicMetres: Double { cellSizeMetres*cellSizeMetres*cellSizeMetres }
    public init(dimensions: SIMD3<Int> = SIMD3(repeating: 16), minimumCornerMetres: SIMD3<Double>, cellSizeMetres: Double) throws {
        guard [dimensions.x,dimensions.y,dimensions.z].allSatisfy({ (1...512).contains($0) }),
              dimensions.x*dimensions.y*dimensions.z <= 262144,
              fluidFinite(minimumCornerMetres), cellSizeMetres.isFinite, cellSizeMetres > 0,
              (cellSizeMetres*cellSizeMetres*cellSizeMetres).isFinite,
              cellSizeMetres*cellSizeMetres*cellSizeMetres > 0 else { throw FluidError.invalidGeometry }
        let extent = SIMD3(Double(dimensions.x),Double(dimensions.y),Double(dimensions.z))*cellSizeMetres
        guard fluidFinite(minimumCornerMetres+extent),
              (0..<3).allSatisfy({ minimumCornerMetres[$0]+cellSizeMetres > minimumCornerMetres[$0] }) else { throw FluidError.invalidGeometry }
        self.dimensions = dimensions; self.minimumCornerMetres = minimumCornerMetres; self.cellSizeMetres = cellSizeMetres
    }
    package func index(_ q: SIMD3<Int>) -> Int { q.x+dimensions.x*(q.y+dimensions.y*q.z) }
    package func centre(_ q: SIMD3<Int>) -> SIMD3<Double> {
        minimumCornerMetres+SIMD3(Double(q.x)+0.5,Double(q.y)+0.5,Double(q.z)+0.5)*cellSizeMetres
    }
    public func cellCentreMetres(index: Int) throws -> SIMD3<Double> {
        guard (0..<cellCount).contains(index) else { throw FluidError.invalidParameter }
        return centre(SIMD3(index%dimensions.x,(index/dimensions.x)%dimensions.y,index/(dimensions.x*dimensions.y)))
    }
    public func cellIndex(at positionMetres: SIMD3<Double>) throws -> Int {
        guard fluidFinite(positionMetres) else { throw FluidError.invalidParameter }
        let q = (positionMetres-minimumCornerMetres)/cellSizeMetres, n = dimensions
        guard fluidFinite(q), q.x >= 0, q.y >= 0, q.z >= 0, q.x <= Double(n.x), q.y <= Double(n.y), q.z <= Double(n.z) else {
            throw FluidError.outsideCoverage
        }
        return index(SIMD3(min(n.x-1,Int(q.x)),min(n.y-1,Int(q.y)),min(n.z-1,Int(q.z))))
    }
}

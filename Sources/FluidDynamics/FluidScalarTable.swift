// SPDX-License-Identifier: Apache-2.0
import Foundation

/// One strictly ordered interpolation coordinate, with explicit name and units.
public struct FluidTableAxis: Sendable, Codable, Equatable {
    public let name: String
    public let units: String
    public let coordinates: [Double]
    public init(name: String, units: String, coordinates: [Double]) throws {
        guard !name.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty, !units.isEmpty,
              !coordinates.isEmpty, coordinates.allSatisfy(\.isFinite),
              zip(coordinates,coordinates.dropFirst()).allSatisfy({ $0 < $1 && ($1-$0).isFinite }) else {
            throw FluidError.malformedData
        }
        self.name = name; self.units = units; self.coordinates = coordinates
    }
    private enum CodingKeys: String, CodingKey { case name, units, coordinates }
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        try self.init(name: c.decode(String.self,forKey: .name),units: c.decode(String.self,forKey: .units),
                      coordinates: c.decode([Double].self,forKey: .coordinates))
    }
}

/// Checked rectangular scalar data. Last axis varies fastest; interpolation is
/// linear in the supplied coordinates (use an explicit log coordinate if intended).
/// No extrapolation, automatic blending of regimes or inferred scientific validity.
public struct FluidScalarTable: Sendable, Codable {
    public let axes: [FluidTableAxis]
    public let values: [Double]
    public let valueUnits: String
    public let source: String
    public let version: String
    public let assumptions: String
    public init(axes: [FluidTableAxis], values: [Double], valueUnits: String,
                source: String, version: String, assumptions: String) throws {
        guard (1...8).contains(axes.count), Set(axes.map(\.name)).count == axes.count,
              [valueUnits,source,version,assumptions].allSatisfy({ !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }),
              values.allSatisfy(\.isFinite) else { throw FluidError.malformedData }
        var count = 1
        for axis in axes {
            guard axis.coordinates.count <= 1_000_000/count else { throw FluidError.malformedData }
            count *= axis.coordinates.count
        }
        guard values.count == count else { throw FluidError.incompleteGrid }
        self.axes = axes; self.values = values; self.valueUnits = valueUnits
        self.source = source; self.version = version; self.assumptions = assumptions
    }
    /// O(sum(axis lengths) + 2^dimensions), at most eight dimensions.
    public func value(at coordinates: [Double]) throws -> Double {
        guard coordinates.count == axes.count, coordinates.allSatisfy(\.isFinite) else {
            throw FluidError.invalidParameter
        }
        var corners: [(index: Int, weight: Double)] = [(0,1)]
        for (axis,x) in zip(axes,coordinates) {
            let grid = axis.coordinates
            guard x >= grid[0], x <= grid[grid.count-1] else { throw FluidError.outsideCoverage }
            let lower = grid.lastIndex(where: { $0 <= x })!
            let upper = min(lower+1,grid.count-1)
            let weight = lower == upper ? 0 : (x-grid[lower])/(grid[upper]-grid[lower])
            guard weight.isFinite else { throw FluidError.nonFiniteResult }
            var next: [(Int,Double)] = []
            for corner in corners {
                if weight < 1 { next.append((corner.index*grid.count+lower,corner.weight*(1-weight))) }
                if weight > 0 { next.append((corner.index*grid.count+upper,corner.weight*weight)) }
            }
            corners = next
        }
        let result = corners.reduce(0.0) { $0+values[$1.index]*$1.weight }
        guard result.isFinite else { throw FluidError.nonFiniteResult }
        return result
    }
    private enum CodingKeys: String, CodingKey { case axes, values, valueUnits, source, version, assumptions }
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        try self.init(axes: c.decode([FluidTableAxis].self,forKey: .axes),values: c.decode([Double].self,forKey: .values),
            valueUnits: c.decode(String.self,forKey: .valueUnits),source: c.decode(String.self,forKey: .source),
            version: c.decode(String.self,forKey: .version),assumptions: c.decode(String.self,forKey: .assumptions))
    }
}

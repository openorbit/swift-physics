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
import FluidDynamics
import Foundation
import Yams

/// Explicit boundary behavior; no automatic physical extrapolation beyond measured samples.
public enum CoefficientBoundaryPolicy: String, Codable, Sendable { case reject, clamp }

/// One measured/computed angle sample. Angles are radians and moment reference is dataset-defined.
public struct AirfoilSample: Sendable, Codable {
    public let angleOfAttackRadians: Double
    public let coefficients: AirfoilCoefficients
    public init(angleOfAttackRadians: Double, coefficients: AirfoilCoefficients) {
        self.angleOfAttackRadians = angleOfAttackRadians; self.coefficients = coefficients
    }
}

/// Validated section polar at one Reynolds/Mach condition, preserving irregular angle spacing.
public struct AirfoilPolar: Sendable, Codable {
    public let reynoldsNumber: Double
    public let machNumber: Double
    public let samples: [AirfoilSample]
    public var angleRangeRadians: ClosedRange<Double> { samples.first!.angleOfAttackRadians...samples.last!.angleOfAttackRadians }
    public init(reynoldsNumber: Double, machNumber: Double = 0, samples: [AirfoilSample]) throws {
        guard reynoldsNumber.isFinite, reynoldsNumber > 0, machNumber.isFinite, machNumber >= 0,
              !samples.isEmpty, samples.allSatisfy({ $0.angleOfAttackRadians.isFinite && abs($0.angleOfAttackRadians) <= .pi && $0.coefficients.isFinite }),
              zip(samples,samples.dropFirst()).allSatisfy({ $0.angleOfAttackRadians < $1.angleOfAttackRadians }) else {
            throw FluidError.invalidCoefficients
        }
        self.reynoldsNumber = reynoldsNumber; self.machNumber = machNumber; self.samples = samples
    }
    public func coefficients(at angle: Double, boundaryPolicy: CoefficientBoundaryPolicy = .reject) throws -> AirfoilCoefficients {
        let a = try boundedCoefficientCoordinate(angle, range: angleRangeRadians, policy: boundaryPolicy)
        let lower = samples.lastIndex(where: { $0.angleOfAttackRadians <= a })!
        guard lower+1 < samples.count else { return samples[lower].coefficients }
        let first = samples[lower], next = samples[lower+1]
        let weight = (a-first.angleOfAttackRadians)/(next.angleOfAttackRadians-first.angleOfAttackRadians)
        return try blendCoefficients(first.coefficients,next.coefficients,weight: weight)
    }
    private enum CodingKeys: String, CodingKey { case reynoldsNumber, machNumber, samples }
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        try self.init(reynoldsNumber: c.decode(Double.self,forKey: .reynoldsNumber),
                      machNumber: c.decode(Double.self,forKey: .machNumber),samples: c.decode([AirfoilSample].self,forKey: .samples))
    }
}

/// Checked section model interface, independent of motion integration and vehicle systems.
public protocol AirfoilCoefficientModel: Sendable {
    func coefficients(angleOfAttackRadians: Double, reynoldsNumber: Double, machNumber: Double) throws -> AirfoilCoefficients
}

/// Rectangular Reynolds/Mach grid of polars, with linear interpolation in each coordinate.
/// Angle grids may differ; all contributing polars must cover a query unless clamping is selected.
public struct TabulatedAirfoil: AirfoilCoefficientModel, Sendable, Codable {
    public let polars: [AirfoilPolar]
    public let boundaryPolicy: CoefficientBoundaryPolicy
    public let source: String
    public let version: String
    private let reynoldsGrid: [Double]
    private let machGrid: [Double]
    private let grid: [AirfoilPolar]
    public init(polars: [AirfoilPolar], boundaryPolicy: CoefficientBoundaryPolicy = .reject,
                source: String, version: String) throws {
        guard !polars.isEmpty, !source.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              !version.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else { throw FluidError.invalidParameter }
        let re = Array(Set(polars.map(\.reynoldsNumber))).sorted(), ma = Array(Set(polars.map(\.machNumber))).sorted()
        var ordered = [AirfoilPolar]()
        for m in ma {
            for r in re {
                let matching = polars.filter { $0.reynoldsNumber == r && $0.machNumber == m }
                guard matching.count <= 1 else { throw FluidError.duplicateCondition }
                guard let polar = matching.first else { throw FluidError.incompleteGrid }
                ordered.append(polar)
            }
        }
        self.polars = polars; self.boundaryPolicy = boundaryPolicy; self.source = source; self.version = version
        reynoldsGrid = re; machGrid = ma; grid = ordered
    }
    public func coefficients(angleOfAttackRadians: Double, reynoldsNumber: Double, machNumber: Double) throws -> AirfoilCoefficients {
        guard reynoldsNumber.isFinite, reynoldsNumber > 0, machNumber.isFinite, machNumber >= 0 else {
            throw FluidError.invalidParameter
        }
        let re = try bracket(reynoldsNumber,grid: reynoldsGrid), ma = try bracket(machNumber,grid: machGrid)
        func evaluate(_ r: Int, _ m: Int) throws -> AirfoilCoefficients {
            try grid[m*reynoldsGrid.count+r].coefficients(at: angleOfAttackRadians,boundaryPolicy: boundaryPolicy)
        }
        func atMach(_ m: Int) throws -> AirfoilCoefficients {
            let first = try evaluate(re.0,m)
            guard re.0 != re.1 else { return first }
            return try blendCoefficients(first,evaluate(re.1,m),weight: re.2)
        }
        let first = try atMach(ma.0)
        guard ma.0 != ma.1 else { return first }
        return try blendCoefficients(first,atMach(ma.1),weight: ma.2)
    }
    private func bracket(_ value: Double, grid: [Double]) throws -> (Int,Int,Double) {
        let v = try boundedCoefficientCoordinate(value,range: grid.first!...grid.last!,policy: boundaryPolicy)
        let lower = grid.lastIndex(where: { $0 <= v })!
        guard lower+1 < grid.count, v != grid[lower] else { return (lower,lower,0) }
        return (lower,lower+1,(v-grid[lower])/(grid[lower+1]-grid[lower]))
    }
    private enum CodingKeys: String, CodingKey { case polars, boundaryPolicy, source, version }
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        try self.init(polars: c.decode([AirfoilPolar].self,forKey: .polars),
                      boundaryPolicy: c.decode(CoefficientBoundaryPolicy.self,forKey: .boundaryPolicy),
                      source: c.decode(String.self,forKey: .source),version: c.decode(String.self,forKey: .version))
    }
    public func encode(to encoder: Encoder) throws {
        var c = encoder.container(keyedBy: CodingKeys.self)
        try c.encode(polars,forKey: .polars); try c.encode(boundaryPolicy,forKey: .boundaryPolicy)
        try c.encode(source,forKey: .source); try c.encode(version,forKey: .version)
    }
}

/// Loader for the existing XFOIL-style YAML: degrees, CL, CD, CDp, CM.
/// Every `data` entry is retained, with optional `mach-number` (default zero).
public enum AirfoilTableLoader {
    public static func yaml(_ text: String, boundaryPolicy: CoefficientBoundaryPolicy = .reject,
                            source: String = "Caller-supplied airfoil YAML", version: String = "1") throws -> TabulatedAirfoil {
        let file: YAMLAirfoilFile
        do { file = try YAMLDecoder().decode(YAMLAirfoilFile.self,from: text) }
        catch { throw FluidError.malformedData }
        guard !file.data.isEmpty else { throw FluidError.malformedData }
        var polars = [AirfoilPolar]()
        for entry in file.data {
            let re = entry.reynoldsNumber, mach = entry.machNumber ?? 0, rows = entry.coefficients
            guard !rows.isEmpty else { throw FluidError.malformedData }
            var samples = [AirfoilSample]()
            for row in rows {
                guard row.count == 5, row.allSatisfy(\.isFinite), (-180.0...180.0).contains(row[0]) else {
                    throw FluidError.malformedData
                }
                samples.append(AirfoilSample(angleOfAttackRadians: row[0] * .pi/180,
                    coefficients: AirfoilCoefficients(lift: row[1],drag: row[2],moment: row[4])))
            }
            polars.append(try AirfoilPolar(reynoldsNumber: re,machNumber: mach,samples: samples))
        }
        return try TabulatedAirfoil(polars: polars,boundaryPolicy: boundaryPolicy,source: source,version: version)
    }
    public static func yamlFile(_ url: URL, boundaryPolicy: CoefficientBoundaryPolicy = .reject) throws -> TabulatedAirfoil {
        try yaml(String(contentsOf: url,encoding: .utf8),boundaryPolicy: boundaryPolicy,source: url.absoluteString)
    }
}

internal func boundedCoefficientCoordinate(_ value: Double, range: ClosedRange<Double>, policy: CoefficientBoundaryPolicy) throws -> Double {
    guard value.isFinite else { throw FluidError.invalidParameter }
    guard policy == .clamp || range.contains(value) else { throw FluidError.outsideCoverage }
    return min(range.upperBound,max(range.lowerBound,value))
}
internal func blendCoefficients(_ first: AirfoilCoefficients, _ last: AirfoilCoefficients, weight: Double) throws -> AirfoilCoefficients {
    let result = AirfoilCoefficients(lift: first.lift*(1-weight)+last.lift*weight,
        drag: first.drag*(1-weight)+last.drag*weight,moment: first.moment*(1-weight)+last.moment*weight)
    guard result.isFinite else { throw FluidError.nonFiniteResult }
    return result
}

private struct YAMLAirfoilFile: Decodable { let data: [YAMLAirfoilEntry] }
private struct YAMLAirfoilEntry: Decodable {
    let reynoldsNumber: Double
    let machNumber: Double?
    let coefficients: [[Double]]
    private enum CodingKeys: String, CodingKey {
        case reynoldsNumber = "reynolds-number", machNumber = "mach-number", coefficients
    }
}

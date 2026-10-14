// SPDX-License-Identifier: Apache-2.0
import Foundation

public struct DragCoefficientQuery: Sendable {
    public let flow: AerodynamicFlow
    public let referenceLengthMetres: Double
    public let axisAngleRadians: Double?
    public init(flow: AerodynamicFlow, referenceLengthMetres: Double, axisAngleRadians: Double? = nil) throws {
        guard referenceLengthMetres.isFinite, referenceLengthMetres > 0,
              axisAngleRadians.map({ $0.isFinite && (0...Double.pi).contains($0) }) ?? true else {
            throw AerodynamicError.invalidParameter
        }
        self.flow = flow; self.referenceLengthMetres = referenceLengthMetres; self.axisAngleRadians = axisAngleRadians
    }
    public func reynoldsNumber() throws -> Double? { try flow.reynoldsNumber(referenceLengthMetres: referenceLengthMetres) }
}
public protocol DragCoefficientModel: Sendable {
    func dragCoefficient(for query: DragCoefficientQuery) throws -> Double
}
public struct DragCoefficientSample: Sendable {
    public let machNumber: Double
    public let reynoldsNumber: Double
    public let dragCoefficient: Double
    public init(machNumber: Double, reynoldsNumber: Double, dragCoefficient: Double) {
        self.machNumber = machNumber; self.reynoldsNumber = reynoldsNumber; self.dragCoefficient = dragCoefficient
    }
}

/// Rectangular Mach/Reynolds grid with linear interpolation and reject-only coverage.
/// Orientation-invariant when maximumAxisAngleRadians is nil. No default real canopy
/// or balloon dataset is implied; area definitions must match the supplied coefficients.
public struct TabulatedDragCoefficients: DragCoefficientModel {
    public let source: String
    public let version: String
    public let maximumAxisAngleRadians: Double?
    public let machNumbers: [Double]
    public let reynoldsNumbers: [Double]
    private let values: [[Double]]
    public init(samples: [DragCoefficientSample], source: String, version: String,
                maximumAxisAngleRadians: Double? = nil) throws {
        guard !samples.isEmpty, [source,version].allSatisfy({ !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }),
              maximumAxisAngleRadians.map({ $0.isFinite && (0...Double.pi).contains($0) }) ?? true,
              samples.allSatisfy({ $0.machNumber.isFinite && $0.machNumber >= 0 && $0.reynoldsNumber.isFinite
                  && $0.reynoldsNumber > 0 && $0.dragCoefficient.isFinite && $0.dragCoefficient >= 0 }) else {
            throw AerodynamicError.invalidCoefficients
        }
        machNumbers = Array(Set(samples.map(\.machNumber))).sorted()
        reynoldsNumbers = Array(Set(samples.map(\.reynoldsNumber))).sorted()
        var grid = [[Double]]()
        for mach in machNumbers {
            var row = [Double]()
            for re in reynoldsNumbers {
                let matching = samples.filter { $0.machNumber == mach && $0.reynoldsNumber == re }
                guard matching.count <= 1 else { throw AerodynamicError.duplicateCondition }
                guard let sample = matching.first else { throw AerodynamicError.incompleteGrid }
                row.append(sample.dragCoefficient)
            }
            grid.append(row)
        }
        self.source = source; self.version = version; self.maximumAxisAngleRadians = maximumAxisAngleRadians; values = grid
    }
    private func bracket(_ value: Double, grid: [Double]) throws -> (Int,Int,Double) {
        guard value.isFinite, value >= grid[0], value <= grid.last! else { throw AerodynamicError.outsideCoverage }
        if grid.count == 1 { return (0,0,0) }
        let upper = grid.firstIndex { $0 >= value }!
        if upper == 0 { return (0,0,0) }
        let lower = upper-1
        return (lower,upper,(value-grid[lower])/(grid[upper]-grid[lower]))
    }
    public func dragCoefficient(for query: DragCoefficientQuery) throws -> Double {
        guard let mach = query.flow.machNumber, let re = try query.reynoldsNumber() else { throw AerodynamicError.missingFlowProperty }
        if let maximumAxisAngleRadians {
            guard let angle = query.axisAngleRadians else { throw AerodynamicError.missingFlowProperty }
            guard angle <= maximumAxisAngleRadians else { throw AerodynamicError.outsideCoverage }
        }
        let (m0,m1,mt) = try bracket(mach, grid: machNumbers), (r0,r1,rt) = try bracket(re, grid: reynoldsNumbers)
        let low = values[m0][r0]*(1-rt)+values[m0][r1]*rt
        let high = values[m1][r0]*(1-rt)+values[m1][r1]*rt
        return low*(1-mt)+high*mt
    }
}

/// Steady drag opposing atmosphere-relative motion. Supplied area already represents
/// current shape, reefing/inflation or wake corrections; this evaluator evolves none.
public struct DragAreaAerodynamics: Sendable {
    public let areaSquareMetres: Double
    public let referenceLengthMetres: Double
    public let forceApplicationPointMetres: SIMD3<Double>
    public let axis: SIMD3<Double>?
    public let model: any DragCoefficientModel
    public init(areaSquareMetres: Double, referenceLengthMetres: Double,
                forceApplicationPointMetres: SIMD3<Double> = .zero, axis: SIMD3<Double>? = nil,
                model: any DragCoefficientModel) throws {
        guard areaSquareMetres.isFinite, areaSquareMetres >= 0, referenceLengthMetres.isFinite,
              referenceLengthMetres > 0, aeroFinite(forceApplicationPointMetres) else { throw AerodynamicError.invalidGeometry }
        if let axis {
            let length = aeroNorm(axis)
            guard aeroFinite(axis), length.isFinite, length > 0 else { throw AerodynamicError.invalidGeometry }
            self.axis = axis/length
        } else { self.axis = nil }
        self.areaSquareMetres = areaSquareMetres; self.referenceLengthMetres = referenceLengthMetres
        self.forceApplicationPointMetres = forceApplicationPointMetres; self.model = model
    }
    public func evaluate(flow: AerodynamicFlow, momentReferencePointMetres: SIMD3<Double> = .zero) throws -> AerodynamicLoads {
        var force = SIMD3<Double>.zero
        if flow.dynamicPressurePascals > 0 && areaSquareMetres > 0 {
            let direction = flow.relativeVelocityMetresPerSecond/flow.speedMetresPerSecond
            let angle = axis.map { acos(max(-1,min(1,($0*direction).sum()))) }
            let coefficient = try model.dragCoefficient(for: DragCoefficientQuery(flow: flow,
                referenceLengthMetres: referenceLengthMetres, axisAngleRadians: angle))
            guard coefficient.isFinite, coefficient >= 0 else { throw AerodynamicError.invalidCoefficients }
            force = -direction*((flow.dynamicPressurePascals*areaSquareMetres)*coefficient)
        }
        return try AerodynamicLoads(forceNewtons: force, momentNewtonMetres: .zero,
            referencePointMetres: forceApplicationPointMetres).referred(to: momentReferencePointMetres)
    }
}

/// Current-canopy steady drag. No opening loads, tether transmission or deployment.
public struct ParachuteAerodynamics: Sendable {
    public let drag: DragAreaAerodynamics
    public init(drag: DragAreaAerodynamics) { self.drag = drag }
    public func evaluate(flow: AerodynamicFlow, momentReferencePointMetres: SIMD3<Double> = .zero) throws -> AerodynamicLoads {
        try drag.evaluate(flow: flow, momentReferencePointMetres: momentReferencePointMetres)
    }
}

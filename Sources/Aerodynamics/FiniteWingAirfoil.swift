// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Foundation

public enum FiniteWingError: Error, Equatable { case unbracketedSolution, iterationLimit }

public struct FiniteWingSolution: Sendable {
    public let coefficients: AirfoilCoefficients
    public let effectiveAngleOfAttackRadians: Double
    public let inducedAngleRadians: Double
    public let inducedDragCoefficient: Double
    public let angleResidualRadians: Double
    public let iterations: Int
}

/// Equivalent-wing lifting-line approximation for steady, attached, low-speed flow.
/// Solves alpha = alphaEffective + CL(alphaEffective)/(pi e AR), then adds
/// CDi = CL²/(pi e AR). The input polar must contain profile drag only.
/// This is a uniform-wing approximation, not a spanwise wake/interference solver.
public struct FiniteWingAirfoil: AirfoilCoefficientModel {
    public let airfoil: any AirfoilCoefficientModel
    public let aspectRatio: Double
    public let spanEfficiency: Double
    public let effectiveAngleRangeRadians: ClosedRange<Double>
    public let operatingAngleRangeRadians: ClosedRange<Double>
    public let angleToleranceRadians: Double
    public let maximumIterations: Int
    public let maximumMachNumber: Double
    private let inducedFactor: Double

    /// The caller declares an attached-flow effective-angle bracket covered by the
    /// polar. This approximation restricts AR>=4, Mach<=0.3 and incidence to ±0.35 rad.
    /// A calibrated span efficiency is explicit; it is not inferred from segmentation.
    public init(airfoil: any AirfoilCoefficientModel, aspectRatio: Double, spanEfficiency: Double,
                effectiveAngleRangeRadians: ClosedRange<Double>,
                operatingAngleRangeRadians: ClosedRange<Double>,
                maximumMachNumber: Double = 0.3, angleToleranceRadians: Double = 1e-10,
                maximumIterations: Int = 64) throws {
        guard aspectRatio.isFinite, aspectRatio >= 4, spanEfficiency.isFinite,
              spanEfficiency > 0, spanEfficiency <= 1,
              maximumMachNumber.isFinite, maximumMachNumber >= 0, maximumMachNumber <= 0.3,
              angleToleranceRadians.isFinite, angleToleranceRadians > 0,
              maximumIterations > 0, maximumIterations <= 256,
              [effectiveAngleRangeRadians.lowerBound,effectiveAngleRangeRadians.upperBound,
               operatingAngleRangeRadians.lowerBound,operatingAngleRangeRadians.upperBound].allSatisfy({ $0.isFinite && abs($0) <= 0.35 }),
              effectiveAngleRangeRadians.lowerBound < effectiveAngleRangeRadians.upperBound,
              operatingAngleRangeRadians.lowerBound < operatingAngleRangeRadians.upperBound else {
            throw FluidError.invalidParameter
        }
        self.airfoil = airfoil; self.aspectRatio = aspectRatio; self.spanEfficiency = spanEfficiency
        self.effectiveAngleRangeRadians = effectiveAngleRangeRadians; self.operatingAngleRangeRadians = operatingAngleRangeRadians
        self.maximumMachNumber = maximumMachNumber; self.angleToleranceRadians = angleToleranceRadians
        self.maximumIterations = maximumIterations
        inducedFactor = (1/aspectRatio)/spanEfficiency/Double.pi
        guard inducedFactor.isFinite, inducedFactor > 0 else { throw FluidError.invalidParameter }
    }
    public func coefficients(angleOfAttackRadians: Double, reynoldsNumber: Double, machNumber: Double) throws -> AirfoilCoefficients {
        try solve(angleOfAttackRadians: angleOfAttackRadians,reynoldsNumber: reynoldsNumber,machNumber: machNumber).coefficients
    }
    /// Bisection requires a bracketed root, not a promise of a unique solution.
    /// Use a monotone attached-flow polar; post-stall branches are outside this model.
    public func solve(angleOfAttackRadians alpha: Double, reynoldsNumber: Double, machNumber: Double) throws -> FiniteWingSolution {
        guard alpha.isFinite, reynoldsNumber.isFinite, reynoldsNumber > 0, machNumber.isFinite, machNumber >= 0 else {
            throw FluidError.invalidParameter
        }
        guard operatingAngleRangeRadians.contains(alpha), machNumber <= maximumMachNumber else {
            throw FluidError.outsideCoverage
        }
        func evaluate(_ effective: Double) throws -> (AirfoilCoefficients,Double) {
            let coefficients = try airfoil.coefficients(angleOfAttackRadians: effective,reynoldsNumber: reynoldsNumber,machNumber: machNumber)
            guard coefficients.isFinite, coefficients.drag >= 0 else { throw FluidError.invalidCoefficients }
            let residual = effective+coefficients.lift*inducedFactor-alpha
            guard residual.isFinite else { throw FluidError.nonFiniteResult }
            return (coefficients,residual)
        }
        func result(_ effective: Double, _ value: (AirfoilCoefficients,Double), _ iterations: Int) throws -> FiniteWingSolution {
            let inducedAngle = value.0.lift*inducedFactor, inducedDrag = value.0.lift*inducedAngle
            let coefficients = AirfoilCoefficients(lift: value.0.lift,drag: value.0.drag+inducedDrag,moment: value.0.moment)
            guard inducedDrag.isFinite, coefficients.isFinite else { throw FluidError.nonFiniteResult }
            return FiniteWingSolution(coefficients: coefficients,effectiveAngleOfAttackRadians: effective,
                inducedAngleRadians: inducedAngle,inducedDragCoefficient: inducedDrag,
                angleResidualRadians: value.1,iterations: iterations)
        }
        var lower = effectiveAngleRangeRadians.lowerBound, upper = effectiveAngleRangeRadians.upperBound
        var lowerValue = try evaluate(lower)
        let upperValue = try evaluate(upper)
        if abs(lowerValue.1) <= angleToleranceRadians { return try result(lower,lowerValue,0) }
        if abs(upperValue.1) <= angleToleranceRadians { return try result(upper,upperValue,0) }
        guard (lowerValue.1 < 0) != (upperValue.1 < 0) else { throw FiniteWingError.unbracketedSolution }
        for iteration in 1...maximumIterations {
            let effective = lower+(upper-lower)*0.5
            let value = try evaluate(effective)
            if abs(value.1) <= angleToleranceRadians { return try result(effective,value,iteration) }
            if (lowerValue.1 < 0) == (value.1 < 0) { lower = effective; lowerValue = value }
            else { upper = effective }
        }
        throw FiniteWingError.iterationLimit
    }
}

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

/// Supplied linear attached-flow polar at one span station; profile drag excludes
/// induced drag. Positive geometric incidence raises the leading edge.
public struct LiftingLineSection: Sendable {
    public let chordMetres: Double
    public let liftSlopePerRadian: Double
    public let zeroLiftAngleRadians: Double
    public let geometricIncidenceRadians: Double
    public let profileDragCoefficient: Double
    public let pitchingMomentCoefficient: Double
    public let effectiveAngleRangeRadians: ClosedRange<Double>
    public let reynoldsRange: ClosedRange<Double>
    public init(chordMetres: Double, liftSlopePerRadian: Double, zeroLiftAngleRadians: Double = 0,
                geometricIncidenceRadians: Double = 0, profileDragCoefficient: Double = 0,
                pitchingMomentCoefficient: Double = 0, effectiveAngleRangeRadians: ClosedRange<Double>,
                reynoldsRange: ClosedRange<Double>) throws {
        guard chordMetres.isFinite, chordMetres > 0, liftSlopePerRadian.isFinite,
              (0.01...20).contains(liftSlopePerRadian),
              [zeroLiftAngleRadians,geometricIncidenceRadians].allSatisfy({ $0.isFinite && abs($0) <= 0.35 }),
              profileDragCoefficient.isFinite, profileDragCoefficient >= 0, pitchingMomentCoefficient.isFinite,
              [effectiveAngleRangeRadians.lowerBound,effectiveAngleRangeRadians.upperBound]
                .allSatisfy({ $0.isFinite && abs($0) <= 0.35 }),
              reynoldsRange.lowerBound.isFinite, reynoldsRange.lowerBound > 0, reynoldsRange.upperBound.isFinite else {
            throw FluidError.invalidParameter
        }
        self.chordMetres = chordMetres; self.liftSlopePerRadian = liftSlopePerRadian; self.zeroLiftAngleRadians = zeroLiftAngleRadians
        self.geometricIncidenceRadians = geometricIncidenceRadians; self.profileDragCoefficient = profileDragCoefficient
        self.pitchingMomentCoefficient = pitchingMomentCoefficient; self.effectiveAngleRangeRadians = effectiveAngleRangeRadians
        self.reynoldsRange = reynoldsRange
    }
}
public struct LiftingLineStationResult: Sendable {
    public let spanPositionMetres: Double
    public let circulationSquareMetresPerSecond: Double
    public let sectionLiftCoefficient: Double
    public let inducedAngleRadians: Double
    public let effectiveAngleRadians: Double
}
public struct LiftingLineResult: Sendable {
    public let stations: [LiftingLineStationResult]
    public let fourierCoefficients: [Double]
    public let liftCoefficient: Double
    public let inducedDragCoefficient: Double
    public let profileDragCoefficient: Double
    /// Undefined at numerically zero total lift, including antisymmetric loading.
    public let spanEfficiency: Double?
    public let maximumEquationResidualRadians: Double
    public let loads: FluidLoads
}

/// Classical linear lifting line for one straight, unswept planar wing along body Y.
/// Uniform incident flow only, zero sideslip, Mach≤0.3, AR≥4, small attached angles.
/// Geometry matrix is factored once (O(N³)); queries cost O(N²). No wake interaction
/// with other surfaces, nonuniform incident speed, body rotation or stall solver.
public struct LiftingLineWing: Sendable {
    public let spanMetres: Double
    public let referenceAreaSquareMetres: Double
    public let momentReferenceLineCentreMetres: SIMD3<Double>
    public let sections: [LiftingLineSection]
    public let source: String
    public let version: String
    public var aspectRatio: Double { (spanMetres/referenceAreaSquareMetres)*spanMetres }
    private let theta: [Double]
    private let widths: [Double]
    private let matrix: [Double]
    private let factor: LiftingLineFactor
    /// sectionAtSpanFraction receives −cos(theta), from left (−1) to right (+1).
    /// Reference area is the dataset's normalization area; it must suit the geometry.
    public init(spanMetres: Double, referenceAreaSquareMetres: Double, stationCount: Int = 8,
                momentReferenceLineCentreMetres: SIMD3<Double> = .zero, source: String, version: String,
                sectionAtSpanFraction: (Double) throws -> LiftingLineSection) throws {
        guard spanMetres.isFinite, spanMetres > 0, referenceAreaSquareMetres.isFinite, referenceAreaSquareMetres > 0,
              (2...64).contains(stationCount), fluidFinite(momentReferenceLineCentreMetres),
              [source,version].allSatisfy({ !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }) else { throw FluidError.invalidGeometry }
        let ar = (spanMetres/referenceAreaSquareMetres)*spanMetres
        guard ar.isFinite, ar >= 4 else { throw FluidError.outsideCoverage }
        let angles = (1...stationCount).map { Double($0)*Double.pi/Double(stationCount+1) }
        let local = try angles.map { try sectionAtSpanFraction(-cos($0)) }
        guard local.allSatisfy({ (1e-6...1).contains($0.chordMetres/spanMetres) }) else { throw FluidError.outsideCoverage }
        var entries = [Double](), weights = [Double]()
        for (j,t) in angles.enumerated() {
            let section = local[j]
            for n in 1...stationCount {
                entries.append(sin(Double(n)*t)*(4*spanMetres/(section.liftSlopePerRadian*section.chordMetres)+Double(n)/sin(t)))
            }
            // Fejer second-rule weights over x=−cos(theta), then dy=span/2 dx.
            var sum = 0.0
            for k in stride(from: 1,through: stationCount,by: 2) { sum += sin(Double(k)*t)/Double(k) }
            weights.append(spanMetres*2*sin(t)*sum/Double(stationCount+1))
        }
        guard entries.allSatisfy(\.isFinite), weights.allSatisfy({ $0.isFinite && $0 > 0 }) else { throw FluidError.nonFiniteResult }
        factor = try LiftingLineFactor(matrix: entries,count: stationCount)
        self.spanMetres = spanMetres; self.referenceAreaSquareMetres = referenceAreaSquareMetres
        self.momentReferenceLineCentreMetres = momentReferenceLineCentreMetres; sections = local; theta = angles
        matrix = entries; widths = weights; self.source = source; self.version = version
    }
    public func evaluate(flow: FluidFlow, momentReferencePointMetres: SIMD3<Double> = .zero) throws -> LiftingLineResult {
        guard let mach = flow.machNumber, let alpha = flow.angleOfAttackRadians,
              let beta = flow.sideslipRadians else { throw FluidError.missingFlowProperty }
        guard flow.speedMetresPerSecond >= 1, mach <= 0.3, abs(alpha) <= 0.35, abs(beta) <= 1e-12 else { throw FluidError.outsideCoverage }
        for section in sections {
            guard let re = try flow.reynoldsNumber(referenceLengthMetres: section.chordMetres) else { throw FluidError.missingFlowProperty }
            guard section.reynoldsRange.contains(re), abs(alpha+section.geometricIncidenceRadians) <= 0.35 else { throw FluidError.outsideCoverage }
        }
        let rhs = sections.map { alpha+$0.geometricIncidenceRadians-$0.zeroLiftAngleRadians }
        let a = try factor.solve(rhs)
        let n = sections.count, speed = flow.speedMetresPerSecond, q = flow.dynamicPressurePascals
        var stations = [LiftingLineStationResult](), profileDrag = 0.0, dragSpanMoment = 0.0, pitch = 0.0, residual = 0.0
        for j in 0..<n {
            let t = theta[j], section = sections[j], y = -spanMetres*cos(t)/2
            var circulationTerm = 0.0, induced = 0.0, equation = 0.0
            for k in 0..<n {
                let harmonic = sin(Double(k+1)*t)*a[k]
                circulationTerm += harmonic; induced += Double(k+1)*harmonic/sin(t)
                equation += matrix[j*n+k]*a[k]
            }
            let effective = alpha+section.geometricIncidenceRadians-induced
            guard section.effectiveAngleRangeRadians.contains(effective) else { throw FluidError.outsideCoverage }
            let cl = 4*spanMetres*circulationTerm/section.chordMetres
            let liftPerSpan = 4*q*spanMetres*circulationTerm
            let profilePerSpan = q*section.chordMetres*section.profileDragCoefficient
            profileDrag += profilePerSpan*widths[j]
            dragSpanMoment += y*(profilePerSpan+liftPerSpan*induced)*widths[j]
            pitch += q*section.chordMetres*section.chordMetres*section.pitchingMomentCoefficient*widths[j]
            residual = max(residual,abs(equation-rhs[j]))
            stations.append(LiftingLineStationResult(spanPositionMetres: y,
                circulationSquareMetresPerSecond: 2*spanMetres*speed*circulationTerm,
                sectionLiftCoefficient: cl,inducedAngleRadians: induced,effectiveAngleRadians: effective))
        }
        guard residual <= 1e-10 else { throw FluidError.nonFiniteResult }
        let sum = a.enumerated().reduce(0.0) { $0+Double($1.offset+1)*$1.element*$1.element }
        let cl = Double.pi*aspectRatio*a[0], cdi = Double.pi*aspectRatio*sum
        let lift = q*referenceAreaSquareMetres*cl, inducedDrag = q*referenceAreaSquareMetres*cdi
        let direction = flow.relativeVelocityMetresPerSecond/speed, liftDirection = SIMD3(direction.z,0,-direction.x)
        let force = liftDirection*lift-direction*(profileDrag+inducedDrag)
        // Exact Fourier integral ∫y L'(y)dy, with y=−span/2 cos(theta).
        let liftSpanMoment = -Double.pi*q*spanMetres*spanMetres*spanMetres*a[1]/4
        let moment = SIMD3(-direction.x,0,-direction.z)*liftSpanMoment
            + SIMD3(-direction.z,0,direction.x)*dragSpanMoment+SIMD3(0,pitch,0)
        let loads = try FluidLoads(forceNewtons: force,momentNewtonMetres: moment,
            referencePointMetres: momentReferenceLineCentreMetres).referred(to: momentReferencePointMetres)
        let cdProfile = sections.enumerated().reduce(0.0) { $0+$1.element.chordMetres*$1.element.profileDragCoefficient*widths[$1.offset] }/referenceAreaSquareMetres
        return LiftingLineResult(stations: stations,fourierCoefficients: a,liftCoefficient: cl,inducedDragCoefficient: cdi,
            profileDragCoefficient: cdProfile,spanEfficiency: sum > 0 && abs(a[0]) > 1e-14 ? a[0]*a[0]/sum : nil,
            maximumEquationResidualRadians: residual,loads: loads)
    }
}

private struct LiftingLineFactor: Sendable {
    let count: Int
    let lu: [Double]
    let swaps: [Int]
    init(matrix: [Double], count: Int) throws {
        var m = matrix, pivots = [Int]()
        let scale = matrix.map { abs($0) }.max()!
        for k in 0..<count {
            let row = (k..<count).max { abs(m[$0*count+k]) < abs(m[$1*count+k]) }!
            guard abs(m[row*count+k]) > scale*1e-12 else { throw FluidError.invalidGeometry }
            pivots.append(row)
            if row != k { for j in 0..<count { m.swapAt(k*count+j,row*count+j) } }
            if k+1 < count {
                for i in k+1..<count {
                    m[i*count+k] /= m[k*count+k]
                    for j in k+1..<count { m[i*count+j] -= m[i*count+k]*m[k*count+j] }
                }
            }
        }
        self.count = count; lu = m; swaps = pivots
    }
    func solve(_ rhs: [Double]) throws -> [Double] {
        var result = rhs
        for k in 0..<count { if swaps[k] != k { result.swapAt(k,swaps[k]) } }
        for i in 0..<count { for j in 0..<i { result[i] -= lu[i*count+j]*result[j] } }
        for i in stride(from: count-1,through: 0,by: -1) {
            if i+1 < count { for j in i+1..<count { result[i] -= lu[i*count+j]*result[j] } }
            result[i] /= lu[i*count+i]
        }
        guard result.allSatisfy(\.isFinite) else { throw FluidError.nonFiniteResult }
        return result
    }
}

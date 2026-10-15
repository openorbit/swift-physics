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

/// Converged fixed-wake solution. Circulation is m²/s and induced air velocity is
/// the medium's physical velocity, so subtract it from body-relative incident flow.
public struct NonlinearLiftingLineResult: Sendable {
    public let circulationSquareMetresPerSecond: [Double]
    public let inducedAirVelocityMetresPerSecond: [SIMD3<Double>]
    public let segments: [WingSegmentResult]
    public let loads: FluidLoads
    public let iterations: Int
    public let maximumNormalizedCirculationResidual: Double
}

/// Nonlinear discrete horseshoe lifting line for multiple unswept lifting surfaces.
/// Uses supplied section polars and nonuniform incident flow, with fixed straight
/// wakes and a finite vortex core. No free-wake rollup, compressibility or separated
/// 3-D flow prediction. Nonlinear polar input does not establish post-stall accuracy.
/// Setup O(N²); solve O(iterations*N²), N≤256. Geometry is aerodynamic, not assembly state.
public struct NonlinearLiftingLine: Sendable {
    public let segments: [WingSegment]
    public let wakeDirection: SIMD3<Double>
    public let wakeLengthMetres: Double
    public let vortexCoreRadiusMetres: Double
    public let relaxation: Double
    public let tolerance: Double
    public let maximumIterations: Int
    private let influence: [[SIMD3<Double>]]
    public init(segments: [WingSegment], wakeDirection: SIMD3<Double> = SIMD3(-1,0,0),
                wakeLengthMetres: Double, vortexCoreRadiusMetres: Double,
                relaxation: Double = 0.1, tolerance: Double = 1e-8, maximumIterations: Int = 1000) throws {
        let directionLength = fluidNorm(wakeDirection)
        guard !segments.isEmpty, segments.count <= 256, fluidFinite(wakeDirection), directionLength.isFinite, directionLength > 0,
              wakeLengthMetres.isFinite, wakeLengthMetres > 0, vortexCoreRadiusMetres.isFinite, vortexCoreRadiusMetres > 0,
              relaxation.isFinite, relaxation > 0, relaxation <= 0.5,
              tolerance.isFinite, tolerance > 0, tolerance <= 1e-3,
              (1...10000).contains(maximumIterations) else { throw FluidError.invalidParameter }
        let direction = wakeDirection/directionLength
        let minimumWidth = segments.map { $0.areaSquareMetres/$0.chordMetres }.min()!
        guard minimumWidth.isFinite, minimumWidth > 0, vortexCoreRadiusMetres <= minimumWidth*0.02,
              segments.allSatisfy({ abs(($0.axes.span*direction).sum()) < 1e-10
                  && -($0.axes.forward*direction).sum() >= cos(0.35) }) else { throw FluidError.invalidGeometry }
        // Duplicate control points are ambiguous even with core regularization.
        for i in segments.indices {
            for j in 0..<i {
                guard fluidNorm(segments[i].positionMetres-segments[j].positionMetres) > 2*vortexCoreRadiusMetres else {
                    throw FluidError.invalidGeometry
                }
            }
        }
        self.segments = segments; self.wakeDirection = direction; self.wakeLengthMetres = wakeLengthMetres
        self.vortexCoreRadiusMetres = vortexCoreRadiusMetres; self.relaxation = relaxation
        self.tolerance = tolerance; self.maximumIterations = maximumIterations
        self.influence = try segments.map { target in
            try segments.map { source in
                try Self.horseshoe(at: target.positionMetres,segment: source,direction: direction,
                    length: wakeLengthMetres,core: vortexCoreRadiusMetres)
            }
        }
    }
    /// Unit-circulation finite straight segment, regularized Biot–Savart kernel.
    private static func vortex(at p: SIMD3<Double>, from a: SIMD3<Double>, to b: SIMD3<Double>, core: Double) throws -> SIMD3<Double> {
        let r1 = p-a, r2 = p-b, line = b-a
        let l1 = fluidNorm(r1), l2 = fluidNorm(r2), lineLength = fluidNorm(line)
        guard l1.isFinite, l2.isFinite, lineLength.isFinite, lineLength > 0,
              l1 > core, l2 > core else { throw FluidError.invalidGeometry }
        let cross = fluidCross(r1,r2)
        let denominator = (cross*cross).sum()+core*core*lineLength*lineLength
        let velocity = cross*((line*(r1/l1-r2/l2)).sum()/(4*Double.pi*denominator))
        guard fluidFinite(velocity) else { throw FluidError.nonFiniteResult }
        return velocity
    }
    private static func horseshoe(at p: SIMD3<Double>, segment: WingSegment, direction: SIMD3<Double>,
                                  length: Double, core: Double) throws -> SIMD3<Double> {
        let halfSpan = segment.axes.span*(0.5*segment.areaSquareMetres/segment.chordMetres)
        let a = segment.positionMetres-halfSpan, b = segment.positionMetres+halfSpan
        return try vortex(at: p,from: a+direction*length,to: a,core: core)
            + vortex(at: p,from: a,to: b,core: core)
            + vortex(at: p,from: b,to: b+direction*length,core: core)
    }
    /// Physical induced medium velocity for external bodies/tails. This includes
    /// bound and trailing vortices; callers must avoid the finite-core/end-point region.
    public func inducedAirVelocity(at positionMetres: SIMD3<Double>, circulation: [Double]) throws -> SIMD3<Double> {
        guard fluidFinite(positionMetres), circulation.count == segments.count,
              circulation.allSatisfy(\.isFinite) else { throw FluidError.invalidParameter }
        var velocity = SIMD3<Double>.zero
        for i in segments.indices {
            velocity += try Self.horseshoe(at: positionMetres,segment: segments[i],direction: wakeDirection,
                length: wakeLengthMetres,core: vortexCoreRadiusMetres)*circulation[i]
        }
        guard fluidFinite(velocity) else { throw FluidError.nonFiniteResult }
        return velocity
    }
    public func evaluate(momentReferencePointMetres: SIMD3<Double> = .zero,
                         flowAtSegment: (WingSegment) throws -> FluidFlow) throws -> NonlinearLiftingLineResult {
        let incident = try segments.map(flowAtSegment)
        for (segment,flow) in zip(segments,incident) {
            guard let mach = flow.machNumber else { throw FluidError.missingFlowProperty }
            let speed = flow.speedMetresPerSecond
            guard speed >= 1, mach <= 0.3, flow.densityKilogramsPerCubicMetre > 0,
                  abs((flow.relativeVelocityMetresPerSecond*segment.axes.span).sum())/speed <= 0.1,
                  -(flow.relativeVelocityMetresPerSecond*wakeDirection).sum()/speed >= cos(0.35) else {
                throw FluidError.outsideCoverage
            }
        }
        var gamma = [Double](repeating: 0,count: segments.count)
        for iteration in 1...maximumIterations {
            var results: [WingSegmentResult] = [], velocities: [SIMD3<Double>] = [], targets: [Double] = []
            var residual = 0.0
            for i in segments.indices {
                let induced = zip(influence[i],gamma).reduce(SIMD3<Double>.zero) { $0+$1.0*$1.1 }
                let local = try incident[i].replacingRelativeVelocity(incident[i].relativeVelocityMetresPerSecond-induced)
                guard let mach = local.machNumber, mach <= 0.3 else { throw FluidError.outsideCoverage }
                let result = try segments[i].evaluate(bodyFlow: local)
                guard let alpha = result.sectionFlow.angleOfAttackRadians, abs(alpha) <= 0.35,
                      let coefficients = result.coefficients, coefficients.drag >= 0,
                      abs(coefficients.lift) <= 2 else { throw FluidError.outsideCoverage }
                let target = 0.5*result.sectionFlow.speedMetresPerSecond*segments[i].chordMetres*coefficients.lift
                let scale = 0.5*incident[i].speedMetresPerSecond*segments[i].chordMetres
                let error = abs(target-gamma[i])/scale
                guard target.isFinite, error.isFinite else { throw FluidError.nonFiniteResult }
                residual = max(residual,error); targets.append(target); results.append(result); velocities.append(induced)
            }
            if residual <= tolerance {
                var loads = try FluidLoads(forceNewtons: .zero,momentNewtonMetres: .zero,referencePointMetres: momentReferencePointMetres)
                for result in results { loads = try loads.adding(result.loads) }
                return NonlinearLiftingLineResult(circulationSquareMetresPerSecond: gamma,
                    inducedAirVelocityMetresPerSecond: velocities,segments: results,loads: loads,
                    iterations: iteration,maximumNormalizedCirculationResidual: residual)
            }
            for i in gamma.indices { gamma[i] += relaxation*(targets[i]-gamma[i]) }
        }
        throw FluidError.nonConvergence
    }
}

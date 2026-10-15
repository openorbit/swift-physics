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

/// Orthonormal section axes expressed in body axes. Forward points toward the
/// leading edge, span points section-right, down = forward × span.
public struct WingSectionAxes: Sendable {
    public let forward: SIMD3<Double>
    public let span: SIMD3<Double>
    public var down: SIMD3<Double> { fluidCross(forward,span) }
    public init(forward: SIMD3<Double> = SIMD3(1,0,0), span: SIMD3<Double> = SIMD3(0,1,0)) throws {
        guard fluidFinite(forward), fluidFinite(span), fluidNorm(forward) > 0, fluidNorm(span) > 0,
              fluidNorm(forward).isFinite, fluidNorm(span).isFinite else { throw FluidError.invalidGeometry }
        let x = forward/fluidNorm(forward), y = span/fluidNorm(span)
        guard abs((x*y).sum()) < 1e-12 else { throw FluidError.invalidGeometry }
        self.forward = x
        let orthogonal = y-x*(x*y).sum()
        self.span = orthogonal/fluidNorm(orthogonal)
    }
    /// Positive geometric incidence raises the leading edge in X-forward/Z-down axes.
    public init(incidenceRadians: Double) throws {
        guard incidenceRadians.isFinite else { throw FluidError.invalidGeometry }
        try self.init(forward: SIMD3(cos(incidenceRadians),0,-sin(incidenceRadians)))
    }
    public func sectionVector(fromBody vector: SIMD3<Double>) -> SIMD3<Double> {
        SIMD3((vector*forward).sum(),(vector*span).sum(),(vector*down).sum())
    }
    public func bodyVector(fromSection vector: SIMD3<Double>) -> SIMD3<Double> {
        forward*vector.x+span*vector.y+down*vector.z
    }
}

/// Inputs to a section coefficient model; section flow excludes spanwise velocity.
/// Deflection is caller-supplied geometry/state, with sign defined by the provider.
public struct WingSectionQuery: Sendable {
    public let flow: FluidFlow
    public let chordMetres: Double
    public let controlDeflectionRadians: Double
    public var angleOfAttackRadians: Double? { flow.angleOfAttackRadians }
    public var machNumber: Double? { flow.machNumber }
    public func reynoldsNumber() throws -> Double? { try flow.reynoldsNumber(referenceLengthMetres: chordMetres) }
    public init(flow: FluidFlow, chordMetres: Double, controlDeflectionRadians: Double = 0) throws {
        guard chordMetres.isFinite, chordMetres > 0, controlDeflectionRadians.isFinite else { throw FluidError.invalidParameter }
        self.flow = flow; self.chordMetres = chordMetres; self.controlDeflectionRadians = controlDeflectionRadians
    }
}
public protocol WingSectionCoefficientModel: Sendable {
    /// CL/CD and CM about the segment's supplied moment-reference position.
    /// Positive lift is section-up; positive moment rotates about section +span.
    func coefficients(for query: WingSectionQuery) throws -> AirfoilCoefficients
}

/// Uses a checked polar. Control deflection requires a separate calibrated provider;
/// it is never silently converted into an incidence change or ignored.
public struct AirfoilWingSectionModel: WingSectionCoefficientModel {
    public let airfoil: any AirfoilCoefficientModel
    public init(airfoil: any AirfoilCoefficientModel) { self.airfoil = airfoil }
    public func coefficients(for query: WingSectionQuery) throws -> AirfoilCoefficients {
        guard query.controlDeflectionRadians == 0 else { throw FluidError.unsupportedQuery }
        guard let angle = query.angleOfAttackRadians, let mach = query.machNumber,
              let re = try query.reynoldsNumber() else { throw FluidError.missingFlowProperty }
        return try airfoil.coefficients(angleOfAttackRadians: angle,reynoldsNumber: re,machNumber: mach)
    }
}

/// One aerodynamic strip, not a vehicle-system component. Coefficient CM and
/// physical strip position must refer to the same moment reference (e.g. quarter chord).
public struct WingSegment: Sendable {
    public let positionMetres: SIMD3<Double>
    public let axes: WingSectionAxes
    public let chordMetres: Double
    public let areaSquareMetres: Double
    public let controlDeflectionRadians: Double
    public let model: any WingSectionCoefficientModel
    public init(positionMetres: SIMD3<Double>, axes: WingSectionAxes, chordMetres: Double,
                areaSquareMetres: Double, controlDeflectionRadians: Double = 0,
                model: any WingSectionCoefficientModel) throws {
        guard fluidFinite(positionMetres), chordMetres.isFinite, chordMetres > 0,
              areaSquareMetres.isFinite, areaSquareMetres > 0, controlDeflectionRadians.isFinite else {
            throw FluidError.invalidGeometry
        }
        self.positionMetres = positionMetres; self.axes = axes; self.chordMetres = chordMetres
        self.areaSquareMetres = areaSquareMetres; self.controlDeflectionRadians = controlDeflectionRadians; self.model = model
    }
    public func evaluate(bodyFlow: FluidFlow) throws -> WingSegmentResult {
        let v = axes.sectionVector(fromBody: bodyFlow.relativeVelocityMetresPerSecond)
        let flow = try FluidFlow(relativeVelocityMetresPerSecond: SIMD3(v.x,0,v.z),
            densityKilogramsPerCubicMetre: bodyFlow.densityKilogramsPerCubicMetre,
            pressurePascals: bodyFlow.pressurePascals,temperatureKelvin: bodyFlow.temperatureKelvin,
            speedOfSoundMetresPerSecond: bodyFlow.speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: bodyFlow.dynamicViscosityPascalSeconds,meanFreePathMetres: bodyFlow.meanFreePathMetres)
        // Zero chord-plane flow or zero density contributes no aerodynamic load;
        // do not query a coefficient model at undefined incidence or Reynolds zero.
        guard flow.dynamicPressurePascals > 0 else {
            return WingSegmentResult(sectionFlow: flow,coefficients: nil,
                loads: try FluidLoads(forceNewtons: .zero,momentNewtonMetres: .zero,referencePointMetres: positionMetres))
        }
        let query = try WingSectionQuery(flow: flow,chordMetres: chordMetres,controlDeflectionRadians: controlDeflectionRadians)
        let c = try model.coefficients(for: query)
        guard c.isFinite else { throw FluidError.invalidCoefficients }
        let direction = flow.relativeVelocityMetresPerSecond/flow.speedMetresPerSecond
        let liftDirection = fluidCross(SIMD3(0,1,0),direction)
        let scale = flow.dynamicPressurePascals*areaSquareMetres
        guard scale.isFinite else { throw FluidError.nonFiniteResult }
        let force = axes.bodyVector(fromSection: (liftDirection*c.lift-direction*c.drag)*scale)
        let moment = axes.span*(scale*c.moment*chordMetres)
        return WingSegmentResult(sectionFlow: flow,coefficients: c,
            loads: try FluidLoads(forceNewtons: force,momentNewtonMetres: moment,referencePointMetres: positionMetres))
    }
}
public struct WingSegmentResult: Sendable {
    public let sectionFlow: FluidFlow
    /// nil for inactive strips (zero chord-plane dynamic pressure).
    public let coefficients: AirfoilCoefficients?
    public let loads: FluidLoads
}
public struct SegmentedWingResult: Sendable {
    public let segments: [WingSegmentResult]
    public let loads: FluidLoads
}

/// Independent-strip evaluator. No induced drag, downwash, wake or finite-wing
/// correction is inferred. Supply corrected local flow or a suitable section model.
public struct SegmentedWing: Sendable {
    public let segments: [WingSegment]
    public init(segments: [WingSegment]) throws {
        guard !segments.isEmpty else { throw FluidError.invalidGeometry }
        self.segments = segments
    }
    /// Uniform rectangular geometry centred on its supplied section-moment-reference line.
    /// Positions are body axes; incidence rotates each section, not the spanwise positions.
    public static func rectangular(spanMetres: Double, chordMetres: Double, segmentCount: Int = 8,
                centreMetres: SIMD3<Double> = .zero, incidenceRadians: Double = 0,
                model: any WingSectionCoefficientModel) throws -> Self {
        guard spanMetres.isFinite, spanMetres > 0, chordMetres.isFinite, chordMetres > 0,
              segmentCount > 0, segmentCount <= 65536, fluidFinite(centreMetres) else { throw FluidError.invalidGeometry }
        let width = spanMetres/Double(segmentCount), axes = try WingSectionAxes(incidenceRadians: incidenceRadians)
        return try Self(segments: (0..<segmentCount).map {
            try WingSegment(positionMetres: centreMetres+SIMD3(0,-spanMetres/2+(Double($0)+0.5)*width,0),
                axes: axes,chordMetres: chordMetres,areaSquareMetres: chordMetres*width,model: model)
        })
    }
    /// Each flow is sampled at its strip position in body axes, enabling external
    /// atmosphere gradients and prescribed induced-flow corrections without shared state.
    public func evaluate(momentReferencePointMetres: SIMD3<Double> = .zero,
                         flowAtSegment: (WingSegment) throws -> FluidFlow) throws -> SegmentedWingResult {
        var total = try FluidLoads(forceNewtons: .zero,momentNewtonMetres: .zero,referencePointMetres: momentReferencePointMetres)
        var results = [WingSegmentResult]()
        for segment in segments {
            let result = try segment.evaluate(bodyFlow: flowAtSegment(segment))
            total = try total.adding(result.loads); results.append(result)
        }
        return SegmentedWingResult(segments: results,loads: total)
    }
    /// Uniform or rigidly rotating medium convenience. Flow's linear velocity is measured at the stated
    /// velocity reference point; both rates use the same reference frame and body axes.
    public func evaluate(flow: FluidFlow, angularVelocityRadiansPerSecond: SIMD3<Double> = .zero,
                mediumAngularVelocityRadiansPerSecond: SIMD3<Double> = .zero,
                velocityReferencePointMetres: SIMD3<Double> = .zero,
                momentReferencePointMetres: SIMD3<Double> = .zero) throws -> SegmentedWingResult {
        guard fluidFinite(angularVelocityRadiansPerSecond), fluidFinite(velocityReferencePointMetres) else { throw FluidError.invalidFlow }
        return try evaluate(momentReferencePointMetres: momentReferencePointMetres) {
            try flow.atLocalPoint(offsetMetres: $0.positionMetres-velocityReferencePointMetres,
                                 angularVelocityRadiansPerSecond: angularVelocityRadiansPerSecond,
                                 mediumAngularVelocityRadiansPerSecond: mediumAngularVelocityRadiansPerSecond)
        }
    }
}

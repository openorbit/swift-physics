// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Orthonormal section axes expressed in body axes. Forward points toward the
/// leading edge, span points section-right, down = forward × span.
public struct WingSectionAxes: Sendable {
    public let forward: SIMD3<Double>
    public let span: SIMD3<Double>
    public var down: SIMD3<Double> { aeroCross(forward,span) }
    public init(forward: SIMD3<Double> = SIMD3(1,0,0), span: SIMD3<Double> = SIMD3(0,1,0)) throws {
        guard aeroFinite(forward), aeroFinite(span), aeroNorm(forward) > 0, aeroNorm(span) > 0,
              aeroNorm(forward).isFinite, aeroNorm(span).isFinite else { throw AerodynamicError.invalidGeometry }
        let x = forward/aeroNorm(forward), y = span/aeroNorm(span)
        guard abs((x*y).sum()) < 1e-12 else { throw AerodynamicError.invalidGeometry }
        self.forward = x
        let orthogonal = y-x*(x*y).sum()
        self.span = orthogonal/aeroNorm(orthogonal)
    }
    /// Positive geometric incidence raises the leading edge in X-forward/Z-down axes.
    public init(incidenceRadians: Double) throws {
        guard incidenceRadians.isFinite else { throw AerodynamicError.invalidGeometry }
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
    public let flow: AerodynamicFlow
    public let chordMetres: Double
    public let controlDeflectionRadians: Double
    public var angleOfAttackRadians: Double? { flow.angleOfAttackRadians }
    public var machNumber: Double? { flow.machNumber }
    public func reynoldsNumber() throws -> Double? { try flow.reynoldsNumber(referenceLengthMetres: chordMetres) }
    public init(flow: AerodynamicFlow, chordMetres: Double, controlDeflectionRadians: Double = 0) throws {
        guard chordMetres.isFinite, chordMetres > 0, controlDeflectionRadians.isFinite else { throw AerodynamicError.invalidParameter }
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
        guard query.controlDeflectionRadians == 0 else { throw AerodynamicError.unsupportedQuery }
        guard let angle = query.angleOfAttackRadians, let mach = query.machNumber,
              let re = try query.reynoldsNumber() else { throw AerodynamicError.missingFlowProperty }
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
        guard aeroFinite(positionMetres), chordMetres.isFinite, chordMetres > 0,
              areaSquareMetres.isFinite, areaSquareMetres > 0, controlDeflectionRadians.isFinite else {
            throw AerodynamicError.invalidGeometry
        }
        self.positionMetres = positionMetres; self.axes = axes; self.chordMetres = chordMetres
        self.areaSquareMetres = areaSquareMetres; self.controlDeflectionRadians = controlDeflectionRadians; self.model = model
    }
    public func evaluate(bodyFlow: AerodynamicFlow) throws -> WingSegmentResult {
        let v = axes.sectionVector(fromBody: bodyFlow.relativeVelocityMetresPerSecond)
        let flow = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(v.x,0,v.z),
            densityKilogramsPerCubicMetre: bodyFlow.densityKilogramsPerCubicMetre,
            pressurePascals: bodyFlow.pressurePascals,temperatureKelvin: bodyFlow.temperatureKelvin,
            speedOfSoundMetresPerSecond: bodyFlow.speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: bodyFlow.dynamicViscosityPascalSeconds,meanFreePathMetres: bodyFlow.meanFreePathMetres)
        // Zero chord-plane flow or zero density contributes no aerodynamic load;
        // do not query a coefficient model at undefined incidence or Reynolds zero.
        guard flow.dynamicPressurePascals > 0 else {
            return WingSegmentResult(sectionFlow: flow,coefficients: nil,
                loads: try AerodynamicLoads(forceNewtons: .zero,momentNewtonMetres: .zero,referencePointMetres: positionMetres))
        }
        let query = try WingSectionQuery(flow: flow,chordMetres: chordMetres,controlDeflectionRadians: controlDeflectionRadians)
        let c = try model.coefficients(for: query)
        guard c.isFinite else { throw AerodynamicError.invalidCoefficients }
        let direction = flow.relativeVelocityMetresPerSecond/flow.speedMetresPerSecond
        let liftDirection = aeroCross(SIMD3(0,1,0),direction)
        let scale = flow.dynamicPressurePascals*areaSquareMetres
        guard scale.isFinite else { throw AerodynamicError.nonFiniteResult }
        let force = axes.bodyVector(fromSection: (liftDirection*c.lift-direction*c.drag)*scale)
        let moment = axes.span*(scale*c.moment*chordMetres)
        return WingSegmentResult(sectionFlow: flow,coefficients: c,
            loads: try AerodynamicLoads(forceNewtons: force,momentNewtonMetres: moment,referencePointMetres: positionMetres))
    }
}
public struct WingSegmentResult: Sendable {
    public let sectionFlow: AerodynamicFlow
    /// nil for inactive strips (zero chord-plane dynamic pressure).
    public let coefficients: AirfoilCoefficients?
    public let loads: AerodynamicLoads
}
public struct SegmentedWingResult: Sendable {
    public let segments: [WingSegmentResult]
    public let loads: AerodynamicLoads
}

/// Independent-strip evaluator. No induced drag, downwash, wake or finite-wing
/// correction is inferred. Supply corrected local flow or a suitable section model.
public struct SegmentedWing: Sendable {
    public let segments: [WingSegment]
    public init(segments: [WingSegment]) throws {
        guard !segments.isEmpty else { throw AerodynamicError.invalidGeometry }
        self.segments = segments
    }
    /// Uniform rectangular geometry centred on its supplied section-moment-reference line.
    /// Positions are body axes; incidence rotates each section, not the spanwise positions.
    public static func rectangular(spanMetres: Double, chordMetres: Double, segmentCount: Int = 8,
                centreMetres: SIMD3<Double> = .zero, incidenceRadians: Double = 0,
                model: any WingSectionCoefficientModel) throws -> Self {
        guard spanMetres.isFinite, spanMetres > 0, chordMetres.isFinite, chordMetres > 0,
              segmentCount > 0, segmentCount <= 65536, aeroFinite(centreMetres) else { throw AerodynamicError.invalidGeometry }
        let width = spanMetres/Double(segmentCount), axes = try WingSectionAxes(incidenceRadians: incidenceRadians)
        return try Self(segments: (0..<segmentCount).map {
            try WingSegment(positionMetres: centreMetres+SIMD3(0,-spanMetres/2+(Double($0)+0.5)*width,0),
                axes: axes,chordMetres: chordMetres,areaSquareMetres: chordMetres*width,model: model)
        })
    }
    /// Each flow is sampled at its strip position in body axes, enabling external
    /// atmosphere gradients and prescribed induced-flow corrections without shared state.
    public func evaluate(momentReferencePointMetres: SIMD3<Double> = .zero,
                         flowAtSegment: (WingSegment) throws -> AerodynamicFlow) throws -> SegmentedWingResult {
        var total = try AerodynamicLoads(forceNewtons: .zero,momentNewtonMetres: .zero,referencePointMetres: momentReferencePointMetres)
        var results = [WingSegmentResult]()
        for segment in segments {
            let result = try segment.evaluate(bodyFlow: flowAtSegment(segment))
            total = try total.adding(result.loads); results.append(result)
        }
        return SegmentedWingResult(segments: results,loads: total)
    }
    /// Uniform medium convenience. Flow's linear velocity is measured at the stated
    /// velocity reference point; ω and positions use the same body axes and origin.
    public func evaluate(flow: AerodynamicFlow, angularVelocityRadiansPerSecond: SIMD3<Double> = .zero,
                velocityReferencePointMetres: SIMD3<Double> = .zero,
                momentReferencePointMetres: SIMD3<Double> = .zero) throws -> SegmentedWingResult {
        guard aeroFinite(angularVelocityRadiansPerSecond), aeroFinite(velocityReferencePointMetres) else { throw AerodynamicError.invalidFlow }
        return try evaluate(momentReferencePointMetres: momentReferencePointMetres) {
            try flow.atLocalPoint(offsetMetres: $0.positionMetres-velocityReferencePointMetres,
                                 angularVelocityRadiansPerSecond: angularVelocityRadiansPerSecond)
        }
    }
}

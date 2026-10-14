// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Aerodynamics

@Suite struct WingAndBodyTests {
    private struct FixedSection: WingSectionCoefficientModel {
        let value: AirfoilCoefficients
        func coefficients(for query: WingSectionQuery) throws -> AirfoilCoefficients { value }
    }
    private struct ThrowingSection: WingSectionCoefficientModel {
        func coefficients(for query: WingSectionQuery) throws -> AirfoilCoefficients { throw AerodynamicError.outsideCoverage }
    }
    private struct ControlledSection: WingSectionCoefficientModel {
        func coefficients(for query: WingSectionQuery) throws -> AirfoilCoefficients {
            AirfoilCoefficients(lift: query.controlDeflectionRadians,drag:0,moment:0)
        }
    }
    private struct FixedPolar: AirfoilCoefficientModel {
        func coefficients(angleOfAttackRadians: Double,reynoldsNumber: Double,machNumber: Double) throws -> AirfoilCoefficients {
            AirfoilCoefficients(lift:1,drag:0.1,moment:0)
        }
    }
    private func flow(_ v: SIMD3<Double> = SIMD3(10,0,0)) throws -> AerodynamicFlow {
        try AerodynamicFlow(relativeVelocityMetresPerSecond:v,densityKilogramsPerCubicMetre:2,
            speedOfSoundMetresPerSecond:100,dynamicViscosityPascalSeconds:0.02)
    }
    @Test(arguments:[1,8,64]) func uniformWingMatchesAnalyticAreaAndMoment(segmentCount:Int) throws {
        let wing = try SegmentedWing.rectangular(spanMetres:8,chordMetres:2,segmentCount:segmentCount,
            centreMetres:SIMD3(3,0,0),model:FixedSection(value:AirfoilCoefficients(lift:1,drag:0.1,moment:0.2)))
        let result = try wing.evaluate(flow:flow())
        // q=100 Pa, S=16 m²; CM about segment reference gives 640 N m.
        #expect(abs(result.loads.forceNewtons.x + 160) < 1e-10)
        #expect(abs(result.loads.forceNewtons.z + 1600) < 1e-10)
        #expect(abs(result.loads.momentNewtonMetres.y - (640+4800)) < 1e-10)
        #expect(abs(result.loads.momentNewtonMetres.x) < 1e-10)
        #expect(abs(result.loads.momentNewtonMetres.z) < 1e-10)
        #expect(result.segments.count == segmentCount)
    }
    @Test func rotationalLocalFlowAndMidpointConvergence() throws {
        var errors = [Double]()
        for count in [8,16,32] {
            let wing = try SegmentedWing.rectangular(spanMetres:8,chordMetres:2,segmentCount:count,
                model:FixedSection(value:AirfoilCoefficients(lift:0,drag:1,moment:0)))
            let result = try wing.evaluate(flow:flow(),angularVelocityRadiansPerSecond:SIMD3(0,0,1))
            // Independently integrate q(y)=rho/2*(10-y)^2 along an 8 m span.
            let exactDrag = 2*(100*8+512.0/12)
            errors.append(abs(-result.loads.forceNewtons.x-exactDrag))
            #expect(result.segments.first!.sectionFlow.speedMetresPerSecond > result.segments.last!.sectionFlow.speedMetresPerSecond)
            #expect(result.loads.momentNewtonMetres.z < 0)
        }
        #expect(abs(errors[0]/errors[1]-4) < 1e-10)
        #expect(abs(errors[1]/errors[2]-4) < 1e-10)
        let wing = try SegmentedWing.rectangular(spanMetres:8,chordMetres:1,
            model:FixedSection(value:AirfoilCoefficients(lift:1,drag:0,moment:0)))
        let result = try wing.evaluate(flow:flow(),angularVelocityRadiansPerSecond:SIMD3(1,0,0))
        #expect(result.segments.first!.sectionFlow.angleOfAttackRadians! < 0)
        #expect(result.segments.last!.sectionFlow.angleOfAttackRadians! > 0)
    }
    @Test func twistAxesAndLocalSampling() throws {
        let axes = try WingSectionAxes(incidenceRadians:0.1)
        let segment = try WingSegment(positionMetres:SIMD3(0,2,0),axes:axes,chordMetres:1,areaSquareMetres:1,
            model:FixedSection(value:AirfoilCoefficients(lift:1,drag:0,moment:0)))
        let result = try segment.evaluate(bodyFlow:flow())
        #expect(abs(result.sectionFlow.angleOfAttackRadians!-0.1) < 1e-14)
        #expect(abs(result.loads.forceNewtons.x) < 1e-12)
        #expect(abs(result.loads.forceNewtons.z+100) < 1e-12)
        let wing = try SegmentedWing(segments:[segment])
        let sampled = try wing.evaluate { _ in try flow(SIMD3(20,0,0)) }
        #expect(abs(sampled.loads.forceNewtons.z+400) < 1e-12)
        #expect(abs(sampled.loads.momentNewtonMetres.x+800) < 1e-12)
        let vertical = try WingSectionAxes(forward:SIMD3(1,0,0),span:SIMD3(0,0,-1))
        #expect(vertical.bodyVector(fromSection:SIMD3(0,0,-1)) == SIMD3(0,-1,0))
        #expect(throws:AerodynamicError.invalidGeometry) { try WingSectionAxes(forward:SIMD3(1,0,0),span:SIMD3(1,1,0)) }
    }
    @Test func zeroFlowControlsAndCoverage() throws {
        let axes = try WingSectionAxes()
        let inactive = try WingSegment(positionMetres:.zero,axes:axes,chordMetres:1,areaSquareMetres:1,model:ThrowingSection())
        #expect(try inactive.evaluate(bodyFlow:flow(.zero)).coefficients == nil)
        #expect(try inactive.evaluate(bodyFlow:flow(SIMD3(0,10,0))).loads.forceNewtons == .zero)
        #expect(throws:AerodynamicError.outsideCoverage) { try inactive.evaluate(bodyFlow:flow()) }
        let controlled = try WingSegment(positionMetres:.zero,axes:axes,chordMetres:1,areaSquareMetres:1,
            controlDeflectionRadians:0.3,model:ControlledSection())
        #expect(try controlled.evaluate(bodyFlow:flow()).loads.forceNewtons.z == -30)
        let polarModel = AirfoilWingSectionModel(airfoil:FixedPolar())
        let unsupported = try WingSegment(positionMetres:.zero,axes:axes,chordMetres:1,areaSquareMetres:1,
            controlDeflectionRadians:0.3,model:polarModel)
        #expect(throws:AerodynamicError.unsupportedQuery) { try unsupported.evaluate(bodyFlow:flow()) }
        let missing = try AerodynamicFlow(relativeVelocityMetresPerSecond:SIMD3(10,0,0),densityKilogramsPerCubicMetre:2)
        let ordinary = try WingSegment(positionMetres:.zero,axes:axes,chordMetres:1,areaSquareMetres:1,model:polarModel)
        #expect(throws:AerodynamicError.missingFlowProperty) { try ordinary.evaluate(bodyFlow:missing) }
    }
    private func samples() throws -> [BodyCoefficientSample] {
        var records = [BodyCoefficientSample]()
        for alpha in [-0.2,0.2] { for beta in [-0.1,0.1] { for mach in [0.1,0.5] { for re in [100.0,500.0] {
            let f = alpha+2*beta+3*mach+re/100
            records.append(BodyCoefficientSample(angleOfAttackRadians:alpha,sideslipRadians:beta,machNumber:mach,reynoldsNumber:re,
                coefficients:try AerodynamicCoefficients(force:SIMD3(-f,2*f,-3*f),moment:SIMD3(f,2*f,3*f))))
        } } } }
        return records
    }
    @Test func fourDimensionalBodyGridAndValidatedLoading() throws {
        let model = try TabulatedBodyCoefficients(samples:samples(),source:"Synthetic affine grid",version:"1")
        let value = try model.coefficients(angleOfAttackRadians:0.05,sideslipRadians:0.025,machNumber:0.3,reynoldsNumber:300)
        #expect(abs(value.force.x+4) < 1e-12)
        #expect(abs(value.moment.z-12) < 1e-12)
        let geometry = try AerodynamicReferenceGeometry(areaSquareMetres:2,momentLengthsMetres:SIMD3(repeating:1))
        let body = try WholeBodyAerodynamics(geometry:geometry,reynoldsReferenceLengthMetres:1,model:model)
        let localFlow = try AerodynamicFlow(relativeVelocityMetresPerSecond:SIMD3(10,0,0),
            densityKilogramsPerCubicMetre:0.6,speedOfSoundMetresPerSecond:100,dynamicViscosityPascalSeconds:0.02)
        let loads = try body.evaluate(flow:localFlow)
        // alpha=beta=0, Mach=.1, Re=300 gives Cx=-3.3 and q*S=60.
        #expect(abs(loads.forceNewtons.x+198) < 1e-10)
        #expect(try body.evaluate(flow:flow(.zero)).forceNewtons == .zero)
        #expect(throws:AerodynamicError.outsideCoverage) { try body.evaluate(flow:flow()) }
        let clamped = try TabulatedBodyCoefficients(samples:model.samples,boundaryPolicy:.clamp,source:"Test",version:"1")
        let boundary = try clamped.coefficients(angleOfAttackRadians:0.3,sideslipRadians:0.2,machNumber:0.6,reynoldsNumber:600)
        #expect(abs(boundary.force.x+6.9) < 1e-12)

        let decoded = try BodyCoefficientTableLoader.json(JSONEncoder().encode(model))
        #expect(try decoded.coefficients(angleOfAttackRadians:0.05,sideslipRadians:0.025,machNumber:0.3,reynoldsNumber:300).force == value.force)
        #expect(throws:AerodynamicError.outsideCoverage) {
            try model.coefficients(angleOfAttackRadians:0.3,sideslipRadians:0,machNumber:0.3,reynoldsNumber:300)
        }
        #expect(throws:AerodynamicError.incompleteGrid) {
            try TabulatedBodyCoefficients(samples:Array(samples().dropLast()),source:"Test",version:"1")
        }
        let duplicate = try samples()
        #expect(throws:AerodynamicError.duplicateCondition) {
            try TabulatedBodyCoefficients(samples:duplicate+[duplicate[0]],source:"Test",version:"1")
        }
    }
    @Test func windTransformAndCapsuleLoads() throws {
        let wind = try WindForceCoefficients(drag:1,lift:0.5)
        let coefficients = try wind.bodyCoefficients(angleOfAttackRadians:0,sideslipRadians:0,momentBody:SIMD3(0,0.1,0))
        let model = try ConstantBodyCoefficients(value:coefficients,source:"Synthetic lifting capsule",version:"1")
        let geometry = try AerodynamicReferenceGeometry(areaSquareMetres:2,momentLengthsMetres:SIMD3(repeating:3),
            momentReferencePointMetres:SIMD3(2,0,0))
        let capsule = try WholeBodyAerodynamics(geometry:geometry,reynoldsReferenceLengthMetres:3,model:model)
        let loads = try capsule.evaluate(flow:flow(),momentReferencePointMetres:.zero)
        #expect(loads.forceNewtons == SIMD3(-200,0,-100))
        #expect(loads.momentNewtonMetres == SIMD3(0,260,0))
        #expect(try capsule.evaluate(flow:flow(.zero)).forceNewtons == .zero)
        let dragOnly = try WindForceCoefficients(drag:2,lift:0)
        let alpha = 0.3, beta = 0.2
        let rotated = try dragOnly.bodyCoefficients(angleOfAttackRadians:alpha,sideslipRadians:beta)
        let direction = SIMD3(cos(alpha)*cos(beta),sin(beta),sin(alpha)*cos(beta))
        #expect(abs((rotated.force*direction).sum()+2) < 1e-14)
        let query = try BodyAerodynamicQuery(flow:flow(),reynoldsReferenceLengthMetres:2,
            momentLengthsMetres:SIMD3(4,2,4),angularVelocityRadiansPerSecond:SIMD3(1,2,3))
        #expect(query.nondimensionalAngularRates == SIMD3(0.2,0.2,0.6))
    }
}

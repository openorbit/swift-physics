// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Testing
import Aerodynamics

@Suite struct BodyRateTests {
    func flow(speed: Double = 10) throws -> FluidFlow {
        try FluidFlow(relativeVelocityMetresPerSecond: SIMD3(speed,0,0),densityKilogramsPerCubicMetre: 2,
            speedOfSoundMetresPerSecond: 100,dynamicViscosityPascalSeconds: 0.02)
    }
    func model(reference: BodyAerodynamicRateReference = .inertialReference) throws -> LinearBodyRateAerodynamics {
        let zero = try FluidCoefficients(force: .zero,moment: .zero)
        let derivatives = BodyAerodynamicRateDerivatives(
            rollRate: try FluidCoefficients(force: SIMD3(0,1,0),moment: SIMD3(-2,0,0)),
            pitchRate: try FluidCoefficients(force: SIMD3(0,0,1),moment: SIMD3(0,-3,0)),
            yawRate: try FluidCoefficients(force: SIMD3(0,1,0),moment: SIMD3(0,0,-4)))
        return try LinearBodyRateAerodynamics(staticModel: ConstantBodyCoefficients(value: zero,source: "Synthetic static",version: "1"),
            derivativeModel: ConstantBodyRateDerivatives(value: derivatives,source: "Synthetic derivatives",version: "1",
                machRange: 0...1,reynoldsRange: 100...10000,incidenceRangeRadians: -0.2...0.2,sideslipRangeRadians: -0.2...0.2,
                momentLengthsMetres: SIMD3(4,2,4),reynoldsReferenceLengthMetres: 2),rateReference: reference,
            minimumSpeedMetresPerSecond: 1,maximumAbsoluteNondimensionalRates: SIMD3(repeating: 0.2))
    }
    @Test func analyticDampingLoadsOpposeRotation() throws {
        let geometry = try FluidReferenceGeometry(areaSquareMetres: 2,momentLengthsMetres: SIMD3(4,2,4))
        let body = try WholeBodyAerodynamics(geometry: geometry,reynoldsReferenceLengthMetres: 2,model: model())
        let rates = SIMD3<Double>(0.5,1,0.5) // All normalized rates = 0.1.
        let loads = try body.evaluate(flow: flow(),angularVelocityRadiansPerSecond: rates)
        #expect(loads.forceNewtons == SIMD3(0,40,20))
        #expect(abs(loads.momentNewtonMetres.x+160) < 1e-12)
        #expect(abs(loads.momentNewtonMetres.y+120) < 1e-12)
        #expect(abs(loads.momentNewtonMetres.z+320) < 1e-12)
        #expect((loads.momentNewtonMetres*rates).sum() < 0) // Negative rotational power.
        let reversed = try body.evaluate(flow: flow(),angularVelocityRadiansPerSecond: -rates)
        #expect(reversed.forceNewtons == -loads.forceNewtons)
        #expect(reversed.momentNewtonMetres == -loads.momentNewtonMetres)
        let fast = try body.evaluate(flow: flow(speed: 20),angularVelocityRadiansPerSecond: rates)
        // Rate coefficient ~1/V and q~V²: dimensional damping scales linearly with V.
        #expect(fast.momentNewtonMetres == 2*loads.momentNewtonMetres)
    }
    @Test func coverageAndRateReferenceAreExplicit() throws {
        let model = try model(reference: .rigidlyRotatingMedium)
        func query(speed: Double = 10, rates: SIMD3<Double> = .zero,
                   reference: BodyAerodynamicRateReference = .rigidlyRotatingMedium,
                   lengths: SIMD3<Double> = SIMD3(4,2,4)) throws -> BodyAerodynamicQuery {
            try BodyAerodynamicQuery(flow: flow(speed: speed),reynoldsReferenceLengthMetres: 2,
                momentLengthsMetres: lengths,angularVelocityRadiansPerSecond: rates,angularRateReference: reference)
        }
        #expect(throws: FluidError.unsupportedQuery) { try model.coefficients(for: query(reference: .inertialReference)) }
        #expect(throws: FluidError.outsideCoverage) { try model.coefficients(for: query(speed: 0.5)) }
        #expect(throws: FluidError.outsideCoverage) { try model.coefficients(for: query(rates: SIMD3(2,0,0))) }
        #expect(throws: FluidError.invalidGeometry) { try model.coefficients(for: query(lengths: SIMD3(repeating: 1))) }
        #expect(try model.coefficients(for: query()).moment == .zero)
    }
    struct ZeroRateStatic: BodyAerodynamicCoefficientModel {
        let source = "Zero-rate query check", version = "1"
        func coefficients(for query: BodyAerodynamicQuery) throws -> FluidCoefficients {
            #expect(query.angularVelocityRadiansPerSecond == .zero)
            return try FluidCoefficients(force: SIMD3(-1,0,-0.5),moment: SIMD3(0,0.2,0))
        }
    }
    @Test func preservesStaticLiftAndRestoringMomentWithoutDoubleCounting() throws {
        let template = try model()
        let wrapper = try LinearBodyRateAerodynamics(staticModel: ZeroRateStatic(),derivativeModel: template.derivativeModel,
            rateReference: .inertialReference,minimumSpeedMetresPerSecond: 1,
            maximumAbsoluteNondimensionalRates: SIMD3(repeating: 0.2))
        let query = try BodyAerodynamicQuery(flow: flow(),reynoldsReferenceLengthMetres: 2,
            momentLengthsMetres: SIMD3(4,2,4),angularVelocityRadiansPerSecond: SIMD3(0,1,0))
        let result = try wrapper.coefficients(for: query)
        #expect(result.force.x == -1)
        #expect(abs(result.force.z+0.4) < 1e-12)
        #expect(abs(result.moment.y+0.1) < 1e-12)
    }
}

// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Aerodynamics

@Suite struct AerodynamicsFoundationTests {
    private func polar(re: Double, mach: Double) throws -> AirfoilPolar {
        // Synthetic affine field has an exact independent interpolation result in all three coordinates.
        try AirfoilPolar(reynoldsNumber: re,machNumber: mach,samples: [-0.2,0.0,0.3].map {
            AirfoilSample(angleOfAttackRadians: $0,coefficients: AirfoilCoefficients(
                lift: 2*$0+re/100+3*mach,drag: 0.1+$0+re/1000+mach,moment: -$0+mach))
        })
    }
    @Test func irregularAnglesReynoldsMachAndCoverage() throws {
        let model = try TabulatedAirfoil(polars: [polar(re:100,mach:0),polar(re:300,mach:0),
            polar(re:100,mach:1),polar(re:300,mach:1)],source:"Synthetic affine field",version:"1")
        let result = try model.coefficients(angleOfAttackRadians:0.15,reynoldsNumber:200,machNumber:0.4)
        #expect(abs(result.lift-3.5) < 1e-12)
        #expect(abs(result.drag-0.85) < 1e-12)
        #expect(abs(result.moment-0.25) < 1e-12)
        #expect(throws: AerodynamicError.outsideCoverage) {
            try model.coefficients(angleOfAttackRadians:0.31,reynoldsNumber:200,machNumber:0.4)
        }
        #expect(throws: AerodynamicError.outsideCoverage) {
            try model.coefficients(angleOfAttackRadians:0,reynoldsNumber:400,machNumber:0.4)
        }
        #expect(throws: AerodynamicError.outsideCoverage) {
            try model.coefficients(angleOfAttackRadians:0,reynoldsNumber:200,machNumber:1.1)
        }
        let clamp = try TabulatedAirfoil(polars:model.polars,boundaryPolicy:.clamp,source:"Test",version:"1")
        #expect(try clamp.coefficients(angleOfAttackRadians:0.5,reynoldsNumber:400,machNumber:2).lift == 6.6)
        let decoded = try JSONDecoder().decode(TabulatedAirfoil.self,from:JSONEncoder().encode(model))
        #expect(try decoded.coefficients(angleOfAttackRadians:0.15,reynoldsNumber:200,machNumber:0.4).lift == result.lift)
    }
    @Test func validatesGridAndDecodedPolars() throws {
        #expect(throws: AerodynamicError.duplicateCondition) {
            try TabulatedAirfoil(polars:[polar(re:100,mach:0),polar(re:100,mach:0)],source:"Test",version:"1")
        }
        #expect(throws: AerodynamicError.incompleteGrid) {
            try TabulatedAirfoil(polars:[polar(re:100,mach:0),polar(re:300,mach:1)],source:"Test",version:"1")
        }
        let bad = Data("{\"reynoldsNumber\":100,\"machNumber\":0,\"samples\":[]}".utf8)
        #expect(throws: AerodynamicError.invalidCoefficients) { try JSONDecoder().decode(AirfoilPolar.self,from:bad) }
    }
    @Test func yamlPreservesSamplesAndUsesAllReynoldsEntries() throws {
        let yaml = """
        data:
        - reynolds-number: 100
          coefficients:
          - [-10, -1, 0.1, 0, -0.1]
          - [3, 1, 0.2, 0, 0.1]
        - reynolds-number: 300
          coefficients:
          - [-10, 1, 0.3, 0, 0.1]
          - [3, 3, 0.4, 0, 0.3]
        """
        let model = try AirfoilTableLoader.yaml(yaml)
        let legacy = try #require(ParametricAirfoil(withString:yaml))
        let angle = -3.5 * .pi/180
        #expect(abs(try model.coefficients(angleOfAttackRadians:angle,reynoldsNumber:200,machNumber:0).lift-1) < 1e-12)
        #expect(abs(legacy.liftCoefficient(at:angle,for:200)-1) < 1e-12)
        #expect(model.polars[0].samples.count == 2)
        for malformed in ["data: []", "data: [{reynolds-number: 100, coefficients: [[0,1]]}]",
            yaml.replacingOccurrences(of:"[3, 1,",with:"[-10, 1,")] {
            #expect(ParametricAirfoil(withString:malformed) == nil)
        }
        let decoded = try JSONDecoder().decode(ParametricAirfoil.self,from:JSONEncoder().encode(legacy))
        #expect(decoded.liftCoefficient(at:angle,for:200) == legacy.liftCoefficient(at:angle,for:200))
    }
    @Test func legacyCodableValidationAndAtmosphereIntegration() throws {
        let values = (-180...180).map { AirfoilCoefficients(lift:Double($0),drag:0.1,moment:0) }
        struct Legacy: Encodable { let reynoldsNumber: Double; let coefficients: [AirfoilCoefficients] }
        for count in [361,362] {
            let encoded = try JSONEncoder().encode(Legacy(reynoldsNumber:1e6,
                coefficients: count == 361 ? values : [values[0]]+values))
            let decoded = try JSONDecoder().decode(ParametricAirfoil.self,from:encoded)
            #expect(abs(decoded.liftCoefficient(at:0.5 * .pi/180,for:1e6)-0.5) < 1e-12)
        }
        let invalid = try JSONEncoder().encode(Legacy(reynoldsNumber:1e6,coefficients:[]))
        #expect(throws:AerodynamicError.invalidCoefficients) { try JSONDecoder().decode(ParametricAirfoil.self,from:invalid) }
        let properties = try ISAAtmosphere().properties(at:0)
        let flow = try AerodynamicFlow(relativeVelocityMetresPerSecond:SIMD3(100,0,0),properties:properties)
        #expect(abs(flow.dynamicPressurePascals-6125) < 0.01)
        #expect(flow.machNumber != nil)
        #expect(try flow.reynoldsNumber(referenceLengthMetres:1) != nil)
        #expect(try flow.knudsenNumber(referenceLengthMetres:1) == nil)
    }

    @Test func nacaCamberAndSurfaceNormal() throws {
        let symmetric = try NACA4(validatingCamberPercentage:0,maximumCamber:0,maximumThicknessPercentage:12)
        #expect(symmetric.meanCamberLine(0) == 0)
        #expect(symmetric.meanCamberSlope(0) == 0)
        let foil = try NACA4(validatingCamberPercentage:2,maximumCamber:4,maximumThicknessPercentage:12)
        #expect(abs(foil.meanCamberLine(0.4)-0.02) < 1e-15)
        #expect(abs(foil.meanCamberSlope(0)-0.1) < 1e-15)
        #expect(foil.meanCamberSlope(0.4) == 0)
        for x in [0.1,0.4,0.7] {
            let (upper,lower) = foil.surfaceCoord(x:x), delta = upper-lower
            let slope = foil.meanCamberSlope(x)
            // Thickness is normal to camber tangent (1,slope), not a vertical offset.
            #expect(abs(delta.x+slope*delta.y) < 1e-14)
            #expect(abs(hypot(delta.x,delta.y)-2*foil.thickness(x)) < 1e-14)
            #expect(abs((upper.y+lower.y)/2-foil.meanCamberLine(x)) < 1e-14)
            let step = 1e-5
            let numericalSlope = (foil.meanCamberLine(x+step)-foil.meanCamberLine(x-step))/(2*step)
            #expect(abs(numericalSlope-slope) < 1e-6)
        }
        #expect(throws:AerodynamicError.invalidGeometry) {
            try NACA4(validatingCamberPercentage:2,maximumCamber:0,maximumThicknessPercentage:12)
        }
        #expect(throws:AerodynamicError.invalidGeometry) { try foil.surfaceCoordinates(atChordFraction:1.01) }
        let cambered = ThinCamberedAirfoil(cl0:0.4)
        #expect(cambered.coefficients(at:0.1,for:1e6).lift == cambered.liftCoefficient(at:0.1,for:1e6))
    }
    @Test func flowNumbersRotationAndMissingProperties() throws {
        let flow = try AerodynamicFlow(relativeVelocityMetresPerSecond:SIMD3(3,0,4),densityKilogramsPerCubicMetre:2,
            speedOfSoundMetresPerSecond:10,dynamicViscosityPascalSeconds:0.5,meanFreePathMetres:0.01)
        #expect(flow.speedMetresPerSecond == 5)
        #expect(flow.dynamicPressurePascals == 25)
        #expect(flow.machNumber == 0.5)
        #expect(try flow.reynoldsNumber(referenceLengthMetres:2) == 40)
        #expect(try flow.knudsenNumber(referenceLengthMetres:2) == 0.005)
        #expect(flow.angleOfAttackRadians == atan2(4.0,3.0))
        let local = try flow.atLocalPoint(offsetMetres:SIMD3(0,2,0),angularVelocityRadiansPerSecond:SIMD3(0,0,1))
        #expect(local.relativeVelocityMetresPerSecond == SIMD3(1,0,4))
        let still = try AerodynamicFlow(relativeVelocityMetresPerSecond:.zero,densityKilogramsPerCubicMetre:0)
        #expect(still.dynamicPressurePascals == 0)
        #expect(still.angleOfAttackRadians == nil)
        #expect(still.sideslipRadians == nil)
        #expect(still.machNumber == nil)
        #expect(try still.reynoldsNumber(referenceLengthMetres:1) == nil)
        #expect(throws:AerodynamicError.invalidFlow) {
            try AerodynamicFlow(relativeVelocityMetresPerSecond:.zero,densityKilogramsPerCubicMetre:-1)
        }
        #expect(throws:AerodynamicError.invalidGeometry) { try flow.reynoldsNumber(referenceLengthMetres:0) }
    }
    @Test func coefficientScalingMomentTransferAndAddition() throws {
        let geometry = try AerodynamicReferenceGeometry(areaSquareMetres:2,momentLengthsMetres:SIMD3(4,3,4),
            momentReferencePointMetres:SIMD3(2,0,0))
        let coefficients = try AerodynamicCoefficients(force:SIMD3(-1,0,-2),moment:SIMD3(0,0.5,0))
        let loads = try coefficients.loads(dynamicPressurePascals:10,geometry:geometry)
        #expect(loads.forceNewtons == SIMD3(-20,0,-40))
        #expect(loads.momentNewtonMetres == SIMD3(0,30,0))
        let shifted = try loads.referred(to:.zero)
        #expect(shifted.momentNewtonMetres == SIMD3(0,110,0))
        #expect(try shifted.referred(to:geometry.momentReferencePointMetres).momentNewtonMetres == loads.momentNewtonMetres)
        #expect(try shifted.adding(loads).momentNewtonMetres == SIMD3(0,220,0))
        #expect(try coefficients.loads(dynamicPressurePascals:0,geometry:geometry).forceNewtons == .zero)
    }
}

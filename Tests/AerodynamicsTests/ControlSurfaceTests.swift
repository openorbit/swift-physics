// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Aerodynamics

@Suite struct ControlSurfaceTests {
    func setting(_ deflection: Double, angles: [Double] = [-0.2,0.2]) throws -> ControlSurfaceSetting {
        var polars = [AirfoilPolar]()
        for mach in [0.0,1.0] { for re in [1000.0,10000.0] {
            polars.append(try AirfoilPolar(reynoldsNumber: re,machNumber: mach,samples: angles.map {
                AirfoilSample(angleOfAttackRadians: $0,coefficients: AirfoilCoefficients(
                    lift: 2*$0+3*deflection+re/10000+mach,drag: 0.2+0.1*deflection+0.01*$0,
                    moment: -0.5*$0-deflection+0.01*mach))
            }))
        } }
        return try ControlSurfaceSetting(deflectionRadians: deflection,
            airfoil: TabulatedAirfoil(polars: polars,source: "Synthetic affine polar",version: "1"))
    }
    func table(_ settings: [ControlSurfaceSetting], policy: CoefficientBoundaryPolicy = .reject) throws -> TabulatedControlSurface {
        try TabulatedControlSurface(settings: settings,deflectionBoundaryPolicy: policy,
            source: "Synthetic control surface",version: "1",deflectionConvention: "Positive trailing edge down",
            momentReferenceDescription: "Quarter chord")
    }
    func query(alpha: Double = 0.1, deflection: Double = 0.05) throws -> WingSectionQuery {
        try WingSectionQuery(flow: AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(50*cos(alpha),0,50*sin(alpha)),
            densityKilogramsPerCubicMetre: 1,speedOfSoundMetresPerSecond: 100,dynamicViscosityPascalSeconds: 0.01),
            chordMetres: 1,controlDeflectionRadians: deflection)
    }
    @Test func fourCoordinateInterpolationAndValidatedJSON() throws {
        let model = try table([setting(0.2),setting(-0.1),setting(0)])
        let value = try model.coefficients(for: query())
        #expect(abs(value.lift-1.35) < 1e-12)
        #expect(abs(value.drag-0.206) < 1e-12)
        #expect(abs(value.moment+0.095) < 1e-12)
        let decoded = try ControlSurfaceTableLoader.load(data: JSONEncoder().encode(model))
        #expect(try decoded.coefficients(for: query()).lift == value.lift)
        #expect(decoded.settings.map(\.deflectionRadians) == [-0.1,0,0.2])
        var json = try #require(JSONSerialization.jsonObject(with: JSONEncoder().encode(model)) as? [String: Any])
        json["settings"] = []
        #expect(throws: AerodynamicError.invalidParameter) {
            try ControlSurfaceTableLoader.load(data: JSONSerialization.data(withJSONObject: json))
        }
    }
    @Test func exactSettingsAndCoveragePolicies() throws {
        let settings = try [setting(0),setting(0.2,angles: [-0.05,0.05])]
        let model = try table(settings)
        #expect(try model.coefficients(for: query(deflection: 0)).lift > 0)
        #expect(throws: AerodynamicError.outsideCoverage) { try model.coefficients(for: query(deflection: 0.1)) }
        #expect(throws: AerodynamicError.outsideCoverage) { try model.coefficients(for: query(deflection: -0.01)) }
        let clamped = try table(settings,policy: .clamp)
        #expect(try clamped.coefficients(for: query(deflection: -0.01)).lift == model.coefficients(for: query(deflection: 0)).lift)
        #expect(throws: AerodynamicError.duplicateCondition) { try table([setting(0),setting(0)]) }
        #expect(throws: AerodynamicError.invalidParameter) { try ControlSurfaceSetting(deflectionRadians: .nan,airfoil: setting(0).airfoil) }
    }
    @Test func suppliedPostStallSamplesArePreserved() throws {
        let samples = [
            AirfoilSample(angleOfAttackRadians: 0,coefficients: AirfoilCoefficients(lift: 0,drag: 0.02,moment: 0)),
            AirfoilSample(angleOfAttackRadians: 0.2,coefficients: AirfoilCoefficients(lift: 1,drag: 0.1,moment: -0.1)),
            AirfoilSample(angleOfAttackRadians: 0.4,coefficients: AirfoilCoefficients(lift: 0.6,drag: 0.4,moment: -0.2))]
        let polars = try [0.0,1.0].flatMap { mach in
            try [1000.0,10000.0].map { try AirfoilPolar(reynoldsNumber: $0,machNumber: mach,samples: samples) }
        }
        let model = try table([ControlSurfaceSetting(deflectionRadians: 0,
            airfoil: TabulatedAirfoil(polars: polars,source: "Synthetic post-stall shape",version: "1"))])
        let coefficient = try model.coefficients(for: query(alpha: 0.3,deflection: 0))
        #expect(abs(coefficient.lift-0.8) < 1e-12)
        #expect(abs(coefficient.drag-0.25) < 1e-12)
        #expect(throws: AerodynamicError.outsideCoverage) { try model.coefficients(for: query(alpha: 0.5,deflection: 0)) }
    }
    @Test func differentialWingAndVerticalTailLoads() throws {
        let model = try table([setting(-0.1),setting(0.1)])
        let axes = try WingSectionAxes()
        let wing = try SegmentedWing(segments: [-1.0,1.0].map { y in
            try WingSegment(positionMetres: SIMD3(0,y,0),axes: axes,chordMetres: 1,areaSquareMetres: 1,
                controlDeflectionRadians: 0.1*y,model: model)
        })
        let flow = try query(alpha: 0).flow
        let loads = try wing.evaluate(flow: flow).loads
        #expect(abs(loads.momentNewtonMetres.x+750) < 1e-10)
        #expect(abs(loads.forceNewtons.z+2500) < 1e-10)
        let tail = try WingSegment(positionMetres: SIMD3(-3,0,0),
            axes: WingSectionAxes(forward: SIMD3(1,0,0),span: SIMD3(0,0,1)),
            chordMetres: 1,areaSquareMetres: 1,controlDeflectionRadians: 0.1,model: model)
        let tailLoads = try tail.evaluate(bodyFlow: flow).loads.referred(to: .zero)
        #expect(abs(tailLoads.forceNewtons.y-1625) < 1e-10)
        #expect(abs(tailLoads.momentNewtonMetres.z+4993.75) < 1e-10)
    }
}

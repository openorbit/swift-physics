// SPDX-License-Identifier: Apache-2.0
import Testing
import Aerodynamics

@Suite struct ParachuteBalloonTests {
    func coefficients(maximumAngle: Double? = nil) throws -> TabulatedDragCoefficients {
        try TabulatedDragCoefficients(samples: [0.0,1.0].flatMap { mach in
            [100.0,10000.0].map { DragCoefficientSample(machNumber: mach,reynoldsNumber: $0,dragCoefficient: 2) }
        }, source: "Synthetic drag reference",version: "1",maximumAxisAngleRadians: maximumAngle)
    }
    func flow(_ v: SIMD3<Double> = SIMD3(10,0,0), density: Double = 2) throws -> AerodynamicFlow {
        try AerodynamicFlow(relativeVelocityMetresPerSecond: v,densityKilogramsPerCubicMetre: density,
            speedOfSoundMetresPerSecond: 100,dynamicViscosityPascalSeconds: 0.02)
    }
    @Test func steadyParachuteDragDirectionAndMoment() throws {
        let drag = try DragAreaAerodynamics(areaSquareMetres: 3,referenceLengthMetres: 1,
            forceApplicationPointMetres: SIMD3(0,2,0),axis: SIMD3(2,0,0),model: coefficients(maximumAngle: 0.2))
        let chute = ParachuteAerodynamics(drag: drag)
        let loads = try chute.evaluate(flow: flow())
        #expect(loads.forceNewtons == SIMD3(-600,0,0))
        #expect(loads.momentNewtonMetres == SIMD3(0,0,1200))
        let shifted = try chute.evaluate(flow: flow(),momentReferencePointMetres: SIMD3(0,2,0))
        #expect(shifted.momentNewtonMetres == .zero)
        #expect(throws: AerodynamicError.outsideCoverage) { try chute.evaluate(flow: flow(SIMD3(-10,0,0))) }
        // Coefficients are orientation-invariant when no angle coverage is declared.
        let oblique = try DragAreaAerodynamics(areaSquareMetres: 2,referenceLengthMetres: 1,model: coefficients())
        #expect(try oblique.evaluate(flow: flow(SIMD3(3,4,0))).forceNewtons == SIMD3(-60,-80,0))
    }
    @Test func terminalDescentBalancesExternalWeight() throws {
        let chute = ParachuteAerodynamics(drag: try DragAreaAerodynamics(areaSquareMetres: 2,
            referenceLengthMetres: 1,model: coefficients()))
        // m=20 kg, g=10, rho=1, Cd=2, A=2 -> terminal speed 10 m/s.
        let loads = try chute.evaluate(flow: flow(SIMD3(0,0,10),density: 1))
        #expect(loads.forceNewtons+SIMD3(0,0,200) == .zero)
    }
    @Test func buoyancyAlternativeFieldsAndBalloonDrag() throws {
        let buoyancy = try Buoyancy(displacedVolumeCubicMetres: 10,centreOfBuoyancyMetres: SIMD3(1,0,0))
        let gravity = BuoyancyField.hydrostaticGravityMetresPerSecondSquared(SIMD3(0,0,10))
        let pressure = BuoyancyField.pressureGradientPascalsPerMetre(SIMD3(0,0,20))
        let a = try buoyancy.evaluate(densityKilogramsPerCubicMetre: 2,field: gravity)
        let b = try buoyancy.evaluate(densityKilogramsPerCubicMetre: 2,field: pressure)
        #expect(a.forceNewtons == SIMD3(0,0,-200))
        #expect(a.momentNewtonMetres == SIMD3(0,200,0))
        #expect(a.forceNewtons == b.forceNewtons)
        #expect(a.momentNewtonMetres == b.momentNewtonMetres)
        // Total external mass=displaced air mass (20 kg) gives equilibrium at rest.
        #expect(a.forceNewtons+SIMD3(0,0,200) == .zero)
        let balloon = BalloonAerodynamics(buoyancy: buoyancy,drag: try DragAreaAerodynamics(
            areaSquareMetres: 3,referenceLengthMetres: 1,forceApplicationPointMetres: SIMD3(0,2,0),model: coefficients()))
        let total = try balloon.evaluate(flow: flow(),field: gravity)
        #expect(total.forceNewtons == SIMD3(-600,0,-200))
        #expect(total.momentNewtonMetres == SIMD3(0,200,1200))
    }
    @Test func dragInterpolationAndSuppliedStateLimits() throws {
        let table = try TabulatedDragCoefficients(samples: [0.0,1.0].flatMap { mach in
            [1000.0,10000.0].map { DragCoefficientSample(machNumber: mach,reynoldsNumber: $0,
                dragCoefficient: 0.1+2*mach+$0/10000) }
        },source: "Synthetic affine coefficients",version: "1")
        let local = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(50,0,0),
            densityKilogramsPerCubicMetre: 1,speedOfSoundMetresPerSecond: 100,dynamicViscosityPascalSeconds: 0.01)
        #expect(abs((try table.dragCoefficient(for: DragCoefficientQuery(flow: local,referenceLengthMetres: 1)))-1.6) < 1e-12)
        let rest = try AerodynamicFlow(relativeVelocityMetresPerSecond: .zero,densityKilogramsPerCubicMetre: 2)
        let drag = try DragAreaAerodynamics(areaSquareMetres: 1,referenceLengthMetres: 1,model: table)
        #expect(try drag.evaluate(flow: rest).forceNewtons == .zero)
        let collapsed = try DragAreaAerodynamics(areaSquareMetres: 0,referenceLengthMetres: 1,model: table)
        #expect(try collapsed.evaluate(flow: flow()).forceNewtons == .zero)
        #expect(throws: AerodynamicError.missingFlowProperty) {
            try drag.evaluate(flow: AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(10,0,0),densityKilogramsPerCubicMetre: 1))
        }
        #expect(throws: AerodynamicError.invalidGeometry) { try Buoyancy(displacedVolumeCubicMetres: -1) }
        #expect(throws: AerodynamicError.invalidGeometry) {
            try DragAreaAerodynamics(areaSquareMetres: 1,referenceLengthMetres: 1,axis: .zero,model: table)
        }
        #expect(throws: AerodynamicError.incompleteGrid) {
            try TabulatedDragCoefficients(samples: [DragCoefficientSample(machNumber: 0,reynoldsNumber: 1,dragCoefficient: 1),
                DragCoefficientSample(machNumber: 1,reynoldsNumber: 2,dragCoefficient: 1)],source: "Incomplete",version: "1")
        }
    }
}

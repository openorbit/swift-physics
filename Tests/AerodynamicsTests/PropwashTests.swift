// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Aerodynamics

@Suite struct PropwashTests {
    func disk(speed: Double = 0, torque: Double = 10, axis: SIMD3<Double> = SIMD3(1,0,0)) throws -> ActuatorDiskPropwash {
        try ActuatorDiskPropwash(centreMetres: .zero,downstreamAxis: axis,radiusMetres: 1,
            thrustNewtons: 200,airTorqueNewtonMetres: torque,densityKilogramsPerCubicMetre: 2,
            upstreamAxialSpeedMetresPerSecond: speed)
    }
    @Test func hoverAndForwardFlightMomentumReferences() throws {
        for speed in [0.0,10,100] {
            let rotor = try disk(speed: speed)
            let vi = rotor.inducedSpeedAtDiskMetresPerSecond
            #expect(abs(200-2*2*Double.pi*vi*(speed+vi)) < 1e-10)
            #expect(abs(try rotor.velocityMetresPerSecond(at: .zero).x-vi) < 1e-12)
            let far = try rotor.velocityMetresPerSecond(at: SIMD3(1e6,0,0)).x
            #expect(abs(far/(2*vi)-1) < 1e-12)
            let radius = try rotor.wakeRadiusMetres(downstreamDistanceMetres: 1e6)
            #expect(abs(2*Double.pi*radius*radius*(speed+far)-rotor.massFlowKilogramsPerSecond) < 1e-10)
            #expect(abs(rotor.massFlowKilogramsPerSecond*far-200) < 1e-9)
            #expect(abs(rotor.idealAxialPowerWatts-200*(speed+vi)) < 1e-10)
        }
    }
    @Test func swirlCarriesTorqueAndRotatesWithDisk() throws {
        let rotor = try disk(speed: 10)
        let z = 2.0, radius = try rotor.wakeRadiusMetres(downstreamDistanceMetres: z)
        // Radial midpoint quadrature of angular-momentum flux through the wake.
        let count = 512, dr = radius/Double(count), axialSpeed = 10+(try rotor.velocityMetresPerSecond(at: SIMD3(z,0,0))).x
        var torque = 0.0
        for i in 0..<count {
            let r = (Double(i)+0.5)*dr
            let v = try rotor.velocityMetresPerSecond(at: SIMD3(z,r,0))
            torque += 2*axialSpeed*r*v.z*2*Double.pi*r*dr
        }
        #expect(abs(torque-10) < 2e-5)
        let reversed = try disk(speed: 10,torque: -10)
        let a = try rotor.velocityMetresPerSecond(at: SIMD3(2,0.2,0))
        let b = try reversed.velocityMetresPerSecond(at: SIMD3(2,0.2,0))
        #expect(a.x == b.x && a.z == -b.z)
        let rotated = try disk(speed: 10,axis: SIMD3(0,0,1)).velocityMetresPerSecond(at: SIMD3(0.2,0,2))
        #expect(abs(rotated.z-a.x)+abs(rotated.y-a.z) < 1e-12)
        let loads = try rotor.reactionLoads(momentReferencePointMetres: SIMD3(0,1,0))
        #expect(loads.forceNewtons == SIMD3(-200,0,0))
        #expect(loads.momentNewtonMetres == SIMD3(-10,0,-200))
    }
    @Test func localWingResponseAndVelocityComposition() throws {
        struct Section: WingSectionCoefficientModel {
            func coefficients(for query: WingSectionQuery) throws -> AirfoilCoefficients {
                AirfoilCoefficients(lift: 1,drag: 0.1,moment: 0)
            }
        }
        let rotor = try disk(speed: 10,axis: SIMD3(-1,0,0))
        let field = try SuperposedAirVelocityField(backgroundMetresPerSecond: SIMD3(-10,0,0),increments: [rotor])
        let segment = try WingSegment(positionMetres: SIMD3(-2,0,0),axes: WingSectionAxes(),chordMetres: 1,
            areaSquareMetres: 1,model: Section())
        let base = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(10,0,0),densityKilogramsPerCubicMetre: 2)
        let local = try base.replacingRelativeVelocity(-field.velocityMetresPerSecond(at: segment.positionMetres))
        #expect(try segment.evaluate(bodyFlow: local).loads.forceNewtons.z < segment.evaluate(bodyFlow: base).loads.forceNewtons.z)
        #expect(try rotor.velocityMetresPerSecond(at: SIMD3(1,0,0)) == .zero)
        #expect(try rotor.velocityMetresPerSecond(at: SIMD3(-1,2,0)) == .zero)
    }
    @Test func volumeSourceForceAndTorqueIntegrals() throws {
        let rotor = try disk(), source = try ActuatorDiskAccelerationField(disk: rotor,thicknessMetres: 0.5)
        // Continuous cylinder quadrature, independent radial integration of torque.
        var force = 0.0, torque = 0.0
        let count = 1024, dr = 1/Double(count)
        for i in 0..<count {
            let r = (Double(i)+0.5)*dr
            let acceleration = try source.accelerationMetresPerSecondSquared(at: SIMD3(0,r,0))
            let mass = 2*2*Double.pi*r*dr*0.5
            force += mass*acceleration.x; torque += mass*r*acceleration.z
        }
        #expect(abs(force-200) < 1e-10)
        #expect(abs(torque-10) < 1e-5)
        #expect(try source.accelerationMetresPerSecondSquared(at: SIMD3(2,0,0)) == .zero)
    }
    @Test func zeroThrustAndInvalidInputs() throws {
        let off = try ActuatorDiskPropwash(centreMetres: .zero,downstreamAxis: SIMD3(1,0,0),radiusMetres: 1,
            thrustNewtons: 0,densityKilogramsPerCubicMetre: 1,upstreamAxialSpeedMetresPerSecond: 0)
        #expect(try off.velocityMetresPerSecond(at: .zero) == .zero)
        #expect(throws: AerodynamicError.invalidParameter) { try disk(axis: .zero) }
        #expect(throws: AerodynamicError.outsideCoverage) { try disk().wakeRadiusMetres(downstreamDistanceMetres: -1) }
        #expect(throws: AerodynamicError.invalidParameter) { try disk().velocityMetresPerSecond(at: SIMD3(.nan,0,0)) }
    }
}

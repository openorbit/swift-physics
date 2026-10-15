// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct ParkingOrbitTransferTests {
    func gm(_ value: Double = 1) throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: value, convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Explicit analytic reference GM", version: "1"))
    }
    func transfer(_ speed: Double = 1) throws -> ParkingOrbitTransfer {
        try ParkingOrbitTransfer(parkingRadiusMetres: 1, bodyRadiusMetres: 0.5,
            minimumPeriapsisAltitudeMetres: 0.5, hyperbolicExcessSpeedMetresPerSecond: speed,
            gravitationalParameter: gm())
    }
    func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: t)) }
    @Test func independentHighPrecisionScalarReferences() throws {
        // Values calculated independently with 55-digit Decimal arithmetic.
        let rows: [(Double,Double,Double,Double,Double,Double,Double)] = [
            (398600441800000,6671000,3000,7729.8918473558424098,11335.892375267005691,3606.0005279111632816,1.1506245194533048307),
            (42828375214000,3590000,2500,3453.9701113880700447,5487.2414800812469787,2033.2713686931769340,1.5238933274467412767)]
        for (mu,r,v,vc,vp,dv,e) in rows {
            let t = try ParkingOrbitTransfer(parkingRadiusMetres: r, bodyRadiusMetres: r-100000,
                hyperbolicExcessSpeedMetresPerSecond: v, gravitationalParameter: gm(mu))
            for (actual,expected) in [(t.circularSpeedMetresPerSecond,vc),
                (t.periapsisSpeedMetresPerSecond,vp),(t.deltaVelocityMetresPerSecond,dv),(t.hyperbolicEccentricity,e)] {
                #expect(abs(actual-expected) < abs(expected)*3e-15)
            }
            #expect(t.characteristicEnergySquareMetresPerSquareSecond == v*v)
        }
        let parabolic = try transfer(0)
        #expect(parabolic.hyperbolicEccentricity == 1)
        #expect(abs(parabolic.deltaVelocityMetresPerSecond-(sqrt(2)-1)) < 1e-15)
        #expect(parabolic.asymptoteAngleRadians == 0)
    }
    @Test(arguments: [ParkingOrbitManeuverKind.escape, .capture])
    func orientedBurnsMatchEnergyAndIndependentAsymptotes(kind: ParkingOrbitManeuverKind) throws {
        let frame = Frame(name: "Local planet frame", pseudoInertial: true)
        let radial = SIMD3<Double>(1,2,2)/3, tangent = SIMD3<Double>(2,-1,0)/sqrt(5)
        let normal = SIMD3<Double>(2,4,-5)/(3*sqrt(5))
        let excess = radial*(kind == .escape ? -0.5 : 0.5)+tangent*(sqrt(3)/2)
        let t = try transfer()
        #expect(abs(t.asymptoteAngleRadians-Double.pi/3) < 1e-15)
        let plan = try t.plan(kind: kind, periapsisEpoch: .j2000, frame: frame, centralBody: .earth,
            hyperbolicExcessVelocityMetresPerSecond: excess, orbitalNormal: normal*3)
        #expect(routeNorm(plan.periapsisRadialDirection-radial) < 1e-14)
        #expect(routeNorm(plan.periapsisProgradeDirection-tangent) < 1e-14)
        let expectedImpulse = tangent*((kind == .escape ? 1 : -1)*(sqrt(3)-1))
        #expect(routeNorm(plan.maneuver.deltaVelocityMetresPerSecond-expectedImpulse) < 1e-14)
        #expect(routeNorm(plan.postBurnState.velocity-plan.preBurnState.velocity-expectedImpulse) < 1e-14)
        let hyperbolic = kind == .escape ? plan.postBurnState : plan.preBurnState
        let circular = kind == .escape ? plan.preBurnState : plan.postBurnState
        func energy(_ s: CartesianState) -> Double { (s.velocity*s.velocity).sum()/2-1/routeNorm(s.position) }
        #expect(abs(energy(hyperbolic)-0.5) < 1e-14)
        #expect(abs(energy(circular)+0.5) < 1e-14)
        let far = try plan.trajectory.state(at: epoch(kind == .escape ? 10000 : -10000))
        #expect(routeNorm(far.velocity-excess) < 1e-3)
        let parkingEpoch = try epoch(kind == .escape ? -2*Double.pi : 2*Double.pi)
        let parking = try plan.trajectory.state(at: parkingEpoch)
        #expect(routeNorm(parking.position-radial) < 1e-8)
        #expect(routeNorm(parking.velocity-tangent) < 1e-8)
        let before = try plan.trajectory.state(at: .j2000, side: .before)
        let after = try plan.trajectory.state(at: .j2000)
        #expect(routeNorm(before.velocity-plan.preBurnState.velocity) < 1e-14)
        #expect(routeNorm(after.velocity-plan.postBurnState.velocity) < 1e-14)
    }
    @Test func safetyAndDirectionValidation() throws {
        let frame = Frame(name: "Local", pseudoInertial: true)
        #expect(throws: ParkingOrbitTransferError.unsafePeriapsis) {
            try ParkingOrbitTransfer(parkingRadiusMetres: 1, bodyRadiusMetres: 0.5,
                minimumPeriapsisAltitudeMetres: 0.6, hyperbolicExcessSpeedMetresPerSecond: 1, gravitationalParameter: gm())
        }
        #expect(throws: ParkingOrbitTransferError.unsafePeriapsis) {
            try ParkingOrbitTransfer(parkingRadiusMetres: 0.5, bodyRadiusMetres: 0.5,
                hyperbolicExcessSpeedMetresPerSecond: 1, gravitationalParameter: gm())
        }
        #expect(throws: ParkingOrbitTransferError.invalidExcessSpeed) { try transfer(-1) }
        #expect(throws: ParkingOrbitTransferError.zeroExcessDirection) {
            try transfer(0).plan(kind: .escape, periapsisEpoch: .j2000, frame: frame, centralBody: .earth,
                hyperbolicExcessVelocityMetresPerSecond: .zero, orbitalNormal: SIMD3(0,0,1))
        }
        #expect(throws: ParkingOrbitTransferError.excessSpeedMismatch) {
            try transfer().plan(kind: .escape, periapsisEpoch: .j2000, frame: frame, centralBody: .earth,
                hyperbolicExcessVelocityMetresPerSecond: SIMD3(2,0,0), orbitalNormal: SIMD3(0,0,1))
        }
        for normal in [SIMD3<Double>.zero, SIMD3(1,0,0), SIMD3(.nan,0,1)] {
            #expect(throws: ParkingOrbitTransferError.invalidDirections) {
                try transfer().plan(kind: .capture, periapsisEpoch: .j2000, frame: frame, centralBody: .earth,
                    hyperbolicExcessVelocityMetresPerSecond: SIMD3(1,0,0), orbitalNormal: normal)
            }
        }
        #expect(throws: PropagationError.nonInertialFrame) {
            try transfer().plan(kind: .escape, periapsisEpoch: .j2000,
                frame: Frame(name: "Rotating", pseudoInertial: false), centralBody: .earth,
                hyperbolicExcessVelocityMetresPerSecond: SIMD3(1,0,0), orbitalNormal: SIMD3(0,0,1))
        }
    }
}

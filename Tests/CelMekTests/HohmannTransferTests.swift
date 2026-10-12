// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct HohmannTransferTests {
    private func gm(_ value: Double = 3.986004418e14) throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: value, convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Explicit test GM",version: "1"))
    }
    private func norm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
    private struct Reference: Decodable {
        let r1,r2,gm,a,e,dv1,dv2,total,seconds: Double
    }
    @Test func independentHighPrecisionReferences() throws {
        let url = try #require(Bundle.module.url(forResource: "hohmann-reference", withExtension: "json", subdirectory: "Fixtures"))
        let rows = try JSONDecoder().decode([Reference].self,from: Data(contentsOf: url))
        #expect(rows.count == 7)
        for row in rows {
            let t = try HohmannTransfer(departureRadiusMetres: row.r1,arrivalRadiusMetres: row.r2,gravitationalParameter: gm(row.gm))
            for (actual,expected) in [(t.semiMajorAxisMetres,row.a),(t.eccentricity,row.e),
                (t.departureDeltaVelocityMetresPerSecond,row.dv1),(t.arrivalDeltaVelocityMetresPerSecond,row.dv2),
                (t.totalDeltaVelocityMetresPerSecond,row.total),(t.duration.seconds,row.seconds)] {
                #expect(abs(actual-expected) <= max(abs(expected)*2e-14,1e-25))
            }
        }
    }
    @Test(arguments: [false,true]) func executablePlansReachCircularOrbits(inward: Bool) throws {
        let r1 = inward ? 42164e3 : 7e6, r2 = inward ? 7e6 : 42164e3
        let transfer = try HohmannTransfer(departureRadiusMetres: r1,arrivalRadiusMetres: r2,gravitationalParameter: gm())
        let frame = Frame(name: "Earth-centred inertial",pseudoInertial: true)
        // Inclined plane and arbitrary start angle; no assumed global XY orbit.
        let radial = SIMD3<Double>(1,2,2)/3, tangent = SIMD3<Double>(2,-1,0)/sqrt(5)
        let plan = try transfer.plan(departureEpoch: .j2000,frame: frame,centralBody: .earth,
            radialDirection: radial*5,progradeDirection: tangent*10)
        #expect((transfer.departureDeltaVelocityMetresPerSecond < 0) == inward)
        #expect((transfer.arrivalDeltaVelocityMetresPerSecond < 0) == inward)
        let before = try plan.trajectory.state(at: .j2000,side: .before)
        #expect(norm(before.position-radial*r1) < 1e-7)
        #expect(norm(before.velocity-tangent*transfer.departureCircularSpeedMetresPerSecond) < 1e-10)
        let arriving = try plan.trajectory.state(at: plan.arrivalEpoch,side: .before)
        let after = try plan.trajectory.state(at: plan.arrivalEpoch)
        #expect(norm(after.position+radial*r2) < 0.002)
        #expect(norm(after.velocity+tangent*transfer.arrivalCircularSpeedMetresPerSecond) < 1e-6)
        #expect(arriving.position == after.position)
        #expect(norm(after.velocity-arriving.velocity-plan.arrivalBurn.deltaVelocityMetresPerSecond) < 1e-10)
        #expect(abs((after.position*after.velocity).sum())/(r2*norm(after.velocity)) < 1e-10)
        let future = try plan.arrivalEpoch.advanced(by: Duration(seconds: 1000))
        #expect(abs(norm(try plan.trajectory.state(at: future).position)-r2) < 0.002)
        let epochs = [.j2000,plan.arrivalEpoch,future]
        let sampled = try TrajectorySampler.sample(plan.trajectory,at: epochs)
        #expect(sampled.map(\.epoch) == epochs)
        for (epoch,state) in zip(epochs.reversed(),sampled.reversed()) {
            #expect(try plan.trajectory.state(at: epoch) == state)
        }
    }
    @Test func equalRadiiAndValidation() throws {
        let transfer = try HohmannTransfer(departureRadiusMetres: 7e6,arrivalRadiusMetres: 7e6,gravitationalParameter: gm())
        #expect(transfer.totalDeltaVelocityMetresPerSecond == 0 && transfer.eccentricity == 0)
        let frame = Frame(name: "inertial",pseudoInertial: true)
        let plan = try transfer.plan(departureEpoch: .j2000,frame: frame,centralBody: .earth,
            radialDirection: SIMD3(1,0,0),progradeDirection: SIMD3(0,-1,0))
        #expect(norm(try plan.trajectory.state(at: plan.arrivalEpoch).position-SIMD3(-7e6,0,0)) < 0.001)
        for radius in [0.0,-1,.nan,.infinity] {
            #expect(throws: HohmannTransferError.invalidRadius) {
                try HohmannTransfer(departureRadiusMetres: radius,arrivalRadiusMetres: 7e6,gravitationalParameter: gm())
            }
            #expect(throws: HohmannTransferError.invalidRadius) {
                try HohmannTransfer(departureRadiusMetres: 7e6,arrivalRadiusMetres: radius,gravitationalParameter: gm())
            }
        }
        let coordinateGM = try GravitationalParameter(cubicMetresPerSecondSquared: 3.986004418e14,
            convention: .tdbCompatible, provenance: PhysicalDataProvenance(source: "coordinate GM",version: "1"))
        #expect(throws: PropagationError.incompatibleGMConvention) {
            try HohmannTransfer(departureRadiusMetres: 7e6,arrivalRadiusMetres: 8e6,gravitationalParameter: coordinateGM)
        }
        #expect(throws: HohmannTransferError.numericalRangeExceeded) {
            try HohmannTransfer(departureRadiusMetres: Double.greatestFiniteMagnitude,
                arrivalRadiusMetres: Double.greatestFiniteMagnitude,gravitationalParameter: gm(1))
        }
        for direction in [SIMD3<Double>.zero,SIMD3(1,0,0),SIMD3(.nan,1,0)] {
            #expect(throws: HohmannTransferError.invalidDirections) {
                try transfer.plan(departureEpoch: .j2000,frame: frame,centralBody: .earth,
                    radialDirection: SIMD3(1,0,0),progradeDirection: direction)
            }
        }
        #expect(throws: PropagationError.nonInertialFrame) {
            try transfer.plan(departureEpoch: .j2000,frame: Frame(name: "rotating",pseudoInertial: false),centralBody: .earth,
                radialDirection: SIMD3(1,0,0),progradeDirection: SIMD3(0,1,0))
        }
    }
}

// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct ExtendedCircularTransferTests {
    func gm() throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: 1, convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Independent normalized references", version: "1"))
    }
    @Test func planeChangeAndCombinedImpulseHaveCorrectVectorsAndStableCosts() throws {
        let pure = try PlaneChangeImpulse(beforeSpeedMetresPerSecond: 5, afterSpeedMetresPerSecond: 5, angleRadians: .pi/3)
        #expect(abs(pure.deltaVelocityMetresPerSecond-5) < 1e-14)
        let combined = try PlaneChangeImpulse(beforeSpeedMetresPerSecond: 3, afterSpeedMetresPerSecond: 4, angleRadians: .pi/2)
        #expect(abs(combined.deltaVelocityMetresPerSecond-5) < 1e-14)
        let frame = Frame(name: "Plane change", pseudoInertial: true)
        for sense in [PlaneChangeSense.positive, .negative] {
            let burn = try combined.maneuver(epoch: .j2000, frame: frame,
                radialDirection: SIMD3(1,0,0), progradeDirection: SIMD3(0,1,0), sense: sense)
            #expect(routeNorm(burn.deltaVelocityMetresPerSecond-SIMD3(0,-3,sense == .positive ? 4 : -4)) < 1e-14)
            #expect(abs(routeNorm(burn.deltaVelocityMetresPerSecond)-combined.deltaVelocityMetresPerSecond) < 1e-14)
        }
        let reversal = try PlaneChangeImpulse(beforeSpeedMetresPerSecond: 3, afterSpeedMetresPerSecond: 4, angleRadians: .pi)
        #expect(abs(reversal.deltaVelocityMetresPerSecond-7) < 1e-14)
        let tiny = try PlaneChangeImpulse(beforeSpeedMetresPerSecond: 1e150, afterSpeedMetresPerSecond: 1e150, angleRadians: 1e-150)
        #expect(abs(tiny.deltaVelocityMetresPerSecond-1) < 1e-14)
        #expect(throws: PlaneChangeError.invalidAngle) {
            try PlaneChangeImpulse(beforeSpeedMetresPerSecond: 1, afterSpeedMetresPerSecond: 1, angleRadians: -0.1)
        }
        #expect(throws: PlaneChangeError.invalidDirections) {
            try pure.maneuver(epoch: .j2000, frame: frame, radialDirection: SIMD3(1,0,0), progradeDirection: SIMD3(1,0,0))
        }
    }
    @Test func independentBiEllipticScalarReferencesAndNearEqualRadii() throws {
        // 55-digit Decimal arithmetic from vis-viva and half-period equations.
        let rows: [(Double,Double,Double,Double,Double,Double,Double)] = [
            (1,20,0.4106912317171965376,0.0230976782991903796,-0.0779045468277846530,3165.1840264250561955,3624.4231896744530571),
            (20,1,0.0779045468277846530,-0.0230976782991903796,-0.4106912317171965376,3624.4231896744530571,3165.1840264250561955)]
        for (r1,r2,d1,d2,d3,t1,t2) in rows {
            let t = try BiEllipticTransfer(departureRadiusMetres: r1, arrivalRadiusMetres: r2,
                intermediateRadiusMetres: 200, gravitationalParameter: gm())
            for (actual,expected) in [(t.departureDeltaVelocityMetresPerSecond,d1),
                (t.intermediateSpeedChangeMetresPerSecond,d2),(t.arrivalDeltaVelocityMetresPerSecond,d3),
                (t.totalDeltaVelocityMetresPerSecond,0.5116934568441715702),
                (t.firstLegDuration.seconds,t1),(t.secondLegDuration.seconds,t2)] {
                #expect(abs(actual-expected) < max(1e-16,abs(expected)*3e-15))
            }
        }
        let near = try BiEllipticTransfer(departureRadiusMetres: 1, arrivalRadiusMetres: 1.nextUp,
            intermediateRadiusMetres: 4, gravitationalParameter: gm())
        #expect(near.intermediateSpeedChangeMetresPerSecond > 0)
        let equal = try BiEllipticTransfer(departureRadiusMetres: 1, arrivalRadiusMetres: 1,
            intermediateRadiusMetres: 4, gravitationalParameter: gm())
        #expect(equal.intermediateSpeedChangeMetresPerSecond == 0)
        #expect(abs(equal.totalDeltaVelocityMetresPerSecond-0.5298221281347034656) < 1e-15)
    }
    @Test(arguments: [0.0, Double.pi/3]) func threeBurnPlansReachExpectedOrbit(angle: Double) throws {
        let frame = Frame(name: "Bi-elliptic", pseudoInertial: true)
        let radial = SIMD3<Double>(1,2,2)/3, tangent = SIMD3<Double>(2,-1,0)/sqrt(5)
        for (r1,r2) in [(1.0,20.0),(20.0,1.0)] {
            let transfer = try BiEllipticTransfer(departureRadiusMetres: r1, arrivalRadiusMetres: r2,
                intermediateRadiusMetres: 200, gravitationalParameter: gm(), planeChangeRadians: angle)
            let plan = try transfer.plan(departureEpoch: .j2000, frame: frame, centralBody: .sun,
                radialDirection: radial, progradeDirection: tangent)
            let apoBefore = try plan.trajectory.state(at: plan.intermediateEpoch, side: .before)
            let apoAfter = try plan.trajectory.state(at: plan.intermediateEpoch)
            #expect(routeNorm(apoBefore.position+radial*200) < 1e-6)
            #expect(apoBefore.position == apoAfter.position)
            #expect(routeNorm(apoAfter.velocity-apoBefore.velocity-plan.intermediateBurn.deltaVelocityMetresPerSecond) < 1e-14)
            let actual = try plan.trajectory.state(at: plan.arrivalEpoch)
            #expect(routeNorm(actual.position-plan.expectedArrivalState.position) < 1e-6)
            #expect(routeNorm(actual.velocity-plan.expectedArrivalState.velocity) < 1e-7)
            #expect(abs((actual.position*actual.velocity).sum()) < 1e-6)
            let spent = [plan.departureBurn,plan.intermediateBurn,plan.arrivalBurn]
                .reduce(0.0) { $0+routeNorm($1.deltaVelocityMetresPerSecond) }
            #expect(abs(spent-transfer.totalDeltaVelocityMetresPerSecond) < 1e-14)
            #expect(plan.arrivalEpoch == (try Epoch.j2000.advanced(by: transfer.duration)))
        }
    }
    @Test func callerObjectivesAndLimitsSelectDifferentStrategies() throws {
        func comparison(_ objective: CircularTransferObjective, time: CelMek.Duration? = nil, delta: Double? = nil) throws -> CircularTransferComparison {
            try CircularTransferComparison(departureRadiusMetres: 1, arrivalRadiusMetres: 20,
                intermediateRadiiMetres: [200,40,200], gravitationalParameter: gm(), objective: objective,
                maximumDuration: time, maximumDeltaVelocityMetresPerSecond: delta)
        }
        let fuel = try comparison(.minimumDeltaVelocity)
        #expect(fuel.rankedCandidates.count == 3)
        if case .biElliptic(let transfer) = fuel.best { #expect(transfer.intermediateRadiusMetres == 200) }
        else { Issue.record("Expected the lower-delta-v bi-elliptic candidate") }
        let fast = try comparison(.minimumDuration)
        if case .hohmann = fast.best {} else { Issue.record("Expected the faster Hohmann candidate") }
        let limit = try HohmannTransfer(departureRadiusMetres: 1, arrivalRadiusMetres: 20, gravitationalParameter: gm()).duration
        let timed = try comparison(.minimumDeltaVelocity, time: limit)
        #expect(timed.rankedCandidates.count == 1 && timed.rejectedCandidates.count == 2)
        let impossible = try comparison(.minimumDeltaVelocity, delta: 0)
        #expect(impossible.best == nil && impossible.rejectedCandidates.count == 3)
        #expect(throws: BiEllipticTransferError.invalidIntermediateRadius) {
            try BiEllipticTransfer(departureRadiusMetres: 1, arrivalRadiusMetres: 20,
                intermediateRadiusMetres: 20, gravitationalParameter: gm())
        }
    }
}

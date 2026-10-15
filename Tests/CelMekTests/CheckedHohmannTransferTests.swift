// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct CheckedHohmannTransferTests {
    let frame = Frame(name: "Checked transfers", pseudoInertial: true)
    let flight = 12.418235332245126543
    func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: t)) }
    func gm() throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: 1, convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Analytic normalized test", version: "1"))
    }
    func orbit(_ radius: Double, phase: Double = 0) throws -> CircularHohmannOrbit {
        try CircularHohmannOrbit(referenceEpoch: .j2000, frame: frame, centralBody: .sun,
            radiusMetres: radius, radialDirection: SIMD3(cos(phase),sin(phase),0),
            progradeDirection: SIMD3(-sin(phase),cos(phase),0), gravitationalParameter: gm())
    }
    struct ChangingProvider: StateProvider {
        let model: CircularHohmannOrbit
        let changeAfter: Epoch
        var frame: Frame { model.frame }
        func state(at epoch: Epoch) throws -> CartesianState {
            let state = try model.state(at: epoch)
            guard epoch > changeAfter else { return state }
            let c = cos(0.01), s = sin(0.01)
            func rotate(_ v: SIMD3<Double>) -> SIMD3<Double> { SIMD3(c*v.x-s*v.y,s*v.x+c*v.y,v.z) }
            return try CartesianState(epoch: epoch, frame: frame, position: rotate(state.position),
                velocity: rotate(state.velocity), acceleration: state.acceleration.map(rotate))
        }
    }
    func provider(_ model: any StateProvider, start: Double = 0, end: Double = 40,
                  origin: BodyIdentifier = .sun) throws -> MovingWaypointRouteTarget {
        try MovingWaypointRouteTarget(provider: model, origin: origin, coverage: epoch(start)...epoch(end))
    }
    @Test func checkedStateProjectionPreservesActualDepartureAndBurnAccounting() throws {
        let actual = try CartesianState(epoch: .j2000, frame: frame, position: SIMD3(1,0,0),
            velocity: SIMD3(1e-9,1+5e-9,0))
        let target = try orbit(4, phase: .pi-flight/8).state(at: .j2000)
        let result = try CheckedHohmannTransferPlanner.plan(departureState: actual, targetState: target,
            centralBody: .sun, gravitationalParameter: gm(), positionToleranceMetres: 1e-8,
            velocityToleranceMetresPerSecond: 1e-8)
        #expect(result.rendezvousSatisfied)
        let before = try result.plan.trajectory.state(at: .j2000, side: .before)
        #expect(before.position == actual.position && before.velocity == actual.velocity)
        #expect(result.departureVelocityProjectionErrorMetresPerSecond > 5e-9)
        let burn = result.plan.departureBurn.deltaVelocityMetresPerSecond
        #expect(abs(burn.x+1e-9) < 1e-16)
        let total = routeNorm(burn)+routeNorm(result.plan.arrivalBurn.deltaVelocityMetresPerSecond)
        #expect(result.totalDeltaVelocityMetresPerSecond == total)
        #expect(abs(total-result.plan.transfer.totalDeltaVelocityMetresPerSecond) > 1e-9)
        let unphased = try CheckedHohmannTransferPlanner.plan(departureState: orbit(1).state(at: .j2000),
            targetState: orbit(4).state(at: .j2000), centralBody: .sun, gravitationalParameter: gm())
        #expect(!unphased.rendezvousSatisfied && unphased.positionErrorMetres > 1)
    }
    @Test func providerCandidatesUseActualArrivalAndRespectCoverageAndBudget() throws {
        let departure = try provider(orbit(1))
        let model = try orbit(4, phase: .pi-flight/8+0.875*2.75)
        let target = try provider(model)
        let result = try CheckedHohmannTransferPlanner.search(departureProvider: departure, targetProvider: target,
            centralBody: .sun, gravitationalParameter: gm(),
            candidateDepartureEpochs: [epoch(2.75),epoch(-1),epoch(0),epoch(2.75)])
        #expect(result.status == .feasible && result.evaluatedCandidates == 3)
        #expect(result.failedAttempts[0].failure == .departureOutsideCoverage)
        #expect(result.failedAttempts[1].failure == .rendezvousResidual)
        let plan = try #require(result.plan)
        #expect(plan.plan.departureBurn.epoch == (try epoch(2.75)))
        #expect(plan.targetArrivalState == (try model.state(at: plan.plan.arrivalEpoch)))
        let changed = try provider(ChangingProvider(model: model, changeAfter: epoch(8)))
        let rejected = try CheckedHohmannTransferPlanner.search(departureProvider: departure, targetProvider: changed,
            centralBody: .sun, gravitationalParameter: gm(), candidateDepartureEpochs: [epoch(2.75)])
        #expect(rejected.status == .candidateSetExhausted)
        #expect(rejected.failedAttempts[0].failure == .rendezvousResidual)
        #expect(try #require(rejected.failedAttempts[0].plan).positionErrorMetres > 0.01)
        let limited = try CheckedHohmannTransferPlanner.search(departureProvider: departure, targetProvider: target,
            centralBody: .sun, gravitationalParameter: gm(), candidateDepartureEpochs: [epoch(0),epoch(2.75)],
            maximumCandidateEvaluations: 1)
        #expect(limited.status == .evaluationLimitReached)
        let uncovered = try CheckedHohmannTransferPlanner.search(departureProvider: departure,
            targetProvider: provider(model, end: 10), centralBody: .sun, gravitationalParameter: gm(),
            candidateDepartureEpochs: [epoch(2.75)])
        #expect(uncovered.failedAttempts[0].failure == .targetOutsideCoverage)
    }
    @Test func relaxedPlaneToleranceDoesNotGuaranteeRendezvous() throws {
        let phase = Double.pi-flight/8, tilt = 1e-5
        let target = try CartesianState(epoch: .j2000, frame: frame,
            position: SIMD3(4*cos(phase),4*sin(phase)*cos(tilt),4*sin(phase)*sin(tilt)),
            velocity: SIMD3(-0.5*sin(phase),0.5*cos(phase)*cos(tilt),0.5*cos(phase)*sin(tilt)))
        #expect(throws: HohmannRendezvousError.nonCoplanarOrbits) {
            try CheckedHohmannTransferPlanner.plan(departureState: orbit(1).state(at: .j2000),
                targetState: target, centralBody: .sun, gravitationalParameter: gm())
        }
        let result = try CheckedHohmannTransferPlanner.plan(departureState: orbit(1).state(at: .j2000),
            targetState: target, centralBody: .sun, gravitationalParameter: gm(),
            validationOptions: CircularHohmannValidationOptions(maximumPlaneNormalDifference: 2e-5),
            positionToleranceMetres: 1e-8, velocityToleranceMetresPerSecond: 1e-8)
        #expect(!result.rendezvousSatisfied)
        #expect(result.velocityErrorMetresPerSecond > 1e-6)
    }
    @Test func unsupportedStatesAndOriginsFailExplicitly() throws {
        let slow = try CartesianState(epoch: .j2000, frame: frame, position: SIMD3(1,0,0), velocity: SIMD3(0,0.8,0))
        #expect(throws: CircularHohmannStateError.nonCircularSpeed) {
            try CircularHohmannOrbit(referenceState: slow, centralBody: .sun, gravitationalParameter: gm())
        }
        let radial = try CartesianState(epoch: .j2000, frame: frame, position: SIMD3(1,0,0), velocity: SIMD3(0.1,1,0))
        #expect(throws: CircularHohmannStateError.nonTangentialVelocity) {
            try CircularHohmannOrbit(referenceState: radial, centralBody: .sun, gravitationalParameter: gm())
        }
        #expect(throws: CircularHohmannStateError.invalidOptions) { try CircularHohmannValidationOptions(relativeSpeedTolerance: .nan) }
        #expect(throws: CircularHohmannStateError.epochMismatch) {
            try CheckedHohmannTransferPlanner.plan(departureState: orbit(1).state(at: epoch(0)),
                targetState: orbit(4).state(at: epoch(1)), centralBody: .sun, gravitationalParameter: gm())
        }
        #expect(throws: HohmannRendezvousError.incompatibleCentralBody) {
            try CircularHohmannOrbit(provider: provider(orbit(1), origin: .earth), at: .j2000,
                centralBody: .sun, gravitationalParameter: gm())
        }
        let coincident = try CheckedHohmannTransferPlanner.search(departureProvider: provider(orbit(1), end: 0),
            targetProvider: provider(orbit(1), end: 0), centralBody: .sun, gravitationalParameter: gm(),
            candidateDepartureEpochs: [.j2000])
        #expect(coincident.status == .alreadyRendezvoused && coincident.plan == nil)
        let separated = try CheckedHohmannTransferPlanner.search(departureProvider: provider(orbit(1)),
            targetProvider: provider(orbit(1, phase: 1)), centralBody: .sun, gravitationalParameter: gm(),
            candidateDepartureEpochs: [.j2000])
        #expect(separated.failedAttempts[0].failure == .sameOrbitPhasingRequired)
    }
}

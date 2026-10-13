// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Math
@testable import CelMek

@Suite struct LambertTransferPlannerTests {
    let frame = Frame(name: "Lambert transfer", pseudoInertial: true)
    func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: t)) }
    func gm() throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: 1, convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Analytic circular transfer", version: "1"))
    }
    struct Linear: StateProvider {
        let frame: Frame
        let p, v: SIMD3<Double>
        func state(at epoch: Epoch) throws -> CartesianState {
            try CartesianState(epoch: epoch, frame: frame,
                position: p + v*epoch.duration(since: .j2000).seconds, velocity: v)
        }
    }
    func provider(_ p: SIMD3<Double>, _ v: SIMD3<Double>, origin: BodyIdentifier = .sun,
        coverage: ClosedRange<Epoch>? = nil) throws -> MovingWaypointRouteTarget {
        try MovingWaypointRouteTarget(provider: Linear(frame: frame,p: p,v: v), origin: origin,
            coverage: coverage ?? epoch(0)...epoch(20))
    }
    func plan(_ condition: LambertArrivalCondition) throws -> LambertTransferPlan {
        let time = Double.pi/2, targetVelocity = SIMD3<Double>(-0.5,0.2,0)
        return try LambertTransferPlanner.plan(departure: provider(SIMD3(1,0,0),SIMD3(0,0.5,0)),
            arrival: provider(SIMD3(0,1,0)-targetVelocity*time,targetVelocity),
            departureEpoch: .j2000, arrivalEpoch: epoch(time), centralBody: .sun,
            gravitationalParameter: gm(), referencePole: SIMD3(0,0,1), direction: .positive,
            arrivalCondition: condition, lambertOptions: LambertOptions(absoluteTimeToleranceSeconds: 1e-13))
    }
    @Test func movingEndpointsAndIndependentVelocityDifferences() throws {
        let intercept = try plan(.intercept), flyby = try plan(.flyby), rendezvous = try plan(.velocityMatchedRendezvous)
        for result in [intercept,flyby,rendezvous] {
            #expect(routeNorm(result.arrivalTargetState.position-SIMD3(0,1,0)) < 1e-12)
            #expect(routeNorm(result.departureBurn.deltaVelocityMetresPerSecond-SIMD3(0,0.5,0)) < 1e-10)
            #expect(routeNorm(result.incomingRelativeVelocityMetresPerSecond-SIMD3(-0.5,-0.2,0)) < 1e-10)
            #expect(result.positionErrorMetres < 1e-10)
            #expect(result.arrivalConditionSatisfied)
            let beforeDeparture = try result.trajectory.state(at: .j2000,side: .before)
            let afterDeparture = try result.trajectory.state(at: .j2000)
            #expect(beforeDeparture == result.departureReferenceState)
            #expect(routeNorm(afterDeparture.velocity-beforeDeparture.velocity-result.departureBurn.deltaVelocityMetresPerSecond) < 1e-12)
        }
        #expect(intercept.arrivalBurn == nil && flyby.arrivalBurn == nil)
        #expect(intercept.maneuvers.count == 1 && flyby.maneuvers.count == 1)
        #expect(intercept.velocityErrorMetresPerSecond > 0.5)
        #expect(flyby.postArrivalState.velocity == flyby.preArrivalState.velocity)
        #expect(abs(flyby.totalDeltaVelocityMetresPerSecond-0.5) < 1e-10)
        let burn = try #require(rendezvous.arrivalBurn)
        #expect(routeNorm(burn.deltaVelocityMetresPerSecond-SIMD3(0.5,0.2,0)) < 1e-10)
        #expect(abs(rendezvous.totalDeltaVelocityMetresPerSecond-(0.5+sqrt(0.29))) < 1e-10)
        #expect(rendezvous.maneuvers.count == 2)
        let before = try rendezvous.trajectory.state(at: rendezvous.arrivalTargetState.epoch,side: .before)
        let after = try rendezvous.trajectory.state(at: rendezvous.arrivalTargetState.epoch)
        #expect(before.position == after.position)
        #expect(routeNorm(after.velocity-before.velocity-burn.deltaVelocityMetresPerSecond) < 1e-12)
        #expect(rendezvous.velocityErrorMetresPerSecond < 1e-12)
    }
    @Test func stationaryTargetAndNoHistoryExtrapolation() throws {
        let depart = try CartesianState(epoch: .j2000,frame: frame,position: SIMD3(1,0,0),velocity: .zero)
        let arrival = try CartesianState(epoch: epoch(.pi/2),frame: frame,position: SIMD3(0,1,0),velocity: .zero)
        let plan = try LambertTransferPlanner.plan(departureState: depart, arrivalTargetState: arrival,
            centralBody: .sun,gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),direction: .positive,
            arrivalCondition: .velocityMatchedRendezvous)
        #expect(plan.postArrivalState.velocity == .zero)
        #expect(plan.arrivalConditionSatisfied)
        #expect(throws: LambertTransferPlanningError.outsideFlight) { try plan.trajectory.state(at: epoch(-1)) }
        #expect(throws: LambertTransferPlanningError.outsideFlight) { try plan.trajectory.state(at: epoch(3)) }
    }
    @Test func residualDoesNotResetArrivalPosition() throws {
        let departure = try CartesianState(epoch: .j2000,frame: frame,position: SIMD3(1,0,0),velocity: .zero)
        let target = try CartesianState(epoch: epoch(2),frame: frame,position: SIMD3(0,1,0),velocity: SIMD3(0.2,0.3,0))
        let result = try LambertTransferPlanner.plan(departureState: departure,arrivalTargetState: target,
            centralBody: .sun,gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),direction: .positive,
            arrivalCondition: .velocityMatchedRendezvous,positionToleranceMetres: 1e-25,
            lambertOptions: LambertOptions(absoluteTimeToleranceSeconds: 1e-6,relativeTimeTolerance: 0))
        #expect(result.positionErrorMetres > 1e-25)
        #expect(!result.arrivalConditionSatisfied)
        #expect(result.preArrivalState.position == result.postArrivalState.position)
        #expect(result.postArrivalState.position != target.position)
        #expect(result.velocityErrorMetresPerSecond < 1e-12)
    }
    @Test func multiRevolutionBranchIsRetained() throws {
        let time = 2 * Double.pi + Double.pi/2
        for branch in [LambertMultiRevolutionBranch.lowerUniversalParameter,.upperUniversalParameter] {
            let departure = try CartesianState(epoch: .j2000,frame: frame,position: SIMD3(1,0,0),velocity: SIMD3(0,1,0))
            let arrival = try CartesianState(epoch: epoch(time),frame: frame,position: SIMD3(0,1,0),velocity: SIMD3(-1,0,0))
            let plan = try LambertTransferPlanner.plan(departureState: departure,arrivalTargetState: arrival,
                centralBody: .sun,gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),direction: .positive,
                arrivalCondition: .flyby,revolutions: 1,multiRevolutionBranch: branch)
            #expect(plan.solution.multiRevolutionBranch == branch)
            #expect(plan.solution.revolutions == 1)
            #expect(plan.positionErrorMetres < 1e-8)
            #expect(plan.arrivalBurn == nil)
        }
    }
    @Test func ephemerisAdapterUsesMetadataAndArrivalEpoch() throws {
        func metadata(_ body: BodyIdentifier) throws -> EphemerisMetadata {
            try EphemerisMetadata(body: body,origin: .sun,frame: frame,coverage: epoch(0)...epoch(20),
                coverageBasis: .callerSelected,modelTimeScale: .tt,velocityMethod: .analytic,
                accelerationMethod: .unavailable,provenance: PhysicalDataProvenance(source: "Test linear ephemeris",version: "1"))
        }
        let time = Double.pi/2, velocity = SIMD3<Double>(-0.5,0.2,0)
        let departure = try StateProviderEphemeris(provider: Linear(frame: frame,p: SIMD3(1,0,0),v: SIMD3(0,0.5,0)),metadata: metadata(.earth))
        let target = try StateProviderEphemeris(provider: Linear(frame: frame,p: SIMD3(0,1,0)-velocity*time,v: velocity),metadata: metadata(.mars))
        let fromEphemerides = try LambertTransferPlanner.plan(departure: departure,arrival: target,
            departureEpoch: .j2000,arrivalEpoch: epoch(time),centralBody: .sun,gravitationalParameter: gm(),
            referencePole: SIMD3(0,0,1),direction: .positive,arrivalCondition: .flyby,
            lambertOptions: LambertOptions(absoluteTimeToleranceSeconds: 1e-13))
        let reference = try plan(.flyby)
        #expect(fromEphemerides.arrivalTargetState == reference.arrivalTargetState)
        #expect(fromEphemerides.departureReferenceState == reference.departureReferenceState)
        #expect(fromEphemerides.departureBurn.deltaVelocityMetresPerSecond == reference.departureBurn.deltaVelocityMetresPerSecond)
    }
    @Test func providerContractFailures() throws {
        let departure = try provider(SIMD3(1,0,0),.zero), target = try provider(SIMD3(0,1,0),.zero,origin: .earth)
        #expect(throws: LambertTransferPlanningError.originMismatch) {
            try LambertTransferPlanner.plan(departure: departure,arrival: target,departureEpoch: .j2000,arrivalEpoch: epoch(1),
                centralBody: .sun,gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),direction: .positive,arrivalCondition: .intercept)
        }
        let short = try provider(SIMD3(0,1,0),.zero,coverage: epoch(0)...epoch(0.5))
        #expect(throws: MovingWaypointTargetError.outsideCoverage) {
            try LambertTransferPlanner.plan(departure: departure,arrival: short,departureEpoch: .j2000,arrivalEpoch: epoch(1),
                centralBody: .sun,gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),direction: .positive,arrivalCondition: .flyby)
        }
        let otherFrame = Frame(name: frame.name,pseudoInertial: true)
        let wrong = try MovingWaypointRouteTarget(provider: Linear(frame: otherFrame,p: SIMD3(0,1,0),v: .zero),origin: .sun,
            coverage: epoch(0)...epoch(20))
        #expect(throws: SpacecraftDynamicsError.frameMismatch) {
            try LambertTransferPlanner.plan(departure: departure,arrival: wrong,departureEpoch: .j2000,arrivalEpoch: epoch(1),
                centralBody: .sun,gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),direction: .positive,arrivalCondition: .flyby)
        }
    }
}

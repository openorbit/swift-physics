// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct HohmannRendezvousTests {
    let frame = Frame(name: "Circular rendezvous", pseudoInertial: true)
    let radial = SIMD3<Double>(1,2,2)/3
    let tangent = SIMD3<Double>(2,-1,0)/sqrt(5)
    // Independent 55-digit Decimal calculations: a=2.5, GM=1,
    // half-period = pi*sqrt(a^3), synodic period = 2*pi/(1-1/8).
    let flight = 12.418235332245126543041938337900438906
    let synodic = 7.180783208205241687914613447496006592
    func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: t)) }
    func gm(_ value: Double = 1) throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: value, convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Dimensionless analytic test", version: "1"))
    }
    func orbit(_ radius: Double, phase: Double = 0, reference: Double = 0, reversed: Bool = false,
               body: BodyIdentifier = .sun) throws -> CircularHohmannOrbit {
        let t = reversed ? -tangent : tangent
        return try CircularHohmannOrbit(referenceEpoch: epoch(reference), frame: frame, centralBody: body,
            radiusMetres: radius, radialDirection: radial*cos(phase)+t*sin(phase),
            progradeDirection: t*cos(phase)-radial*sin(phase), gravitationalParameter: gm())
    }
    func pair(inward: Bool, launch: Double = 2.75) throws -> (CircularHohmannOrbit, CircularHohmannOrbit) {
        let r1 = inward ? 4.0 : 1.0, r2 = inward ? 1.0 : 4.0
        let n1 = inward ? 0.125 : 1.0, n2 = inward ? 1.0 : 0.125
        // Target has a different reference epoch. Phase required at launch:
        // thetaTarget(launch)-thetaDeparture(launch) = pi-n2*flight.
        let phaseAtReference = .pi-n2*flight+(n1-n2)*launch+n2*3
        return try (orbit(r1), orbit(r2, phase: phaseAtReference, reference: 3))
    }
    @Test(arguments: [false,true]) func independentlyPhasedTransfersReachMovingTarget(inward: Bool) throws {
        let (departure,target) = try pair(inward: inward)
        let result = try HohmannRendezvousPlanner.search(departureOrbit: departure, targetOrbit: target,
            departureWindow: epoch(0)...epoch(6), positionToleranceMetres: 1e-8,
            velocityToleranceMetresPerSecond: 1e-8)
        #expect(result.status == .feasible)
        let launch = try #require(result.departureEpoch)
        #expect(abs(try launch.duration(since: .j2000).seconds-2.75) < 1e-12)
        let expectedPhase = inward ? 3.2897279357038396493 : 1.5893132370591524206
        #expect(abs(result.requiredPhaseAngleRadians-expectedPhase) < 2e-15)
        #expect(abs(result.relativeAngularRateRadiansPerSecond-(inward ? 0.875 : -0.875)) < 1e-15)
        let plan = try #require(result.plan)
        #expect(abs(plan.transfer.duration.seconds-flight) < 1e-14)
        let actual = try plan.trajectory.state(at: plan.arrivalEpoch)
        let expected = try target.state(at: plan.arrivalEpoch)
        #expect(routeNorm(actual.position-expected.position) < 1e-8)
        #expect(routeNorm(actual.velocity-expected.velocity) < 1e-8)
        #expect(result.targetArrivalState == expected)
        #expect(abs(result.attempts[0].phaseErrorRadians) < 1e-14)
        // Before ignition the plan remains on the actual departure circular orbit.
        let before = try plan.trajectory.state(at: launch, side: .before)
        let circular = try departure.state(at: launch)
        #expect(routeNorm(before.position-circular.position) < 1e-14)
        #expect(routeNorm(before.velocity-circular.velocity) < 1e-14)
    }
    @Test func inclusiveWindowsAndSynodicRepetition() throws {
        let (departure,target) = try pair(inward: false)
        let exact = try HohmannRendezvousPlanner.search(departureOrbit: departure, targetOrbit: target,
            departureWindow: epoch(2.75)...epoch(2.75))
        #expect(exact.status == .feasible); #expect(exact.departureEpoch == (try epoch(2.75)))
        let next = try HohmannRendezvousPlanner.search(departureOrbit: departure, targetOrbit: target,
            departureWindow: epoch(3)...epoch(11))
        #expect(abs(try #require(next.departureEpoch).duration(since: .j2000).seconds-(2.75+synodic)) < 1e-12)
        let negative = try HohmannRendezvousPlanner.search(departureOrbit: departure, targetOrbit: target,
            departureWindow: epoch(-5)...epoch(0))
        #expect(abs(try #require(negative.departureEpoch).duration(since: .j2000).seconds-(2.75-synodic)) < 1e-12)
        let missed = try HohmannRendezvousPlanner.search(departureOrbit: departure, targetOrbit: target,
            departureWindow: epoch(2.75001)...epoch(2.8))
        #expect(missed.status == .departureWindowExhausted)
        #expect(missed.attempts.isEmpty && missed.plan == nil)
    }
    @Test func equalRadiiAreNotArtificialRendezvous() throws {
        let departure = try orbit(1)
        let same = try HohmannRendezvousPlanner.search(departureOrbit: departure,
            targetOrbit: orbit(1, phase: 2, reference: 2), departureWindow: epoch(1)...epoch(100))
        #expect(same.status == .alreadyRendezvoused); #expect(same.plan == nil)
        #expect(same.departureEpoch == (try epoch(1)))
        let separated = try HohmannRendezvousPlanner.search(departureOrbit: departure,
            targetOrbit: orbit(1, phase: 1), departureWindow: epoch(0)...epoch(100))
        #expect(separated.status == .sameOrbitPhasingRequired)
        #expect(separated.plan == nil && separated.departureEpoch == nil)
        #expect(separated.relativeAngularRateRadiansPerSecond == 0)
    }
    @Test func residualFailuresAndBudgetsRemainExplicit() throws {
        let (departure,target) = try pair(inward: false)
        let bounded = try HohmannRendezvousPlanner.search(departureOrbit: departure, targetOrbit: target,
            departureWindow: epoch(0)...epoch(15), maximumOpportunityEvaluations: 1,
            positionToleranceMetres: 1e-30, velocityToleranceMetresPerSecond: 1e-30)
        #expect(bounded.status == .evaluationLimitReached)
        #expect(bounded.attempts.count == 1 && bounded.plan == nil)
        #expect(bounded.attempts[0].positionErrorMetres > 1e-30 ||
            bounded.attempts[0].velocityErrorMetresPerSecond > 1e-30)
        let exhausted = try HohmannRendezvousPlanner.search(departureOrbit: departure, targetOrbit: target,
            departureWindow: epoch(0)...epoch(6), maximumOpportunityEvaluations: 1,
            positionToleranceMetres: 1e-30, velocityToleranceMetresPerSecond: 1e-30)
        #expect(exhausted.status == .departureWindowExhausted)
        #expect(throws: HohmannRendezvousError.invalidBudget) {
            try HohmannRendezvousPlanner.search(departureOrbit: departure, targetOrbit: target,
                departureWindow: epoch(0)...epoch(1), maximumOpportunityEvaluations: 0)
        }
        #expect(throws: HohmannRendezvousError.invalidTolerance) {
            try HohmannRendezvousPlanner.search(departureOrbit: departure, targetOrbit: target,
                departureWindow: epoch(0)...epoch(1), positionToleranceMetres: .nan)
        }
    }
    @Test func modelCompatibilityAndCircularDerivatives() throws {
        let departure = try orbit(1)
        #expect(throws: HohmannRendezvousError.incompatibleCentralBody) {
            try HohmannRendezvousPlanner.search(departureOrbit: departure,
                targetOrbit: orbit(4, body: .earth), departureWindow: epoch(0)...epoch(10))
        }
        #expect(throws: HohmannRendezvousError.opposingOrbitalDirections) {
            try HohmannRendezvousPlanner.search(departureOrbit: departure,
                targetOrbit: orbit(4, reversed: true), departureWindow: epoch(0)...epoch(10))
        }
        let tilted = try CircularHohmannOrbit(referenceEpoch: .j2000, frame: frame, centralBody: .sun,
            radiusMetres: 4, radialDirection: radial, progradeDirection: SIMD3(-4,-4,6)/sqrt(68),
            gravitationalParameter: gm())
        #expect(throws: HohmannRendezvousError.nonCoplanarOrbits) {
            try HohmannRendezvousPlanner.search(departureOrbit: departure, targetOrbit: tilted,
                departureWindow: epoch(0)...epoch(10))
        }
        let otherGM = try CircularHohmannOrbit(referenceEpoch: .j2000, frame: frame, centralBody: .sun,
            radiusMetres: 4, radialDirection: radial, progradeDirection: tangent, gravitationalParameter: gm(2))
        #expect(throws: HohmannRendezvousError.incompatibleGM) {
            try HohmannRendezvousPlanner.search(departureOrbit: departure, targetOrbit: otherGM,
                departureWindow: epoch(0)...epoch(10))
        }
        let otherFrame = try CircularHohmannOrbit(referenceEpoch: .j2000,
            frame: Frame(name: "Circular rendezvous", pseudoInertial: true), centralBody: .sun,
            radiusMetres: 4, radialDirection: radial, progradeDirection: tangent, gravitationalParameter: gm())
        #expect(throws: SpacecraftDynamicsError.frameMismatch) {
            try HohmannRendezvousPlanner.search(departureOrbit: departure, targetOrbit: otherFrame,
                departureWindow: epoch(0)...epoch(10))
        }
        let at = try departure.state(at: epoch(0.7)), h = 1e-4
        let plus = try departure.state(at: epoch(0.7+h)), minus = try departure.state(at: epoch(0.7-h))
        #expect(routeNorm((plus.position-minus.position)/(2*h)-at.velocity) < 1e-8)
        #expect(routeNorm((plus.velocity-minus.velocity)/(2*h)-(try #require(at.acceleration))) < 1e-8)
        let near = try HohmannRendezvousPlanner.search(departureOrbit: departure,
            targetOrbit: orbit(1+1e-10), departureWindow: epoch(0)...epoch(0))
        #expect(near.relativeAngularRateRadiansPerSecond < 0)
        #expect(abs(near.relativeAngularRateRadiansPerSecond) < 2e-10)
    }
}

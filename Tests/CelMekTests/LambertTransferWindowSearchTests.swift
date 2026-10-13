// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct LambertTransferWindowSearchTests {
    let frame = Frame(name: "Transfer windows",pseudoInertial: true)
    func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: t)) }
    func gm(_ value: Double = 1) throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: value,convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Explicit test model",version: "1"))
    }
    struct Linear: StateProvider {
        let frame: Frame
        let p,v: SIMD3<Double>
        func state(at epoch: Epoch) throws -> CartesianState {
            try CartesianState(epoch: epoch,frame: frame,position: p+v*epoch.duration(since: .j2000).seconds,velocity: v)
        }
    }
    func sources() throws -> (MovingWaypointRouteTarget,MovingWaypointRouteTarget) {
        try (MovingWaypointRouteTarget(provider: Linear(frame: frame,p: SIMD3(1,0,0),v: SIMD3(0,0.5,0)),origin: .sun,
                coverage: epoch(0)...epoch(10)),
            MovingWaypointRouteTarget(provider: Linear(frame: frame,p: SIMD3(0,1,0),v: .zero),origin: .sun,
                coverage: epoch(0)...epoch(10)))
    }
    func search(budget: Int = 1000,objective: LambertWindowObjective = .minimumDeltaVelocity,
        cost: Double? = nil,speed: Double? = nil) throws -> LambertWindowSearchResult {
        let (departure,arrival) = try sources()
        let positive = try LambertWindowBranch(direction: .positive), negative = try LambertWindowBranch(direction: .negative)
        return try LambertTransferWindowSearch.search(departure: departure,arrival: arrival,centralBody: .sun,
            gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),departureEpochs: [epoch(0.5),.j2000,.j2000],
            flightTimes: [CelMek.Duration(seconds: .pi/2),CelMek.Duration(seconds: .pi/4),CelMek.Duration(seconds: .pi/2)],
            branches: [positive,negative,positive],arrivalCondition: .velocityMatchedRendezvous,objective: objective,
            maximumCandidateEvaluations: budget,maximumDeltaVelocityMetresPerSecond: cost,
            maximumArrivalRelativeSpeedMetresPerSecond: speed)
    }
    @Test func gridDeduplicationRankingAndBudgets() throws {
        let all = try search()
        #expect(all.totalCandidates == 8 && all.evaluatedCandidates == 8)
        #expect(all.status == .gridExhausted)
        #expect(all.rows.allSatisfy { $0.feasible })
        #expect(all.rows.first?.candidate.departureEpoch == .j2000)
        let costs = all.rankedFeasibleRows.compactMap(\.totalDeltaVelocityMetresPerSecond)
        #expect(costs == costs.sorted())
        let fastest = try search(objective: .minimumFlightTime)
        let minimumTime = try CelMek.Duration(seconds: .pi/4)
        #expect(fastest.best?.candidate.flightTime == minimumTime)
        let limited = try search(budget: 1)
        #expect(limited.status == .evaluationLimitReached && limited.evaluatedCandidates == 1)
        #expect(limited.best != nil && limited.totalCandidates == 8)
        #expect(try search(budget: 8).status == .gridExhausted)
        for row in all.rows {
            let plan = try #require(row.plan)
            let endpoint = try plan.trajectory.state(at: plan.arrivalTargetState.epoch)
            #expect(routeNorm(endpoint.position-plan.arrivalTargetState.position) < 1e-7)
            #expect(routeNorm(endpoint.velocity-plan.arrivalTargetState.velocity) < 1e-7)
            #expect(row.arrivalRelativeSpeedMetresPerSecond! > 0)
        }
    }
    @Test func limitsKeepRejectedPlansAndIncomingSpeed() throws {
        let expensive = try search(cost: 0)
        #expect(expensive.best == nil)
        #expect(expensive.rows.allSatisfy { $0.failure == .deltaVelocityLimit && $0.plan != nil })
        let fast = try search(speed: 0)
        #expect(fast.best == nil)
        #expect(fast.rows.allSatisfy { $0.failure == .arrivalRelativeSpeedLimit && $0.plan!.velocityErrorMetresPerSecond < 1e-7 })
        let reference = try search(budget: 1)
        let inclusive = try search(budget: 1,cost: reference.best!.totalDeltaVelocityMetresPerSecond!,
            speed: reference.best!.arrivalRelativeSpeedMetresPerSecond!)
        #expect(inclusive.best != nil)
    }
    @Test func coverageGeometryAndMultiRevolutionFailures() throws {
        let (departure,arrival) = try sources()
        let zero = try LambertWindowBranch(direction: .positive)
        let coverage = try LambertTransferWindowSearch.search(departure: departure,arrival: arrival,centralBody: .sun,
            gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),departureEpochs: [epoch(-1),epoch(9)],
            flightTimes: [CelMek.Duration(seconds: 2)],branches: [zero],arrivalCondition: .flyby)
        #expect(coverage.rows[0].failure == .departureOutsideCoverage)
        #expect(coverage.rows[1].failure == .arrivalOutsideCoverage)
        let collinear = try MovingWaypointRouteTarget(provider: Linear(frame: frame,p: SIMD3(-1,0,0),v: .zero),origin: .sun,
            coverage: epoch(0)...epoch(10))
        let failed = try LambertTransferWindowSearch.search(departure: departure,arrival: collinear,centralBody: .sun,
            gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),departureEpochs: [.j2000],
            flightTimes: [CelMek.Duration(seconds: 2)],branches: [zero],arrivalCondition: .intercept)
        #expect(failed.rows[0].failure == .lambert(.degenerateGeometry))
        let lower = try LambertWindowBranch(direction: .positive,revolutions: 1,multiRevolutionBranch: .lowerUniversalParameter)
        let upper = try LambertWindowBranch(direction: .positive,revolutions: 1,multiRevolutionBranch: .upperUniversalParameter)
        let multi = try LambertTransferWindowSearch.search(departure: departure,arrival: arrival,centralBody: .sun,
            gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),departureEpochs: [.j2000],
            flightTimes: [CelMek.Duration(seconds: 1),CelMek.Duration(seconds: 8)],branches: [lower,upper],arrivalCondition: .flyby)
        #expect(multi.rows.count == 4 && multi.rankedFeasibleRows.count == 2)
        if case .lambert(.belowMinimumFlightTime) = multi.rows[0].failure {} else { Issue.record("Expected minimum-time failure") }
        #expect(multi.rows[2].candidate.branch.multiRevolutionBranch == .lowerUniversalParameter)
        #expect(multi.rows[3].candidate.branch.multiRevolutionBranch == .upperUniversalParameter)
    }
    struct Broken: StateProvider {
        let frame: Frame
        func state(at epoch: Epoch) throws -> CartesianState { throw LambertError.invalidInput }
    }
    @Test func configurationAndProviderErrorsPropagate() throws {
        let (departure,arrival) = try sources()
        #expect(throws: LambertWindowSearchError.invalidBudget) { try search(budget: 0) }
        #expect(throws: LambertWindowSearchError.invalidConstraints) { try search(cost: -.infinity) }
        let broken = try MovingWaypointRouteTarget(provider: Broken(frame: frame),origin: .sun,coverage: epoch(0)...epoch(10))
        #expect(throws: LambertError.invalidInput) {
            try LambertTransferWindowSearch.search(departure: broken,arrival: arrival,centralBody: .sun,
                gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),departureEpochs: [.j2000],
                flightTimes: [CelMek.Duration(seconds: 1)],branches: [LambertWindowBranch(direction: .positive)],arrivalCondition: .flyby)
        }
        #expect(throws: LambertWindowSearchError.invalidGrid) {
            try LambertTransferWindowSearch.search(departure: departure,arrival: arrival,centralBody: .sun,
                gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),departureEpochs: Array(repeating: .j2000,count: 400),
                flightTimes: Array(repeating: CelMek.Duration(seconds: 1),count: 300),
                branches: [LambertWindowBranch(direction: .positive)],arrivalCondition: .flyby)
        }
    }
    @Test func epochRangeAndResolutionFailuresAreRows() throws {
        let (departure,arrival) = try sources()
        let branch = try LambertWindowBranch(direction: .positive)
        let maximum = Epoch(offsetFromJ2000: try CelMek.Duration(wholeSeconds: CelMek.Duration.maximumWholeSeconds))
        let range = try LambertTransferWindowSearch.search(departure: departure,arrival: arrival,centralBody: .sun,
            gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),departureEpochs: [maximum],
            flightTimes: [CelMek.Duration(seconds: 2)],branches: [branch],arrivalCondition: .flyby)
        #expect(range.rows[0].failure == .timeArithmetic(.outOfRange))
        #expect(range.rows[0].arrivalEpoch == nil)
        let resolution = try LambertTransferWindowSearch.search(departure: departure,arrival: arrival,centralBody: .sun,
            gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),departureEpochs: [epoch(0.5)],
            flightTimes: [CelMek.Duration(seconds: 1e-30)],branches: [branch],arrivalCondition: .flyby)
        #expect(resolution.rows[0].failure == .arrivalEpochResolutionExceeded)
        #expect(resolution.best == nil && resolution.status == .gridExhausted)
        #expect(throws: LambertWindowSearchError.invalidPole) {
            try LambertTransferWindowSearch.search(departure: departure,arrival: arrival,centralBody: .sun,
                gravitationalParameter: gm(),referencePole: .zero,departureEpochs: [maximum],
                flightTimes: [CelMek.Duration(seconds: 2)],branches: [branch],arrivalCondition: .flyby)
        }
    }
    @Test func interplanetaryKeplerEphemerisGrid() throws {
        // Ideal heliocentric Earth/Mars-sized circular models, not a dated mission prediction.
        let mu = 1.32712440018e20, au = CartesianUnits.metresPerAstronomicalUnit
        let parameter = try gm(mu), coverage = try epoch(0)...epoch(600*86400)
        func body(_ id: BodyIdentifier,_ radius: Double,_ phase: Double) throws -> KeplerEphemeris {
            let state = try CartesianState(epoch: .j2000,frame: frame,
                position: SIMD3(cos(phase),sin(phase),0)*radius,
                velocity: SIMD3(-sin(phase),cos(phase),0)*sqrt(mu/radius))
            return try KeplerEphemeris(body: id,propagator: TwoBodyPropagator(initialState: state,centralBody: .sun,
                gravitationalParameter: parameter),coverage: coverage,
                provenance: PhysicalDataProvenance(source: "Ideal circular demonstration",version: "1"))
        }
        let departure = try body(.earth,au,0), arrival = try body(.marsBarycentre,1.524*au,0.7)
        let result = try LambertTransferWindowSearch.search(departure: departure,arrival: arrival,centralBody: .sun,
            gravitationalParameter: parameter,referencePole: SIMD3(0,0,1),departureEpochs: [epoch(0),epoch(30*86400)],
            flightTimes: [CelMek.Duration(seconds: 180*86400),CelMek.Duration(seconds: 240*86400)],
            branches: [LambertWindowBranch(direction: .positive)],arrivalCondition: .velocityMatchedRendezvous,
            positionToleranceMetres: 100,velocityToleranceMetresPerSecond: 1e-4)
        #expect(result.status == .gridExhausted && result.rows.count == 4)
        #expect(result.rankedFeasibleRows.count == 4)
        for row in result.rows {
            let plan = try #require(row.plan)
            let sampled = try arrival.state(at: plan.arrivalTargetState.epoch)
            #expect(sampled == plan.arrivalTargetState)
            #expect(row.totalDeltaVelocityMetresPerSecond! > 0)
            #expect(row.departureRelativeSpeedMetresPerSecond! > 0)
            #expect(row.arrivalRelativeSpeedMetresPerSecond! > 0)
        }
    }
}

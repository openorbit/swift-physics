// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct MilestoneSixPlanningTests {
    let frame = Frame(name: "M6 planning", pseudoInertial: true)
    func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: t)) }
    struct LinearTarget: StateProvider {
        let frame: Frame
        let p: SIMD3<Double>
        let v: SIMD3<Double>
        func state(at epoch: Epoch) throws -> CartesianState {
            try CartesianState(epoch: epoch, frame: frame,
                position: p + v*epoch.duration(since: .j2000).seconds, velocity: v, acceleration: .zero)
        }
    }
    func target(p: SIMD3<Double>, v: SIMD3<Double> = .zero) throws -> MovingWaypointRouteTarget {
        try MovingWaypointRouteTarget(provider: LinearTarget(frame: frame, p: p, v: v),
            origin: .sun, coverage: epoch(0)...epoch(100))
    }
    func mass() throws -> SpacecraftMassProperties {
        try SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50)
    }
    func background() throws -> SpacecraftDynamics {
        try SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun))
    }
    func limits() throws -> WaypointRouteEngineLimits {
        try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 1,
            maximumThrustNewtons: 1000, effectiveExhaustVelocityMetresPerSecond: 100)
    }
    @Test func relativeOffsetsHaveIndependentRotatingDerivatives() throws {
        let tree = FrameTree()
        let translated = try tree.makeFrame(parent: tree.root, name: "Unrelated origin",
            provider: FixedTransformProvider(translation: SIMD3(1000,-2000,500)))
        let axes = try tree.makeFrame(parent: translated, name: "Offset axes",
            provider: UniformRotationProvider(axis: SIMD3(0,0,1), rate: 0.5))
        let native = try MovingWaypointRouteTarget(provider:
            LinearTarget(frame: tree.root, p: SIMD3(10,0,0), v: SIMD3(1,0,0)),
            origin: .sun, coverage: epoch(0)...epoch(10))
        let shifted = try native.offset(positionMetres: SIMD3(2,0,0), in: axes, using: tree)
        let t = 2.0, theta = 0.5*t
        let state = try shifted.state(at: epoch(t))
        #expect(routeNorm(state.position-SIMD3(12+2*cos(theta),2*sin(theta),0)) < 1e-12)
        #expect(routeNorm(state.velocity-SIMD3(1-sin(theta),cos(theta),0)) < 1e-12)
        #expect(routeNorm(try #require(state.acceleration)-SIMD3(-0.5*cos(theta),-0.5*sin(theta),0)) < 1e-12)
        let fixed = try native.offset(positionMetres: SIMD3(0,3,0)).state(at: epoch(t))
        #expect(fixed.position == SIMD3(12,3,0)); #expect(fixed.velocity == SIMD3(1,0,0))
        #expect(shifted.origin == .sun); #expect(shifted.coverage == native.coverage)
        #expect(throws: CartesianStateError.nonFinitePosition) { try native.offset(positionMetres: SIMD3(.nan,0,0)) }
        #expect(throws: MovingWaypointTargetError.outsideCoverage) { try shifted.state(at: epoch(11)) }
    }
    @Test func relativePassageRetainsFlyThroughAndMatchesOffsetRendezvous() throws {
        let shifted = try target(p: SIMD3(1,0,0), v: SIMD3(0.1,0,0)).offset(positionMetres: SIMD3(0,1,0))
        let result = try WaypointRouteTimingPlanner.search(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: [BrachistochroneWaypoint(positionMetres: .zero),
                BrachistochroneWaypoint(positionMetres: SIMD3(-99,0,0)),
                BrachistochroneWaypoint(positionMetres: SIMD3(6,0,0))],
            legAccelerationsMetresPerSecondSquared: [1,1], candidateArrivalEpochs: [[epoch(10)],[epoch(20)]],
            initialMass: mass(), backgroundDynamics: background(), legEngineLimits: [limits(),limits()],
            movingWaypoints: [1: MovingWaypointPassage(provider: shifted)])
        let route = try #require(result.route)
        #expect(route.waypoints[1].positionMetres == SIMD3(2,1,0))
        #expect(route.waypointVelocitiesMetresPerSecond[1].x > 0)
        #expect(try shifted.waypoint(at: epoch(10)).velocityMetresPerSecond == SIMD3(0.1,0,0))
    }
    @Test func multipleSeedsShareBudgetAndKeepEarliestFeasibleResult() throws {
        func run(_ budget: Int) throws -> WaypointTimingMultiStartResult {
            try WaypointRouteTimingPlanner.optimizeTimingFromSeeds(departureEpoch: .j2000, frame: frame, origin: .sun,
                waypoints: [BrachistochroneWaypoint(positionMetres: .zero),
                    BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0))],
                legAccelerationsMetresPerSecondSquared: [1], arrivalWindows: [epoch(4)...epoch(12)],
                seedArrivalEpochs: [[epoch(4)],[epoch(10)],[epoch(10)]], maximumTotalCandidateEvaluations: budget,
                optimizationOptions: WaypointTimingOptimizationOptions(initialStepSeconds: 2, minimumStepSeconds: 0.01,
                    maximumIterations: 1), initialMass: mass(), backgroundDynamics: background(),
                legEngineLimits: [limits()])
        }
        let result = try run(10)
        #expect(result.status == .seedsExhausted); #expect(result.starts.count == 2)
        #expect(result.starts[0].status == .seedRejected)
        #expect(result.evaluatedCandidates == 3)
        #expect(try #require(result.route).arrivalEpoch == epoch(8))
        let bounded = try run(2)
        #expect(bounded.status == .evaluationLimitReached)
        #expect(bounded.evaluatedCandidates == 2)
        #expect(try #require(bounded.route).arrivalEpoch == epoch(10))
    }
    @Test func multipleSeedsExploreSeparatedObstacleTimingBasins() throws {
        let obstacle = try MovingSphericalExclusionRegion(identifier: "crossing",
            centreProvider: target(p: SIMD3(5,-5,0), v: SIMD3(0,1,0)), radiusMetres: 0.3,
            maximumCentreSpeedMetresPerSecond: 1)
        let result = try WaypointRouteTimingPlanner.optimizeTimingFromSeeds(departureEpoch: .j2000,
            frame: frame, origin: .sun,
            waypoints: [BrachistochroneWaypoint(positionMetres: .zero),
                BrachistochroneWaypoint(positionMetres: SIMD3(10,0,0))],
            legAccelerationsMetresPerSecondSquared: [1], arrivalWindows: [epoch(7)...epoch(20)],
            seedArrivalEpochs: [[epoch(20)],[epoch(8)]], maximumTotalCandidateEvaluations: 400,
            optimizationOptions: WaypointTimingOptimizationOptions(initialStepSeconds: 1, minimumStepSeconds: 0.02),
            initialMass: mass(), backgroundDynamics: background(), legEngineLimits: [limits()],
            movingObstacles: [obstacle])
        #expect(result.starts.count == 2)
        #expect(try #require(result.starts[0].route).arrivalEpoch > epoch(10))
        #expect(try #require(result.starts[1].route).arrivalEpoch < epoch(9))
        #expect(result.route?.arrivalEpoch == result.starts[1].route?.arrivalEpoch)
    }
    @Test func heatEnvelopeCarriesAcrossPassageAndCooldown() throws {
        let v = SIMD3<Double>(1,0,0)
        let route = try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: [BrachistochroneWaypoint(positionMetres: .zero, velocityMetresPerSecond: v),
                BrachistochroneWaypoint(positionMetres: SIMD3(3,0,0), velocityMetresPerSecond: v),
                BrachistochroneWaypoint(positionMetres: SIMD3(6,0,0), velocityMetresPerSecond: v)],
            legAccelerationsMetresPerSecondSquared: [1,1],
            legDurations: [CelMek.Duration(seconds: 3), CelMek.Duration(seconds: 3)])
        let burn = try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 1,
            maximumThrustNewtons: 2, effectiveExhaustVelocityMetresPerSecond: 100)
        let coast = try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 0,
            maximumThrustNewtons: 0, effectiveExhaustVelocityMetresPerSecond: 100)
        let phase = try WaypointRouteEnginePhase(startEpoch: epoch(3), endEpoch: epoch(6), limits: coast)
        func execute(cap: Double) throws -> WaypointRoutePropagationResult {
            try WaypointRoutePropagator(route: route, initialMass: mass(), backgroundDynamics: background(),
                legEngineLimits: [burn,burn], timedEnginePhases: [phase],
                physicalConstraints: WaypointRoutePhysicalConstraints(thermalEnvelope:
                    WaypointThermalEnvelope(maximumHeatJoules: cap, heatingWattsPerNewton: 1,
                        coolingRatePerSecond: log(2)/3))).propagation(at: epoch(6))
        }
        let result = try execute(cap: 5)
        let peak = try #require(result.segments[0].finalHeatUpperBoundJoules)
        #expect(abs(peak-3/log(2)) < 1e-11)
        #expect(result.segments[1].initialHeatUpperBoundJoules == peak)
        #expect(abs(try #require(result.finalHeatUpperBoundJoules)-peak/2) < 1e-11)
        #expect(result.segments[0].finalState.mass == result.segments[1].initialState.mass)
        #expect(result.state.mass == (try mass()))
        #expect(throws: WaypointRoutePropagationError.thermalBoundExceeded(leg: 0)) { try execute(cap: 4) }
        let linear = try WaypointThermalEnvelope(initialHeatJoules: 1, maximumHeatJoules: 10,
            heatingWattsPerNewton: 2, coolingRatePerSecond: 0)
        #expect(abs(try linear.advance(heat: 1, duration: 2, thrustCeiling: 1, leg: 0)-5) < 1e-12)
    }
    @Test func physicalConstraintsReachAdaptiveAndMultiStartSearches() throws {
        let points = try [BrachistochroneWaypoint(positionMetres: .zero),
            BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0))]
        let constraints = try WaypointRoutePhysicalConstraints(thermalEnvelope:
            WaypointThermalEnvelope(maximumHeatJoules: 1, heatingWattsPerNewton: 1, coolingRatePerSecond: 0))
        let adaptive = try WaypointRouteTimingPlanner.adaptiveSearch(departureEpoch: .j2000,
            frame: frame, origin: .sun, waypoints: points, legAccelerationsMetresPerSecondSquared: [1],
            arrivalWindows: [epoch(10)...epoch(10)],
            refinementOptions: WaypointTimingRefinementOptions(initialSpacingSeconds: 1, minimumSpacingSeconds: 1),
            initialMass: mass(), backgroundDynamics: background(), legEngineLimits: [limits()],
            physicalConstraints: constraints)
        #expect(adaptive.route == nil)
        #expect(adaptive.failedAttempts.first?.failure == .engine(.thermalBoundExceeded(leg: 0)))
        let multiple = try WaypointRouteTimingPlanner.optimizeTimingFromSeeds(departureEpoch: .j2000,
            frame: frame, origin: .sun, waypoints: points, legAccelerationsMetresPerSecondSquared: [1],
            arrivalWindows: [epoch(10)...epoch(10)], seedArrivalEpochs: [[epoch(10)]],
            optimizationOptions: WaypointTimingOptimizationOptions(initialStepSeconds: 1, minimumStepSeconds: 1),
            initialMass: mass(), backgroundDynamics: background(), legEngineLimits: [limits()],
            physicalConstraints: constraints)
        #expect(multiple.route == nil)
        #expect(multiple.starts[0].status == .seedRejected)
        #expect(multiple.starts[0].failedAttempts.first?.failure == .engine(.thermalBoundExceeded(leg: 0)))
    }
    @Test func attitudeConeChecksActualRequiredThrustAndSearchDiagnostics() throws {
        let points = try [BrachistochroneWaypoint(positionMetres: .zero),
            BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0))]
        let fixed = try WaypointThrustCone(frame: frame, halfAngleRadians: 0.1) { _ in SIMD3(1,0,0) }
        let failed = try WaypointRouteTimingPlanner.search(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points, legAccelerationsMetresPerSecondSquared: [1], candidateArrivalEpochs: [[epoch(10)]],
            initialMass: mass(), backgroundDynamics: background(), legEngineLimits: [limits()],
            physicalConstraints: WaypointRoutePhysicalConstraints(thrustCone: fixed))
        #expect(failed.status == .candidateSetExhausted)
        #expect(failed.failedAttempts.first?.failure == .engine(.attitudeLimitExceeded(leg: 0)))
        let rotating = try WaypointThrustCone(frame: frame, halfAngleRadians: 0) { epoch in
            try epoch.duration(since: .j2000).seconds <= 5 ? SIMD3(1,0,0) : SIMD3(-1,0,0)
        }
        let route = try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points, legAccelerationsMetresPerSecondSquared: [1],
            legDurations: [CelMek.Duration(seconds: 10)])
        let result = try WaypointRoutePropagator(route: route, initialMass: mass(), backgroundDynamics: background(),
            legEngineLimits: [limits()], physicalConstraints:
                WaypointRoutePhysicalConstraints(thrustCone: rotating)).propagation(at: epoch(10))
        #expect(result.checkpoints[0].positionErrorMetres < 1e-8)
    }
}

// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct EncounterSequenceTests {
    let frame = Frame(name: "Sequence central inertial",pseudoInertial: true)
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: seconds)) }
    func gm(_ value: Double = 1) throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: value,convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Independent analytic fixture",version: "1"))
    }
    func target(_ p: SIMD3<Double>, _ v: SIMD3<Double> = .zero, origin: BodyIdentifier = .sun) throws -> MovingWaypointRouteTarget {
        try MovingWaypointRouteTarget(provider: LambertTransferPlannerTests.Linear(frame: frame,p: p,v: v),origin: origin,coverage: epoch(0)...epoch(30))
    }
    func mass() throws -> SpacecraftMassProperties { try SpacecraftMassProperties(dryMassKilograms: 50,propellantMassKilograms: 50) }
    @Test func directTransferSearchAndNumericalConvergence() throws {
        let end = try epoch(.pi/2), context = try ForceEvaluationContext(frame: frame,origin: .sun)
        let initial = try SpacecraftState(cartesian: CartesianState(epoch: .j2000,frame: frame,position: SIMD3(1,0,0),velocity: SIMD3(0,0.5,0)),mass: mass())
        let arrival = try target(SIMD3(0,1,0))
        let visit = try EncounterSequenceVisit(identifier: "Mars passage",target: arrival,arrivalWindow: end...end)
        let branch = try LambertWindowBranch(direction: .positive)
        let good = EncounterSequenceCandidate(epochs: [end],branches: [branch])
        let bad = try EncounterSequenceCandidate(epochs: [epoch(1)],branches: [branch])
        let search = try EncounterSequencePlanner.search(initialState: initial,visits: [visit],candidates: [bad,good],
            context: context,gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),effectiveExhaustVelocityMetresPerSecond: 100)
        #expect(search.status == .candidatesExhausted)
        #expect(search.rows[0].failure?.contains("outsideWindow") == true)
        let seed = try #require(search.rankedSeeds.first)
        #expect(abs(seed.totalDeltaVelocityMetresPerSecond-0.5) < 1e-10)
        #expect(abs(seed.propellantUsedKilograms-100 * -expm1(-0.005)) < 1e-10)
        let result = try NumericalEncounterSequence.refine(initialState: initial,arrivalEpoch: end,nodes: seed.maneuvers,
            requirements: seed.requirements,context: context,backgroundModels: [CentralGravity(body: .sun,gravitationalParameter: gm())],maximumEvaluations: 10,numericalOptions: .preset(for: .full),
            eventOptions: PropagationEventOptions(maximumCheckIntervalSeconds: 0.1))
        #expect(result.status == .converged)
        #expect(result.toleranceReport?.consistent == true)
        #expect(result.best.encounters[0].relativeSpeedMetresPerSecond > 0.99)
        #expect(abs(result.best.arrivalState.mass.totalMassKilograms-seed.finalMass.totalMassKilograms) < 1e-10)
        let limited = try EncounterSequencePlanner.search(initialState: initial,visits: [visit],candidates: [good,bad],
            context: context,gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),effectiveExhaustVelocityMetresPerSecond: 100,maximumEvaluations: 1)
        #expect(limited.status == .evaluationLimitReached && limited.rows.count == 1)
    }
    @Test func finiteLocalAssistWithIndependentInteriorPassage() throws {
        let tree = try FrameTree(root: frame)
        let local = try tree.makeFrame(parent: frame,name: "Encounter-centred",pseudoInertial: true,provider: FixedTransformProvider())
        let peri = try epoch(10), context = try ForceEvaluationContext(frame: frame,origin: .earth)
        let request = try FlybyEncounterRequest(body: .earth,periapsisEpoch: peri,frame: local,
            gravitationalParameter: gm(),bodyRadiusMetres: 0.1,incomingExcessVelocityMetresPerSecond: SIMD3(1,0,0),outgoingExcessVelocityMetresPerSecond: SIMD3(0,1,0))
        let flyby = try UnpoweredFlyby.solve(request)
        func physical(_ epoch: Epoch) throws -> CartesianState {
            try tree.transform(from: local,to: frame,at: epoch).applying(to: flyby.trajectory.state(at: epoch),destination: frame)
        }
        let initial = SpacecraftState(cartesian: try physical(.j2000),mass: try mass())
        let body = try target(.zero,origin: .earth), end = try epoch(12), interiorEpoch = try epoch(10.5)
        let final = try target(physical(end).position,origin: .earth)
        let visits = try [EncounterSequenceVisit(identifier: "Earth assist",target: body,arrivalWindow: epoch(9.8)...epoch(10.2),
            condition: .unpoweredFlyby(request: request,incomingSeconds: 1,outgoingSeconds: 1,tree: tree)),
            EncounterSequenceVisit(identifier: "Departure passage",target: final,arrivalWindow: end...end)]
        let interior = try NumericalEncounterRequirement(identifier: "Independent physical waypoint",target: target(physical(interiorEpoch).position,origin: .earth),
            epoch: interiorEpoch,window: epoch(10.4)...epoch(10.6),positionToleranceMetres: 1e-4)
        let branch = try LambertWindowBranch(direction: .positive)
        let seed = try EncounterSequencePlanner.seed(initialState: initial,visits: visits,
            candidate: EncounterSequenceCandidate(epochs: [peri,end],branches: [branch,branch]),context: context,
            gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),effectiveExhaustVelocityMetresPerSecond: 100,additionalPassages: [interior])
        let unsafe = try FlybyEncounterRequest(body: .earth,periapsisEpoch: peri,frame: local,gravitationalParameter: gm(),bodyRadiusMetres: 0.5,
            incomingExcessVelocityMetresPerSecond: SIMD3(1,0,0),outgoingExcessVelocityMetresPerSecond: SIMD3(0,1,0))
        let rejectedVisits = try [EncounterSequenceVisit(identifier: "Unsafe assist",target: body,arrivalWindow: visits[0].arrivalWindow,
            condition: .unpoweredFlyby(request: unsafe,incomingSeconds: 1,outgoingSeconds: 1,tree: tree)),visits[1]]
        let rejected = try EncounterSequencePlanner.search(initialState: initial,visits: rejectedVisits,candidates: [seed.candidate],context: context,
            gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),effectiveExhaustVelocityMetresPerSecond: 100)
        #expect(rejected.rows[0].failure?.contains("unsafePeriapsis") == true)
        #expect(rejected.rankedSeeds.isEmpty)
        #expect(seed.arcs.count == 3)
        #expect(seed.requirements.count == 3)
        #expect(seed.arcs[1].coverage.contains(interiorEpoch))
        let result = try NumericalEncounterSequence.refine(initialState: initial,arrivalEpoch: end,nodes: seed.maneuvers,
            requirements: seed.requirements,context: context,backgroundModels: [CentralGravity(body: .earth,gravitationalParameter: gm())],
            maximumEvaluations: 10,numericalOptions: .preset(for: .full),
            eventOptions: PropagationEventOptions(maximumCheckIntervalSeconds: 0.05,timeToleranceSeconds: 1e-7))
        #expect(result.status == .converged)
        let encounter = result.best.encounters[0]
        #expect(abs(try encounter.closestApproachEpoch.duration(since: peri).seconds) < 1e-5)
        #expect(abs(encounter.closestApproachDistanceMetres-(sqrt(2)-1)) < 1e-5)
        #expect(encounter.relativeSpeedMetresPerSecond > 2)
        #expect(result.best.encounters[1].relativeSpeedMetresPerSecond > 1)
        #expect(abs(result.best.arrivalState.mass.totalMassKilograms-initial.mass.totalMassKilograms) < 1e-7)
        for (a,b) in zip(seed.arcs,seed.arcs.dropFirst()) {
            #expect(routeNorm(try a.state(at: a.coverage.upperBound).position-b.state(at: b.coverage.lowerBound).position) < 1e-10)
        }
    }
    @Test func refinementBurnsAndIndependentEpochs() throws {
        let context = try ForceEvaluationContext(frame: frame,origin: .sun), end = try epoch(10)
        let initial = try SpacecraftState(cartesian: CartesianState(epoch: .j2000,frame: frame,position: SIMD3(1,0,0),velocity: SIMD3(0,1,0)),mass: mass())
        let impulse = try ImpulsiveManeuver(epoch: epoch(2),frame: frame,deltaVelocityMetresPerSecond: SIMD3(0.8,0,0))
        let node = try DeepSpaceManeuverNode(window: epoch(1)...epoch(3),command: .impulse(impulse,effectiveExhaustVelocityMetresPerSecond: 100))
        let requirement = try NumericalEncounterRequirement(identifier: "End",target: target(SIMD3(9,10,0)),epoch: end,window: end...end,positionToleranceMetres: 1e-6)
        let result = try NumericalEncounterSequence.refine(initialState: initial,arrivalEpoch: end,nodes: [node],requirements: [requirement],
            context: context,backgroundModels: [],maximumEvaluations: 100,initialTimeStepSeconds: 0.1,initialDeltaVelocityStepMetresPerSecond: 0.2)
        #expect(result.status == .converged)
        #expect(result.best.normalizedResidual < result.initial.normalizedResidual)
        let exhausted = try NumericalEncounterSequence.refine(initialState: initial,arrivalEpoch: end,nodes: [node],requirements: [requirement],
            context: context,backgroundModels: [],maximumEvaluations: 1)
        #expect(exhausted.status == .evaluationLimitReached)
        let burn = try FiniteBurn(start: epoch(2),end: epoch(4),thrust: FiniteThrust(frame: frame,thrustNewtons: 10,
            effectiveExhaustVelocityMetresPerSecond: 100,direction: SIMD3(1,0,0)))
        let exactNode = try DeepSpaceManeuverNode(window: epoch(1)...epoch(5),command: .finiteBurn(burn))
        let exactTrajectory = try NumericalEncounterTrajectory(initialState: initial,arrivalEpoch: end,nodes: [exactNode],context: context,backgroundModels: [])
        let exact = try exactTrajectory.state(at: end)
        let endpoint = try NumericalEncounterRequirement(identifier: "Finite endpoint",target: target(exact.position),epoch: end,window: end...end,positionToleranceMetres: 1e-5)
        let shorter = try FiniteBurn(start: epoch(2),end: epoch(3.6),thrust: burn.thrust)
        let finite = try NumericalEncounterSequence.refine(initialState: initial,arrivalEpoch: end,
            nodes: [DeepSpaceManeuverNode(window: exactNode.window,command: .finiteBurn(shorter))],requirements: [endpoint],context: context,
            backgroundModels: [],maximumEvaluations: 200,initialTimeStepSeconds: 0.1)
        #expect(finite.best.normalizedResidual < finite.initial.normalizedResidual)
        #expect(finite.best.nodes.allSatisfy { $0.window.contains(sequenceInterval($0).lowerBound) && $0.window.contains(sequenceInterval($0).upperBound) })
        #expect(abs(try exactTrajectory.spacecraftState(at: epoch(3)).mass.totalMassKilograms-99.9) < 1e-8)
    }
    @Test func independentPassageEpochAndPhysicalConstraints() throws {
        let context = try ForceEvaluationContext(frame: frame,origin: .sun), end = try epoch(3)
        let initial = try SpacecraftState(cartesian: CartesianState(epoch: .j2000,frame: frame,position: SIMD3(1,0,0),velocity: SIMD3(1,0,0)),mass: mass())
        let passage = try NumericalEncounterRequirement(identifier: "Interior",target: target(SIMD3(3,0,0)),epoch: epoch(1),window: epoch(0.5)...epoch(2.5),positionToleranceMetres: 1e-6)
        let timed = try NumericalEncounterSequence.refine(initialState: initial,arrivalEpoch: end,nodes: [],requirements: [passage],context: context,
            backgroundModels: [],maximumEvaluations: 20,initialTimeStepSeconds: 1,eventOptions: PropagationEventOptions(maximumCheckIntervalSeconds: 0.1))
        #expect(timed.status == .converged)
        #expect(timed.best.requirements[0].epoch == (try epoch(2)))
        let engine = try FiniteThrust(frame: frame,thrustNewtons: 10,effectiveExhaustVelocityMetresPerSecond: 100,direction: SIMD3(1,0,0))
        let node = try DeepSpaceManeuverNode(window: epoch(0.5)...epoch(2),command: .finiteBurn(FiniteBurn(start: epoch(1),end: epoch(2),thrust: engine)))
        let endpoint = try NumericalEncounterRequirement(identifier: "Endpoint",target: target(SIMD3(4,0,0)),epoch: end,window: end...end)
        let cone = try WaypointThrustCone(frame: frame,halfAngleRadians: 0.1) { _ in SIMD3(1,0,0) }
        let thermal = try WaypointThermalEnvelope(maximumHeatJoules: 20,heatingWattsPerNewton: 1,coolingRatePerSecond: 0)
        let evaluated = try NumericalEncounterSequence.evaluate(initialState: initial,arrivalEpoch: end,nodes: [node],requirements: [endpoint],context: context,
            backgroundModels: [],constraints: NumericalSequenceConstraints(physical: WaypointRoutePhysicalConstraints(thrustCone: cone,thermalEnvelope: thermal),maximumThrustNewtons: 10))
        #expect(abs(evaluated.finalHeatJoules!-10) < 1e-10)
        let wrongCone = try WaypointThrustCone(frame: frame,halfAngleRadians: 0.1) { _ in SIMD3(0,1,0) }
        #expect(throws: WaypointRoutePropagationError.attitudeLimitExceeded(leg: 0)) {
            try NumericalEncounterSequence.evaluate(initialState: initial,arrivalEpoch: end,nodes: [node],requirements: [endpoint],context: context,
                backgroundModels: [],constraints: NumericalSequenceConstraints(physical: WaypointRoutePhysicalConstraints(thrustCone: wrongCone)))
        }
        #expect(throws: NumericalSequenceError.constraintViolation("engine ceiling")) {
            try NumericalEncounterSequence.evaluate(initialState: initial,arrivalEpoch: end,nodes: [node],requirements: [endpoint],context: context,
                backgroundModels: [],constraints: NumericalSequenceConstraints(maximumThrustNewtons: 5))
        }
        let hot = try WaypointThermalEnvelope(maximumHeatJoules: 5,heatingWattsPerNewton: 1,coolingRatePerSecond: 0)
        #expect(throws: WaypointRoutePropagationError.thermalBoundExceeded(leg: 1)) {
            try NumericalEncounterSequence.evaluate(initialState: initial,arrivalEpoch: end,nodes: [node],requirements: [endpoint],context: context,
                backgroundModels: [],constraints: NumericalSequenceConstraints(physical: WaypointRoutePhysicalConstraints(thermalEnvelope: hot)))
        }
        #expect(throws: NumericalSequenceError.insufficientReserve) {
            try NumericalEncounterSequence.evaluate(initialState: initial,arrivalEpoch: end,nodes: [node],requirements: [endpoint],context: context,
                backgroundModels: [],constraints: NumericalSequenceConstraints(minimumRemainingPropellantKilograms: 50))
        }
        let bounded = try NumericalEncounterSequence.refine(initialState: initial,arrivalEpoch: end,nodes: [],requirements: [passage],context: context,
            backgroundModels: [],maximumEvaluations: 10,maximumProposalAttempts: 1,initialTimeStepSeconds: 100)
        #expect(bounded.status == .proposalLimitReached)
    }
    @Test func interplanetaryScaleSeedUnderMultiBodyGravity() throws {
        let au = 149_597_870_700.0, mu = 1.32712440018e20, days = 200*86400.0
        let end = try epoch(days), coverage = Epoch.j2000...end
        let earthSpeed = sqrt(mu/au), marsRadius = 1.524*au
        let marsPosition = SIMD3<Double>(0,marsRadius,0), marsVelocity = SIMD3<Double>(-sqrt(mu/marsRadius),0,0)
        let departure = try SpacecraftState(cartesian: CartesianState(epoch: .j2000,frame: frame,position: SIMD3(au,0,0),velocity: SIMD3(0,earthSpeed,0)),mass: mass())
        let mars = try MovingWaypointRouteTarget(provider: LambertTransferPlannerTests.Linear(frame: frame,p: marsPosition-marsVelocity*days,v: marsVelocity),origin: .sun,coverage: coverage)
        let moonMetadata = try EphemerisMetadata(body: .jupiter,origin: .sun,frame: frame,coverage: coverage,
            coverageBasis: .callerSelected,modelTimeScale: .tt,velocityMethod: .analytic,accelerationMethod: .unavailable,
            provenance: PhysicalDataProvenance(source: "Synthetic distant perturber",version: "1"))
        let jupiter = try StateProviderEphemeris(provider: LambertTransferPlannerTests.Linear(frame: frame,p: SIMD3(1000*au,0,0),v: .zero),metadata: moonMetadata)
        let context = try ForceEvaluationContext(frame: frame,origin: .sun,ephemerides: EphemerisCollection(providers: [jupiter]))
        let visit = try EncounterSequenceVisit(identifier: "Mars",target: mars,arrivalWindow: end...end)
        let seed = try EncounterSequencePlanner.seed(initialState: departure,visits: [visit],candidate: EncounterSequenceCandidate(epochs: [end],branches: [LambertWindowBranch(direction: .positive)]),
            context: context,gravitationalParameter: gm(mu),referencePole: SIMD3(0,0,1),effectiveExhaustVelocityMetresPerSecond: 50_000,positionToleranceMetres: 100)
        let options = try NumericalPropagationOptions(absolutePositionToleranceMetres: 0.001,absoluteVelocityToleranceMetresPerSecond: 1e-8,
            relativeTolerance: 1e-13,initialStepSeconds: 10,maximumStepSeconds: 3600)
        let requirement = try NumericalEncounterRequirement(identifier: "Mars",target: mars,epoch: end,window: end...end,positionToleranceMetres: 100)
        let result = try NumericalEncounterSequence.refine(initialState: departure,arrivalEpoch: end,nodes: seed.maneuvers,requirements: [requirement],context: context,
            backgroundModels: [CentralGravity(body: .sun,gravitationalParameter: gm(mu)),ThirdBodyGravity(body: .jupiter,centralBody: .sun,gravitationalParameter: gm(1e14))],
            maximumEvaluations: 5,numericalOptions: options,convergencePositionMetres: 1,convergenceVelocityMetresPerSecond: 1e-5)
        #expect(result.status == .converged)
        #expect(result.best.encounters[0].positionErrorMetres < 100)
        #expect(result.best.encounters[0].relativeSpeedMetresPerSecond > 1000)
        #expect(result.best.propellantUsedKilograms > 0)
    }
    @Test func seedDsmCandidatesFiniteConversionAndVerifiedRanking() throws {
        let end = try epoch(.pi/2), context = try ForceEvaluationContext(frame: frame,origin: .sun)
        let initial = try SpacecraftState(cartesian: CartesianState(epoch: .j2000,frame: frame,position: SIMD3(1,0,0),velocity: SIMD3(0,0.5,0)),mass: mass())
        let visit = try EncounterSequenceVisit(identifier: "Arrival",target: target(SIMD3(0,1,0)),arrivalWindow: end...end)
        let branch = try LambertWindowBranch(direction: .positive)
        let dsm = try DeepSpaceManeuverNode(window: epoch(0.4)...epoch(0.6),command: .impulse(
            ImpulsiveManeuver(epoch: epoch(0.5),frame: frame,deltaVelocityMetresPerSecond: SIMD3(0.1,0,0)),effectiveExhaustVelocityMetresPerSecond: 100))
        let grid = try EncounterSequencePlanner.candidates(epochOptions: [[end,epoch(1)]],branchOptions: [[branch]],deepSpaceOptions: [[],[dsm]])
        #expect(grid.count == 4)
        #expect(throws: EncounterSequenceError.invalidInput) {
            try EncounterSequencePlanner.candidates(epochOptions: [[end,epoch(1)]],branchOptions: [[branch]],deepSpaceOptions: [[],[dsm]],maximumCandidates: 3)
        }
        let candidates = [EncounterSequenceCandidate(epochs: [end],branches: [branch],deepSpaceManeuvers: [dsm]),
            EncounterSequenceCandidate(epochs: [end],branches: [branch])]
        let searched = try EncounterSequencePlanner.search(initialState: initial,visits: [visit],candidates: candidates,context: context,
            gravitationalParameter: gm(),referencePole: SIMD3(0,0,1),effectiveExhaustVelocityMetresPerSecond: 100,objective: .propellant)
        #expect(searched.rankedSeeds.count == 2)
        let seed = try #require(searched.rankedSeeds.first)
        #expect(seed.maneuvers.count == 1)
        #expect(searched.rankedSeeds[1].propellantUsedKilograms > seed.propellantUsedKilograms)
        let finite = try NumericalEncounterSequence.finiteSeed(seed,burnWindows: [.j2000...epoch(0.4)],thrustNewtons: 1000)
        let propagated = try NumericalEncounterTrajectory(initialState: initial,arrivalEpoch: end,nodes: finite,context: context,
            backgroundModels: [CentralGravity(body: .sun,gravitationalParameter: gm())]).spacecraftState(at: end)
        #expect(abs(propagated.mass.totalMassKilograms-seed.finalMass.totalMassKilograms) < 1e-8)
        let verified = try NumericalEncounterSequence.search(seeds: [seed,seed],context: context,
            backgroundModels: [CentralGravity(body: .sun,gravitationalParameter: gm())],maximumCandidates: 1,
            maximumEvaluationsPerCandidate: 5,numericalOptions: .preset(for: .full))
        #expect(verified.status == .evaluationLimitReached)
        #expect(verified.rows[0].status == .converged)
        #expect(verified.rankedConvergedRows.count == 1)
        let reserve = try NumericalEncounterSequence.search(seeds: [seed],context: context,
            backgroundModels: [CentralGravity(body: .sun,gravitationalParameter: gm())],
            constraints: NumericalSequenceConstraints(minimumRemainingPropellantKilograms: 50))
        #expect(reserve.rows[0].status == .infeasible)
        #expect(reserve.rankedConvergedRows.isEmpty)
        let budget = try NumericalEncounterSequence.search(seeds: [seed],context: context,backgroundModels: [],maximumStateQueriesPerEvaluation: 1)
        #expect(budget.rows[0].status == .workLimitReached)
    }
    @Test func clearanceResolvedUnresolvedAndRejected() throws {
        let initial = try SpacecraftState(cartesian: CartesianState(epoch: .j2000,frame: frame,position: SIMD3(0,2,0),velocity: SIMD3(1,0,0)),mass: mass())
        let context = try ForceEvaluationContext(frame: frame,origin: .sun), end = try epoch(2)
        let requirement = try NumericalEncounterRequirement(identifier: "Pass",target: target(SIMD3(2,2,0)),epoch: end,window: end...end)
        let centre = try target(.zero)
        func evaluate(_ radius: Double, _ bound: Double?) throws -> NumericalSequenceEvaluation {
            try NumericalEncounterSequence.evaluate(initialState: initial,arrivalEpoch: end,nodes: [],requirements: [requirement],
                context: context,backgroundModels: [],constraints: NumericalSequenceConstraints(clearances:
                    [NumericalSequenceClearance(identifier: "Sun clearance",centre: centre,minimumDistanceMetres: radius,relativeSpeedBoundMetresPerSecond: bound)]),
                eventOptions: PropagationEventOptions(maximumCheckIntervalSeconds: 0.1))
        }
        #expect(try evaluate(1,1).clearance[0].status == .satisfied)
        #expect(try evaluate(1,nil).clearance[0].status == .unresolved)
        #expect(try evaluate(3,1).clearance[0].status == .violated)
        #expect(throws: NumericalSequenceError.queryLimitReached) {
            try NumericalEncounterSequence.evaluate(initialState: initial,arrivalEpoch: end,nodes: [],requirements: [requirement],
                context: context,backgroundModels: [],maximumStateQueries: 1)
        }
        let unresolved = try NumericalEncounterSequence.refine(initialState: initial,arrivalEpoch: end,nodes: [],requirements: [requirement],context: context,
            backgroundModels: [],constraints: NumericalSequenceConstraints(clearances: [NumericalSequenceClearance(identifier: "Unbounded",centre: centre,minimumDistanceMetres: 1)]))
        #expect(unresolved.status == .unresolvedConstraints)
    }
}

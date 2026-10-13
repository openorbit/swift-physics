// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct EncounterPlanningTests {
    let fixture = UnpoweredFlybyTests()
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: seconds)) }
    @Test func bPlaneSignsAndRoundTrip() throws {
        let flyby = try UnpoweredFlyby.solve(fixture.request())
        let b = try FlybyBPlane(flyby: flyby,referencePole: SIMD3(0,0,1))
        // p=mu=vInf=1 implies impact parameter 1, B=-Y.
        #expect(routeNorm(b.s-SIMD3(1,0,0)) < 1e-14)
        #expect(routeNorm(b.t-SIMD3(0,-1,0)) < 1e-14)
        #expect(abs(b.bTMetres-1) < 1e-14)
        #expect(abs(b.bRMetres) < 1e-14)
        let recovered = try b.solve(fixture.request())
        #expect(routeNorm(recovered.periapsisState.position-flyby.periapsisState.position) < 1e-13)
        #expect(routeNorm(recovered.outgoingExcessVelocityMetresPerSecond-SIMD3(0,1,0)) < 1e-13)
        #expect(throws: EncounterPlanningError.degenerateAxes) {
            try FlybyBPlane(incomingVelocity: SIMD3(1,0,0),referencePole: SIMD3(1,0,0),bTMetres: 1,bRMetres: 0)
        }
        #expect(throws: FlybyError.unsafePeriapsis) {
            try FlybyBPlane(incomingVelocity: SIMD3(1,0,0),referencePole: SIMD3(0,0,1),bTMetres: 0.01,bRMetres: 0).solve(fixture.request())
        }
    }
    @Test func poweredUnequalSpeedsFuelAndAsymptotes() throws {
        // rp=mu=1, incoming e=2, outgoing e=5; common radial +X.
        let incoming = SIMD3<Double>(0.5,sqrt(0.75),0)
        let outgoing = SIMD3<Double>(-0.4,2*sqrt(0.96),0)
        let request = try fixture.request(incoming: incoming,outgoing: outgoing)
        let mass = try SpacecraftMassProperties(dryMassKilograms: 50,propellantMassKilograms: 50)
        let result = try PoweredFlyby.solve(request,periapsisRadiusMetres: 1,radialDirection: SIMD3(1,0,0),
            mass: mass,effectiveExhaustVelocityMetresPerSecond: 100)
        let dv = sqrt(6.0)-sqrt(3.0)
        #expect(abs(result.maneuver.deltaVelocityMetresPerSecond.y-dv) < 1e-13)
        #expect(abs(result.propellantUsedKilograms-100 * -expm1(-dv/100)) < 1e-12)
        #expect(result.preBurnState.cartesian.position == result.postBurnState.cartesian.position)
        for (trajectory,seconds,expected) in [(result.incomingTrajectory,-1e6,incoming),(result.outgoingTrajectory,1e6,outgoing)] {
            let state = try trajectory.state(at: epoch(seconds))
            #expect(routeNorm(state.velocity-expected) < 1e-5)
        }
        #expect(throws: EncounterPlanningError.invalidGeometry) {
            try PoweredFlyby.solve(request,periapsisRadiusMetres: 1,radialDirection: SIMD3(0,1,0),mass: mass,effectiveExhaustVelocityMetresPerSecond: 100)
        }
        #expect(throws: EncounterPlanningError.insufficientPropellant) {
            try PoweredFlyby.solve(request,periapsisRadiusMetres: 1,radialDirection: SIMD3(1,0,0),
                mass: SpacecraftMassProperties(dryMassKilograms: 50),effectiveExhaustVelocityMetresPerSecond: 100)
        }
    }
    @Test func poweredPlaneChangeAndUnpoweredLimit() throws {
        let mass = try SpacecraftMassProperties(dryMassKilograms: 50,propellantMassKilograms: 50)
        let flyby = try UnpoweredFlyby.solve(fixture.request())
        let zero = try PoweredFlyby.solve(flyby.request,periapsisRadiusMetres: flyby.periapsisRadiusMetres,
            radialDirection: flyby.periapsisRadialDirection,mass: mass,effectiveExhaustVelocityMetresPerSecond: 100)
        #expect(routeNorm(zero.maneuver.deltaVelocityMetresPerSecond) < 1e-13)
        let incoming = SIMD3<Double>(0.5,sqrt(0.75),0), outgoing = SIMD3<Double>(-0.5,0,sqrt(0.75))
        let turn = try PoweredFlyby.solve(fixture.request(incoming: incoming,outgoing: outgoing),
            periapsisRadiusMetres: 1,radialDirection: SIMD3(1,0,0),mass: mass,effectiveExhaustVelocityMetresPerSecond: 100)
        #expect(routeNorm(turn.maneuver.deltaVelocityMetresPerSecond-SIMD3(0,-sqrt(3),sqrt(3))) < 1e-13)
        #expect(throws: FlybyError.unsafePeriapsis) {
            try PoweredFlyby.solve(fixture.request(),periapsisRadiusMetres: 0.1,radialDirection: SIMD3(1,0,0),
                mass: mass,effectiveExhaustVelocityMetresPerSecond: 100)
        }
    }
    @Test func mixedNodesCarryVelocityAndMass() throws {
        let f = fixture.frame, context = try ForceEvaluationContext(frame: f,origin: .sun)
        let seed = try SpacecraftState(cartesian: CartesianState(epoch: .j2000,frame: f,position: SIMD3(1,0,0),velocity: SIMD3(0,1,0)),
            mass: SpacecraftMassProperties(dryMassKilograms: 50,propellantMassKilograms: 50))
        let impulse = try ImpulsiveManeuver(epoch: epoch(2),frame: f,deltaVelocityMetresPerSecond: SIMD3(1,0,0))
        let burn = try FiniteBurn(start: epoch(4),end: epoch(6),thrust: FiniteThrust(frame: f,thrustNewtons: 10,
            effectiveExhaustVelocityMetresPerSecond: 100,direction: SIMD3(1,0,0)))
        let nodes = try [DeepSpaceManeuverNode(window: epoch(1)...epoch(3),command: .impulse(impulse,effectiveExhaustVelocityMetresPerSecond: 100)),
            DeepSpaceManeuverNode(window: epoch(4)...epoch(6),command: .finiteBurn(burn))]
        let finiteSeed = try nodes[0].finiteBurnSeed(initialState: seed,thrustNewtons: 100)
        if case .finiteBurn(let converted) = finiteSeed.command {
            #expect(finiteSeed.window.contains(converted.start))
            #expect(finiteSeed.window.contains(converted.end))
            #expect(abs(try converted.end.duration(since: converted.start).seconds-100 * -expm1(-0.01)) < 1e-12)
        } else { Issue.record("Expected finite seed") }
        let result = try DeepSpaceManeuverNode.execute(nodes.reversed(),initialState: seed,arrivalEpoch: epoch(10),context: context)
        let m = 100*exp(-0.01), expectedVX = 1+100*log(m/(m-0.2))
        #expect(abs(result.mass.totalMassKilograms-(m-0.2)) < 1e-10)
        #expect(abs(result.cartesian.velocity.x-expectedVX) < 1e-8)
        #expect(abs(result.cartesian.position.y-10) < 1e-8)
        #expect(throws: EncounterPlanningError.invalidWindow) {
            try DeepSpaceManeuverNode(window: epoch(5)...epoch(7),command: .finiteBurn(burn))
        }
    }
    @Test func refinementRespectsEveryWindow() throws {
        let f = FiniteTransferPlannerTests(), context = try ForceEvaluationContext(frame: f.frame,origin: .sun)
        let burn = try FiniteBurn(start: f.epoch(2),end: f.epoch(4),thrust: FiniteThrust(frame: f.frame,thrustNewtons: 10,
            effectiveExhaustVelocityMetresPerSecond: 100,direction: SIMD3(1,0,0)))
        let schedule = try ScheduledThrust(frame: f.frame,burns: [burn])
        let target = try CartesianState(epoch: f.epoch(10),frame: f.frame,position: SIMD3(20,10,0),velocity: SIMD3(2,1,0))
        let window = try f.epoch(2)...f.epoch(4)
        let result = try DeepSpaceManeuverNode.refine([DeepSpaceManeuverNode(window: window,command: .finiteBurn(burn))],
            initialState: f.seed(),target: target,context: context,maximumEvaluations: 20)
        #expect(result.rejectedTrials > 0)
        #expect(result.best.schedule.burns.allSatisfy { window.contains($0.start) && window.contains($0.end) })
        #expect(throws: FiniteTransferPlanningError.burnsDoNotFit) {
            try FiniteTransferPlanner.refine(initialState: f.seed(),initialSchedule: schedule,target: target,
                context: context,burnWindows: [f.epoch(3)...f.epoch(4)])
        }
    }
    @Test func patchTranslationRotationRoundTrip() throws {
        let tree = FrameTree(), root = tree.root
        let translated = try tree.makeFrame(parent: root,name: "Body origin",pseudoInertial: true,
            provider: FixedTransformProvider(translation: SIMD3(100,200,300)))
        let rotating = try tree.makeFrame(parent: translated,name: "Rotating body",provider:
            UniformRotationProvider(axis: SIMD3(0,0,1),rate: 0.01))
        let source = try ForceEvaluationContext(frame: root,origin: .sun)
        let destination = try ForceEvaluationContext(frame: rotating,origin: .earth)
        let state = try SpacecraftState(cartesian: CartesianState(epoch: epoch(10),frame: root,
            position: SIMD3(1000,2000,3000),velocity: SIMD3(10,20,30)),mass: SpacecraftMassProperties(dryMassKilograms: 50,propellantMassKilograms: 40))
        let local = try EncounterPatchBoundary(epoch: state.epoch,source: source,destination: destination).handover(state,using: tree)
        let restored = try EncounterPatchBoundary(epoch: state.epoch,source: destination,destination: source).handover(local,using: tree)
        #expect(routeNorm(restored.cartesian.position-state.cartesian.position) < 1e-10)
        #expect(routeNorm(restored.cartesian.velocity-state.cartesian.velocity) < 1e-12)
        #expect(restored.mass == state.mass)
        #expect(local.epoch == state.epoch)
        #expect(routeNorm(local.cartesian.velocity-state.cartesian.velocity) > 1)
    }
}

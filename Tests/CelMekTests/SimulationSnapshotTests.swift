// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import CelMek
import Math

@Suite struct SimulationSnapshotTests {
    private let frame = Frame(name: "save inertial",pseudoInertial: true)
    private func model() throws -> TwoBodyPropagator {
        let gm = try GravitationalParameter(cubicMetresPerSecondSquared: 3.986004418e14,convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "save test",version: "1",reference: "test reference"))
        let state = try CartesianState(epoch: .j2000,frame: frame,position: SIMD3(7e6,0,0),velocity: SIMD3(0,7500,0))
        return try TwoBodyPropagator(initialState: state,centralBody: .custom(namespace: "game",name: "Earth"),
            gravitationalParameter: gm,options: TwoBodyPropagationOptions(absoluteTimeToleranceSeconds: 1e-8,
            relativeTimeTolerance: 1e-12,maximumIterations: 100,maximumBracketExpansions: 90))
    }
    @Test func roundTripPreservesInputsStateAndClockContinuation() throws {
        let model = try model()
        var clock = try SimulationClock(epoch: .j2000,rate: 0.3)
        try clock.advance(realDuration: CelMek.Duration(seconds: 0.7))
        let mass = try SpacecraftMassProperties(dryMassKilograms: 100,propellantMassKilograms: 50)
        let state = try model.state(at: clock.epoch)
        let savedState = try SavedCartesianState(id: "ship",frameID: "earth-inertial",state: state,mass: mass)
        let savedInput = try SavedTwoBodyInput(id: "coast",frameID: "earth-inertial",propagator: model)
        let snapshot = try SimulationSnapshot(clock: clock,states: [savedState],twoBodyInputs: [savedInput],
            dependencies: [SimulationDependency(id: "earth-inertial",kind: "frame",version: "1",parameters: ["origin":"game:Earth"])])
        let data = try JSONEncoder().encode(snapshot)
        let decoded = try JSONDecoder().decode(SimulationSnapshot.self,from: data)
        let restored = try decoded.states[0].restore { id in
            #expect(id == "earth-inertial"); return frame
        }
        #expect(restored == state)
        #expect(restored.acceleration == state.acceleration)
        #expect(try decoded.states[0].restoreMass() == mass)
        let restoredModel = try decoded.twoBodyInputs[0].restore { _ in frame }
        #expect(restoredModel.options == model.options)
        #expect(restoredModel.gravitationalParameter == model.gravitationalParameter)
        #expect(restoredModel.centralBody == model.centralBody)
        let later = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: 900))
        #expect(try restoredModel.state(at: later) == model.state(at: later))
        var restoredClock = try #require(decoded.clock).restore()
        for _ in 0..<20 {
            let step = try CelMek.Duration(seconds: 0.11)
            try clock.advance(realDuration: step); try restoredClock.advance(realDuration: step)
            #expect(clock.epoch == restoredClock.epoch)
        }
        clock.setPaused(true)
        var paused = try clock.saved().restore()
        try paused.advance(realDuration: CelMek.Duration(seconds: 10))
        #expect(paused.epoch == clock.epoch && paused.isPaused)
    }
    @Test func preciseTimeAndMissingAcceleration() throws {
        let epoch = Epoch(offsetFromJ2000: try CelMek.Duration(wholeSeconds: 8_000_000_000_000,fractionalSecond: 0.125))
        let state = try CartesianState(epoch: epoch,frame: frame,position: .zero,velocity: .zero)
        let snapshot = try SimulationSnapshot(states: [SavedCartesianState(id: "point",frameID: "frame",state: state)])
        let data = try JSONEncoder().encode(snapshot)
        #expect(String(decoding: data,as: UTF8.self).contains("\"8000000000000\""))
        let decoded = try JSONDecoder().decode(SimulationSnapshot.self,from: data)
        #expect(try decoded.states[0].restore { _ in frame } == state)
        #expect(try decoded.states[0].restoreMass() == nil)
    }
    private func corrupted(_ change: (inout [String:Any]) -> Void) throws -> Data {
        let model = try model()
        let state = try SavedCartesianState(id: "ship",frameID: "frame",state: model.initialState)
        let snapshot = try SimulationSnapshot(clock: SimulationClock(epoch: .j2000,rate: 2),states: [state],
            twoBodyInputs: [SavedTwoBodyInput(id: "coast",frameID: "frame",propagator: model)])
        var json = try #require(JSONSerialization.jsonObject(with: JSONEncoder().encode(snapshot)) as? [String:Any])
        change(&json)
        return try JSONSerialization.data(withJSONObject: json)
    }
    @Test func rejectsUnknownVersionsDuplicateIDsAndMalformedPayloads() throws {
        let changes: [(inout [String:Any]) -> Void] = [
            { $0["schemaVersion"] = 2 },
            { $0["format"] = "other" },
            { json in let states = json["states"] as! [Any]; json["states"] = states+states },
            { json in var c = json["clock"] as! [String:Any]; c["rate"] = 3; c["epoch"] = ["offsetFromJ2000":["wholeSeconds":"1","fractionalSecond":0]]; json["clock"] = c },
            { json in var s = (json["states"] as! [[String:Any]])[0]; s["frameID"] = " "; json["states"] = [s] },
            { json in var s = (json["states"] as! [[String:Any]])[0]; s["epoch"] = ["offsetFromJ2000":["wholeSeconds":"00","fractionalSecond":0]]; json["states"] = [s] },
            { json in var s = (json["states"] as! [[String:Any]])[0]; s["dryMassKilograms"] = 100; json["states"] = [s] },
            { json in var s = (json["twoBodyInputs"] as! [[String:Any]])[0]; s["maximumIterations"] = 0; json["twoBodyInputs"] = [s] },
            { json in var s = (json["twoBodyInputs"] as! [[String:Any]])[0]; s["gmConvention"] = "unknown"; json["twoBodyInputs"] = [s] }
        ]
        for change in changes {
            let data = try corrupted(change)
            #expect(throws: (any Error).self) { try JSONDecoder().decode(SimulationSnapshot.self,from: data) }
        }
    }
    @Test func versionOneFixtureRemainsReadable() throws {
        let url = try #require(Bundle.module.url(forResource: "simulation-save-v1",withExtension: "json",subdirectory: "Fixtures"))
        let snapshot = try JSONDecoder().decode(SimulationSnapshot.self,from: Data(contentsOf: url))
        var clock = try #require(snapshot.clock).restore()
        #expect(clock.epoch.offsetFromJ2000 == (try CelMek.Duration(seconds: 2.5)))
        try clock.advance(realDuration: CelMek.Duration(seconds: 0.25))
        #expect(clock.epoch.offsetFromJ2000 == (try CelMek.Duration(seconds: 3)))
        let state = try snapshot.states[0].restore { _ in frame }
        #expect(state.epoch.offsetFromJ2000.wholeSeconds == 8_000_000_000_000)
        #expect(state.epoch.offsetFromJ2000.fractionalSecond == 0.125)
        #expect(state.acceleration == nil)
        #expect(try snapshot.states[0].restoreMass()?.totalMassKilograms == 125)
        let coast = try snapshot.twoBodyInputs[0].restore { _ in frame }
        #expect(coast.centralBody == .earth)
        #expect(coast.gravitationalParameter.cubicMetresPerSecondSquared == 3.986004418e14)
        #expect(coast.options.maximumIterations == 100)
        #expect(snapshot.dependencies[0].parameters["origin"] == "naif:399")
    }
    @Test func callerBindingErrorsAndPhysicalFrameValidationPropagate() throws {
        enum BindingError: Error { case missing }
        let input = try SavedTwoBodyInput(id: "coast",frameID: "frame",propagator: model())
        #expect(throws: BindingError.missing) { try input.restore { _ in throw BindingError.missing } }
        #expect(throws: PropagationError.nonInertialFrame) {
            try input.restore { _ in Frame(name: "rotating",pseudoInertial: false) }
        }
    }
}

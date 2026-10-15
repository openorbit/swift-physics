import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct ManeuverAndSamplingTests {
    private func norm(_ x: SIMD3<Double>) -> Double { hypot(hypot(x.x, x.y), x.z) }
    private func gm() throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: 3.986004418e14, convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Test", version: "1", reference: "game model"))
    }
    private func base(_ frame: Frame) throws -> TwoBodyPropagator {
        let initial = try CartesianState(epoch: .j2000, frame: frame,
            position: SIMD3(7e6,0,0), velocity: SIMD3(0,7546.053290107542,0))
        return try TwoBodyPropagator(initialState: initial, centralBody: .earth, gravitationalParameter: gm())
    }
    private func epoch(_ seconds: Double) throws -> Epoch { try .j2000.advanced(by: Duration(seconds: seconds)) }

    @Test func burnsAreSortedAndRightContinuous() throws {
        let tree = FrameTree(), model = try base(tree.root)
        let first = try ImpulsiveManeuver(epoch: epoch(100), frame: tree.root, deltaVelocityMetresPerSecond: SIMD3(1,100,3))
        let second = try ImpulsiveManeuver(epoch: epoch(500), frame: tree.root, deltaVelocityMetresPerSecond: SIMD3(-30,5,0))
        let trajectory = try ManeuveredTwoBodyPropagator(base: model, maneuvers: [second,first])
        #expect(trajectory.resolvedManeuvers.map(\.maneuver.epoch) == [first.epoch, second.epoch])
        for record in trajectory.resolvedManeuvers {
            let before = try trajectory.state(at: record.maneuver.epoch, side: .before)
            let after = try trajectory.state(at: record.maneuver.epoch)
            #expect(before.position == after.position)
            #expect(norm(after.velocity-before.velocity-record.deltaVelocityInPropagationFrame) < 1e-12)
            #expect(before == record.preBurnState)
            #expect(after == record.postBurnState)
        }
        let earlier = try epoch(-1000)
        #expect(try trajectory.state(at: earlier) == model.state(at: earlier))
        let dates = try [1000.0, 100, -100, 500, 1000].map { try epoch($0) }
        let states = try dates.map { try trajectory.state(at: $0) }
        for (date, state) in zip(dates.reversed(), states.reversed()) { #expect(try trajectory.state(at: date) == state) }
        let empty = try ManeuveredTwoBodyPropagator(base: model, maneuvers: [])
        #expect(try empty.state(at: epoch(1000)) == model.state(at: epoch(1000)))
    }

    @Test func rotatingTranslatedBurnFrameUsesOnlyRotation() throws {
        let tree = FrameTree(), model = try base(tree.root)
        struct MovingAxes: TransformProvider {
            func transform(at epoch: Epoch) throws -> FrameTransform {
                try FrameTransform(epoch: epoch, rotation: Quaternion(angle: -.pi/2, axis: SIMD3(0,0,1)),
                    translation: SIMD3(1e9,2e9,3e9), velocity: SIMD3(100,200,300),
                    acceleration: SIMD3(1,2,3), angularVelocity: SIMD3(0,0,-0.1))
            }
        }
        let axes = try tree.makeFrame(parent: tree.root, name: "Rotating and translated", provider: MovingAxes())
        let burn = try ImpulsiveManeuver(epoch: epoch(100), frame: axes, deltaVelocityMetresPerSecond: SIMD3(10,0,0))
        let trajectory = try ManeuveredTwoBodyPropagator(base: model, maneuvers: [burn], using: tree)
        let record = trajectory.resolvedManeuvers[0]
        #expect(norm(record.deltaVelocityInPropagationFrame-SIMD3(0.0,10,0)) < 1e-12)
        #expect(record.preBurnState.position == record.postBurnState.position)
        #expect(norm(record.postBurnState.velocity-record.preBurnState.velocity-SIMD3(0.0,10,0)) < 1e-12)
        #expect(throws: ManeuverError.frameTreeRequired) { try ManeuveredTwoBodyPropagator(base: model, maneuvers: [burn]) }
    }

    @Test func maneuverValidationAndReferenceEpochBurn() throws {
        let tree = FrameTree(), model = try base(tree.root)
        let burn = try ImpulsiveManeuver(epoch: .j2000, frame: tree.root, deltaVelocityMetresPerSecond: SIMD3(0,10,0))
        let trajectory = try ManeuveredTwoBodyPropagator(base: model, maneuvers: [burn])
        #expect(try trajectory.state(at: .j2000, side: .before) == model.state(at: .j2000))
        #expect(try trajectory.state(at: .j2000).velocity == model.initialState.velocity + SIMD3(0.0,10,0))
        #expect(throws: ManeuverError.duplicateEpoch) { try ManeuveredTwoBodyPropagator(base: model, maneuvers: [burn,burn]) }
        let past = try ImpulsiveManeuver(epoch: epoch(-1), frame: tree.root, deltaVelocityMetresPerSecond: .zero)
        #expect(throws: ManeuverError.beforeReferenceEpoch) { try ManeuveredTwoBodyPropagator(base: model, maneuvers: [past]) }
        #expect(throws: ManeuverError.nonFiniteDeltaVelocity) {
            try ImpulsiveManeuver(epoch: .j2000, frame: tree.root, deltaVelocityMetresPerSecond: SIMD3(.nan,0,0))
        }
        // A burn that removes all tangential motion creates an unsupported radial segment.
        let radial = try ImpulsiveManeuver(epoch: .j2000, frame: tree.root,
                                           deltaVelocityMetresPerSecond: -model.initialState.velocity)
        #expect(throws: PropagationError.radialTrajectory) { try ManeuveredTwoBodyPropagator(base: model, maneuvers: [radial]) }
    }

    @Test func samplingOrderEndpointsAndBounds() throws {
        let tree = FrameTree(), model = try base(tree.root)
        let dates = try [100.0, 0, -100, 100].map { try epoch($0) }
        let samples = try TrajectorySampler.sample(model, at: dates)
        #expect(samples.map(\.epoch) == dates)
        #expect(samples[0] == samples[3])
        let end = try epoch(10)
        let uniform = try TrajectorySampler.sample(model, from: .j2000, through: end, step: Duration(seconds: 3))
        #expect(try uniform.map { try $0.epoch.duration(since: .j2000).seconds } == [0,3,6,9,10])
        let reversed = try TrajectorySampler.sample(model, from: end, through: .j2000, step: Duration(seconds: 3))
        #expect(try reversed.map { try $0.epoch.duration(since: .j2000).seconds } == [10,7,4,1,0])
        let fractional = try TrajectorySampler.sample(model, from: .j2000, through: epoch(1),
                                                        step: Duration(seconds: 0.1))
        #expect(fractional.count == 11)
        #expect(try fractional.last?.epoch == epoch(1))
        #expect(try TrajectorySampler.sample(model, at: []).isEmpty)
        #expect(try TrajectorySampler.sample(model, from: end, through: end, step: Duration(seconds: 1)).count == 1)
        #expect(throws: TrajectorySamplingError.sampleLimitExceeded) {
            try TrajectorySampler.sample(model, from: .j2000, through: end, step: Duration(seconds: 1), maximumSamples: 5)
        }
        #expect(throws: TrajectorySamplingError.invalidStep) {
            try TrajectorySampler.sample(model, from: .j2000, through: end, step: .zero)
        }
        #expect(throws: TrajectorySamplingError.invalidSampleLimit) {
            try TrajectorySampler.sample(model, at: dates, maximumSamples: 0)
        }
        #expect(throws: TrajectorySamplingError.sampleLimitExceeded) {
            try TrajectorySampler.sample(model, at: dates, maximumSamples: 2)
        }
    }

    @Test func samplingValidatesProviderContractsBeforeConversion() throws {
        final class Bad: StateProvider {
            let frame: Frame
            let other: Frame
            var calls = 0
            var wrongEpoch = false
            init(frame: Frame, other: Frame) { self.frame = frame; self.other = other }
            func state(at epoch: Epoch) throws -> CartesianState {
                calls += 1
                return try CartesianState(epoch: wrongEpoch ? epoch.advanced(by: Duration(seconds: 1)) : epoch,
                    frame: wrongEpoch ? frame : other, position: .zero, velocity: .zero)
            }
        }
        let tree = FrameTree(), other = FrameTree(), bad = Bad(frame: tree.root, other: other.root)
        #expect(throws: FrameQueryError.providerFrameMismatch) { try TrajectorySampler.sample(bad, at: [.j2000]) }
        bad.wrongEpoch = true
        #expect(throws: FrameTransformError.epochMismatch) { try TrajectorySampler.sample(bad, at: [.j2000]) }
        let before = bad.calls
        #expect(throws: TrajectorySamplingError.sampleLimitExceeded) {
            try TrajectorySampler.sample(bad, from: .j2000, through: epoch(10), step: Duration(seconds: 1), maximumSamples: 2)
        }
        #expect(bad.calls == before)
        #expect(throws: TrajectorySamplingError.frameTreeRequired) {
            try TrajectorySampler.sample(bad, at: [.j2000], in: other.root)
        }
    }

    @Test func ephemerisAndSpacecraftGameCheckpoint() throws {
        let tree = try CelestialFrames.makeTree()
        let coverage = try epoch(-1000)...epoch(2000)
        let ephemeris = try VSOP87Ephemeris(body: earth, coverage: coverage, tdbModel: FairheadTDBModel(), accuracy: .low)
        let bodies = try EphemerisCollection(providers: [ephemeris])
        let centre = try FrameConvertedStateProvider(ephemeris, in: CelestialFrames.icrf, using: tree)
        let inertial = try tree.makeGCRF(earth: centre, accelerationPolicy: .assumeZero)
        let spin = try FixedAxisOrientationProvider(referenceFrame: CelestialFrames.icrf, axis: SIMD3(0,0,1), rate: 7.292115e-5)
        let rotating = try tree.makeBodyFixedFrame(name: "Game Earth rotation", centre: centre,
                                                   orientation: spin, accelerationPolicy: .assumeZero)
        let spacecraft = try base(inertial)
        let dates = try [0.0, 600, -300, 0, 1200].map { try epoch($0) }
        let samples = try TrajectorySampler.sample(spacecraft, at: dates, in: rotating, using: tree)
        #expect(samples[0] == samples[3]) // Pause/resume and time warp don't advance stored state.
        for (date, sample) in zip(dates, samples) {
            let planet = try bodies.state(of: .earth, at: date)
            #expect(planet.epoch == date)
            #expect(planet.frame === CelestialFrames.vsop87E)
            #expect(abs(norm(sample.position)-7e6) < 0.001)
            let original = try spacecraft.state(at: date)
            let roundTrip = try tree.convert(sample, to: inertial)
            #expect(norm(roundTrip.position-original.position) < 0.001)
            #expect(norm(roundTrip.velocity-original.velocity) < 1e-7)
        }
    }
}

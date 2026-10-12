// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct BrachistochroneTests {
    private func epoch(_ t: Double) throws -> Epoch { try .j2000.advanced(by: Duration(seconds: t)) }
    private func norm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
    private func vector(_ v: [Double]) -> SIMD3<Double> { SIMD3(v[0],v[1],v[2]) }
    private func trajectory(frame: Frame = Frame(name: "inertial",pseudoInertial: true)) throws -> IdealBrachistochroneTrajectory {
        try IdealBrachistochroneTrajectory(departureEpoch: .j2000,frame: frame,
            departurePositionMetres: .zero,arrivalPositionMetres: SIMD3(100,0,0),
            maximumAccelerationMetresPerSecondSquared: 4)
    }
    private struct Reference: Decodable {
        struct Sample: Decodable { let seconds: Double; let position,velocity,acceleration: [Double] }
        let start,end: [Double]
        let acceleration,distance,seconds,peak: Double
        let samples: [Sample]
    }
    @Test func independentIntegratedReferences() throws {
        let url = try #require(Bundle.module.url(forResource: "brachistochrone-reference",withExtension: "json",subdirectory: "Fixtures"))
        let rows = try JSONDecoder().decode([Reference].self,from: Data(contentsOf: url))
        #expect(rows.count == 4)
        for row in rows {
            let model = try IdealBrachistochroneTrajectory(departureEpoch: .j2000,
                frame: Frame(name: "test inertial",pseudoInertial: true),
                departurePositionMetres: vector(row.start),arrivalPositionMetres: vector(row.end),
                maximumAccelerationMetresPerSecondSquared: row.acceleration)
            #expect(abs(model.duration.seconds-row.seconds) <= max(1e-12,row.seconds*1e-14))
            #expect(abs(model.distanceMetres-row.distance) <= max(1e-12,row.distance*1e-14))
            #expect(abs(model.peakSpeedMetresPerSecond-row.peak) <= max(1e-12,row.peak*1e-14))
            for (index,sample) in row.samples.enumerated() {
                let date = index == 3 ? model.switchEpoch : index == row.samples.count-1 ? model.arrivalEpoch : try epoch(sample.seconds)
                let state = try model.state(at: date)
                #expect(norm(state.position-vector(sample.position)) <= max(1e-10,norm(vector(sample.position))*2e-15))
                #expect(norm(state.velocity-vector(sample.velocity)) <= max(1e-10,row.peak*1e-14))
                #expect(norm(try #require(state.acceleration)-vector(sample.acceleration)) < 1e-12)
                #expect(state.epoch == date && state.frame === model.frame)
            }
        }
    }
    @Test func exactSwitchEndpointsAndCoverage() throws {
        let model = try trajectory()
        #expect(model.duration.seconds == 10)
        #expect(model.switchEpoch == (try epoch(5)))
        #expect(model.peakSpeedMetresPerSecond == 20)
        #expect(model.totalDeltaVelocityMetresPerSecond == 40)
        #expect(model.accelerationDirection == SIMD3(1,0,0) && model.brakingDirection == SIMD3(-1,0,0))
        let start = try model.state(at: .j2000), middle = try model.state(at: model.switchEpoch), end = try model.state(at: model.arrivalEpoch)
        #expect(start.position == .zero && start.velocity == .zero && start.acceleration == SIMD3(4,0,0))
        #expect(middle.position == SIMD3(50,0,0) && middle.velocity == SIMD3(20,0,0) && middle.acceleration == SIMD3(-4,0,0))
        #expect(end.position == SIMD3(100,0,0) && end.velocity == .zero && end.acceleration == .zero)
        #expect(try model.phase(at: .j2000) == .accelerating)
        #expect(try model.phase(at: model.switchEpoch) == .braking)
        #expect(try model.phase(at: model.arrivalEpoch) == .arrived)
        let before = try model.state(at: epoch(5-1e-8)), after = try model.state(at: epoch(5+1e-8))
        #expect(norm(before.velocity-after.velocity) < 1e-12)
        #expect(before.acceleration == SIMD3(4,0,0) && after.acceleration == SIMD3(-4,0,0))
        #expect(norm(before.position-after.position) < 1e-6)
        for date in [try epoch(-1e-8),try epoch(10+1e-8)] {
            #expect(throws: BrachistochroneError.outsideCoverage) { try model.state(at: date) }
            #expect(throws: BrachistochroneError.outsideCoverage) { try model.phase(at: date) }
        }
    }
    @Test func derivativesSamplingAndRotatingFrames() throws {
        let tree = FrameTree(), model = try trajectory(frame: tree.root)
        let rotating = try tree.makeFrame(parent: tree.root,name: "rotating view",
            provider: UniformRotationProvider(axis: SIMD3(0,0,1),rate: 0.1))
        for t in [2.5,7.5] {
            let date = try epoch(t), h = 0.001
            let state = try model.state(at: date), before = try model.state(at: epoch(t-h)), after = try model.state(at: epoch(t+h))
            #expect(norm((after.position-before.position)/(2*h)-state.velocity) < 1e-9)
            #expect(norm((after.velocity-before.velocity)/(2*h)-(try #require(state.acceleration))) < 1e-9)
            let converted = try model.state(at: date,in: rotating,using: tree)
            let b = try model.state(at: epoch(t-h),in: rotating,using: tree)
            let a = try model.state(at: epoch(t+h),in: rotating,using: tree)
            #expect(norm((a.position-b.position)/(2*h)-converted.velocity) < 1e-6)
            #expect(norm((a.velocity-b.velocity)/(2*h)-(try #require(converted.acceleration))) < 1e-6)
        }
        let forward = try TrajectorySampler.sample(model,from: .j2000,through: model.arrivalEpoch,step: Duration(seconds: 2))
        let backward = try TrajectorySampler.sample(model,from: model.arrivalEpoch,through: .j2000,step: Duration(seconds: 2))
        #expect(backward == Array(forward.reversed()))
        let dates = [model.arrivalEpoch,model.switchEpoch,.j2000,model.switchEpoch]
        let samples = try TrajectorySampler.sample(model,at: dates,in: rotating,using: tree)
        #expect(samples.map(\.epoch) == dates && samples[1] == samples[3])
    }
    @Test func coincidentEndpointsAndInvalidInputs() throws {
        let frame = Frame(name: "inertial",pseudoInertial: true)
        let point = SIMD3<Double>(1,2,3)
        let zero = try IdealBrachistochroneTrajectory(departureEpoch: .j2000,frame: frame,
            departurePositionMetres: point,arrivalPositionMetres: point,maximumAccelerationMetresPerSecondSquared: 1)
        #expect(zero.duration == .zero && zero.switchEpoch == .j2000 && zero.arrivalEpoch == .j2000)
        #expect(zero.peakSpeedMetresPerSecond == 0 && zero.accelerationDirection == .zero)
        #expect(try zero.state(at: .j2000).position == point)
        #expect(try zero.state(at: .j2000).acceleration == .zero)
        for value in [0.0,-1,.nan,.infinity] {
            #expect(throws: BrachistochroneError.invalidAcceleration) {
                try IdealBrachistochroneTrajectory(departureEpoch: .j2000,frame: frame,
                    departurePositionMetres: .zero,arrivalPositionMetres: point,maximumAccelerationMetresPerSecondSquared: value)
            }
        }
        #expect(throws: BrachistochroneError.nonInertialFrame) {
            try trajectory(frame: Frame(name: "rotating",pseudoInertial: false))
        }
        #expect(throws: BrachistochroneError.nonFiniteEndpoint) {
            try IdealBrachistochroneTrajectory(departureEpoch: .j2000,frame: frame,
                departurePositionMetres: SIMD3(.nan,0,0),arrivalPositionMetres: point,maximumAccelerationMetresPerSecondSquared: 1)
        }
        #expect(throws: BrachistochroneError.numericalRangeExceeded) {
            try IdealBrachistochroneTrajectory(departureEpoch: .j2000,frame: frame,
                departurePositionMetres: SIMD3(-Double.greatestFiniteMagnitude,0,0),
                arrivalPositionMetres: SIMD3(Double.greatestFiniteMagnitude,0,0),maximumAccelerationMetresPerSecondSquared: 1)
        }
    }
    @Test func scaledProductsAndFineRemainingTime() throws {
        let frame = Frame(name: "inertial",pseudoInertial: true)
        let huge = try IdealBrachistochroneTrajectory(departureEpoch: .j2000,frame: frame,
            departurePositionMetres: .zero,arrivalPositionMetres: SIMD3(1e200,0,0),
            maximumAccelerationMetresPerSecondSquared: 1e200)
        #expect(huge.duration.seconds == 2)
        #expect(huge.peakSpeedMetresPerSecond.isFinite)
        #expect(try huge.state(at: huge.switchEpoch).position.x == 5e199)
        let long = try IdealBrachistochroneTrajectory(departureEpoch: .j2000,frame: frame,
            departurePositionMetres: .zero,arrivalPositionMetres: SIMD3(1e18,0,0),
            maximumAccelerationMetresPerSecondSquared: 1)
        let date = try long.arrivalEpoch.advanced(by: Duration(seconds: -1e-9))
        #expect(try date.duration(since: .j2000).seconds == long.duration.seconds)
        let remaining = try long.arrivalEpoch.duration(since: date).seconds
        let state = try long.state(at: date)
        #expect(remaining > 0 && state.velocity.x > 0)
        #expect(abs(state.velocity.x-remaining) < 1e-20)
        #expect(state.acceleration == SIMD3(-1,0,0))
        #expect(throws: BrachistochroneError.numericalRangeExceeded) {
            try IdealBrachistochroneTrajectory(departureEpoch: epoch(0.5),frame: frame,
                departurePositionMetres: .zero,arrivalPositionMetres: SIMD3(1e-300,0,0),
                maximumAccelerationMetresPerSecondSquared: 1e300)
        }
    }

}

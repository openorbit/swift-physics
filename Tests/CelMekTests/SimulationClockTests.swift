// SPDX-License-Identifier: Apache-2.0
import Testing
import CelMek

@Suite struct SimulationClockTests {
    private func duration(_ seconds: Double) throws -> CelMek.Duration { try CelMek.Duration(seconds: seconds) }
    private func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: duration(seconds)) }
    @Test func pauseReverseSeekAndIndependentPrediction() throws {
        var clock = try SimulationClock(epoch: .j2000,rate: 60,coverage: epoch(-100)...epoch(100))
        #expect(try clock.advance(realDuration: duration(1)) == (try epoch(60)))
        clock.setPaused(true)
        #expect(try clock.advance(realDuration: duration(1000)) == (try epoch(60)))
        #expect(try clock.prediction(after: duration(20)) == (try epoch(80)))
        #expect(clock.epoch == (try epoch(60)))
        try clock.setRate(-10); clock.setPaused(false)
        #expect(try clock.advance(realDuration: duration(2)) == (try epoch(40)))
        try clock.seek(to: epoch(-20))
        #expect(try clock.advance(realDuration: duration(1)) == (try epoch(-30)))
        try clock.setRate(0)
        #expect(try clock.advance(realDuration: duration(1000)) == (try epoch(-30)))
        var copy = clock; try copy.seek(to: .j2000)
        #expect(clock.epoch == (try epoch(-30)) && copy.epoch == .j2000)
    }
    @Test func repeatedSubsecondUpdatesPreserveDistantEpochFractions() throws {
        let start = Epoch(offsetFromJ2000: try CelMek.Duration(wholeSeconds: 8_000_000_000_000,fractionalSecond: 0.125))
        var many = try SimulationClock(epoch: start,rate: 86400), one = many
        for _ in 0..<1000 { try many.advance(realDuration: duration(0.001)) }
        try one.advance(realDuration: duration(1))
        let difference = try many.epoch.duration(since: one.epoch).seconds
        #expect(abs(difference) < 1e-7)
        #expect(one.epoch.offsetFromJ2000.fractionalSecond == 0.125)
        try many.setRate(0.5)
        let before = many.epoch
        try many.advance(realDuration: duration(0.25))
        #expect(try many.epoch.duration(since: before).seconds == 0.125)
    }
    @Test func failuresAreTransactionalAndBoundsAreInclusive() throws {
        var clock = try SimulationClock(epoch: .j2000,rate: 2,coverage: epoch(-2)...epoch(2))
        #expect(throws: SimulationClockError.outsideCoverage) { try clock.advance(realDuration: duration(2)) }
        #expect(clock.epoch == .j2000)
        #expect(try clock.advance(realDuration: duration(1)) == (try epoch(2)))
        #expect(throws: SimulationClockError.outsideCoverage) { try clock.prediction(after: duration(0.1)) }
        #expect(throws: SimulationClockError.outsideCoverage) { try clock.seek(to: epoch(3)) }
        #expect(throws: SimulationClockError.invalidRate) { try clock.setRate(.nan) }
        #expect(clock.rate == 2)
        #expect(clock.epoch == (try epoch(2)))
        clock.setPaused(true)
        #expect(throws: SimulationClockError.invalidRealDuration) { try clock.advance(realDuration: duration(-1)) }
        #expect(throws: SimulationClockError.invalidRate) { try SimulationClock(epoch: .j2000,rate: .infinity) }
        #expect(throws: SimulationClockError.outsideCoverage) { try SimulationClock(epoch: .j2000,coverage: epoch(1)...epoch(2)) }
        var overflow = try SimulationClock(epoch: .j2000,rate: Double.greatestFiniteMagnitude)
        #expect(throws: (any Error).self) { try overflow.advance(realDuration: duration(2)) }
        #expect(overflow.epoch == .j2000)
    }
}

//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2026 Mattias Holm
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

public enum SimulationClockError: Error, Equatable, Sendable {
    case invalidRate, invalidRealDuration, outsideCoverage
}

/// Caller-owned atomic simulation time. No wall-clock reads or provider queries.
/// Positive real elapsed durations can drive forward, paused or reverse playback.
/// A value-type clock is copied independently; serialize mutations in its owner.
public struct SimulationClock: Sendable {
    public private(set) var epoch: Epoch
    /// Simulation seconds per real second; zero freezes and negative rates reverse.
    public private(set) var rate: Double
    public private(set) var isPaused: Bool
    public let coverage: ClosedRange<Epoch>?
    private var anchor: Epoch
    private var realElapsed: Duration = .zero

    public init(epoch: Epoch, rate: Double = 1, isPaused: Bool = false,
                coverage: ClosedRange<Epoch>? = nil) throws {
        guard rate.isFinite else { throw SimulationClockError.invalidRate }
        guard coverage?.contains(epoch) ?? true else { throw SimulationClockError.outsideCoverage }
        self.epoch = epoch; self.anchor = epoch; self.rate = rate
        self.isPaused = isPaused; self.coverage = coverage
    }
    /// Reanchors at the current instant; no jump when changing speed or direction.
    public mutating func setRate(_ rate: Double) throws {
        guard rate.isFinite else { throw SimulationClockError.invalidRate }
        if rate != self.rate { reanchor(); self.rate = rate }
    }
    public mutating func setPaused(_ paused: Bool) {
        if paused != isPaused { reanchor(); isPaused = paused }
    }
    public mutating func seek(to epoch: Epoch) throws {
        try validate(epoch); self.epoch = epoch; reanchor()
    }
    /// Caller supplies an elapsed monotonic duration, not an absolute wall-clock timestamp.
    /// Failure leaves clock state unchanged. Bounds throw; there is no silent clamp.
    @discardableResult public mutating func advance(realDuration: Duration) throws -> Epoch {
        guard realDuration >= .zero else { throw SimulationClockError.invalidRealDuration }
        guard !isPaused, rate != 0 else { return epoch }
        let elapsed = try realElapsed.adding(realDuration)
        // Anchor the whole interval to prevent repeated fractional epoch additions drifting.
        let whole = try Duration(seconds: Double(elapsed.wholeSeconds)*rate)
        let fraction = try Duration(seconds: elapsed.fractionalSecond*rate)
        let next = try anchor.advanced(by: whole.adding(fraction))
        try validate(next)
        realElapsed = elapsed; epoch = next
        return next
    }
    /// Offset is in simulation SI seconds, independent of playback rate/pause state.
    /// Prediction leaves epoch and the accumulated real-time anchor unchanged.
    public func prediction(after offset: Duration) throws -> Epoch {
        let predicted = try epoch.advanced(by: offset)
        try validate(predicted)
        return predicted
    }
    private mutating func reanchor() { anchor = epoch; realElapsed = .zero }
    private func validate(_ epoch: Epoch) throws {
        guard coverage?.contains(epoch) ?? true else { throw SimulationClockError.outsideCoverage }
    }
}

extension SimulationClock {
    /// Includes the elapsed anchor to preserve continuation after a checkpoint.
    public func saved() -> SavedSimulationClock {
        SavedSimulationClock(epoch: SavedEpoch(epoch),anchor: SavedEpoch(anchor),realElapsed: SavedDuration(realElapsed),
            rate: rate,isPaused: isPaused,coverageStart: coverage.map { SavedEpoch($0.lowerBound) },
            coverageEnd: coverage.map { SavedEpoch($0.upperBound) })
    }
    internal static func restoring(_ saved: SavedSimulationClock) throws -> SimulationClock {
        let coverage: ClosedRange<Epoch>?
        switch (saved.coverageStart,saved.coverageEnd) {
        case (nil,nil): coverage = nil
        case let (start?,end?):
            let lower = try start.restore(), upper = try end.restore()
            guard lower <= upper else { throw SimulationSnapshotError.invalidClock }
            coverage = lower...upper
        default: throw SimulationSnapshotError.invalidClock
        }
        let elapsed = try saved.realElapsed.restore()
        guard elapsed >= .zero, !(saved.isPaused || saved.rate == 0) || elapsed == .zero else {
            throw SimulationSnapshotError.invalidClock
        }
        var clock = try SimulationClock(epoch: saved.anchor.restore(),rate: saved.rate,isPaused: saved.isPaused,coverage: coverage)
        try clock.advance(realDuration: elapsed)
        guard clock.epoch == (try saved.epoch.restore()) else { throw SimulationSnapshotError.invalidClock }
        return clock
    }
}

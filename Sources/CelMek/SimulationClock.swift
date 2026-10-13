// SPDX-License-Identifier: Apache-2.0

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

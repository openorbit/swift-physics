// SPDX-License-Identifier: Apache-2.0
import Math

public enum TrajectorySamplingError: Error, Equatable {
    case invalidSampleLimit, sampleLimitExceeded, invalidStep, frameTreeRequired
}

/// Bounded, discrete samples for rendering or prediction. No interpolation or
/// hidden integration steps: each returned value is queried at its own epoch.
public enum TrajectorySampler {
    /// Preserves input order and duplicate epochs. Provider contracts are checked
    /// before optional frame conversion; errors propagate without a partial result.
    public static func sample(_ provider: any StateProvider, at epochs: [Epoch],
                              in destination: Frame? = nil, using tree: FrameTree? = nil,
                              maximumSamples: Int = 100_000) throws -> [CartesianState] {
        guard maximumSamples > 0 else { throw TrajectorySamplingError.invalidSampleLimit }
        guard epochs.count <= maximumSamples else { throw TrajectorySamplingError.sampleLimitExceeded }
        let destination = destination ?? provider.frame
        if destination !== provider.frame && tree == nil { throw TrajectorySamplingError.frameTreeRequired }
        return try epochs.map { epoch in
            let state = try provider.state(at: epoch)
            guard state.epoch == epoch else { throw FrameTransformError.epochMismatch }
            guard state.frame === provider.frame else { throw FrameQueryError.providerFrameMismatch }
            if destination === provider.frame { return state }
            return try tree!.convert(state, to: destination)
        }
    }

    /// Includes both endpoints, with a possibly shorter final interval. A
    /// positive step is used in either time direction. Scheduling is validated
    /// in full before providers are queried; the sample count is bounded.
    public static func sample(_ provider: any StateProvider, from start: Epoch, through end: Epoch,
                              step: Duration, in destination: Frame? = nil, using tree: FrameTree? = nil,
                              maximumSamples: Int = 100_000) throws -> [CartesianState] {
        guard maximumSamples > 0 else { throw TrajectorySamplingError.invalidSampleLimit }
        guard step > .zero else { throw TrajectorySamplingError.invalidStep }
        let forward = start <= end
        var epochs = [start], current = start
        while current != end {
            guard epochs.count < maximumSamples else { throw TrajectorySamplingError.sampleLimitExceeded }
            // Avoid advancing outside the supported epoch range at the final interval.
            let remaining = forward ? try end.duration(since: current) : try current.duration(since: end)
            let next: Epoch
            if remaining <= step { next = end }
            else {
                // Anchor every interval to start; repeated fractional additions
                // would drift and could create a spurious near-duplicate endpoint.
                let (whole, overflow) = step.wholeSeconds.multipliedReportingOverflow(by: Int64(epochs.count))
                guard !overflow else { throw TimeArithmeticError.outOfRange }
                let offset = try Duration(wholeSeconds: whole,
                                          fractionalSecond: step.fractionalSecond * Double(epochs.count))
                let candidate = try start.advanced(by: forward ? offset : offset.negated())
                next = forward ? min(candidate, end) : max(candidate, end)
            }
            guard forward ? next > current : next < current else { throw TrajectorySamplingError.invalidStep }
            epochs.append(next); current = next
        }
        return try sample(provider, at: epochs, in: destination, using: tree, maximumSamples: maximumSamples)
    }
}

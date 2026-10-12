// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

/// Mass-aware trajectory queries, including immutable-seed numerical replay.
public protocol SpacecraftStateProvider: StateProvider {
    func spacecraftState(at epoch: Epoch) throws -> SpacecraftState
}
extension NumericalSpacecraftPropagator: SpacecraftStateProvider {}
extension ScheduledThrustPropagator: SpacecraftStateProvider {}

public enum EventCrossingDirection: Sendable { case any, increasing, decreasing }
public enum PropagationEventAction: Sendable { case record, stop }
public enum PropagationEventError: Error, Equatable {
    case invalidOptions, invalidDetector, nonFiniteValue, evaluationLimitExceeded
    case eventLimitExceeded, localizationLimitExceeded, epochResolutionExceeded
}

/// A continuous scalar zero-crossing function. Pure evaluations may be repeated
/// and out of order. Sign direction is defined with increasing physical time.
public protocol PropagationEventDetector {
    var identifier: String { get }
    var direction: EventCrossingDirection { get }
    var action: PropagationEventAction { get }
    func value(at state: SpacecraftState) throws -> Double
}

public struct EpochEvent: PropagationEventDetector {
    public let identifier: String
    public let epoch: Epoch
    public let direction: EventCrossingDirection = .any
    public let action: PropagationEventAction
    public init(identifier: String, epoch: Epoch, action: PropagationEventAction = .stop) throws {
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else { throw PropagationEventError.invalidDetector }
        self.identifier = identifier; self.epoch = epoch; self.action = action
    }
    public func value(at state: SpacecraftState) throws -> Double { try state.epoch.duration(since: epoch).seconds }
}

/// Sphere centred at a fixed location in explicit trajectory axes. The caller
/// selects the physical origin and body radius. Moving shapes use custom detectors.
public struct SphericalImpactEvent: PropagationEventDetector {
    public let identifier: String
    public let frame: Frame
    public let centreMetres: SIMD3<Double>
    public let radiusMetres: Double
    public let direction: EventCrossingDirection
    public let action: PropagationEventAction
    public init(identifier: String, frame: Frame, centreMetres: SIMD3<Double> = .zero,
                radiusMetres: Double, direction: EventCrossingDirection = .decreasing,
                action: PropagationEventAction = .stop) throws {
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              radiusMetres.isFinite, radiusMetres > 0,
              [centreMetres.x,centreMetres.y,centreMetres.z].allSatisfy(\.isFinite) else {
            throw PropagationEventError.invalidDetector
        }
        self.identifier = identifier; self.frame = frame; self.centreMetres = centreMetres
        self.radiusMetres = radiusMetres; self.direction = direction; self.action = action
    }
    public func value(at state: SpacecraftState) throws -> Double {
        guard state.frame === frame else { throw SpacecraftDynamicsError.frameMismatch }
        let p = state.cartesian.position-centreMetres
        return hypot(hypot(p.x,p.y),p.z)-radiusMetres
    }
}

public struct PropagationEventOptions: Equatable, Sendable {
    public let maximumCheckIntervalSeconds: Double
    public let timeToleranceSeconds: Double
    public let maximumLocalizationIterations: Int
    public let maximumStateEvaluations: Int
    public let maximumEvents: Int
    public init(maximumCheckIntervalSeconds: Double = 60, timeToleranceSeconds: Double = 1e-6,
                maximumLocalizationIterations: Int = 80, maximumStateEvaluations: Int = 100_000,
                maximumEvents: Int = 10_000) throws {
        guard maximumCheckIntervalSeconds.isFinite, maximumCheckIntervalSeconds > 0,
              timeToleranceSeconds.isFinite, timeToleranceSeconds > 0,
              maximumLocalizationIterations > 0, maximumStateEvaluations > 0, maximumEvents > 0 else {
            throw PropagationEventError.invalidOptions
        }
        self.maximumCheckIntervalSeconds = maximumCheckIntervalSeconds
        self.timeToleranceSeconds = timeToleranceSeconds
        self.maximumLocalizationIterations = maximumLocalizationIterations
        self.maximumStateEvaluations = maximumStateEvaluations; self.maximumEvents = maximumEvents
    }
}

public struct PropagationEventOccurrence {
    public let identifier: String
    public let action: PropagationEventAction
    public let state: SpacecraftState
    /// Width of the final sign bracket; zero for an exactly sampled root.
    public let bracketDurationSeconds: Double
}
public struct EventPropagationResult {
    public let state: SpacecraftState
    public let events: [PropagationEventOccurrence]
    public let reachedTarget: Bool
    public let stateEvaluations: Int
}

/// Bounded zero-crossing search over a mass-aware trajectory. Localization queries
/// the actual provider (no linear state interpolation). Underlying provider errors
/// propagate. This searches smooth trajectory segments; it does not apply impulses,
/// change forces or localize a root beyond the provider's valid domain.
public struct EventDetectedPropagation {
    public let trajectory: any SpacecraftStateProvider
    public let detectors: [any PropagationEventDetector]
    public let options: PropagationEventOptions
    public init(trajectory: any SpacecraftStateProvider, detectors: [any PropagationEventDetector],
                options: PropagationEventOptions = try! PropagationEventOptions()) throws {
        let identifiers = detectors.map { $0.identifier }
        guard identifiers.allSatisfy({ !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }),
              Set(identifiers).count == identifiers.count else { throw PropagationEventError.invalidDetector }
        self.trajectory = trajectory; self.detectors = detectors; self.options = options
    }

    public func propagation(from start: Epoch, to target: Epoch, includeInitialEvents: Bool = true) throws -> EventPropagationResult {
        var evaluations = 0, occurrences: [PropagationEventOccurrence] = []
        func query(_ epoch: Epoch) throws -> SpacecraftState {
            guard evaluations < options.maximumStateEvaluations else { throw PropagationEventError.evaluationLimitExceeded }
            evaluations += 1
            let state = try trajectory.spacecraftState(at: epoch)
            guard state.frame === trajectory.frame else { throw SpacecraftDynamicsError.frameMismatch }
            guard state.epoch == epoch else { throw SpacecraftDynamicsError.epochMismatch }
            return state
        }
        func value(_ detector: any PropagationEventDetector, _ state: SpacecraftState) throws -> Double {
            let result = try detector.value(at: state)
            guard result.isFinite else { throw PropagationEventError.nonFiniteValue }
            return result
        }
        func append(_ occurrence: PropagationEventOccurrence) throws {
            guard occurrences.count < options.maximumEvents else { throw PropagationEventError.eventLimitExceeded }
            occurrences.append(occurrence)
        }
        var current = try query(start)
        var values = try detectors.map { try value($0,current) }
        if includeInitialEvents {
            for i in detectors.indices where values[i] == 0 {
                try append(PropagationEventOccurrence(identifier: detectors[i].identifier, action: detectors[i].action,
                    state: current, bracketDurationSeconds: 0))
            }
            if occurrences.contains(where: { $0.action == .stop }) {
                return EventPropagationResult(state: current, events: occurrences, reachedTarget: start == target, stateEvaluations: evaluations)
            }
        }
        let forward = target >= start
        while current.epoch != target {
            let remaining = abs(try target.duration(since: current.epoch).seconds)
            let end = remaining <= options.maximumCheckIntervalSeconds ? target :
                try current.epoch.advanced(by: Duration(seconds: forward ? options.maximumCheckIntervalSeconds : -options.maximumCheckIntervalSeconds))
            guard forward ? end > current.epoch : end < current.epoch else { throw PropagationEventError.epochResolutionExceeded }
            let next = try query(end)
            let nextValues = try detectors.map { try value($0,next) }
            var found: [(Int,PropagationEventOccurrence)] = []
            for i in detectors.indices {
                let a = values[i], b = nextValues[i], detector = detectors[i]
                // A root at the left endpoint was already handled. Tangencies and
                // multiple crossings within one check interval are not guaranteed.
                guard a != 0, b == 0 || (a < 0) != (b < 0) else { continue }
                let increasing = forward ? b > a : a > b
                guard detector.direction == .any || (detector.direction == .increasing) == increasing else { continue }
                var left = current, right = next, leftValue = a
                var width = abs(try right.epoch.duration(since: left.epoch).seconds)
                var root = next
                if let timed = detector as? EpochEvent {
                    root = try query(timed.epoch); width = 0
                } else if b == 0 { width = 0 } else {
                    var iterations = 0
                    while width > options.timeToleranceSeconds {
                        guard iterations < options.maximumLocalizationIterations else { throw PropagationEventError.localizationLimitExceeded }
                        iterations += 1
                        let date = try left.epoch.advanced(by: Duration(seconds: right.epoch.duration(since: left.epoch).seconds/2))
                        guard date != left.epoch, date != right.epoch else { throw PropagationEventError.epochResolutionExceeded }
                        let middle = try query(date), midValue = try value(detector,middle)
                        if midValue == 0 { left = middle; right = middle; width = 0; break }
                        if (midValue < 0) == (leftValue < 0) { left = middle; leftValue = midValue } else { right = middle }
                        width = abs(try right.epoch.duration(since: left.epoch).seconds)
                    }
                    root = abs(try value(detector,left)) <= abs(try value(detector,right)) ? left : right
                }
                found.append((i,PropagationEventOccurrence(identifier: detector.identifier, action: detector.action,
                    state: root, bracketDurationSeconds: width)))
            }
            found.sort {
                if $0.1.state.epoch == $1.1.state.epoch { return $0.0 < $1.0 }
                return forward ? $0.1.state.epoch < $1.1.state.epoch : $0.1.state.epoch > $1.1.state.epoch
            }
            for (_,event) in found {
                try append(event)
                if event.action == .stop {
                    return EventPropagationResult(state: event.state, events: occurrences,
                        reachedTarget: event.state.epoch == target, stateEvaluations: evaluations)
                }
            }
            current = next; values = nextValues
        }
        return EventPropagationResult(state: current, events: occurrences, reachedTarget: true, stateEvaluations: evaluations)
    }
}

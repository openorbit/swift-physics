// SPDX-License-Identifier: Apache-2.0
import Foundation

public enum CircularTransferObjective: Equatable, Sendable { case minimumDeltaVelocity, minimumDuration }
public enum CircularTransferStrategy {
    case hohmann(HohmannTransfer), biElliptic(BiEllipticTransfer)
    public var duration: Duration {
        switch self { case .hohmann(let t): return t.duration; case .biElliptic(let t): return t.duration }
    }
    public var totalDeltaVelocityMetresPerSecond: Double {
        switch self {
        case .hohmann(let t): return t.totalDeltaVelocityMetresPerSecond
        case .biElliptic(let t): return t.totalDeltaVelocityMetresPerSecond
        }
    }
}
public struct CircularTransferComparison {
    public let rankedCandidates: [CircularTransferStrategy]
    public let rejectedCandidates: [CircularTransferStrategy]
    public var best: CircularTransferStrategy? { rankedCandidates.first }
    /// Coplanar ideal-model comparison only. Bounds filter duration and delta-v;
    /// they do not certify phasing, clearance or perturbed trajectory feasibility.
    public init(departureRadiusMetres: Double, arrivalRadiusMetres: Double,
                intermediateRadiiMetres: [Double], gravitationalParameter: GravitationalParameter,
                objective: CircularTransferObjective, maximumDuration: Duration? = nil,
                maximumDeltaVelocityMetresPerSecond: Double? = nil) throws {
        guard intermediateRadiiMetres.count <= 10_000,
              maximumDuration == nil || maximumDuration! > .zero,
              maximumDeltaVelocityMetresPerSecond == nil ||
                (maximumDeltaVelocityMetresPerSecond!.isFinite && maximumDeltaVelocityMetresPerSecond! >= 0) else {
            throw BiEllipticTransferError.invalidComparison
        }
        var candidates: [CircularTransferStrategy] = [.hohmann(try HohmannTransfer(
            departureRadiusMetres: departureRadiusMetres, arrivalRadiusMetres: arrivalRadiusMetres,
            gravitationalParameter: gravitationalParameter))]
        for radius in Array(Set(intermediateRadiiMetres)).sorted() {
            candidates.append(.biElliptic(try BiEllipticTransfer(departureRadiusMetres: departureRadiusMetres,
                arrivalRadiusMetres: arrivalRadiusMetres, intermediateRadiusMetres: radius,
                gravitationalParameter: gravitationalParameter)))
        }
        var allowed: [CircularTransferStrategy] = [], rejected: [CircularTransferStrategy] = []
        for candidate in candidates {
            if (maximumDuration == nil || candidate.duration <= maximumDuration!) &&
                (maximumDeltaVelocityMetresPerSecond == nil ||
                    candidate.totalDeltaVelocityMetresPerSecond <= maximumDeltaVelocityMetresPerSecond!) {
                allowed.append(candidate)
            } else { rejected.append(candidate) }
        }
        rankedCandidates = allowed.enumerated().sorted { a,b in
            let ad = a.element.totalDeltaVelocityMetresPerSecond, bd = b.element.totalDeltaVelocityMetresPerSecond
            let at = a.element.duration, bt = b.element.duration
            switch objective {
            case .minimumDeltaVelocity:
                if ad != bd { return ad < bd }; if at != bt { return at < bt }
            case .minimumDuration:
                if at != bt { return at < bt }; if ad != bd { return ad < bd }
            }
            return a.offset < b.offset
        }.map(\.element)
        rejectedCandidates = rejected
    }
}

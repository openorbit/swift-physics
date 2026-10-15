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
import Foundation
import Math

public struct DE441RecordCacheStatistics: Equatable, Sendable {
    public let cachedRecords: Int
    public let cacheHits: Int
    /// Successful record reads, including constructor validation. Clearing the
    /// cache does not reset these cumulative counters.
    public let diskReads: Int
}

extension DE441Body {
    public var bodyIdentifier: BodyIdentifier {
        switch self {
        case .mercury: return .mercury
        case .venus: return .venus
        case .earth: return .earth
        case .marsBarycentre: return .marsBarycentre
        case .jupiterBarycentre: return .jupiterBarycentre
        case .saturnBarycentre: return .saturnBarycentre
        case .uranusBarycentre: return .uranusBarycentre
        case .neptuneBarycentre: return .neptuneBarycentre
        case .plutoBarycentre: return .plutoBarycentre
        case .moon: return .moon
        case .sun: return .sun
        case .solarSystemBarycentre: return .solarSystemBarycentre
        case .earthMoonBarycentre: return .earthMoonBarycentre
        }
    }
}

extension DE441Ephemeris {
    public func vector(target: DE441Body, center: DE441Body, at date: JulianDateComponents) throws -> EphemerisVector {
        try vector(target: target.rawValue, center: center.rawValue, at: date)
    }
}

/// Common SI adapter for one file-backed target/centre pair. The supplied TDB
/// model defines coordinate-time conversion and bounds usable epoch coverage.
public struct DE441StateProvider: AccuracyAwareEphemerisProvider {
    public let metadata: EphemerisMetadata
    public let accuracy: AccuracyLevel
    public let effectiveAccuracy: AccuracyLevel = .full
    public var availableAccuracyLevels: Set<AccuracyLevel> { [.full] }
    public let target: DE441Body
    public let center: DE441Body
    public let reader: DE441Ephemeris
    private let tdbModel: any TDBTimeModel

    /// For a non-SSB centre, frame is mandatory and must represent that centre
    /// with ICRF-aligned nonrotating axes. No frame relabelling/conversion is inferred.
    public init(reader: DE441Ephemeris, target: DE441Body, center: DE441Body = .solarSystemBarycentre,
                frame: Frame? = nil, coverage: ClosedRange<Epoch>, tdbModel: any TDBTimeModel,
                accuracy: AccuracyLevel = .full) throws {
        let outputFrame: Frame
        if let frame { outputFrame = frame }
        else {
            guard center == .solarSystemBarycentre else { throw EphemerisError.invalidMetadata }
            outputFrame = CelestialFrames.icrf
        }
        guard outputFrame.pseudoInertial else { throw EphemerisError.invalidMetadata }
        guard center != .solarSystemBarycentre || outputFrame === CelestialFrames.icrf else {
            throw EphemerisError.frameMismatch
        }
        guard coverage.lowerBound >= tdbModel.coverage.lowerBound,
              coverage.upperBound <= tdbModel.coverage.upperBound else { throw EphemerisError.outsideCoverage }
        try reader.validateDate(coverage.lowerBound.julianDate(in: .tdb, using: tdbModel))
        try reader.validateDate(coverage.upperBound.julianDate(in: .tdb, using: tdbModel))
        self.reader = reader; self.target = target; self.center = center
        self.tdbModel = tdbModel; self.accuracy = accuracy
        metadata = try EphemerisMetadata(body: target.bodyIdentifier, origin: center.bodyIdentifier,
            frame: outputFrame, coverage: coverage, coverageBasis: .callerSelected, modelTimeScale: .tdb,
            velocityMethod: .analytic, accelerationMethod: .unavailable,
            provenance: PhysicalDataProvenance(source: "JPL native DE441; \(tdbModel.provenance)", version: "DE441",
                reference: "https://ssd.jpl.nasa.gov/planets/eph_export.html"),
            accuracy: accuracy, effectiveAccuracy: .full)
    }

    public func state(at epoch: Epoch) throws -> CartesianState {
        guard metadata.coverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
        let correction = try epoch.tdbCorrection(using: tdbModel)
        let date = try epoch.tdbDate(using: correction)
        let native = try reader.vector(target: target, center: center, at: date)
        let length = reader.metadata.astronomicalUnitKilometres * 1000
        return try CartesianState(epoch: epoch, frame: metadata.frame, position: native.position * length,
            velocity: native.velocity * (length / 86400 * correction.rate))
    }
}

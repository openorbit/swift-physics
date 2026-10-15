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

public enum SimulationSnapshotError: Error, Equatable {
    case unsupportedVersion, invalidIdentifier, duplicateIdentifier, invalidTime, invalidClock
}

/// Stable wire time: decimal whole seconds avoid JSON's 53-bit integer limitation.
/// Fractions must be canonical, finite and in [0,1).
public struct SavedDuration: Codable, Equatable, Sendable {
    public let wholeSeconds: String
    public let fractionalSecond: Double
    public init(_ duration: Duration) {
        wholeSeconds = String(duration.wholeSeconds); fractionalSecond = duration.fractionalSecond
    }
    public func restore() throws -> Duration {
        guard let whole = Int64(wholeSeconds), String(whole) == wholeSeconds,
              fractionalSecond.isFinite, fractionalSecond >= 0, fractionalSecond < 1 else {
            throw SimulationSnapshotError.invalidTime
        }
        return try Duration(wholeSeconds: whole, fractionalSecond: fractionalSecond)
    }
}

/// Atomic offset from JD 2451545.0 TT. No UTC/TDB conversion occurs on restore.
public struct SavedEpoch: Codable, Equatable, Sendable {
    public let offsetFromJ2000: SavedDuration
    public init(_ epoch: Epoch) { offsetFromJ2000 = SavedDuration(epoch.offsetFromJ2000) }
    public func restore() throws -> Epoch { try Epoch(offsetFromJ2000: offsetFromJ2000.restore()) }
}

public struct SavedVector3: Codable, Equatable, Sendable {
    public let x, y, z: Double
    public init(_ vector: SIMD3<Double>) { x = vector.x; y = vector.y; z = vector.z }
    public var vector: SIMD3<Double> { SIMD3(x,y,z) }
}

/// Application IDs identify physical frame bindings; names and object addresses do not.
public struct SavedCartesianState: Codable, Sendable {
    public let id, frameID: String
    public let epoch: SavedEpoch
    public let positionMetres, velocityMetresPerSecond: SavedVector3
    public let accelerationMetresPerSecondSquared: SavedVector3?
    public let dryMassKilograms, propellantMassKilograms: Double?
    public init(id: String, frameID: String, state: CartesianState, mass: SpacecraftMassProperties? = nil) throws {
        try snapshotID(id); try snapshotID(frameID)
        self.id = id; self.frameID = frameID; epoch = SavedEpoch(state.epoch)
        positionMetres = SavedVector3(state.position); velocityMetresPerSecond = SavedVector3(state.velocity)
        accelerationMetresPerSecondSquared = state.acceleration.map(SavedVector3.init)
        dryMassKilograms = mass?.dryMassKilograms; propellantMassKilograms = mass?.propellantMassKilograms
    }
    public func restore(using resolveFrame: (String) throws -> Frame) throws -> CartesianState {
        try snapshotID(id); try snapshotID(frameID)
        return try CartesianState(epoch: epoch.restore(), frame: resolveFrame(frameID),
            position: positionMetres.vector, velocity: velocityMetresPerSecond.vector,
            acceleration: accelerationMetresPerSecondSquared?.vector)
    }
    public func restoreMass() throws -> SpacecraftMassProperties? {
        switch (dryMassKilograms,propellantMassKilograms) {
        case (nil,nil): return nil
        case let (dry?,fuel?): return try SpacecraftMassProperties(dryMassKilograms: dry,propellantMassKilograms: fuel)
        default: throw SpacecraftDynamicsError.invalidMass
        }
    }
}

public struct SavedSimulationClock: Codable, Sendable {
    public let epoch, anchor: SavedEpoch
    public let realElapsed: SavedDuration
    public let rate: Double
    public let isPaused: Bool
    public let coverageStart, coverageEnd: SavedEpoch?
    public func restore() throws -> SimulationClock { try SimulationClock.restoring(self) }
}

public struct SavedBodyIdentifier: Codable, Sendable {
    public let naif: Int32?
    public let namespace, name: String?
    public init(_ body: BodyIdentifier) {
        switch body {
        case .naif(let code): naif = code; namespace = nil; name = nil
        case .custom(let space,let label): naif = nil; namespace = space; name = label
        }
    }
    public func restore() throws -> BodyIdentifier {
        switch (naif,namespace,name) {
        case let (code?,nil,nil): return .naif(code)
        case let (nil,space?,label?):
            let body = BodyIdentifier.custom(namespace: space,name: label); try body.validate(); return body
        default: throw SimulationSnapshotError.invalidIdentifier
        }
    }
}

/// Explicit two-body inputs, including all solver options; no changing defaults on load.
public struct SavedTwoBodyInput: Codable, Sendable {
    public let id: String
    public let initialState: SavedCartesianState
    public let centralBody: SavedBodyIdentifier
    public let gmCubicMetresPerSecondSquared: Double
    public let gmConvention, gmSource, gmVersion: String
    public let gmReference: String?
    public let absoluteTimeToleranceSeconds, relativeTimeTolerance: Double
    public let maximumIterations, maximumBracketExpansions: Int
    public init(id: String, frameID: String, propagator: TwoBodyPropagator) throws {
        try snapshotID(id); self.id = id
        initialState = try SavedCartesianState(id: id,frameID: frameID,state: propagator.initialState)
        centralBody = SavedBodyIdentifier(propagator.centralBody)
        let gm = propagator.gravitationalParameter
        gmCubicMetresPerSecondSquared = gm.cubicMetresPerSecondSquared
        switch gm.convention {
        case .newtonian: gmConvention = "newtonian"
        case .ttCompatible: gmConvention = "ttCompatible"
        case .tdbCompatible: gmConvention = "tdbCompatible"
        case .tcbCompatible: gmConvention = "tcbCompatible"
        }
        gmSource = gm.provenance.source; gmVersion = gm.provenance.version; gmReference = gm.provenance.reference
        absoluteTimeToleranceSeconds = propagator.options.absoluteTimeToleranceSeconds
        relativeTimeTolerance = propagator.options.relativeTimeTolerance
        maximumIterations = propagator.options.maximumIterations
        maximumBracketExpansions = propagator.options.maximumBracketExpansions
    }
    public func restore(using resolveFrame: (String) throws -> Frame) throws -> TwoBodyPropagator {
        try snapshotID(id)
        let convention: GravitationalParameterConvention
        switch gmConvention {
        case "newtonian": convention = .newtonian
        case "ttCompatible": convention = .ttCompatible
        case "tdbCompatible": convention = .tdbCompatible
        case "tcbCompatible": convention = .tcbCompatible
        default: throw SimulationSnapshotError.unsupportedVersion
        }
        let gm = try GravitationalParameter(cubicMetresPerSecondSquared: gmCubicMetresPerSecondSquared,
            convention: convention,provenance: PhysicalDataProvenance(source: gmSource,version: gmVersion,reference: gmReference))
        return try TwoBodyPropagator(initialState: initialState.restore(using: resolveFrame),centralBody: centralBody.restore(),
            gravitationalParameter: gm,options: TwoBodyPropagationOptions(absoluteTimeToleranceSeconds: absoluteTimeToleranceSeconds,
            relativeTimeTolerance: relativeTimeTolerance,maximumIterations: maximumIterations,maximumBracketExpansions: maximumBracketExpansions))
    }
}

/// Application-owned reconstruction metadata, not an executable provider or asset loader.
public struct SimulationDependency: Codable, Sendable {
    public let id, kind, version: String
    public let parameters: [String:String]
    public init(id: String, kind: String, version: String, parameters: [String:String] = [:]) throws {
        try snapshotID(id); try snapshotID(kind); try snapshotID(version)
        self.id = id; self.kind = kind; self.version = version; self.parameters = parameters
    }
}

/// Version 1 portable SI snapshots. Decode checks structural/version invariants;
/// restoration uses validating runtime constructors and caller-owned frame bindings.
public struct SimulationSnapshot: Codable, Sendable {
    public let schemaVersion: Int
    public let format: String
    public let clock: SavedSimulationClock?
    public let states: [SavedCartesianState]
    public let twoBodyInputs: [SavedTwoBodyInput]
    public let dependencies: [SimulationDependency]
    public init(clock: SimulationClock? = nil, states: [SavedCartesianState] = [],
                twoBodyInputs: [SavedTwoBodyInput] = [], dependencies: [SimulationDependency] = []) throws {
        schemaVersion = 1; format = "celmek.simulation"
        self.clock = clock?.saved(); self.states = states; self.twoBodyInputs = twoBodyInputs; self.dependencies = dependencies
        try validate()
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, format, clock, states, twoBodyInputs, dependencies }
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self,forKey: .schemaVersion)
        format = try c.decode(String.self,forKey: .format)
        guard schemaVersion == 1, format == "celmek.simulation" else { throw SimulationSnapshotError.unsupportedVersion }
        clock = try c.decodeIfPresent(SavedSimulationClock.self,forKey: .clock)
        states = try c.decode([SavedCartesianState].self,forKey: .states)
        twoBodyInputs = try c.decode([SavedTwoBodyInput].self,forKey: .twoBodyInputs)
        dependencies = try c.decode([SimulationDependency].self,forKey: .dependencies)
        try validate()
    }
    private func validate() throws {
        func unique(_ ids: [String]) throws {
            for id in ids { try snapshotID(id) }
            guard Set(ids).count == ids.count else { throw SimulationSnapshotError.duplicateIdentifier }
        }
        try unique(states.map(\.id)); try unique(twoBodyInputs.map(\.id)); try unique(dependencies.map(\.id))
        _ = try clock?.restore()
        // Temporary identities validate numeric payloads without trusting saved frame descriptions.
        let inertial = Frame(name: "snapshot validation",pseudoInertial: true)
        for state in states { _ = try state.restore { _ in inertial }; _ = try state.restoreMass() }
        for input in twoBodyInputs {
            _ = try input.initialState.restoreMass()
            _ = try input.restore { _ in inertial }
        }
        for dependency in dependencies { try snapshotID(dependency.kind); try snapshotID(dependency.version) }
    }
}

private func snapshotID(_ value: String) throws {
    guard !value.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
        throw SimulationSnapshotError.invalidIdentifier
    }
}

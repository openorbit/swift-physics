// SPDX-License-Identifier: Apache-2.0
import Foundation
import CelMek
import Math

// Run in release mode. These are measurements, never CI timing assertions.
struct Timing: Encodable {
    let workload: String
    let operationsPerBatch: Int
    let medianNanosecondsPerOperation, minimumNanosecondsPerOperation, maximumNanosecondsPerOperation: Double
}
struct Fidelity: Encodable {
    let mode: String
    let retainedTerms, originalTerms, comparisonEpochs: Int
    let maximumPositionDifferenceMetres, maximumVelocityDifferenceMetresPerSecond: Double
    let preparationSeconds: Double
}
struct Report: Encodable {
    let schemaVersion: Int
    let operatingSystem, architecture, timingMethod, coverage: String
    let logicalProcessors, repetitions: Int
    let checksum: Double
    let timings: [Timing]
    let fidelity: [Fidelity]
}
func now() -> UInt64 { DispatchTime.now().uptimeNanoseconds }
func magnitude(_ v: SIMD3<Double>) -> Double { sqrt((v*v).sum()) }
var checksum = 0.0, timings: [Timing] = []
let repetitions = 7
func measure(_ name: String, operations: Int, prepare: () throws -> Void = {},
             work: () throws -> Double) throws {
    try prepare(); checksum += try work() // Warmup outside timing.
    var values: [Double] = []
    for _ in 0..<repetitions {
        try prepare()
        let start = now()
        let result = try work()
        let elapsed = now()-start
        checksum += result
        values.append(Double(elapsed)/Double(operations))
    }
    values.sort()
    timings.append(Timing(workload: name, operationsPerBatch: operations,
        medianNanosecondsPerOperation: values[repetitions/2],
        minimumNanosecondsPerOperation: values[0], maximumNanosecondsPerOperation: values.last!))
}
let end = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: 365.25*86400))
let coverage = Epoch.j2000...end
let epochs = try (0..<64).map { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: Double($0)*365.25*86400/63)) }
let tdb = try FairheadTDBModel()
let bodies = [sun,mercury,venus,earth,mars,jupiter,saturn,uranus,neptune]
var models: [[VSOP87Ephemeris]] = [], preparation: [Double] = []
for mode in AccuracyLevel.allCases {
    let start = now()
    models.append(try bodies.map { try VSOP87Ephemeris(body: $0,coverage: coverage,tdbModel: tdb,accuracy: mode) })
    preparation.append(Double(now()-start)*1e-9)
}
var fidelity: [Fidelity] = []
for (index, mode) in AccuracyLevel.allCases.enumerated() {
    let providers = models[index]
    var position = 0.0, velocity = 0.0
    for epoch in epochs {
        for body in bodies.indices {
            let state = try providers[body].state(at: epoch), reference = try models[3][body].state(at: epoch)
            position = max(position,magnitude(state.position-reference.position))
            velocity = max(velocity,magnitude(state.velocity-reference.velocity))
        }
    }
    fidelity.append(Fidelity(mode: mode.rawValue,
        retainedTerms: providers.reduce(0) { $0+$1.coefficientSelection.retainedTermCount },
        originalTerms: providers.reduce(0) { $0+$1.coefficientSelection.originalTermCount }, comparisonEpochs: epochs.count,
        maximumPositionDifferenceMetres: position, maximumVelocityDifferenceMetresPerSecond: velocity,
        preparationSeconds: preparation[index]))
    try measure("nine-body VSOP87E \(mode.rawValue)",operations: epochs.count*providers.count) {
        var sum = 0.0
        for epoch in epochs { for provider in providers { let state = try provider.state(at: epoch); sum += state.position.x+state.velocity.y } }
        return sum
    }
}
let cache = try MemoryEphemerisCache(capacity: 128)
let cached = try CachedEphemerisProvider(provider: models[0][3],cache: cache)
try measure("Earth low cache miss",operations: epochs.count,prepare: { try cached.invalidate() }) {
    try cached.states(at: epochs).reduce(0) { $0+$1.position.x }
}
try measure("Earth low cache hit",operations: epochs.count) {
    try cached.states(at: epochs).reduce(0) { $0+$1.position.x }
}
let frame = Frame(name: "benchmark inertial",pseudoInertial: true)
let gm = try GravitationalParameter(cubicMetresPerSecondSquared: 3.986004418e14,convention: .newtonian,
    provenance: PhysicalDataProvenance(source: "benchmark Newtonian Earth",version: "1"))
let seed = try CartesianState(epoch: .j2000,frame: frame,position: SIMD3(7e6,0,0),velocity: SIMD3(0,sqrt(3.986004418e14/7e6),0))
let coast = try TwoBodyPropagator(initialState: seed,centralBody: .earth,gravitationalParameter: gm)
let orbitEpochs = try (0..<256).map { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: Double($0)*60)) }
try measure("two-body state",operations: orbitEpochs.count) {
    var sum = 0.0
    for epoch in orbitEpochs { let state = try coast.state(at: epoch); sum += state.position.x+state.velocity.y }
    return sum
}
try measure("256-point orbit sampling",operations: orbitEpochs.count) {
    try TrajectorySampler.sample(coast,at: orbitEpochs).reduce(0) { $0+$1.position.x }
}
let tree = try CelestialFrames.makeTree()
let native = try models[0][3].states(at: epochs)
try measure("VSOP axes to ICRF conversion",operations: native.count) {
    var sum = 0.0
    for state in native { sum += try tree.convert(state,to: CelestialFrames.icrf).position.x }
    return sum
}
let origin = try CartesianState(epoch: .j2000,frame: models[0][3].frame,position: native[0].position,velocity: native[0].velocity)
let renderer = try RendererAdapter(origin: origin,metresPerRenderUnit: 1e6)
let renderStates = try (0..<1024).map { index in
    try CartesianState(epoch: .j2000,frame: origin.frame,position: origin.position+SIMD3(Double(index)*1000,2000,-3000),velocity: origin.velocity)
}
try measure("moving-origin Float adapter",operations: renderStates.count) {
    var sum = 0.0
    for state in renderStates { sum += Double(try renderer.adapt(state).position.x) }
    return sum
}
#if arch(arm64)
let architecture = "arm64"
#elseif arch(x86_64)
let architecture = "x86_64"
#else
let architecture = "other"
#endif
let report = Report(schemaVersion: 1,operatingSystem: ProcessInfo.processInfo.operatingSystemVersionString,
    architecture: architecture,timingMethod: "Release required; monotonic uptime; warmup then seven sequential batches; setup and cache invalidation excluded; median/min/max per operation; checksum consumed",
    coverage: "J2000 through 365.25 TT days; 64 equally spaced ephemeris epochs",logicalProcessors: ProcessInfo.processInfo.processorCount,
    repetitions: repetitions,checksum: checksum,timings: timings,fidelity: fidelity)
let encoder = JSONEncoder(); encoder.outputFormatting = [.prettyPrinted,.sortedKeys]
print(String(decoding: try encoder.encode(report),as: UTF8.self))

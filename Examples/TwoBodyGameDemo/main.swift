// SPDX-License-Identifier: Apache-2.0
import Foundation
import CelMek

// swift run --disable-sandbox TwoBodyGameDemo
// The output is SI state data suitable for an application's rendering adapter.
let start = Epoch.j2000
let end = try start.advanced(by: Duration(seconds: 600))
let tree = try CelestialFrames.makeTree()
let planets = try VSOP87Ephemeris(body: earth, coverage: start...end,
                                 tdbModel: FairheadTDBModel(), accuracy: .low)
let centre = try FrameConvertedStateProvider(planets, in: CelestialFrames.icrf, using: tree)
// The prescribed origin lacks acceleration; this policy is chosen explicitly.
let inertial = try tree.makeGCRF(earth: centre, accelerationPolicy: .assumeZero)
let spin = try FixedAxisOrientationProvider(referenceFrame: CelestialFrames.icrf,
                                            axis: SIMD3(0,0,1), rate: 7.292115e-5)
let rotating = try tree.makeBodyFixedFrame(name: "Game Earth spin", centre: centre,
                                          orientation: spin, accelerationPolicy: .assumeZero)
// Explicit Newtonian game model; no relativistic scaling of catalogue GM is inferred.
let gm = try GravitationalParameter(cubicMetresPerSecondSquared: 3.986004418e14,
    convention: .newtonian, provenance: PhysicalDataProvenance(source: "Earth game model",
    version: "1", reference: "ORBITAL_MECHANICS.md"))
let seed = try CartesianState(epoch: start, frame: inertial, position: SIMD3(7e6,0,0),
                               velocity: SIMD3(0,sqrt(gm.cubicMetresPerSecondSquared/7e6),0))
let coast = try TwoBodyPropagator(initialState: seed, centralBody: .earth, gravitationalParameter: gm)
let burn = try ImpulsiveManeuver(epoch: start.advanced(by: Duration(seconds: 300)), frame: inertial,
                                deltaVelocityMetresPerSecond: SIMD3(0,20,0))
let trajectory = try ManeuveredTwoBodyPropagator(base: coast, maneuvers: [burn])
let samples = try TrajectorySampler.sample(trajectory, from: start, through: end,
    step: Duration(seconds: 120), in: rotating, using: tree)
struct Sample: Encodable {
    let secondsFromStart: Double
    let earthBarycentricPositionMetres: [Double]
    let spacecraftEarthRotatingPositionMetres: [Double]
    let spacecraftEarthRotatingVelocityMetresPerSecond: [Double]
}
let output = try samples.map { sample -> Sample in
    let earthState = try planets.state(at: sample.epoch, in: CelestialFrames.icrf, using: tree)
    return Sample(secondsFromStart: try sample.epoch.duration(since: start).seconds,
        earthBarycentricPositionMetres: [earthState.position.x,earthState.position.y,earthState.position.z],
        spacecraftEarthRotatingPositionMetres: [sample.position.x,sample.position.y,sample.position.z],
        spacecraftEarthRotatingVelocityMetresPerSecond: [sample.velocity.x,sample.velocity.y,sample.velocity.z])
}
let encoder = JSONEncoder()
encoder.outputFormatting = [.prettyPrinted, .sortedKeys]
print(String(decoding: try encoder.encode(output), as: UTF8.self))

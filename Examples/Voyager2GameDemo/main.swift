// SPDX-License-Identifier: Apache-2.0
import Foundation
import CelMek
import Math

// An asset-free, bounded date search and patched-conic diagnostic, not a mission reconstruction.
func norm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
func cross(_ a: SIMD3<Double>,_ b: SIMD3<Double>) -> SIMD3<Double> {
    SIMD3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
}
func array(_ v: SIMD3<Double>) -> [Double] { [v.x,v.y,v.z] }
func utc(_ y: Int,_ m: Int,_ d: Int,_ seconds: Double = 43200) throws -> Epoch {
    try Epoch(utc: UTCDate(modifiedJulianDay: GregorianDay(year: y,month: m,day: d).modifiedJulianDay,
        secondOfDay: CelMek.Duration(seconds: seconds)),using: .iersBulletinC72)
}
func label(_ epoch: Epoch) throws -> String {
    let date = try epoch.utc(using: .iersBulletinC72), day = try GregorianDay(modifiedJulianDay: date.modifiedJulianDay)
    return String(format: "%04d-%02d-%02d UTC",day.year,day.month,day.day)
}
let provenance = try PhysicalDataProvenance(source: "Rounded Newtonian demo constants; spherical body model",version: "Voyager2-demo-1")
func gm(_ value: Double) throws -> GravitationalParameter {
    try GravitationalParameter(cubicMetresPerSecondSquared: value,convention: .newtonian,provenance: provenance)
}
let launch = try utc(1977,8,20,14*3600+29*60+44)
let anchors = try [utc(1979,7,9),utc(1981,8,25),utc(1986,1,24),utc(1989,8,25)]
let coverage = try launch.advanced(by: CelMek.Duration(seconds: -86400))...anchors.last!.advanced(by: CelMek.Duration(seconds: 400*86400))
let tdb = try FairheadTDBModel()
let coefficients = [earth,jupiter,saturn,uranus,neptune]
let ephemerides = try coefficients.map { try VSOP87Ephemeris(body: $0,coverage: coverage,tdbModel: tdb,accuracy: .full) }
let sunEphemeris = try VSOP87Ephemeris(body: sun,coverage: coverage,tdbModel: tdb,accuracy: .full)
let heliocentric = Frame(name: "Sun origin; VSOP87E native axes",pseudoInertial: true)
@MainActor
func bodyState(_ body: Int,_ epoch: Epoch) throws -> CartesianState {
    let planet = try ephemerides[body].checkedState(at: epoch), centre = try sunEphemeris.checkedState(at: epoch)
    guard planet.frame === centre.frame else { throw EphemerisError.frameMismatch }
    return try CartesianState(epoch: epoch,frame: heliocentric,position: planet.position-centre.position,velocity: planet.velocity-centre.velocity)
}
let solarGM = try gm(1.32712440018e20), earthGM = try gm(3.986004418e14)
let bodies: [BodyIdentifier] = [.jupiterBarycentre,.saturnBarycentre,.uranusBarycentre,.neptuneBarycentre]
let names = ["Jupiter","Saturn","Uranus","Neptune"]
let mus = try [1.26686534e17,3.7931187e16,5.793939e15,6.836529e15].map(gm)
let radii = [71_492e3,60_268e3,25_559e3,24_764e3]
let minAltitudes = [100_000e3,60_000e3,25_000e3,5_000e3] // Demo limits, not historical closest approaches.
let offsets = [-30.0,0,30]
let dates = try anchors.map { date in try offsets.map { try date.advanced(by: CelMek.Duration(seconds: $0*86400)) } }
let endpoints = try dates.enumerated().map { index,options in try options.map { try bodyState(index+1,$0) } }
let departure = try bodyState(0,launch)
struct LegKey: Hashable { let leg,from,to: Int }
var connectors: [LegKey:LambertSolution] = [:]
for leg in 0..<4 {
    for from in 0..<(leg == 0 ? 1 : 3) { for to in 0..<3 {
        let a = leg == 0 ? departure : endpoints[leg-1][from], b = endpoints[leg][to]
        connectors[LegKey(leg: leg,from: from,to: to)] = try LambertSolver.solve(
            departurePosition: a.position,arrivalPosition: b.position,departureEpoch: a.epoch,arrivalEpoch: b.epoch,
            frame: heliocentric,centralBody: .sun,gravitationalParameter: solarGM,
            referencePole: SIMD3(0,0,1),direction: .positive)
    } }
}
struct FlybyDiagnostic {
    let result: UnpoweredFlybyResult
    let residual: Double
    let requestedTurn: Double
    let speedMismatch: Double
}
func flyby(_ index: Int,_ epoch: Epoch,_ incoming: SIMD3<Double>,_ desired: SIMD3<Double>) throws -> FlybyDiagnostic {
    let speed = norm(incoming), desiredSpeed = norm(desired), u = incoming/speed, w = desired/desiredSpeed
    let angle = acos(max(-1,min(1,(u*w).sum())))
    let maximum = 2*asin(1/(1+(radii[index]+minAltitudes[index])*speed*speed/mus[index].cubicMetresPerSecondSquared))
    let turn = min(angle,maximum*(1-1e-10))
    let perpendicular = w-u*(u*w).sum(), length = norm(perpendicular)
    guard length > 1e-12, turn > 1e-12 else { throw FlybyError.degenerateTurn }
    let outgoing = (u*cos(turn)+perpendicular/length*sin(turn))*speed
    let frame = Frame(name: "\(names[index]) system origin; VSOP87E native axes",pseudoInertial: true)
    let request = try FlybyEncounterRequest(body: bodies[index],periapsisEpoch: epoch,frame: frame,
        gravitationalParameter: mus[index],bodyRadiusMetres: radii[index],minimumPeriapsisAltitudeMetres: minAltitudes[index],
        incomingExcessVelocityMetresPerSecond: incoming,outgoingExcessVelocityMetresPerSecond: outgoing)
    let result = try UnpoweredFlyby.solve(request)
    return FlybyDiagnostic(result: result,residual: norm(result.outgoingExcessVelocityMetresPerSecond-desired),
        requestedTurn: angle,speedMismatch: desiredSpeed-speed)
}
struct Candidate {
    let indices: [Int]
    let legs: [LambertSolution]
    let flybys: [FlybyDiagnostic]
    let escape: ParkingOrbitTransfer
    let score: Double
}
var best: Candidate?
for j in 0..<3 { for s in 0..<3 { for u in 0..<3 { for n in 0..<3 {
    let indices = [j,s,u,n]
    let legs = (0..<4).map { connectors[LegKey(leg: $0,from: $0 == 0 ? 0 : indices[$0-1],to: indices[$0])]! }
    var flybys: [FlybyDiagnostic] = []
    for i in 0..<3 {
        let centre = endpoints[i][indices[i]]
        flybys.append(try flyby(i,centre.epoch,legs[i].arrivalState.velocity-centre.velocity,legs[i+1].departureState.velocity-centre.velocity))
    }
    let escape = try ParkingOrbitTransfer(parkingRadiusMetres: 6_378_137+200_000,bodyRadiusMetres: 6_378_137,
        hyperbolicExcessSpeedMetresPerSecond: norm(legs[0].departureState.velocity-departure.velocity),gravitationalParameter: earthGM)
    let score = escape.deltaVelocityMetresPerSecond+flybys.reduce(0) { $0+$1.residual }
    if best == nil || score < best!.score { best = Candidate(indices: indices,legs: legs,flybys: flybys,escape: escape,score: score) }
} } } }
let selected = best!
let earthFrame = Frame(name: "Earth origin; VSOP87E native axes",pseudoInertial: true)
let excess = selected.legs[0].departureState.velocity-departure.velocity
let normal = cross(excess,SIMD3<Double>(0,0,1)); let unitNormal = normal/norm(normal)
let injection = try selected.escape.plan(kind: .escape,periapsisEpoch: launch,frame: earthFrame,centralBody: .earth,
    hyperbolicExcessVelocityMetresPerSecond: excess,orbitalNormal: unitNormal)
// Neptune terminal local flyby targets a chosen out-of-plane departure; no fifth Lambert leg.
let lastCentre = endpoints[3][selected.indices[3]]
let lastIncoming = selected.legs[3].arrivalState.velocity-lastCentre.velocity
let away = lastIncoming/norm(lastIncoming)
let downward = SIMD3<Double>(0,0,-1)+away*away.z
let terminalTurn = 2*asin(1/(1+(radii[3]+minAltitudes[3])*norm(lastIncoming)*norm(lastIncoming)/mus[3].cubicMetresPerSecondSquared))*(1-1e-10)
let terminalDesired = (away*cos(terminalTurn)+downward/norm(downward)*sin(terminalTurn))*norm(lastIncoming)
let terminal = try flyby(3,lastCentre.epoch,lastIncoming,terminalDesired)
let allFlybys = selected.flybys+[terminal]

struct Sample: Encodable {
    let phase,utc: String
    let secondsFromLaunch: Double
    let heliocentricPositionMetres,heliocentricVelocityMetresPerSecond: [Double]
}
var samples: [Sample] = []
@MainActor
func append(_ state: CartesianState,_ phase: String) throws {
    samples.append(Sample(phase: phase,utc: try label(state.epoch),secondsFromLaunch: try state.epoch.duration(since: launch).seconds,
        heliocentricPositionMetres: array(state.position),heliocentricVelocityMetresPerSecond: array(state.velocity)))
}
// Initial orbit and escape are genuine Earth-centred propagated states, translated for plotting.
for (seconds,phase) in [(-5400.0,"parking orbit"),(-3600,"parking orbit"),(-1800,"parking orbit"),(0,"injection"),(3600,"Earth escape") ] {
    let e = try launch.advanced(by: CelMek.Duration(seconds: seconds)), local = try injection.trajectory.state(at: e), earth = try bodyState(0,e)
    try append(CartesianState(epoch: e,frame: heliocentric,position: earth.position+local.position,velocity: earth.velocity+local.velocity),phase)
}
struct LegReport: Encodable { let from,to,departureUTC,arrivalUTC: String; let days,endpointResidualMetres: Double }
var reports: [LegReport] = []
for i in 0..<4 {
    let leg = selected.legs[i], coast = try leg.trajectory(), duration = try leg.arrivalState.epoch.duration(since: leg.departureState.epoch).seconds
    let epochs = try (0...48).map { k in
        k == 48 ? leg.arrivalState.epoch : try leg.departureState.epoch.advanced(by: CelMek.Duration(seconds: duration*Double(k)/48))
    }
    for state in try TrajectorySampler.sample(coast,at: epochs) { try append(state,"Lambert coast \(i+1)") }
    reports.append(LegReport(from: i == 0 ? "Earth" : names[i-1],to: names[i],departureUTC: try label(leg.departureState.epoch),
        arrivalUTC: try label(leg.arrivalState.epoch),days: duration/86400,
        endpointResidualMetres: norm(try coast.state(at: leg.arrivalState.epoch).position-leg.arrivalState.position)))
}
struct FlybyReport: Encodable {
    let body,utc: String
    let incomingExcessKilometresPerSecond,requestedSpeedMismatchMetresPerSecond,connectorVelocityResidualMetresPerSecond: Double
    let requestedTurnDegrees,achievedTurnDegrees,periapsisAltitudeKilometres: Double
    let referenceDateOffsetDays: Double
}
let encounters = try allFlybys.enumerated().map { i,f in
    FlybyReport(body: names[i],utc: try label(f.result.request.periapsisEpoch),incomingExcessKilometresPerSecond: f.result.excessSpeedMetresPerSecond/1000,
        requestedSpeedMismatchMetresPerSecond: f.speedMismatch,connectorVelocityResidualMetresPerSecond: f.residual,
        requestedTurnDegrees: f.requestedTurn*180 / .pi,achievedTurnDegrees: f.result.turningAngleRadians*180 / .pi,
        periapsisAltitudeKilometres: f.result.periapsisAltitudeMetres/1000,
        referenceDateOffsetDays: try f.result.request.periapsisEpoch.duration(since: anchors[i]).seconds/86400)
}
// Local flyby samples are separate views; they are not silently substituted into Lambert legs.
for (i,f) in allFlybys.enumerated() {
    for seconds in [-21600.0,-10800,0,10800,21600] {
        let e = try f.result.request.periapsisEpoch.advanced(by: CelMek.Duration(seconds: seconds))
        let local = try f.result.trajectory.state(at: e), planet = try bodyState(i+1,e)
        try append(CartesianState(epoch: e,frame: heliocentric,position: planet.position+local.position,velocity: planet.velocity+local.velocity),"local \(names[i]) flyby")
    }
}
// Playback selects a bounded coast explicitly; exact patch epochs use the following leg.
let playbackCoasts = try selected.legs.map { try $0.trajectory() }
@MainActor
func spacecraft(at epoch: Epoch) throws -> CartesianState {
    guard epoch >= launch, epoch <= lastCentre.epoch else { throw LambertTransferPlanningError.outsideFlight }
    for i in 0..<4 where epoch < selected.legs[i].arrivalState.epoch || i == 3 {
        return try playbackCoasts[i].state(at: epoch)
    }
    throw LambertTransferPlanningError.outsideFlight
}
struct Playback: Encodable {
    let mode,utc: String
    let simulationSecondsFromLaunch,simulationSecondsPerRealSecond: Double
    let isPaused: Bool
    let spacecraftPositionMetres: [Double]
    let prescribedPlanetPositionsMetres: [[Double]]
    let spacecraftCentredPlanetPositionsRenderUnits: [[Float]]
    let predictionUTC: String?
    let predictedSpacecraftPositionMetres: [Double]?
}
var clock = try SimulationClock(epoch: launch,rate: 365*86400,coverage: launch...lastCentre.epoch)
var playback: [Playback] = []
@MainActor
func capture(_ mode: String) throws {
    let epoch = clock.epoch, remaining = try lastCentre.epoch.duration(since: epoch)
    let lookAhead = try CelMek.Duration(seconds: 30*86400)
    let future = remaining >= lookAhead ? try clock.prediction(after: lookAhead) : nil
    let renderAdapter = try RendererAdapter(origin: spacecraft(at: epoch), metresPerRenderUnit: 1e9)
    let renderPlanets = try (0..<5).map { index in
        let position = try renderAdapter.adapt(bodyState(index,epoch)).position
        return [position.x,position.y,position.z]
    }
    playback.append(Playback(mode: mode,utc: try label(epoch),simulationSecondsFromLaunch: try epoch.duration(since: launch).seconds,
        simulationSecondsPerRealSecond: clock.rate,isPaused: clock.isPaused,spacecraftPositionMetres: array(try spacecraft(at: epoch).position),
        prescribedPlanetPositionsMetres: try (0..<5).map { array(try bodyState($0,epoch).position) },
        spacecraftCentredPlanetPositionsRenderUnits: renderPlanets,
        predictionUTC: try future.map(label),predictedSpacecraftPositionMetres: try future.map { array(try spacecraft(at: $0).position) }))
}
try capture("start")
for _ in 0..<12 { try clock.advance(realDuration: CelMek.Duration(seconds: 1)); try capture("forward warp") }
clock.setPaused(true); try clock.advance(realDuration: CelMek.Duration(seconds: 100)); try capture("paused")
try clock.setRate(-365*86400); clock.setPaused(false)
try clock.advance(realDuration: CelMek.Duration(seconds: 1)); try capture("reverse warp")
try clock.seek(to: selected.legs[0].arrivalState.epoch); try capture("seek Jupiter patch")

struct Output: Encodable {
    let model,axes,nativeEphemerisAxes,referenceURL: String
    let flightDays,referenceFlightDays,flightTimeDifferenceDays,injectionDeltaVelocityMetresPerSecond,departureC3SquareKilometresPerSquareSecond: Double
    let searchedCandidates: Int
    let isContinuousMissionTrajectory: Bool
    let assumptions: [String]
    let legs: [LegReport]; let flybys: [FlybyReport]; let samples: [Sample]; let playback: [Playback]
}
let output = Output(model: "Voyager 2 inspired bounded-date Lambert / unpowered flyby diagnostic",axes: heliocentric.name,
    nativeEphemerisAxes: ephemerides[0].frame.name,referenceURL: "https://science.nasa.gov/mission/voyager/fact-sheet/",
    flightDays: try lastCentre.epoch.duration(since: launch).seconds/86400,
    referenceFlightDays: try anchors[3].duration(since: launch).seconds/86400,
    flightTimeDifferenceDays: try lastCentre.epoch.duration(since: anchors[3]).seconds/86400,
    injectionDeltaVelocityMetresPerSecond: selected.escape.deltaVelocityMetresPerSecond,
    departureC3SquareKilometresPerSquareSecond: selected.escape.characteristicEnergySquareMetresPerSquareSecond/1e6,
    searchedCandidates: 81,isContinuousMissionTrajectory: false,
    assumptions: ["Encounter reference dates use noon UTC, not exact closest-approach times.",
        "Search uses -30/0/+30 days around each historical date; timing agreement is constrained, not a free prediction.",
        "Full VSOP87E positions; system barycentres approximate outer planet centres for local spherical flybys.",
        "Zero-revolution prograde solar two-body Lambert legs; no planetary perturbations during coasts.",
        "200 km circular parking orbit is illustrative; no launch vehicle ascent or historical fuel model.",
        "Rank by injection delta-v plus connector velocity residuals. Residuals are unresolved asymptotic mismatches, not implemented burns.",
        "Independent local flyby arcs and Earth escape are not finitely patched to coasts; their positions and velocities can differ at transitions.",
        "Neptune departure is a chosen out-of-plane local flyby, not a reconstructed interstellar trajectory."],
    legs: reports,flybys: encounters,samples: samples,playback: playback)
let encoder = JSONEncoder(); encoder.outputFormatting = [.prettyPrinted,.sortedKeys]
print(String(decoding: try encoder.encode(output),as: UTF8.self))

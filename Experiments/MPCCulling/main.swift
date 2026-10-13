// SPDX-License-Identifier: Apache-2.0
import Foundation
import CelMek

func timed<T>(_ name: String, _ body: () throws -> T) rethrows -> T {
    let start = DispatchTime.now().uptimeNanoseconds
    let result = try body()
    print(String(format: "%@: %.6f s", name, Double(DispatchTime.now().uptimeNanoseconds - start) / 1e9))
    return result
}
let args = Array(CommandLine.arguments.dropFirst())
func argument(_ key: String, default fallback: String) -> String {
    guard let i = args.firstIndex(of: key), i + 1 < args.count else { return fallback }
    return args[i + 1]
}
do {
    try checkSyntheticCases()
    let root = URL(fileURLWithPath: argument("--data", default: "Data/MPC"))
    let limit: Int?
    if args.contains("--limit") {
        guard let n = Int(argument("--limit", default: "")), n > 0 else { throw ExperimentError.invalidQuery }; limit = n
    } else { limit = nil }
    guard let halfDays = Double(argument("--half-window-days", default: "1")),
          let cellAU = Double(argument("--cell-au", default: "0.25")),
          halfDays >= 0, halfDays.isFinite, cellAU > 0, cellAU.isFinite else { throw ExperimentError.invalidQuery }
    let records = try timed("Load/prepare asteroids") { try loadAsteroids(root.appendingPathComponent("MPCORB.DAT"), limit: limit) }
        + timed("Load/prepare comets") { try loadComets(root.appendingPathComponent("cometels.json")) }
    guard Set(records.map(\.id)).count == records.count else { throw ExperimentError.malformed("duplicate identities") }
    print("Records: \(records.count); missing element epochs: \(records.filter { $0.elementEpoch == nil }.count)")
    let epoch = try calendarEpoch(2026, 10, 13)
    let index = try timed("Build swept radial/spatial index (includes full evaluation)") {
        try CullingIndex(records: records, epoch: epoch, halfWindowSeconds: halfDays * 86400, cellWidth: cellAU * au)
    }
    print("Cells: \(index.buckets.count); half-window: \(halfDays) days; cell width: \(cellAU) AU")
    // Fixed locations in the index's heliocentric ecliptic axes; these are query probes, not planet states.
    let probes: [(String, SIMD3<Double>, Double)] = [
        ("inner", SIMD3(1, 0, 0) * au, 0.1 * au),
        ("belt", SIMD3(2.7, 0, 0) * au, 0.1 * au),
        ("Jupiter distance", SIMD3(5.2, 0, 0) * au, 0.25 * au),
        ("outer", SIMD3(40, 0, 0) * au, au)]
    for offset in [-halfDays, 0, halfDays] {
        let queryEpoch = try epoch.advanced(by: CelMek.Duration(seconds: offset * 86400))
        for (name, observer, radius) in probes {
            print("Query \(name), offset \(offset) days")
            let radial = try timed("  Radial selection") { try index.candidates(observer: observer, radius: radius, epoch: queryEpoch, spatial: false) }
            let spatial = try timed("  Radial + spatial selection") { try index.candidates(observer: observer, radius: radius, epoch: queryEpoch, spatial: true) }
            let hits = try timed("  Evaluate spatial candidates") { try matches(records, indices: spatial, observer: observer, radius: radius, epoch: queryEpoch) }
            let oracle = try timed("  Full evaluation oracle") { try matches(records, indices: Array(records.indices), observer: observer, radius: radius, epoch: queryEpoch) }
            guard hits == oracle, oracle.isSubset(of: Set(radial)) else { throw ExperimentError.malformed("culling missed bodies") }
            print("  Candidates radial=\(radial.count), spatial=\(spatial.count), exact=\(hits.count); oracle PASS")
        }
    }
    let outside = try epoch.advanced(by: CelMek.Duration(seconds: halfDays * 86400 + 1))
    do {
        _ = try index.candidates(observer: .zero, radius: au, epoch: outside, spatial: true)
        throw ExperimentError.malformed("stale index accepted")
    } catch ExperimentError.outsideIndexWindow { print("Stale-index rejection PASS") }
} catch {
    fputs("MPCCullingExperiment: \(error)\n", stderr)
    exit(1)
}

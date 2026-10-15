// SPDX-License-Identifier: Apache-2.0
import Foundation
import CelMek

func checkSyntheticCases() throws {
    guard try packedEpoch("K2669") == calendarEpoch(2026, 6, 9) else { throw ExperimentError.malformed("packed epoch check") }
    let records = try [0.0, 0.8, 0.999999, 1.0, 1.2].enumerated().map { i, e in
        try OrbitRecord(id: "synthetic:\(i)", epoch: .j2000, elementEpoch: .j2000,
                        q: au * 0.1, e: e, inclination: 0.3, node: 0.4, peri: 0.5, anomaly: 0)
    }
    let index = try CullingIndex(records: records, epoch: .j2000, halfWindowSeconds: 86400, cellWidth: 0.01 * au)
    for seconds in [-86400.0, -43200, 0, 43200, 86400] {
        let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: seconds))
        for record in records {
            let position = try record.orbit.state(at: epoch).position
            // Zero-radius queries exactly on a body exercise cell crossings and closed bounds.
            for radius in [0.0, 0.01 * au, au] {
                let oracle = try matches(records, indices: Array(records.indices), observer: position, radius: radius, epoch: epoch)
                for spatial in [false, true] {
                    let candidates = try index.candidates(observer: position, radius: radius, epoch: epoch, spatial: spatial)
                    guard oracle.isSubset(of: Set(candidates)) else { throw ExperimentError.malformed("synthetic false negative") }
                }
            }
        }
    }
    do {
        _ = try index.candidates(observer: .zero, radius: -1, epoch: .j2000, spatial: true)
        throw ExperimentError.malformed("negative radius accepted")
    } catch ExperimentError.invalidQuery {}
    print("Synthetic circular/eccentric/near-parabolic/parabolic/hyperbolic checks PASS")
}

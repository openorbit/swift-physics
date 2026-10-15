// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Math
@testable import CelMek

private struct DE441References: Decodable {
    struct Row: Decodable {
        let jd: Double
        let target, center, component: Int
        let expected: Double
    }
    let source, sourceURL, sourceSHA256, units: String
    let binaryCoverage: [Double]
    let recordDays, tolerance: Double
    let rows: [Row]
    static func load() throws -> Self {
        let url = try #require(Bundle.module.url(forResource: "de441-testpo", withExtension: "json", subdirectory: "Fixtures"))
        return try JSONDecoder().decode(Self.self, from: Data(contentsOf: url))
    }
}

@Suite struct DE441ReferenceFixtureTests {
    @Test func independentFixtureIncludesEveryComponentAndBoundaryCases() throws {
        let fixture = try DE441References.load()
        #expect(fixture.source.contains("JPL"))
        #expect(fixture.sourceSHA256 == "e60c08ced7741a07d4dfe89b124afb4b1648bd1e3d0eb68dd9a79bf6990bb9c2")
        #expect(fixture.recordDays == 32 && fixture.tolerance == 1e-12)
        #expect(fixture.binaryCoverage == [-3027215.5,7930192.5])
        for target in 1...14 {
            for component in 1...(target == 14 ? 4 : 6) {
                #expect(fixture.rows.contains { $0.target == target && $0.component == component })
            }
        }
        for target in [3,10,13,14] {
            for boundary in [false,true] {
                for derivative in [false,true] {
                    #expect(fixture.rows.contains {
                        $0.target == target && (($0.jd-fixture.binaryCoverage[0]).truncatingRemainder(dividingBy: 32) == 0) == boundary
                        && ($0.component > (target == 14 ? 2 : 3)) == derivative
                    })
                }
            }
        }
        #expect(fixture.rows.allSatisfy { $0.expected.isFinite && $0.jd >= fixture.binaryCoverage[0] && $0.jd <= fixture.binaryCoverage[1] })
    }
}

@Suite(.enabled(if: FileManager.default.fileExists(atPath:
    LocalDE441Assets.directory.appendingPathComponent("linux_m13000p17000.441").path),
    "Optional DE441 binary is not installed"))
struct DE441PinnedReferenceTests {
    @Test func nativeValuesMatchIndependentJPLFixtures() throws {
        let fixture = try DE441References.load()
        let reader = try DE441Ephemeris(url: LocalDE441Assets.directory.appendingPathComponent("linux_m13000p17000.441"))
        for row in fixture.rows {
            let date = try JulianDateComponents(value: row.jd, scale: .tdb)
            let actual: Double
            if row.target == 14 {
                let value = try reader.nutation(at: date)
                actual = [value.longitude,value.obliquity,value.longitudeRate,value.obliquityRate][row.component-1]
            } else {
                let value = try reader.vector(target: row.target, center: row.center, at: date)
                actual = row.component <= 3 ? value.position[row.component-1] : value.velocity[row.component-4]
            }
            #expect(abs(actual-row.expected) <= fixture.tolerance,
                "JPL JD \(row.jd), target \(row.target), centre \(row.center), component \(row.component)")
        }
    }
    @Test func SIAdapterMatchesIndependentRelativeReferences() throws {
        let fixture = try DE441References.load()
        let reader = try DE441Ephemeris(url: LocalDE441Assets.directory.appendingPathComponent("linux_m13000p17000.441"))
        let model = try FairheadTDBModel()
        var checked = 0
        for row in fixture.rows where [3,10,13].contains(row.target)
            && row.jd >= 2433282.5 && row.jd <= 2469807.5 {
            let epoch = try Epoch(julianDate: JulianDateComponents(value: row.jd,scale: .tdb), using: model)
            let provider = try DE441StateProvider(reader: reader, target: DE441Body(rawValue: row.target)!,
                center: DE441Body(rawValue: row.center)!,
                frame: row.center == 12 ? nil : Frame(name: "ICRF-aligned reference origin", pseudoInertial: true),
                coverage: epoch...epoch, tdbModel: model)
            let state = try provider.checkedState(at: epoch)
            let rate = try epoch.tdbCorrection(using: model).rate
            let au = reader.metadata.astronomicalUnitKilometres*1000
            let scale = row.component <= 3 ? au : au/86400*rate
            let actual = row.component <= 3 ? state.position[row.component-1] : state.velocity[row.component-4]
            #expect(abs(actual-row.expected*scale) <= fixture.tolerance*scale,
                "SI JPL JD \(row.jd), target \(row.target), component \(row.component)")
            checked += 1
        }
        #expect(checked >= 18)
    }
}

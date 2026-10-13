// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
#if canImport(CryptoKit)
import CryptoKit
#endif
@testable import CelMek

private struct SpiceOracleReferences: Decodable {
    struct SPK: Decodable { let body,centre: Int32; let t: Double; let state: [Double] }
    struct CK: Decodable { let instrument: Int32; let tick,tolerance,selected: Double; let matrix: [[Double]]; let rate: [Double]? }
    struct PCK: Decodable { let t: Double; let matrix: [[Double]] }
    struct UTC: Decodable { let label: String; let mjd: Int64; let seconds,tt: Double }
    struct SCLK: Decodable { let tick,tdb: Double }
    let schema: Int
    let toolkit,wrapper,basis: String
    let assets: [String:String]
    let spk: [SPK]; let ck: [CK]; let pck: [PCK]; let utc: [UTC]; let sclk: [SCLK]
}
@Suite struct SpiceOracleReferenceTests {
    private func url(_ name: String) throws -> URL {
        try #require(Bundle.module.url(forResource: name,withExtension: nil,subdirectory: "Fixtures/spice-oracle"))
    }
    private func references() throws -> SpiceOracleReferences {
        try JSONDecoder().decode(SpiceOracleReferences.self,from: Data(contentsOf: url("reference.json")))
    }
    private func epoch(_ seconds: Double) throws -> Epoch { Epoch(offsetFromJ2000: try CelMek.Duration(seconds: seconds)) }
    private func matrixMatches(_ rotation: Quaternion,_ expected: [[Double]],tolerance: Double = 2e-13) -> Bool {
        for column in 0..<3 {
            var axis = SIMD3<Double>.zero; axis[column] = 1
            let actual = rotation.matrix*axis
            for row in 0..<3 { if abs(actual[row]-expected[row][column]) > tolerance { return false } }
        }
        return true
    }
    @Test func pinnedAssetIntegrityAndOracleIdentity() throws {
        let refs = try references()
        #if canImport(CryptoKit)
        let referenceBytes = try Data(contentsOf: url("reference.json"))
        #expect(SHA256.hash(data: referenceBytes).map { String(format: "%02x",$0) }.joined() ==
            "666a60aa902cd7bfffb479f8cf9ecc136575b144e54b3c29caf74fbdc17850be")
        #endif
        #expect(refs.schema == 1 && refs.toolkit == "CSPICE_N0067" && refs.wrapper == "8.3.0")
        #expect(refs.assets.count == 6 && refs.basis.contains("binary oracle only"))
        for (name,hash) in refs.assets {
            let bytes = try Data(contentsOf: url(name))
            #expect(!bytes.isEmpty && hash.count == 64)
            #if canImport(CryptoKit)
            #expect(SHA256.hash(data: bytes).map { String(format: "%02x",$0) }.joined() == hash)
            #endif
        }
    }
    @Test func geometricSPKStatesAndSIAdapterAgainstOracle() throws {
        let refs = try references(), set = try SpiceKernelSet(); try set.load(url("states.bsp"))
        let snapshot = set.snapshot(), coverage = try epoch(-20)...epoch(20)
        let helper = SpiceSPKAdapterTests()
        for body: Int32 in [399,499] {
            let provider = try SpiceSPKStateProvider(snapshot: snapshot,body: .naif(body),
                frameBindings: helper.bindings(),coverage: coverage,tdbModel: SpiceSPKAdapterTests.Clock(coverage: coverage),dataProvenance: helper.provenance)
            for row in refs.spk where row.body == body {
                let native = try snapshot.spkState(targetID: body,centreID: row.centre,atSecondsTDB: row.t)
                let state = try provider.checkedState(at: epoch(row.t))
                #expect(native.contributions.count == 2)
                for i in 0..<3 {
                    #expect(abs(native.positionKilometres[i]-row.state[i]) < 2e-12)
                    #expect(abs(native.velocityKilometresPerTDBSecond[i]-row.state[i+3]) < 2e-13)
                    #expect(abs(state.position[i]-1000*row.state[i]) < 2e-9)
                    #expect(abs(state.velocity[i]-1000*row.state[i+3]) < 2e-10)
                }
            }
        }
    }
    @Test func allCKTypesAndOptionalRatesAgainstOracle() throws {
        let refs = try references(), set = try SpiceKernelSet(); try set.load(url("attitudes.bc"))
        let snapshot = set.snapshot()
        for row in refs.ck {
            let p = try snapshot.ckPointing(instrumentID: row.instrument,atEncodedTicks: row.tick,toleranceTicks: row.tolerance)
            #expect(p.encodedTicks == row.selected)
            #expect(matrixMatches(p.referenceToInstrument,row.matrix))
            if let expected = row.rate {
                let actual = try #require(p.angularVelocityReferenceRadiansPerSecond)
                for i in 0..<3 { #expect(abs(actual[i]-expected[i]) < 2e-13) }
            } else { #expect(p.angularVelocityReferenceRadiansPerSecond == nil) }
        }
        #expect(throws: SpiceCKError.coverageGap) { try snapshot.ckPointing(instrumentID: -77003,atEncodedTicks: 40) }
    }
    @Test func binaryPCKOrientationAgainstOracle() throws {
        let refs = try references(), set = try SpiceKernelSet(); try set.load(url("rotation.bpc"))
        let helper = SpiceSPKAdapterTests(), coverage = try epoch(-20)...epoch(20)
        let provider = try SpiceBinaryPCKOrientationProvider(snapshot: set.snapshot(),frameClassID: 3000,
            referenceFrame: CelestialFrames.icrf,frameBindings: helper.bindings(),coverage: coverage,
            tdbModel: SpiceSPKAdapterTests.Clock(coverage: coverage),provenance: helper.provenance)
        for row in refs.pck { #expect(matrixMatches(try provider.orientation(at: epoch(row.t)).rotation,row.matrix)) }
    }
    @Test func leapLabelsAndFractionalClockCorrelationAgainstOracle() throws {
        let refs = try references(), set = try SpiceKernelSet(); try set.load(url("time.tls")); try set.load(url("clock.tsc"))
        let table = try set.snapshot().leapSecondData().table(validUntilModifiedJulianDay: 42000,version: "oracle-fixture")
        for row in refs.utc {
            let label = try UTCDate(modifiedJulianDay: row.mjd,secondOfDay: CelMek.Duration(seconds: row.seconds))
            let instant = try Epoch(utc: label,using: table)
            #expect(abs(instant.offsetFromJ2000.seconds-row.tt) <= 3e-7)
            #expect(try instant.utc(using: table) == label)
        }
        let clock = try SpiceSCLKClock(snapshot: set.snapshot(),spacecraftID: -77,encodedTickCoverage: 0...2000,
            provenance: SpiceSPKAdapterTests().provenance)
        for row in refs.sclk { #expect(abs(try clock.sample(atEncodedTicks: row.tick).parallelSeconds-row.tdb) < 1e-12) }
    }
}

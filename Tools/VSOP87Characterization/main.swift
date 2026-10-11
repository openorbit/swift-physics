// SPDX-License-Identifier: Apache-2.0
import Foundation
import CelMek
import Math

// Reproducible public-API characterization, not a CI wall-clock performance test.
struct Measurement: Codable {
  let body: String, mode: String, effectiveMode: String
  let originalTerms: Int, retainedTerms: Int, sampleCount: Int
  let maximumPositionErrorMetres: Double, rmsPositionErrorMetres: Double
  let maximumVelocityErrorMetresPerSecond: Double, rmsVelocityErrorMetresPerSecond: Double
  let maximumPositionErrorTTJulianDate: Double, maximumVelocityErrorTTJulianDate: Double
  let perAxisPositionBoundMetres: [Double]
  let preparationSeconds: Double, medianQueryNanoseconds: Double, speedupAgainstFull: Double
}
struct Report: Codable {
  let schemaVersion: Int
  let operatingSystem: String
  let processors: Int
  let coverageStartTTJulianDate: Double, coverageEndTTJulianDate: Double, samplingStepDays: Int
  let tdbModel: String, timingDescription: String
  let timingRepetitions: Int, queriesPerRepetition: Int
  let checksum: Double
  let measurements: [Measurement]
}
struct MutableMeasurement {
  let body: String, mode: AccuracyLevel, provider: VSOP87Ephemeris, preparationSeconds: Double
  var positionMaximum = 0.0, positionSumSquares = 0.0, positionEpoch = 0.0
  var velocityMaximum = 0.0, velocitySumSquares = 0.0, velocityEpoch = 0.0
  var count = 0, times = [Double]()
}
func magnitude(_ value: SIMD3<Double>) -> Double { sqrt((value*value).sum()) }
func now() -> Double { Double(DispatchTime.now().uptimeNanoseconds)*1e-9 }

let arguments = CommandLine.arguments
let output = arguments.count > 1 ? arguments[1] : "/tmp/celmek-vsop-characterization.json"
let stepDays = arguments.count > 2 ? Int(arguments[2]) ?? 0 : 8
guard stepDays > 0 else { fatalError("Sampling step must be a positive integer number of days") }
let tdb = try FairheadTDBModel()
let coverage = FairheadTDBModel.referenceCoverage
let bodies: [(String,vsop87_body)] = [("Sun",sun),("Mercury",mercury),("Venus",venus),("Earth",earth),
  ("Mars system",mars),("Jupiter system",jupiter),("Saturn system",saturn),("Uranus system",uranus),("Neptune system",neptune)]
var rows = [MutableMeasurement]()
for (name,body) in bodies {
  for mode in AccuracyLevel.allCases {
    let start = now()
    let provider = try VSOP87Ephemeris(body: body,coverage: coverage,tdbModel: tdb,accuracy: mode)
    rows.append(MutableMeasurement(body: name,mode: mode,provider: provider,preparationSeconds: now()-start))
  }
}
let span = try coverage.upperBound.duration(since: coverage.lowerBound).seconds
var epochs = [Epoch]()
var elapsed = 0.0
while elapsed < span {
  epochs.append(try coverage.lowerBound.advanced(by: CelMek.Duration(seconds: elapsed)))
  elapsed += Double(stepDays)*86400
}
epochs.append(coverage.upperBound)
for (number,epoch) in epochs.enumerated() {
  for bodyIndex in bodies.indices {
    let base = bodyIndex*4, full = try rows[base+3].provider.state(at: epoch)
    for index in base..<(base+4) {
      let actual = index == base+3 ? full : try rows[index].provider.state(at: epoch)
      let positionDifference = actual.position-full.position
      let p = magnitude(positionDifference), v = magnitude(actual.velocity-full.velocity)
      let bounds = rows[index].provider.coefficientSelection.positionErrorBoundMetres
      for axis in 0..<3 {
        guard abs(positionDifference[axis]) <= bounds[axis]+0.1 else {
          fatalError("Position bound exceeded for \(rows[index].body) \(rows[index].mode) at \(epoch)")
        }
      }
      rows[index].positionSumSquares += p*p; rows[index].velocitySumSquares += v*v; rows[index].count += 1
      if p > rows[index].positionMaximum {
        rows[index].positionMaximum = p; rows[index].positionEpoch = try epoch.julianDate(in: .tt).value
      }
      if v > rows[index].velocityMaximum {
        rows[index].velocityMaximum = v; rows[index].velocityEpoch = try epoch.julianDate(in: .tt).value
      }
    }
  }
  if number % 500 == 0 { FileHandle.standardError.write(Data("Compared \(number)/\(epochs.count) epochs\n".utf8)) }
}
let repetitions = 5, queries = 64
var checksum = 0.0
// Warm all modes before timing. Cycle evaluation order to reduce ordering bias.
for row in rows { checksum += try row.provider.state(at: .j2000).position.x }
for repetition in 0..<repetitions {
  for offset in rows.indices {
    let index = (offset+repetition*7)%rows.count
    let provider = rows[index].provider
    let start = now()
    for query in 0..<queries {
      let epoch = epochs[(query*71+repetition*131)%epochs.count]
      let state = try provider.state(at: epoch)
      checksum += state.position.x + state.velocity.y
    }
    rows[index].times.append((now()-start)*1e9/Double(queries))
  }
}
let medianTimes = rows.map { $0.times.sorted()[repetitions/2] }
let measurements = rows.enumerated().map { index,row in
  let selection = row.provider.coefficientSelection
  return Measurement(body: row.body,mode: row.mode.rawValue,effectiveMode: row.provider.effectiveAccuracy.rawValue,
    originalTerms: selection.originalTermCount,retainedTerms: selection.retainedTermCount,sampleCount: row.count,
    maximumPositionErrorMetres: row.positionMaximum,rmsPositionErrorMetres: sqrt(row.positionSumSquares/Double(row.count)),
    maximumVelocityErrorMetresPerSecond: row.velocityMaximum,rmsVelocityErrorMetresPerSecond: sqrt(row.velocitySumSquares/Double(row.count)),
    maximumPositionErrorTTJulianDate: row.positionEpoch,maximumVelocityErrorTTJulianDate: row.velocityEpoch,
    perAxisPositionBoundMetres: [selection.positionErrorBoundMetres.x,selection.positionErrorBoundMetres.y,selection.positionErrorBoundMetres.z],
    preparationSeconds: row.preparationSeconds,medianQueryNanoseconds: medianTimes[index],
    speedupAgainstFull: medianTimes[(index/4)*4+3]/medianTimes[index])
}
let report = Report(schemaVersion: 1,operatingSystem: ProcessInfo.processInfo.operatingSystemVersionString,
  processors: ProcessInfo.processInfo.processorCount,coverageStartTTJulianDate: try coverage.lowerBound.julianDate(in: .tt).value,
  coverageEndTTJulianDate: try coverage.upperBound.julianDate(in: .tt).value,samplingStepDays: stepDays,
  tdbModel: tdb.provenance,timingDescription: "Release public state query; includes Fairhead TDB conversion, allocations and SI conversion; warm providers; median of rotated-order batches; excludes provider preparation",
  timingRepetitions: repetitions,queriesPerRepetition: queries,checksum: checksum,measurements: measurements)
let encoder = JSONEncoder(); encoder.outputFormatting = [.prettyPrinted,.sortedKeys]
try encoder.encode(report).write(to: URL(fileURLWithPath: output))
print("Saved \(measurements.count) measurements at \(epochs.count) epochs to \(output)")

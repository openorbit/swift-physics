// Fairhead evaluator adapted from ERFA 2.0.1: Copyright 2013-2023 NumFOCUS.
// BSD terms and SOFA heritage: Documentation/licenses/ERFA.txt.
import Foundation

public enum TDBModelError: Error, Equatable, Sendable {
  case outsideCoverage, invalidCorrection, epochMismatch, noConvergence
}

/// Offset and analytic coordinate-time derivatives at a TT/atomic instant.
/// rate is dTDB/dTT (approximately 1); acceleration is d²TDB/dTT² in 1/s.
public struct TDBCorrection: Sendable {
  public let epoch: Epoch
  public let tdbMinusTT: Duration
  public let rate: Double
  public let acceleration: Double
  public init(epoch: Epoch, tdbMinusTT: Duration, rate: Double, acceleration: Double) throws {
    guard rate.isFinite, rate > 0, acceleration.isFinite else { throw TDBModelError.invalidCorrection }
    self.epoch = epoch; self.tdbMinusTT = tdbMinusTT; self.rate = rate; self.acceleration = acceleration
  }
}

/// Injectable geocentric TT↔TDB model. Coverage is in atomic Epochs, not TDB labels.
/// Implementations own accuracy and provenance; sample derivatives use atomic seconds.
public protocol TDBTimeModel {
  var coverage: ClosedRange<Epoch> { get }
  var provenance: String { get }
  func correction(at epoch: Epoch) throws -> TDBCorrection
}

extension Epoch {
  /// Validates the provider contract before using offsets or derivatives.
  public func tdbCorrection(using model: any TDBTimeModel) throws -> TDBCorrection {
    guard model.coverage.contains(self) else { throw TDBModelError.outsideCoverage }
    let sample = try model.correction(at: self)
    guard sample.epoch == self else { throw TDBModelError.epochMismatch }
    return sample
  }
  internal func tdbDate(using sample: TDBCorrection) throws -> JulianDateComponents {
    guard sample.epoch == self else { throw TDBModelError.epochMismatch }
    let shifted = Epoch(offsetFromJ2000: try offsetFromJ2000.adding(sample.tdbMinusTT))
    let date = try shifted.julianDate(in: .tt)
    return try JulianDateComponents(wholeDay: date.wholeDay,fractionalDay: date.fractionalDay,scale: .tdb)
  }
  internal static func invertTDB(_ date: JulianDateComponents, using model: any TDBTimeModel) throws -> Epoch {
    let clock = try Epoch(julianDate: JulianDateComponents(wholeDay: date.wholeDay,
      fractionalDay: date.fractionalDay,scale: .tt))
    // TDB may label a covered TT boundary just outside the TT numeric interval.
    var estimate = min(max(clock,model.coverage.lowerBound),model.coverage.upperBound)
    for _ in 0..<12 {
      let sample = try estimate.tdbCorrection(using: model)
      // Newton iteration uses the model's derivative, not an assumed unit rate.
      let predicted = try estimate.offsetFromJ2000.adding(sample.tdbMinusTT)
      let residual = try predicted.subtracting(clock.offsetFromJ2000)
      let boundary = estimate == model.coverage.lowerBound || estimate == model.coverage.upperBound
      if abs(residual.seconds) <= (boundary ? 3e-11 : 1e-12) { return estimate }
      let next = try estimate.advanced(by: Duration(seconds: -residual.seconds/sample.rate))
      if !model.coverage.contains(next) {
        let boundary = min(max(next,model.coverage.lowerBound),model.coverage.upperBound)
        let correction = try boundary.tdbCorrection(using: model)
        let clockAtBoundary = try boundary.offsetFromJ2000.adding(correction.tdbMinusTT)
        // Split JD conversion rounds; accept only a label within 30 ps of a
        // covered boundary, returning that boundary rather than extrapolating.
        if abs(try clockAtBoundary.subtracting(clock.offsetFromJ2000).seconds) <= 3e-11 { return boundary }
        throw TDBModelError.outsideCoverage
      }
      estimate = next
    }
    throw TDBModelError.noConvergence
  }
}

/// The existing two-term approximation with caller-selected, enforced coverage.
/// The selected interval does not create an accuracy guarantee; see TIME.md.
public struct TwoTermTDBModel: TDBTimeModel {
  public let coverage: ClosedRange<Epoch>
  public let provenance = "CelMek two-term geocentric TDB approximation"
  public init(coverage: ClosedRange<Epoch>) { self.coverage = coverage }
  public func correction(at epoch: Epoch) throws -> TDBCorrection {
    guard coverage.contains(epoch) else { throw TDBModelError.outsideCoverage }
    return try TDBCorrection(epoch: epoch,tdbMinusTT: Duration(seconds: Epoch.tdbMinusTT(at: epoch.offsetFromJ2000)),
      rate: epoch.coordinateTimeRate(in: .tdb),acceleration: epoch.coordinateTimeAcceleration(in: .tdb))
  }
}

/// Full 787-term geocentric Fairhead–Bretagnon model with ERFA's planetary-mass
/// adjustments. Native Swift; adapted from ERFA 2.0.1 under its BSD terms.
/// No observer-dependent topocentric terms. TT is used as the series argument,
/// as permitted by the reference implementation. Derivatives are analytic.
public struct FairheadTDBModel: TDBTimeModel {
  public static let referenceCoverage: ClosedRange<Epoch> = {
    let start = try! Epoch(julianDate: 2433282.5,scale: .tt) // 1950-01-01
    let end = try! Epoch(julianDate: 2469807.5,scale: .tt) // 2050-01-01
    return start...end
  }()
  public let coverage: ClosedRange<Epoch>
  public let provenance = "Fairhead–Bretagnon 1990, full geocentric ERFA 2.0.1 adaptation"
  public init(coverage: ClosedRange<Epoch> = Self.referenceCoverage) throws {
    guard coverage.lowerBound >= Self.referenceCoverage.lowerBound,
          coverage.upperBound <= Self.referenceCoverage.upperBound else { throw TDBModelError.outsideCoverage }
    self.coverage = coverage
  }
  public func correction(at epoch: Epoch) throws -> TDBCorrection {
    guard coverage.contains(epoch) else { throw TDBModelError.outsideCoverage }
    let tt = epoch.offsetFromJ2000
    let secondsPerMillennium = 365250.0*86400
    let t = Double(tt.wholeSeconds)/secondsPerMillennium + tt.fractionalSecond/secondsPerMillennium
    // Store value and first/second derivatives with respect to Julian millennia.
    var value = 0.0, first = 0.0, second = 0.0
    let ranges = [0..<474,474..<679,679..<764,764..<784,784..<787]
    for order in (0..<5).reversed() {
      var v = 0.0, d = 0.0, dd = 0.0
      for index in ranges[order].reversed() {
        let row = FairheadTDBSeries.terms[index]
        let amplitude = row[0], frequency = row[1], phase = frequency*t+row[2]
        let harmonic = amplitude*sin(phase)
        v += harmonic; d += amplitude*frequency*cos(phase); dd -= harmonic*frequency*frequency
      }
      // Differentiate Horner's rule, retaining the previous derivatives.
      second = second*t + 2*first + dd
      first = first*t + value + d
      value = value*t + v
    }
    for row in [[0.00065e-6,6069.776754,4.021194],[0.00033e-6,213.299095,5.543132],
                [-0.00196e-6,6208.294251,5.696701],[-0.00173e-6,74.781599,2.435900]] {
      let harmonic = row[0]*sin(row[1]*t+row[2])
      value += harmonic
      first += row[0]*row[1]*cos(row[1]*t+row[2])
      second -= harmonic*row[1]*row[1]
    }
    value += 0.03638e-6*t*t
    first += 2*0.03638e-6*t
    second += 2*0.03638e-6
    return try TDBCorrection(epoch: epoch,tdbMinusTT: Duration(seconds: value),
      rate: 1+first/secondsPerMillennium,acceleration: second/(secondsPerMillennium*secondsPerMillennium))
  }
}

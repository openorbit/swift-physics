// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Coefficient subset prepared once, outside queries. Bounds describe truncation
/// relative to the complete built-in series, not physical ephemeris accuracy.
public struct VSOP87CoefficientSelection {
  public let originalTermCount: Int
  public let retainedTermCount: Int
  /// Sum of absolute omitted amplitudes times the largest |t| power in coverage,
  /// per native Cartesian axis, in metres. No velocity error bound is implied.
  public let positionErrorBoundMetres: SIMD3<Double>
  /// Individual term contribution cutoff in AU over the selected interval.
  public let termCutoffAU: Double
  internal let selectedBody: vsop87_body

  internal init(body: vsop87_body, coverage: ClosedRange<Epoch>, tdbModel: any TDBTimeModel,
                accuracy: AccuracyLevel) throws {
    switch accuracy {
    case .low: termCutoffAU = 1e-6
    case .balanced: termCutoffAU = 1e-8
    case .high: termCutoffAU = 1e-10
    case .full: termCutoffAU = 0
    }
    func argument(_ epoch: Epoch) throws -> Double {
      let date = try epoch.julianDate(in: .tdb,using: tdbModel)
      return (Double(date.wholeDay-2451545)+date.fractionalDay)/365250
    }
    let maximumT = try max(abs(argument(coverage.lowerBound)),abs(argument(coverage.upperBound))) + 1e-12
    var count = 0, retained = 0, bounds = SIMD3<Double>.zero
    let cutoff = termCutoffAU
    func select(_ series: [[(Double,Double,Double)]], axis: Int) throws -> [[(Double,Double,Double)]] {
      var power = 1.0, output = [[(Double,Double,Double)]]()
      for terms in series {
        var selected = [(Double,Double,Double)]()
        for term in terms {
          count += 1
          let bound = abs(term.0)*power
          guard bound.isFinite else { throw TimeArithmeticError.outOfRange }
          if accuracy == .full || bound >= cutoff {
            selected.append(term); retained += 1
          } else {
            // Round the bound upward at each addition to cover floating summation.
            bounds[axis] = (bounds[axis]+bound).nextUp
          }
        }
        output.append(selected); power *= maximumT
      }
      return output
    }
    let x = try select(body.terms_x,axis: 0), y = try select(body.terms_y,axis: 1), z = try select(body.terms_z,axis: 2)
    originalTermCount = count; retainedTermCount = retained
    positionErrorBoundMetres = bounds * CartesianUnits.metresPerAstronomicalUnit
    selectedBody = accuracy == .full ? body : vsop87_body(body: body.body_id,tx: x,ty: y,tz: z)
  }
}

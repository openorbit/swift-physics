// Fairhead evaluator adapted from ERFA 2.0.1: Copyright 2013-2023 NumFOCUS.
// BSD terms and SOFA heritage: Sources/ERFASwift/LICENSE.txt.
import Foundation

/// Geocentric offset and analytic derivatives from the full Fairhead series.
/// Input is split TT seconds from J2000; the output offset is in SI seconds.
package enum FairheadTDBEvaluator {
  package static func correction(wholeSeconds: Int64, fractionalSecond: Double)
    -> (offsetSeconds: Double, rate: Double, acceleration: Double) {
    let secondsPerMillennium = 365250.0*86400
    let t = Double(wholeSeconds)/secondsPerMillennium + fractionalSecond/secondsPerMillennium
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
    return (value, 1+first/secondsPerMillennium,
            second/(secondsPerMillennium*secondsPerMillennium))
  }
}

// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Math
@testable import CelMek

@Suite
struct FastNutationAccuracyTests {
  @Test
  func dailyComparisonWithFullModelFrom1900Through2100() {
    // Empirical regression budgets on this sampled interval, not universal
    // accuracy guarantees or a comparison with modern observed orientation.
    let start = GregorianDate(year: 1900, month: .January, day: 1).toJD()
    let end = GregorianDate(year: 2100, month: .December, day: 31).toJD()
    var maximum = SIMD3<Double>.zero
    var worstEpoch = SIMD3<Double>.zero
    var squareSum = SIMD3<Double>.zero
    var count = 0
    for jd in stride(from: start, through: end, by: 1) {
      let fast = fastNutation(jd: jd)
      let full = nutation(jd: jd)
      let differences = SIMD3(abs(fast.nutationInLongitude - full.nutationInLongitude).asArcSec,
                              abs(fast.nutationInObliquity - full.nutationInObliquity).asArcSec,
                              abs(fastTrueObliquityOfTheEcliptic(jd: jd)
                                  - trueObliquityOfTheEcliptic(jd: jd)).asArcSec)
      squareSum += differences * differences
      for component in 0..<3 {
        if differences[component] > maximum[component] {
          maximum[component] = differences[component]
          worstEpoch[component] = jd
        }
      }
      count += 1
    }
    #expect(count == 73_414)
    let maximumBudgets = SIMD3<Double>(0.35, 0.1, 0.1)
    let rmsBudgets = SIMD3<Double>(0.12, 0.03, 0.03)
    for component in 0..<3 {
      #expect(maximum[component].isFinite && maximum[component] < maximumBudgets[component],
              "Component \(component), worst JD \(worstEpoch[component])")
      #expect(sqrt(squareSum[component] / Double(count)) < rmsBudgets[component])
    }
  }

  @Test(arguments: [2415020.75, 2446895.625, 2451545.0, 2465910.125, 2481630.25, 2488433.9])
  func fractionalEpochsStayWithinApproximationBudget(jd: Double) {
    let fast = fastNutation(jd: jd)
    let full = nutation(jd: jd)
    #expect(abs(fast.nutationInLongitude - full.nutationInLongitude) < 0.35.arcsec)
    #expect(abs(fast.nutationInObliquity - full.nutationInObliquity) < 0.1.arcsec)
    #expect(abs(fastTrueObliquityOfTheEcliptic(jd: jd)
                - trueObliquityOfTheEcliptic(jd: jd)) < 0.1.arcsec)
  }
}

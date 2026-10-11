// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
@testable import CelMek

@Suite
struct VSOP87DerivativeTests {
  @Test
  func independentLinearComponents() {
    let body = vsop87_body(body: .vsop87_earth,
                           tx: [[], [(2, 0, 0)]],
                           ty: [[], [(5, 0, 0)]],
                           tz: [[], [(-3, 0, 0)]])
    let (_, velocity, _, _) = body.position(jd: 2451545)
    for (component, slope) in [2.0, 5, -3].enumerated() {
      #expect(abs(velocity[component] - slope / 365250) < 1e-18)
    }
  }

  @Test(arguments: [-0.5, 0.0, 0.25], [1, 2, 3, 4, 5])
  func polynomialAndPeriodicDerivatives(t: Double, order: Int) {
    // p(t) = 3tⁿ cos(0.4 + 2t), independently differentiated below.
    var terms = [[(Double, Double, Double)]](repeating: [], count: order + 1)
    terms[order] = [(3, 0.4, 2)]
    let body = vsop87_body(body: .vsop87_earth, tx: [], ty: terms, tz: terms)
    let jd = 2451545 + t * 365250
    let (position, velocity, _, _) = body.position(jd: jd)
    let argument = 0.4 + 2 * t
    let expectedPosition = 3 * pow(t, Double(order)) * cos(argument)
    let expectedVelocity = (3 * Double(order) * pow(t, Double(order - 1)) * cos(argument)
                            - 6 * pow(t, Double(order)) * sin(argument)) / 365250
    #expect(position.x == 0)
    #expect(velocity.x == 0)
    for component in 1..<3 {
      #expect(abs(position[component] - expectedPosition) < 1e-14)
      #expect(abs(velocity[component] - expectedVelocity) < 1e-18)
    }
  }

  @Test(arguments: [2451545.0, 2415020.0, 2488070.0])
  func allBodyVelocitiesMatchPositionDifferences(jd: Double) {
    let bodies = [sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune]
    let step = 0.001
    // Use the actual representable JD interval in the denominator.
    let beforeDate = jd - step
    let afterDate = jd + step
    for body in bodies {
      let (before, _, _, _) = body.position(jd: beforeDate)
      let (after, _, _, _) = body.position(jd: afterDate)
      let (_, velocity, _, identifier) = body.position(jd: jd)
      let difference = (after - before) / (afterDate - beforeDate)
      for component in 0..<3 {
        #expect(abs(velocity[component] - difference[component]) < 1e-9,
                "Body \(identifier), component \(component), JD \(jd)")
      }
    }
  }
}

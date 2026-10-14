// SPDX-License-Identifier: Apache-2.0
import Foundation
import Aerodynamics
import Atmospheres

// Synthetic coefficients demonstrate evaluator contracts, not validated vehicle aerodynamics.
let sectionSamples = [-0.2,0.0,0.2].map {
    AirfoilSample(angleOfAttackRadians:$0,coefficients:AirfoilCoefficients(lift:5*$0,drag:0.02,moment:0))
}
var polars = [AirfoilPolar]()
for re in [100_000.0,10_000_000.0] {
    for mach in [0.0,0.8] { polars.append(try AirfoilPolar(reynoldsNumber:re,machNumber:mach,samples:sectionSamples)) }
}
let airfoil = try TabulatedAirfoil(polars:polars,source:"Synthetic demo polars",version:"1")
let wing = try SegmentedWing.rectangular(spanMetres:8,chordMetres:2,incidenceRadians:0.08,
    model:AirfoilWingSectionModel(airfoil:airfoil))
let wingFlow = try AerodynamicFlow(relativeVelocityMetresPerSecond:SIMD3(40,0,0),properties:ISAAtmosphere().properties(at:0))
let wingResult = try wing.evaluate(flow:wingFlow)
print("Wing: \(wingResult.segments.count) strips; force N = \(wingResult.loads.forceNewtons); moment N m = \(wingResult.loads.momentNewtonMetres)")

var bodySamples = [BodyCoefficientSample]()
for alpha in [-0.2,0.0,0.2] { for beta in [-0.1,0.0,0.1] { for mach in [0.0,25.0] { for re in [100_000.0,1_000_000_000.0] {
    let wind = try WindForceCoefficients(drag:1.5,lift:0.5*alpha)
    bodySamples.append(BodyCoefficientSample(angleOfAttackRadians:alpha,sideslipRadians:beta,machNumber:mach,reynoldsNumber:re,
        coefficients:try wind.bodyCoefficients(angleOfAttackRadians:alpha,sideslipRadians:beta,momentBody:SIMD3(0,-0.2*alpha,0))))
} } } }
let capsuleData = try TabulatedBodyCoefficients(samples:bodySamples,source:"Synthetic lifting capsule demo",version:"1")
let geometry = try AerodynamicReferenceGeometry(areaSquareMetres:10,momentLengthsMetres:SIMD3(repeating:3))
let capsule = try WholeBodyAerodynamics(geometry:geometry,reynoldsReferenceLengthMetres:3,model:capsuleData)
let capsuleFlow = try AerodynamicFlow(relativeVelocityMetresPerSecond:SIMD3(6000,0,600),properties:ISAAtmosphere().properties(at:50_000))
let capsuleLoads = try capsule.evaluate(flow:capsuleFlow)
print("Capsule: Mach = \(capsuleFlow.machNumber!); force N = \(capsuleLoads.forceNewtons); moment N m = \(capsuleLoads.momentNewtonMetres)")

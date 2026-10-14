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
let finiteWing = try FiniteWingAirfoil(airfoil:airfoil,aspectRatio:4,spanEfficiency:0.8,
    effectiveAngleRangeRadians:-0.2...0.2,operatingAngleRangeRadians:-0.2...0.2)
let wing = try SegmentedWing.rectangular(spanMetres:8,chordMetres:2,incidenceRadians:0.08,
    model:AirfoilWingSectionModel(airfoil:finiteWing))
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

// Independent NASA TFAWS reference inputs; mean free path and coverage are supplied
// demonstration assumptions, not an upper-atmosphere/transport model.
let heatingCoverage = try StagnationHeatingCoverage(speedMetresPerSecond: 1000...8000,
    densityKilogramsPerCubicMetre: 1e-6...0.01, noseRadiusMetres: 0.1...5, machNumber: 5...40)
let heating = SuttonGravesHeating(parameters: .earthAirColdWall, coverage: heatingCoverage)
let heatingFlow = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(3535,0,0),
    densityKilogramsPerCubicMetre: 3.1459e-4, speedOfSoundMetresPerSecond: 300, meanFreePathMetres: 1e-4)
let heat = try heating.heatFlux(for: StagnationHeatingQuery(flow: heatingFlow,
    gasIdentifier: "earth-air", noseRadiusMetres: 1, surfacePositionMetres: SIMD3(1,0,0)))
print("NASA reference cold-wall convective heat flux: \(heat.wattsPerSquareMetre) W/m²")

let marsHeating = SuttonGravesHeating(parameters: .marsCarbonDioxideNitrogenColdWall,
                                      coverage: heatingCoverage)
let marsHeat = try marsHeating.heatFlux(for: StagnationHeatingQuery(flow: heatingFlow,
    gasIdentifier: marsHeating.parameters.gasIdentifier, noseRadiusMetres: 1, surfacePositionMetres: SIMD3(1,0,0)))
print("Mars mixture, same supplied flow: \(marsHeat.wattsPerSquareMetre) W/m²")

let venusParameters = SuttonGravesMixtureParameters.venusCarbonDioxideNitrogen
let venusHeating = SuttonGravesMixtureHeating(parameters: venusParameters, coverage: heatingCoverage)
let venusHeat = try venusHeating.heatFlux(for: StagnationHeatingQuery(flow: heatingFlow,
    gasIdentifier: venusParameters.gasIdentifier, noseRadiusMetres: 1, surfacePositionMetres: SIMD3(1,0,0),
    gasState: StagnationGasState(pressurePascals: 101325, edgeSpecificEnthalpyJoulesPerKilogram: 10e6)))
print("Venus mixture, supplied shock-layer state: \(venusHeat.wattsPerSquareMetre) W/m²")

let giantMixture = HydrogenHeliumHeatingMixture.neptuneReference
let giantHeating = HydrogenHeliumColdWallHeating(mixture: giantMixture)
let giantFlow = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(19000,0,0),
    densityKilogramsPerCubicMetre: 1e-4, temperatureKelvin: 141.2,
    speedOfSoundMetresPerSecond: 1000, meanFreePathMetres: 1e-5)
let giantResult = try giantHeating.evaluate(StagnationHeatingQuery(flow: giantFlow,
    gasIdentifier: giantMixture.gasIdentifier, noseRadiusMetres: 0.4, surfacePositionMetres: SIMD3(1,0,0),
    wall: AerothermalWallState(temperatureKelvin: 300)))
print("H₂/He reference mixture: \(giantResult.heatFlux.wattsPerSquareMetre) W/m²")

// Supplied-state parachute and balloon examples; coefficients are synthetic.
func demoDragTable(_ coefficient: Double) throws -> TabulatedDragCoefficients {
    try TabulatedDragCoefficients(samples: [0.0,0.5].flatMap { mach in
        [1e4,1e8].map { DragCoefficientSample(machNumber: mach,reynoldsNumber: $0,dragCoefficient: coefficient) }
    },source: "Synthetic steady drag",version: "1")
}
let descentFlow = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(0,0,10),
                                     properties: ISAAtmosphere().properties(at: 0))
let parachute = ParachuteAerodynamics(drag: try DragAreaAerodynamics(areaSquareMetres: 20,
    referenceLengthMetres: 5,forceApplicationPointMetres: SIMD3(0,0,-2),model: demoDragTable(1.5)))
print("Steady parachute force: \(try parachute.evaluate(flow: descentFlow).forceNewtons) N")
let balloon = BalloonAerodynamics(buoyancy: try Buoyancy(displacedVolumeCubicMetres: 100),
    drag: try DragAreaAerodynamics(areaSquareMetres: 5,referenceLengthMetres: 5,model: demoDragTable(0.5)))
let balloonLoads = try balloon.evaluate(flow: descentFlow,
    field: .hydrostaticGravityMetresPerSecondSquared(SIMD3(0,0,9.80665)))
print("Balloon gross buoyancy plus drag (weight external): \(balloonLoads.forceNewtons) N")

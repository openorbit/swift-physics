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

let coupledWing = try LiftingLineWing(spanMetres: 8,referenceAreaSquareMetres: 16,
    stationCount: 16,source: "Synthetic linear rectangular wing",version: "1") { fraction in
    try LiftingLineSection(chordMetres: 2,liftSlopePerRadian: 5,
        geometricIncidenceRadians: 0.08-0.02*fraction,profileDragCoefficient: 0.02,
        effectiveAngleRangeRadians: -0.2...0.2,reynoldsRange: 1e5...1e7)
}
let coupledWingResult = try coupledWing.evaluate(flow: wingFlow)
print("Spanwise coupled wing: CL = \(coupledWingResult.liftCoefficient), CDi = \(coupledWingResult.inducedDragCoefficient); moment \(coupledWingResult.loads.momentNewtonMetres) N m")

let controlSettings = try [-0.1,0.0,0.1].map { deflection in
    let controlledPolars = try polars.map { polar in
        try AirfoilPolar(reynoldsNumber: polar.reynoldsNumber,machNumber: polar.machNumber,
            samples: polar.samples.map { sample in
                AirfoilSample(angleOfAttackRadians: sample.angleOfAttackRadians,coefficients: AirfoilCoefficients(
                    lift: sample.coefficients.lift+2*deflection,
                    drag: sample.coefficients.drag+0.5*deflection*deflection,
                    moment: sample.coefficients.moment-0.3*deflection))
            })
    }
    return try ControlSurfaceSetting(deflectionRadians: deflection,airfoil: TabulatedAirfoil(
        polars: controlledPolars,source: "Synthetic control polar",version: "1"))
}
let controlModel = try TabulatedControlSurface(settings: controlSettings,
    source: "Synthetic tail controls",version: "1",deflectionConvention: "Positive trailing edge toward section down",
    momentReferenceDescription: "Quarter chord")
let tail = try WingSegment(positionMetres: SIMD3(-3,0,0),
    axes: WingSectionAxes(forward: SIMD3(1,0,0),span: SIMD3(0,0,1)),
    chordMetres: 1,areaSquareMetres: 1,controlDeflectionRadians: 0.05,model: controlModel)
let tailLoads = try tail.evaluate(bodyFlow: wingFlow).loads.referred(to: .zero)
print("Synthetic controlled vertical tail: \(tailLoads.forceNewtons) N; \(tailLoads.momentNewtonMetres) N m")

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

let zeroRateColumn = try AerodynamicCoefficients(force: .zero,moment: .zero)
let capsuleRateDerivatives = BodyAerodynamicRateDerivatives(rollRate: zeroRateColumn,
    pitchRate: try AerodynamicCoefficients(force: .zero,moment: SIMD3(0,-10,0)),yawRate: zeroRateColumn)
let capsuleDamping = try ConstantBodyRateDerivatives(value: capsuleRateDerivatives,
    source: "Synthetic capsule damping example",version: "1",machRange: 0...25,
    reynoldsRange: 1e5...1e9,incidenceRangeRadians: -0.2...0.2,sideslipRangeRadians: -0.1...0.1,
    momentLengthsMetres: SIMD3(repeating: 3),reynoldsReferenceLengthMetres: 3)
let dampedCapsule = try WholeBodyAerodynamics(geometry: geometry,reynoldsReferenceLengthMetres: 3,
    model: LinearBodyRateAerodynamics(staticModel: capsuleData,derivativeModel: capsuleDamping,
        rateReference: .inertialReference,minimumSpeedMetresPerSecond: 100,
        maximumAbsoluteNondimensionalRates: SIMD3(repeating: 0.1)))
let dampedLoads = try dampedCapsule.evaluate(flow: capsuleFlow,angularVelocityRadiansPerSecond: SIMD3(0,0.2,0))
print("Synthetic capsule rate-damping increment: \(dampedLoads.momentNewtonMetres-capsuleLoads.momentNewtonMetres) N m")

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

let molecularGas = try FreeMolecularGas(identifier: "Synthetic single-species example",
    specificGasConstantJoulesPerKilogramKelvin: 300,source: "Demonstration assumptions")
let molecularModel = try MaxwellianFreeMolecularTraction(gas: molecularGas,
    surface: GasSurfaceReflection(specularFraction: 0.2,wallTemperatureKelvin: 300))
let molecularPanel = try AerodynamicSurfacePanel(areaSquareMetres: 1,
    centreMetres: SIMD3(1,0,0),outwardNormal: SIMD3(1,0,0))
let molecularFlow = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(7000,0,0),
    densityKilogramsPerCubicMetre: 1e-8,temperatureKelvin: 200,meanFreePathMetres: 20)
let molecularLoads = try SurfaceTractionAerodynamics(panels: [molecularPanel],referenceLengthMetres: 1,
    model: molecularModel).evaluate(flow: molecularFlow).loads
print("Free-molecular example panel force: \(molecularLoads.forceNewtons) N")

// Nonlinear section/wake coupling; do not wrap the section polar in a separate
// finite-wing correction, which would count induced drag twice.
let wakeSegments = try SegmentedWing.rectangular(spanMetres: 8,chordMetres: 2,
    incidenceRadians: 0.08,model: AirfoilWingSectionModel(airfoil: airfoil)).segments
let wakeSolver = try NonlinearLiftingLine(segments: wakeSegments,wakeLengthMetres: 800,vortexCoreRadiusMetres: 0.001)
let wakeResult = try wakeSolver.evaluate { _ in wingFlow }
print("Nonlinear fixed wake: \(wakeResult.iterations) iterations; residual \(wakeResult.maximumNormalizedCirculationResidual); force \(wakeResult.loads.forceNewtons) N")
try runEntryComparison()
try runPropwashExample()
if CommandLine.arguments.contains("--cfd") { try runCFDBenchmark() }

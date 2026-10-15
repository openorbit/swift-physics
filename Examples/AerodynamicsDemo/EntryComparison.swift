// SPDX-License-Identifier: Apache-2.0
import Foundation
import Aerodynamics
import Atmospheres

// Prescribed-attitude, planar entry harness owned by the example, not Aerodynamics.
// Synthetic capsule coefficients and an ideal pressure-only flat wing are illustrative.
struct EntryDemoReport {
    let speed: Double
    let peakDeceleration: Double
    let peakDynamicPressure: Double
    let peakHeatFlux: Double
    let heatLoad: Double
}

func entryComparison(winged: Bool, step: Double) throws -> EntryDemoReport {
    let alpha = winged ? 0.35 : 0.2, mass = winged ? 2000.0 : 1500.0
    let atmosphere = ISAAtmosphere()
    let geometry = try AerodynamicReferenceGeometry(areaSquareMetres: 5,momentLengthsMetres: SIMD3(repeating: 2))
    let capsule = try WholeBodyAerodynamics(geometry: geometry,reynoldsReferenceLengthMetres: 2,
        model: ConstantBodyCoefficients(value: WindForceCoefficients(drag: 1.4,lift: 0.3)
            .bodyCoefficients(angleOfAttackRadians: alpha,sideslipRadians: 0),
            source: "Synthetic entry demonstration",version: "1"))
    let panels = try (0..<8).map { i in
        try AerodynamicSurfacePanel(areaSquareMetres: 1.5,centreMetres: SIMD3(0,Double(i)-3.5,0),outwardNormal: SIMD3(0,0,1))
    }
    let wing = try SurfacePressureAerodynamics(panels: panels,referenceLengthMetres: 2,
        model: ModifiedNewtonianPressure(specificHeatRatio: 1.4,machRange: 5...25))
    let heating = SuttonGravesHeating(parameters: .earthAirColdWall,coverage: try StagnationHeatingCoverage(
        speedMetresPerSecond: 1500...8000,densityKilogramsPerCubicMetre: 1e-6...1,
        noseRadiusMetres: 0.05...1,machNumber: 5...25))
    let leadingEdge = try SweptCylinderStagnationHeating(referenceSphereModel: heating,
        referenceSphereRadiusMetres: 1,coverage: heating.coverage)
    // state = altitude, horizontal speed, downward speed, externally integrated heat per unit area.
    var state = SIMD4<Double>(70000,5000*cos(0.18),5000*sin(0.18),0)
    var peakDeceleration = 0.0, peakQ = 0.0, peakHeat = 0.0
    func derivative(_ s: SIMD4<Double>) throws -> SIMD4<Double> {
        let properties = try atmosphere.properties(at: s.x), speed = hypot(s.y,s.z)
        // Explicit illustrative mean-free-path approximation; not a transport dataset.
        let path = 6.6e-8*1.225/properties.densityKilogramsPerCubicMetre
        let flow = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(speed*cos(alpha),0,speed*sin(alpha)),
            properties: properties,meanFreePathMetres: path)
        guard let mach = flow.machNumber, (5...25).contains(mach) else { throw AerodynamicError.outsideCoverage }
        let loads = try winged ? wing.evaluate(flow: flow).loads : capsule.evaluate(flow: flow)
        let theta = atan2(s.z,s.y)-alpha
        let horizontalForce = cos(theta)*loads.forceNewtons.x-sin(theta)*loads.forceNewtons.z
        let downwardForce = sin(theta)*loads.forceNewtons.x+cos(theta)*loads.forceNewtons.z
        let query = try StagnationHeatingQuery(flow: flow,gasIdentifier: "earth-air",noseRadiusMetres: winged ? 0.1 : 1,
            surfacePositionMetres: .zero)
        let heat = try winged ? leadingEdge.heatFlux(for: query,effectiveSweepRadians: 0.4) : heating.heatFlux(for: query)
        peakDeceleration = max(peakDeceleration,-(horizontalForce*s.y+downwardForce*s.z)/(mass*speed))
        peakQ = max(peakQ,flow.dynamicPressurePascals); peakHeat = max(peakHeat,heat.wattsPerSquareMetre)
        return SIMD4(-s.z,horizontalForce/mass,9.80665+downwardForce/mass,heat.wattsPerSquareMetre)
    }
    let count = Int((25/step).rounded())
    for _ in 0..<count {
        let a = try derivative(state)
        let b = try derivative(state+a*(step/2))
        let c = try derivative(state+b*(step/2))
        let d = try derivative(state+c*step)
        state += (a+2*b+2*c+d)*(step/6)
    }
    return EntryDemoReport(speed: hypot(state.y,state.z),peakDeceleration: peakDeceleration,
        peakDynamicPressure: peakQ,peakHeatFlux: peakHeat,heatLoad: state.w)
}

func runEntryComparison() throws {
    print("Illustrative 25 s entry; prescribed attitude, ISA below 86 km, synthetic capsule / pressure-only wing. Thermal response external.")
    for winged in [false,true] {
        let coarse = try entryComparison(winged: winged,step: 0.2)
        let fine = try entryComparison(winged: winged,step: 0.1)
        print("\(winged ? "Winged" : "Capsule"): final speed \(fine.speed) m/s; peak deceleration \(fine.peakDeceleration) m/s²; peak q \(fine.peakDynamicPressure) Pa; peak local heat \(fine.peakHeatFlux) W/m²; heat load \(fine.heatLoad) J/m²")
        print("Step refinement 0.2→0.1 s: |Δspeed| \(abs(fine.speed-coarse.speed)) m/s; |Δheat load| \(abs(fine.heatLoad-coarse.heatLoad)) J/m². Numerical agreement does not validate the physical models.")
    }
}

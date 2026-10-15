// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Foundation

/// Local boundary-layer-edge flow, not freestream flow across an intervening shock.
/// Distance follows the surface streamline from the thermal/aerodynamic leading edge.
public struct SurfaceHeatingQuery: Sendable {
    public let edgeFlow: FluidFlow
    public let distanceFromLeadingEdgeMetres: Double
    public let surfacePositionMetres: SIMD3<Double>
    public let wall: AerothermalWallState
    public let transport: ThermalTransportProperties
    public init(edgeFlow: FluidFlow, distanceFromLeadingEdgeMetres: Double,
                surfacePositionMetres: SIMD3<Double>, wall: AerothermalWallState,
                transport: ThermalTransportProperties) throws {
        guard distanceFromLeadingEdgeMetres.isFinite, distanceFromLeadingEdgeMetres > 0 else {
            throw FluidError.invalidGeometry
        }
        guard fluidFinite(surfacePositionMetres) else { throw FluidError.invalidGeometry }
        self.edgeFlow = edgeFlow; self.distanceFromLeadingEdgeMetres = distanceFromLeadingEdgeMetres
        self.surfacePositionMetres = surfacePositionMetres; self.wall = wall; self.transport = transport
    }
}

/// Signed convection: positive into the wall, negative cooling into the gas.
/// Unlike incoming-only SurfaceHeatFlux, this retains convective cooling.
public struct SurfaceConvectionResult: Sendable {
    public let wattsPerSquareMetre: Double
    public let surfacePositionMetres: SIMD3<Double>
    public let heatTransferCoefficientWattsPerSquareMetreKelvin: Double
    public let localNusseltNumber: Double
    public let localReynoldsNumber: Double
    public let prandtlNumber: Double
    public init(wattsPerSquareMetre: Double, surfacePositionMetres: SIMD3<Double>,
                heatTransferCoefficientWattsPerSquareMetreKelvin: Double,
                localNusseltNumber: Double, localReynoldsNumber: Double, prandtlNumber: Double) throws {
        guard wattsPerSquareMetre.isFinite, fluidFinite(surfacePositionMetres),
              [heatTransferCoefficientWattsPerSquareMetreKelvin,localNusseltNumber,localReynoldsNumber,prandtlNumber]
                .allSatisfy({ $0.isFinite && $0 > 0 }) else { throw FluidError.nonFiniteResult }
        self.wattsPerSquareMetre = wattsPerSquareMetre; self.surfacePositionMetres = surfacePositionMetres
        self.heatTransferCoefficientWattsPerSquareMetreKelvin = heatTransferCoefficientWattsPerSquareMetreKelvin
        self.localNusseltNumber = localNusseltNumber; self.localReynoldsNumber = localReynoldsNumber
        self.prandtlNumber = prandtlNumber
    }
}

/// Replaceable local surface convection, separate from blunt-body stagnation heating.
public protocol SurfaceConvectiveHeatingModel: Sendable {
    func evaluate(_ query: SurfaceHeatingQuery) throws -> SurfaceConvectionResult
}

/// Isothermal flat plate, attached laminar parallel flow with zero pressure gradient.
/// Nu_x = 0.332 sqrt(Re_x) Pr^(1/3). O(1) per point; no surface mesh or thermal state.
/// Conservative implementation limits: Mach ≤ 0.3, 1000 ≤ Re_x ≤ 500000,
/// 0.6 ≤ Pr ≤ 15, Kn_x ≤ 0.01 and 0.8 ≤ Tw/Te ≤ 1.25.
/// The caller must establish the plate/flow assumptions; these bounds do not detect transition.
public struct LaminarFlatPlateHeating: SurfaceConvectiveHeatingModel {
    public init() {}
    public var source: String {
        "https://archive.nptel.ac.in/content/storage2/courses/112104159/module2/Lec8.pdf"
    }
    public func evaluate(_ query: SurfaceHeatingQuery) throws -> SurfaceConvectionResult {
        try flatPlateConvection(query,coefficient: 0.332,reynoldsExponent: 0.5,reynoldsRange: 1000...500000)
    }
}


/// Fully turbulent isothermal flat plate from its leading edge, not a transition
/// model. Nu_x=0.0296 Re_x^0.8 Pr^(1/3). Same low-Mach/Kn/temperature restrictions
/// as the laminar model, but 500000 ≤ Re_x ≤ 10000000 and 0.6 ≤ Pr ≤ 15.
public struct TurbulentFlatPlateHeating: SurfaceConvectiveHeatingModel {
    public init() {}
    public var source: String { "https://archive.nptel.ac.in/content/storage2/courses/112104159/module2/Lec12.pdf" }
    public func evaluate(_ query: SurfaceHeatingQuery) throws -> SurfaceConvectionResult {
        try flatPlateConvection(query,coefficient: 0.0296,reynoldsExponent: 0.8,reynoldsRange: 500000...10000000)
    }
}

private func flatPlateConvection(_ query: SurfaceHeatingQuery, coefficient: Double,
                            reynoldsExponent: Double, reynoldsRange: ClosedRange<Double>) throws -> SurfaceConvectionResult {
    let flow = query.edgeFlow, properties = query.transport, x = query.distanceFromLeadingEdgeMetres
    guard let temperature = flow.temperatureKelvin, let mach = flow.machNumber,
          let kn = try flow.knudsenNumber(referenceLengthMetres: x) else {
        throw FluidError.missingFlowProperty
    }
    let wallTemperature = query.wall.temperatureKelvin
    let filmTemperature = 0.5*temperature+0.5*wallTemperature
    guard abs(properties.temperatureKelvin-filmTemperature) <= 1e-8*filmTemperature else {
        throw FluidError.unsupportedQuery
    }
    // Transport properties are at film conditions, independently of edge-flow density/viscosity.
    let re = properties.densityKilogramsPerCubicMetre*flow.speedMetresPerSecond*x/properties.dynamicViscosityPascalSeconds
    let pr = properties.prandtlNumber
    guard re.isFinite else { throw FluidError.nonFiniteResult }
    guard mach <= 0.3, reynoldsRange.contains(re), (0.6...15).contains(pr),
          kn <= 0.01, (0.8...1.25).contains(wallTemperature/temperature) else {
        throw FluidError.outsideCoverage
    }
    let nu = coefficient*pow(re,reynoldsExponent)*pow(pr,1.0/3.0)
    let h = nu*properties.thermalConductivityWattsPerMetreKelvin/x
    return try SurfaceConvectionResult(wattsPerSquareMetre: h*(temperature-wallTemperature),
        surfacePositionMetres: query.surfacePositionMetres,
        heatTransferCoefficientWattsPerSquareMetreKelvin: h,
        localNusseltNumber: nu, localReynoldsNumber: re, prandtlNumber: pr)
}

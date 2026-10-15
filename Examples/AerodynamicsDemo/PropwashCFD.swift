// SPDX-License-Identifier: Apache-2.0
import Foundation
import Dispatch
import Aerodynamics

func runPropwashExample() throws {
    let disk = try ActuatorDiskPropwash(centreMetres: .zero,downstreamAxis: SIMD3(-1,0,0),radiusMetres: 1,
        thrustNewtons: 500,airTorqueNewtonMetres: 40,densityKilogramsPerCubicMetre: 1.2,upstreamAxialSpeedMetresPerSecond: 30)
    let air = try SuperposedAirVelocityField(backgroundMetresPerSecond: SIMD3(-30,0,0),increments: [disk])
    print("Actuator disk: induced speed at disk \(disk.inducedSpeedAtDiskMetresPerSecond) m/s; air at tail \(try air.velocityMetresPerSecond(at: SIMD3(-3,0.4,0))) m/s")
}

// Opt-in workload benchmark, not a frame-budget guarantee or drag validation.
func runCFDBenchmark() throws {
    print("CPU local incompressible CFD: fixed 8 m box, free-slip voxel sphere, prescribed rotor volume forcing; three 0.01 s steps. Use -c release for budget measurements.")
    for count in [16,20,32] {
        let config = try LocalCFDConfiguration(dimensions: SIMD3(repeating: count),minimumCornerMetres: SIMD3(repeating: -4),
            cellSizeMetres: 8/Double(count),densityKilogramsPerCubicMetre: 1.2,speedOfSoundMetresPerSecond: 340,meanFreePathMetres: 1e-7,
            kinematicViscositySquareMetresPerSecond: 1.5e-5,boundaryVelocityMetresPerSecond: SIMD3(-10,0,0))
        var grid = try LocalCFDGrid(configuration: config,solidAtCellCentre: { p in (p*p).sum() < 1 })
        _ = try grid.project(timeStepSeconds: 0.01)
        let disk = try ActuatorDiskPropwash(centreMetres: SIMD3(2,0,0),downstreamAxis: SIMD3(-1,0,0),radiusMetres: 0.75,
            thrustNewtons: 80,airTorqueNewtonMetres: 2,densityKilogramsPerCubicMetre: 1.2,upstreamAxialSpeedMetresPerSecond: 10)
        let source = try ActuatorDiskAccelerationField(disk: disk,thicknessMetres: 1)
        var milliseconds: [Double] = [], iterations = 0, divergence = 0.0
        for _ in 0..<3 {
            let start = DispatchTime.now().uptimeNanoseconds
            let report = try grid.advance(timeStepSeconds: 0.01,accelerationMetresPerSecondSquared: source.accelerationMetresPerSecondSquared)
            let elapsed = DispatchTime.now().uptimeNanoseconds-start
            milliseconds.append(Double(elapsed)*1e-6)
            iterations += report.totalPressureIterations; divergence = max(divergence,report.maximumDivergenceAfterPerSecond)
        }
        let pressure = try grid.pressureLoadsOnSolids()
        print("\(count)³ (\(config.cellCount) cells): mean \(milliseconds.reduce(0,+)/3) ms/step, max \(milliseconds.max()!) ms; \(iterations) pressure iterations total; max divergence \(divergence) 1/s; voxel pressure force \(pressure.forceNewtons) N")
    }
}

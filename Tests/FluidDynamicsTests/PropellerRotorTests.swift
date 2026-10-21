import Testing
@testable import FluidDynamics
@Test func coefficientRotorScalesWithDensityRpmAndConservesShaftPower() throws {
    let rotor=try PropellerRotor(radiusMetres:0.6,maximumThrustCoefficient:0.015,profilePowerCoefficient:0.0006,inertiaKilogramMetresSquared:0.002)
    let a=try rotor.sample(angularSpeed:250,collective:0.5,density:0.018)
    let b=try rotor.sample(angularSpeed:500,collective:0.5,density:0.018)
    let c=try rotor.sample(angularSpeed:250,collective:0.5,density:0.036)
    #expect(abs(b.thrustNewtons-4*a.thrustNewtons)<1e-12)
    #expect(abs(c.thrustNewtons-2*a.thrustNewtons)<1e-12)
    #expect(abs(a.shaftPowerWatts-a.loadTorqueNewtonMetres*250)<1e-12)
    let stopped=try rotor.sample(angularSpeed:0,collective:1,density:0.018)
    #expect(stopped.thrustNewtons==0 && stopped.shaftPowerWatts==0)
    #expect(throws:FluidError.self){try rotor.sample(angularSpeed:250,collective:2,density:0.018)}
}

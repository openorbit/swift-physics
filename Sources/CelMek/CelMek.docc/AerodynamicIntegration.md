# Aerodynamic loads in coupled dynamics

Use ``AerodynamicWrench`` to connect an Aerodynamics evaluator to
``CoupledSpacecraftDynamics``. Its closure receives the current coupled state and
context at each integration stage. Return force and moment in the attitude's body
axes, with an explicit moment reference point. The adapter transfers that moment
to the supplied centre of mass, rotates force into reference axes and returns one
force/torque contribution with zero propellant consumption.

For a supplied whole-body evaluator `capsule`, atmosphere provider `atmosphere`
and thermodynamic profile `profile`, an evaluator can follow this pattern:

```swift
let wrench = try AerodynamicWrench(bodyFrame: bodyFrame) { state, context in
    let air = try atmosphere.sample(for: state.spacecraft,
                                    in: context.referenceContext)
    let p = state.spacecraft.cartesian.position
    let altitude = hypot(hypot(p.x, p.y), p.z) - atmosphere.referenceRadiusMetres
    let properties = try profile.properties(at: altitude)
    let flow = try state.aerodynamicFlow(
        atmosphere: air, origin: context.referenceContext.origin,
        properties: properties)
    return try capsule.evaluate(flow: flow,
        angularVelocityRadiansPerSecond: state.attitude.angularVelocityBodyRadiansPerSecond)
}
let dynamics = try CoupledSpacecraftDynamics(context: context,
                                            wrenchModels: [wrench])
```

This example assumes a spherical, body-centred atmosphere and uses its profile's
radius datum. Sample velocity and properties from a consistent
atmosphere model. `SphericalBodyAtmosphere` supplies co-rotation in inertial axes;
the flow helper subtracts this velocity before rotating the difference into body
axes. It does not add another rotation-induced velocity.

For a segmented wing in a rigidly co-rotating atmosphere, pass the atmospheric
angular velocity transformed into body axes as
`mediumAngularVelocityRadiansPerSecond`, alongside the attitude's inertial body
rate. Alternatively, sample air velocity and properties separately at every strip.
Whole-body damping rates must follow the chosen coefficient dataset's convention;
the adapter does not infer that convention or subtract a medium rate automatically.

Load locations and `centreOfMassBodyMetres` share one body-coordinate origin.
Spacecraft Cartesian position and velocity describe the centre of mass. Do not
register another drag force or torque model for the same aerodynamic load.
Closures must be repeatable evaluations: deployment, actuator evolution and
thermal integration belong to external systems models.

Aerodynamics remains independent of CelMek. CelMek depends on Aerodynamics to
provide this integration layer; neither evaluator owns a vehicle assembly.

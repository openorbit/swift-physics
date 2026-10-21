# Driven suspended wheels

`WheelDynamics.swift` supplies `SuspendedWheel`, `WheelState`, `WheelCommand` and
`SuspendedWheelSet` for the existing independent rigid-body mechanical system.
Body +X is forward, +Y left and +Z up. Positive wheel spin is about the steered
+Y axle and rolls forward. Wheel mass and locked-wheel inertia must already be
included in the chassis mass/tensor. Do not add them a second time.

The initial suspension is a massless strut along body -Z. The caller supplies a
local terrain plane, normal and instantaneous surface velocity for each wheel.
The plane determines hub extension; spring compression, normal velocity damping
and a stiffer bottom-out term produce a nonnegative normal load. Unsupported,
overextended or overturned wheels produce no ground force. Wheel travel is
bounded, but bottom-out is compliant rather than a hard geometric constraint.
No normal attraction, grounded-state snap or render-LOD collision is used.

Longitudinal slip compares contact surface speed to wheel spin times effective
rolling radius. Lateral slip comes from the contact point's relative velocity.
Linear creep stiffness produces longitudinal/lateral demands, limited together
by a Coulomb friction circle (mu times normal load). This permits wheelspin,
sliding and loss of traction. Wheel spin responds to drive torque, ground-contact
moment, bounded brakes and a load-dependent rolling-loss torque. Brake torque
opposes rotation or holds the rotor when within capacity. Regularized creep is
not an exact static-friction constraint; slight drift can remain on slopes.

The complete ground force acts at the contact point. The chassis also receives
the negative derivative/transport of relative rotor angular momentum, including
steering-axis transport and the gyroscopic term. Thus motor reaction is present
in flight; spinning an airborne wheel cannot propel the chassis. Rotor speed and
steering angle are real integrated states, not distance-driven visual animation.
Steering is an ideal bounded-rate command (0.7 rad/s, 0.15 s response), not a rack
and motor model. The body tensor stays fixed as suspension/steering geometry
changes; small wheel translations and inertia reorientation are approximated.

`SuspendedWheelSet.step` evaluates body and wheel rates at the same beginning and
midpoint stages. It commits both only after a successful step; a terrain/load
error leaves the previous state intact. Ground and external-load callbacks must
be pure and share the mechanical system's frame/epoch. External gravity and any
other vehicle/cable loads remain separate. Do not apply wheel wrenches again.

Use steps appropriate to suspension and tyre stiffness/rotor inertia. The 0.01 s
API cap is not a stability guarantee. Open Orbit's declared original rover uses
2 ms steps; arbitrary wheels/terrain need their own convergence checks. This is
not a Pacejka tyre, deformable-soil, terramechanics or full articulated rocker-
bogie solver. Sinkage, tread/rock collision, hub masses, wheel damage, motor
power/energy/thermal limits and drivetrain compliance remain future extensions.

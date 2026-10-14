# Attitude and rotational dynamics

All rotations are right handed. Quantities use SI units and elapsed atomic
seconds. Providers and propagators are immutable queries: a prediction restarts
from its seed and does not consume fuel or advance hidden simulation time.

## State and conventions

`AttitudeState` contains an explicit epoch, inertial reference-frame identity,
body-frame identity, a normalized Hamilton `bodyToReference` quaternion, physical
angular velocity in body axes (rad/s), and optional physical angular acceleration
in body axes (rad/s²). Body vectors transform as `reference = q.matrix * body`.
Quaternions q and −q represent the same rotation. Input normalization handles
large finite scales; a zero/non-finite quaternion is rejected.

`AttitudeProvider.checkedAttitude` verifies exact frame identities and epoch.
`AttitudeOrientationProvider` adapts the physical state to the frame system:
its reference-to-body rotation is q conjugate, and its coordinate-map angular
velocity and acceleration are the negatives of the physical body rates.
The adapter requires angular acceleration, so missing derivatives never silently
become zero. `OrientationAttitudeProvider` performs the reverse conversion and
supports existing planetary and local-orbital orientation providers.

`FrameTree.registerSpacecraftFrame` registers the provider's exact body-frame
identity. Construct that frame with the inertial reference as its parent, then
supply the spacecraft centre-state provider and attitude provider. The existing
origin-acceleration policy still applies. As with other frame providers, a provider
must not capture its owning tree strongly. Forces rotate as free vectors; centre
translation and frame-rate terms are used for point-state transforms, not thrust.

## Prescribed attitudes

`InertialHoldAttitudeProvider` keeps a fixed body-to-reference orientation and
supplies zero angular rates. `UniformSpinAttitudeProvider` prescribes constant
physical spin in body axes with analytic quaternion evolution.

`TrackingAttitudeProvider` supports target pointing, velocity tracking (relative
to an optional moving target), and nadir pointing toward a supplied centre or the
reference origin. Body X follows that direction; body Z follows the projected
explicit roll-reference vector, and Y = Z × X. Coincident targets, zero velocity
and parallel roll references fail instead of selecting arbitrary axes. Inputs
must share a physical origin as well as exact reference axes.

Five-point derivatives of the rotation matrix supply angular velocity and
acceleration. Provider coverage at epoch ±2*derivativeStep is required. The step
must resolve the motion; it is an explicit numerical setting, not an accuracy
certificate. Existing `LocalOrbitalOrientationProvider` laws and mounting rotations
can also be adapted through `OrientationAttitudeProvider`.

Prescribed attitude is a kinematic command. It does not certify that available
hardware can realize the angular rates or torques.

## Inertia, torques and control

`RigidBodyInertia` represents a constant centre-of-mass tensor in body axes
(kg m²), including off-diagonal entries. Its parameters are tensor entries, not
unsigned products of inertia. Scaled Cholesky validation requires positive
definiteness; second-moment positive-semidefiniteness enforces the physical
principal-moment triangle inequalities with a small floating-point allowance.
Planar mass distributions are allowed; singular rotational inertia is rejected.
The same factorization solves torque-to-angular-acceleration systems without
forming an inverse. Angular momentum and kinetic energy are exposed for checks.

`TorqueModel` evaluates a pure `TorqueContribution` about the centre of mass in
body axes (N m). `RotationalEvaluationContext` carries the tensor, body identity,
and inertial force context/physical origin. Torque output epochs and identities
are checked. Translation is optional for standalone rotation and mandatory for
models that need position or environmental fields.

`RigidBodyDynamics` uses Euler's equation
`I ωdot = τ − ω × (I ω)` and quaternion kinematics `qdot = 0.5 q * (ω, 0)`.
The tensor remains constant, and internal wheel momentum is not inferred.
`ConstantBodyTorque` is an explicit external torque input. `AttitudePDController`
uses a shortest-rotation quaternion error and body-relative rate error, with a
norm limit on commanded torque. The exact 180° tie is deterministic and invariant
under quaternion sign. It is an ideal external torque controller; hardware
allocation, sensor delays, wheel saturation/desaturation and power demand belong
to additional models.

`GravityGradientTorque` supplies the leading-order point-mass torque
`3 μ/r³ * rHatBody × (I rHatBody)` using a checked moving centre. It assumes the
spacecraft's extent is small relative to its distance from the source. It adds
torque only; add a translational gravity model separately.

`MagnetorquerTorque` evaluates a bounded commanded body dipole (A m²) against an
explicit `MagneticFieldProvider` (tesla): `τBody = dipoleBody × fieldBody`.
The sample's epoch, frame and physical origin are checked. A constant command also
represents a residual dipole. `UniformMagneticField` is a simple prescribed field,
not an inferred Earth magnetic model. Field gradients, hysteresis, electronics
and electrical power are separate concerns.

The equations follow [MIT rigid-body dynamics notes](https://ocw.mit.edu/courses/8-09-classical-mechanics-iii-fall-2014/6fe39e8d5ce4ce746ca256dfea665eda_MIT8_09F14_Chapter_2.pdf),
[NASA gravity-gradient analysis](https://ntrs.nasa.gov/api/citations/19690028658/downloads/19690028658.pdf),
and [NASA magnetic-torque specifications](https://ntrs.nasa.gov/api/citations/19690020961/downloads/19690020961.pdf).
Implementations are native Swift and use no NASA implementation source.

## Rotation and coupled propagation

`RigidBodyAttitudePropagator` integrates rotation alone with an optional prescribed
`SpacecraftStateProvider` for environmental torque inputs. It reports the attitude,
angular acceleration, energy, reference-frame angular momentum and numerical work.
`AttitudePropagationOptions` supplies angle/rate tolerances, step bounds, sampled
rotation-per-step limits and work budgets. Accuracy presets change local tolerances;
explicit settings remain authoritative.

The native Dormand–Prince 5(4) method compares embedded physical quaternion angles
with q/−q equivalence, angular velocities and quaternion norm drift. Stage vector
transforms use a unit quaternion while quaternion derivatives retain the raw RK
state. Accepted states normalize the quaternion without a sign reset. A sampled
angular-speed step limiter avoids advancing too far through rapid rotation. These
are local numerical controls, not a guarantee of global physical accuracy.

`CoupledSpacecraftState` combines translation, mass and attitude at the same epoch
and reference frame. `SpacecraftWrenchModel` returns force in reference axes, torque
in body axes and reservoir mass rate from one physical model. `BodyMountedThruster`
rotates the spacecraft force direction with the actual attitude, contributes
`leverArmBody × forceBody`, and consumes propellant using an explicit effective
exhaust velocity. The direction is force on the spacecraft, not exhaust direction.

`CoupledSpacecraftDynamics` composes translational forces, torque models and
wrenches. `CoupledSpacecraftPropagator` evolves position, velocity, propellant,
quaternion and angular velocity at the same RK stages. This lets an off-centre
engine rotate its own inertial thrust direction during a burn. Separate attitude
and translation provider views support frames, sampling and event queries.

For prescribed attitude, `BodyMountedThruster.forceModel(attitude:)` integrates
translation and fuel through the existing force API, including its dry-boundary
interface. The adapter does not integrate lever-arm torque. Use the coupled wrench
for that; do not add both representations of the same engine to coupled dynamics.

The current rigid-body approximation holds inertia and centre of mass fixed even
when the propellant reservoir changes. Changing mass geometry, internal reaction
wheels, moving appendages and outflow angular-momentum models are not silently
inferred. Propagation expects continuous models; split and reseed with the carried
state at command discontinuities. Coupled fuel exhaustion fails explicitly; no
automatic engine cutoff or torque schedule is inferred. Forward/backward immutable
replay supports predictions and consistency checks, not physical fuel recovery.

## Demonstration and validation

Run the off-centre Earth-orbit engine example:

```sh
swift run --disable-sandbox AttitudeDynamicsDemo
```

It prints SI position, inertial force direction, physical body rate, quaternion and
mass samples. The initial body +X engine force points along inertial +Y; its lever
arm drives negative body-Z rotation and changes the force direction. Its coupled
provider views also register a fully kinematic spacecraft frame.

Swift Testing checks quaternion signs and frame derivatives, target/nadir/velocity
tracking, full and scaled inertia solves, invalid tensors, an analytic torque-free
symmetric top, asymmetric-body energy and inertial momentum conservation, constant
principal-axis torque, bounded control, magnetic/gravity-gradient signs, prescribed
and coupled thrust, independent Simpson integrals of the off-centre engine's
translation, tighter-tolerance convergence, backward replay, contracts, dry firing
and work-budget exhaustion.

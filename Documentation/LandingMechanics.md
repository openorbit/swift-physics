# Compliant landing mechanics

`CelMek/LandingMechanics.swift` provides independent rigid bodies, cable edges,
spherical contact lobes/feet, local terrain planes and adiabatic gas cushions.
This is a compliant engineering model, not a rigid constraint solver or deformable
fabric model. All quantities use SI units.

A `MechanicalSystem` declares one pseudo-inertial `Frame` and an `Epoch`. Body
positions and velocities are in that reference, orientations map body axes into
it, angular velocities and inertia are in body axes, and mount coordinates are
relative to each body's centre of mass. The connection graph is separate from
frame parentage. Multiple cables between the same pair of bodies are allowed.
Rigid assemblies can still use `MassComponent`; detach rigid components before
creating independent bodies and do not count their masses twice.

Cable extension and attachment-point relative speed produce tension. Slack cables
exert no force, even when the endpoints are moving apart. Tension is clamped to
zero rather than allowing compression. Payout subtracts from extension rate and
advances rest length; retrieval can do positive mechanical work. Each load acts
at both mounts with equal and opposite reference forces and body-axis torques.
Cutting changes connectivity only. Breaking commits at accepted substep boundaries;
`cableLoads` is pure and never changes connectivity. A broken cable contributes no
load after the break boundary. Event time uncertainty is up to one substep.

Contacts use a plane supplied by the terrain provider, normal spring support,
compression-only normal damping and speed-regularized Coulomb friction. Forces
act at the lobe's surface point and include angular point velocity. There is no
adhesion. Optional rolling resistance applies a dissipative torque proportional
to normal load and a resistance length, opposing angular motion tangent to the
plane. Its default is zero; a rigid sphere can roll without friction stopping it.
Regularized friction approaches zero with slip speed; this is not an
exact static-friction constraint. Supply local surface velocity explicitly;
for a rotating surface, evaluate a sufficiently local patch per contact. No
terrain queries depend on visual LOD.

Gas pressure follows `P V^gamma = gasInventory`, reduced by ambient pressure.
Pressure support ramps with penetration over half a lobe radius, approximating
the growing contact footprint rather than applying full-area preload at first
touch. This footprint approximation is distinct from the gas compression volume.
Volume has a lower bound; normal contact stiffness supplies the bottom-out guard.
`vent(retaining:)` removes gas inventory monotonically. Creating a cushion at its
initial absolute pressure is an explicit inflation initialization, not a modeled
inflator or gas flow network. Geometry changes, gas vent rates and deployment
controllers belong to the application. Fabric stress, leaks, rupture, temperature
transport and inter-lobe plumbing are not included.

`step` uses explicit midpoint integration and normalizes quaternions. External
loads are evaluated at the beginning and midpoint and must be pure functions of
sample epoch/state. Gravity, aerodynamic loads and controls belong there; do not
also apply the internal cable reactions. A failed step leaves the system unchanged.
Steps must be positive and at most 0.01 s, but that cap is **not** a stability
guarantee for arbitrary stiffness. Choose a smaller step relative to the highest
contact/cable natural frequency and check convergence. Open Orbit uses 0.002 s
for its declared low-mass landing demonstrations. No continuous collision search,
inter-body collision, rigid joints, line/cable terrain collision or long-duration
tether orbital integrator is provided.

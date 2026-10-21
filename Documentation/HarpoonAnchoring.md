# Harpoons and terrain tethers

`HarpoonAnchoring.swift` extends the compliant landing mechanics with explicit
terrain anchoring and mass-conserving projectile launch. This is an engineering
model; supply measured material properties only when they are actually available.

`MechanicalBody.launching` separates a small rigid component whose centre is
specified in the original body's axes. It subtracts component inertia and the
parallel-axis contributions, moves the retained centre of mass, and carries
rotational point velocity into both new bodies. Equal/opposite firing impulses
supply separation energy and produce recoil/mount torque. No mass is created or
counted twice. Component inertia is in the same body axes. Invalid resulting
inertia throws. Use `retainedCentreBodyMetres` to adjust every retained contact
and attachment coordinate. The projectile starts with the parent's orientation
and angular velocity; an impulse through its own centre adds no projectile spin.

`HarpoonAnchor` declares the surface `Frame` identity. Firing is one-shot and
can explicitly misfire. Incoming impacts test incidence and normal kinetic
energy, using a constant penetration-resistance force:

    estimated depth = m * normalSpeed² / (2 * resistance)

Attachment requires depth between the declared minimum and maximum supported
burial depths and an adequate incidence cosine. Too-hard, too-shallow, grazing
and overpenetrating impacts are rejected. This is not a soil mechanics, ice
fracture or fabric/rock interaction solver. Holding strength is a separate input;
it is not inferred from a visual texture, illumination or the penetration depth.
`surfacePointMetres` records the cable exit point on the terrain in surface axes;
penetration depth is recorded separately, not added as rendered surface relief.
The application supplies rejected-projectile rebound/contact behavior.

`TerrainTether` evaluates the existing cable tension law against a prescribed
moving endpoint. Transform the surface-local cable exit point and its velocity
into the inertial integration frame at every stage. Apply the returned force at
the lander's body-local mount and record the equal/opposite reaction on the
prescribed celestial body. Do not integrate a fictitious anchor mass, aggregate
its inertia with the lander, or reparent the lander's physical motion to it.
The full celestial ephemeris is not perturbed by these ground reactions.

`commitLoad` checks breaking and pullout at accepted substep boundaries. Cable
breakage takes precedence if both thresholds are exceeded. A broken or cut cable
leaves a successfully embedded anchor in place. Pullout instead releases the
projectile to independent dynamics while its cable can remain connected. Cutting
or changing an attachment state applies no hidden impulse or velocity reset to
the lander. The application transfers projectile stopping momentum into the
prescribed terrain, rather than removing that momentum from the lander twice.

Winch payout/retraction is commanded mechanical work. Motor energy/force limits,
capstan friction, subsurface stress distribution, temperature-dependent strength,
foot screws and hold-down thrusters are outside this initial model. Event timing
is resolved to the chosen mechanical substep, not continuously certified.

Harpoon and tether/contact APIs belong to `Mechanics`; import it explicitly.

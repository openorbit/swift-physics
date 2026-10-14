# Parachute and balloon loads

These evaluators compute steady physical loads for caller-supplied geometry and
local flow. They do not own a vehicle, deployment sequence, inflation state,
opening shock, gas inventory or thermal evolution.

## Steady drag

``DragAreaAerodynamics`` applies `F = −velocityDirection q Cd area` at an explicit
force application point and transfers the resulting moment to the requested point.
``ParachuteAerodynamics`` uses this contract for a supplied current canopy.
Area must follow the coefficient dataset's definition, such as nominal canopy area
or projected area. Supply shape, reefing and wake effects through the current area
or a calibrated coefficient provider; count each correction once.

``DragCoefficientModel`` permits custom providers. ``TabulatedDragCoefficients``
validates a complete rectangular Mach/Reynolds grid, nonnegative coefficients,
source and version, and interpolates linearly in each coordinate. It rejects queries
outside the grid and requires sound speed and viscosity. A supplied axis points
along positive vehicle-relative-to-air motion at zero incidence; its angle with
that velocity is available to custom providers. The table may restrict that angle
but does not invent orientation dependence inside the permitted range. Without an
angle limit it declares orientation-invariant coefficients. It does not project the
area automatically or generate lift, pitching moment or canopy stability.

```swift
let data = try TabulatedDragCoefficients(samples: [
    DragCoefficientSample(machNumber: 0, reynoldsNumber: 1e4, dragCoefficient: 1.5),
    DragCoefficientSample(machNumber: 0, reynoldsNumber: 1e6, dragCoefficient: 1.5),
    DragCoefficientSample(machNumber: 0.5, reynoldsNumber: 1e4, dragCoefficient: 1.5),
    DragCoefficientSample(machNumber: 0.5, reynoldsNumber: 1e6, dragCoefficient: 1.5)
], source: "Synthetic demonstration coefficients", version: "1")
let drag = try DragAreaAerodynamics(areaSquareMetres: 20, referenceLengthMetres: 5,
    forceApplicationPointMetres: SIMD3(0, 0, -2), model: data)
let parachute = ParachuteAerodynamics(drag: drag)
let loads = try parachute.evaluate(flow: flow)
```

The example coefficients are synthetic, not canopy performance data. Flow is
sampled at the drag resultant location and must already include vehicle/atmosphere
rotation and winds. All positions and vectors use common axes and origin. A canopy
resultant and a force transmitted to a payload attachment are different physical
loads: the external tether/multibody model establishes transmission. This evaluator
does not move a canopy force to an attachment by assumption.

Zero dynamic pressure or zero supplied area produces zero drag without querying
coefficients. Area cannot be negative; reference length stays positive. A zero area
can represent externally supplied collapsed geometry, but this contract provides
no transient inflation loads or opening-load prediction.

## Displaced-volume buoyancy

``Buoyancy`` accepts volume in m³ and centre of buoyancy in metres. Under local
hydrostatic conditions, it returns `F = −ρ volume gravity`. Alternatively, with an
approximately uniform pressure gradient, it returns `F = −volume gradientPressure`.
``BuoyancyField`` selects one description. The gradient includes the full ambient
pressure field and is not added to a second hydrostatic buoyancy force. Density is
used by the gravity form; the gradient form uses the supplied gradient directly.
For significant atmospheric rotation or other departures from static hydrostatic
conditions, supply the actual pressure gradient rather than substituting vehicle
acceleration for gravity. Large bodies spanning significant gradients need a
distributed pressure integration instead of the local uniform approximation.

``BalloonAerodynamics`` adds this pressure force to optional supplied-shape drag,
retaining each force's moment arm. It returns gross buoyancy, without subtracting
lifting-gas, envelope or payload weight. External dynamics must use the full mass
and apply gravity once. Volume, lifting-gas mass, gas heating/expansion, venting and
envelope constraints remain external.

Tests verify terminal parachute descent against external weight, oblique drag,
force moment arms, buoyancy equilibrium, equivalent hydrostatic/gradient inputs,
coefficient interpolation and supplied-state limits. Results can feed CelMek's
existing aerodynamic wrench adapter without adding systems logic here.

References: [NASA parachute recovery forces](https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/rktvrecv.html)
and [NASA balloon design](https://www.grc.nasa.gov/www/k-12/Numbers/Math/Mathematical_Thinking/designing_a_high_altitude.htm).

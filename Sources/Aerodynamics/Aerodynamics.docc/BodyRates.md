# Body angular-rate aerodynamics

``LinearBodyRateAerodynamics`` adds supplied small-rate derivatives to an existing
whole-body coefficient model. It preserves static lift and trim/restoring moments:

`C = Cstatic + Droll pHat + Dpitch qHat + Dyaw rHat`

The normalization is `rateHat[i] = angularVelocity[i] momentLength[i] / (2 speed)`.
Aircraft typically use span/chord/span; capsule datasets must use their stated
reference dimensions. ``BodyAerodynamicRateDerivatives`` holds three columns of
six signed body-axis coefficient derivatives, allowing force, moment and cross-axis
effects. A negative diagonal moment derivative provides damping in the selected
convention; arbitrary derivatives are not assumed to imply a stable vehicle.

``BodyAerodynamicRateDerivativeModel`` permits flow-dependent derivative datasets.
``ConstantBodyRateDerivatives`` supplies prescribed derivatives with explicit
Mach, Reynolds, incidence and sideslip ranges and source/version metadata. It also
checks the query's moment lengths and Reynolds reference length against its dataset.
Area, axes and moment reference point must match the static coefficient dataset and
the supplied geometry; the query cannot independently certify those conventions.
No real capsule or aircraft derivative database is bundled here.

The wrapper requires a minimum speed and per-axis absolute normalized-rate limits.
These are declared applicability limits, not measured accuracy guarantees. It queries
both the static and derivative providers with zero angular rates, then adds one rate
increment. Do not wrap another rate-corrected model or add these loads to a second
provider representing the same effect. Missing properties and unsupported ranges
throw instead of silently extrapolating.

## Rate reference and rotating air

``BodyAerodynamicRateReference`` explicitly distinguishes inertial-reference rates
from rates relative to a rigidly rotating medium. The query tag must match the
wrapper's configured convention. For inertial-reference data, supply the attitude's
inertial angular velocity expressed in body axes. For medium-relative data, subtract
the atmosphere's angular velocity, also expressed in body axes, then tag the query
as `rigidlyRotatingMedium`. The wrapper does not transform axes or infer atmospheric
rotation. The underlying linear flow remains vehicle-relative-to-air in either case.
Do not subtract atmospheric rotation again from an already medium-relative rate.

Pass the convention through ``WholeBodyAerodynamics``:

```swift
let loads = try body.evaluate(flow: flow,
    angularVelocityRadiansPerSecond: vehicleRateBody - atmosphereRateBody,
    angularRateReference: .rigidlyRotatingMedium)
```

Here `body` uses a wrapper configured for that same medium-relative convention.
CelMek's aerodynamic wrench adapter consumes the returned force and moment as usual;
the caller chooses the rate convention required by the dataset.

## Limits and validation

This is a local linear angular-rate model. It does not implement incidence/sideslip
time derivatives, hysteresis, dynamic stall, oscillatory-frequency corrections or
unsteady wakes. Some experimental pitch-damping datasets combine pitch-rate and
incidence-rate derivatives; they cannot be imported as a pure pitch-rate derivative
without a justified conversion.

Direct coefficient queries below the minimum speed throw. Whole-body evaluation
continues to return zero loads at zero dynamic pressure without querying a model.
This does not estimate spin drag when a body rotates in otherwise stationary air;
surface/segment methods with sampled local velocities are needed for that case.

Tests check an analytic three-axis damping torque, negative rotational power for
the supplied stable diagonal example, reversal and velocity scaling, preservation
of static lift/moment, zero-rate provider queries, reference tags and coverage.
The demo's derivative values are synthetic rather than validated capsule data.

Reference: [NASA body-axis dynamic derivative conventions](https://ntrs.nasa.gov/api/citations/20040090578/downloads/20040090578.pdf).

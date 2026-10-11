# Coordinate and API conventions

This contract applies to the epoch-aware `CelMek` states, orientations and frame
transforms. Legacy astronomical routines retain their individually documented units
and coordinate conventions. See [states](STATES.md), [frames](FRAMES.md) and
[time representations](TIME.md) for model-specific details.

## Coordinates, matrices and units

Cartesian coordinates are right-handed: `x × y = z`. Vectors are column vectors;
`Matrix3x3 * SIMD3<Double>` applies the matrix on the left. Matrix subscripts are
`[row, column]`; `Matrix3x3(elements:)` takes nine entries in row order. In `B * A`,
`A` acts first. Frame names do not prescribe a rendering engine's up or forward axis.

| Quantity at the shared API boundary | Unit / convention |
| --- | --- |
| Position and translation | metres |
| Velocity | metres per elapsed atomic second |
| Acceleration | metres per elapsed atomic second squared |
| Angles | radians |
| Angular velocity | radians per elapsed atomic second |
| Angular acceleration | radians per elapsed atomic second squared |
| Epoch arithmetic | continuous SI seconds, stored in split `Duration` |

`CartesianUnits` converts km/seconds or AU/days inputs once at construction. It
cannot convert axes, origins, spherical coordinates or coordinate time scales.
Derivatives from a model evaluated in TDB must include the time-scale chain rule
before becoming shared states or orientations. A nil acceleration means unavailable.

Angle conversions such as `90.0.deg` and `1.0.arcsec` yield radians; `.asDeg` and
`.asArcSec` convert radians outward. These are `Double` conveniences, not distinct
unit types. Angle normalization is explicit: general angle inputs are not implicitly
wrapped to a universal interval. A quaternion's extracted `angle` is in `0...2π`;
model longitudes and rotation phases may use other documented ranges.

## Quaternions and rotation signs

`Math.Quaternion` uses Hamilton multiplication, with vector part `(x, y, z)` and
scalar part `w`. Constructor/component order is **(x, y, z, w)**. For a unit quaternion,
the matrix action agrees with `q * (v.x, v.y, v.z, 0) * q.conj`.

A positive axis-angle rotates a vector by the right-hand rule. For example:

```swift
import Math

let q = Quaternion(angle: .pi / 2, axis: SIMD3<Double>(0, 0, 1))
let y = q.matrix * SIMD3<Double>(1, 0, 0) // approximately (0, 1, 0)
```

`q2 * q1` applies `q1` first, then `q2`. The quaternions `q` and `-q` represent the
same rotation; component equality alone is not a rotation equivalence test.
Conjugation reverses a unit quaternion's rotation. For a general non-unit quaternion,
the algebraic inverse is `q.conj / dot(q, q)`, not just the conjugate.

Use `.identity` for the identity rotation. The legacy `Quaternion()` initializer
creates four zero components, which do not define a rotation. Raw quaternion component
construction does not validate or normalize. `Quaternion(angle:axis:)` normalizes its
axis but uses preconditions for invalid input: angle/axis must be finite, and a zero
axis is permitted only for zero angle. These preconditions are programming contracts,
not catchable errors. General `Math` operations are not uniformly checked APIs.

`FrameTransform` validates and normalizes a finite, nonzero quaternion using scaled
arithmetic. It throws for a zero or nonfinite quaternion. Use this checked boundary
for externally supplied transform data; do not rely on raw `Quaternion.matrix` to
validate invalid or extreme component values.

A rotating frame's axes and a rotated vector have opposite coordinate signs.
`UniformRotationProvider` takes a positive physical axes rate about its parent-frame
axis. At identity initial orientation, child axes turning +90° about parent +z map
parent +x to child −y. Its parent-to-child coordinate mapping has angular velocity
`−rate * axis`, expressed in child coordinates after the initial orientation.

## Transform direction and derivatives

Every registered edge maps **parent coordinates to child coordinates**. For a query,
`from` is the source frame and `to` is the destination:

```text
p_to = R * p_from + t
v_to = R * v_from + ω × (R * p_from) + ṫ
a_to = R * a_from + 2ω × (R * v_from)
     + α × (R * p_from) + ω × (ω × (R * p_from)) + ẗ
```

`R`, translation `t`, its derivatives, angular velocity `ω` and angular acceleration
`α` are expressed in destination coordinates. All derivatives use the same elapsed
atomic seconds as `Epoch`; `d(R*p)/dt = ω × (R*p)` for a source-fixed vector.
Here `ω` describes the coordinate mapping, including the sign convention above.

Translation is the source origin's coordinates in the destination, not the destination
origin's coordinates in the source. For a body centre `c` in parent axes, the mapping
is `R * (p_parent - c)`, so `t = -R*c`. Translation derivatives include this rotation
as well as the centre's motion.

`a.followed(by: b)` applies `a` first and `b` second: rotation is `b.rotation * a.rotation`
and translation is `b.R * a.t + b.t`. Composition and inversion also carry all derivative
terms. A transform belongs to one epoch; composition and application reject mismatches.
Acceleration remains nil when the input state's acceleration is unavailable.

Frame equality is object identity (`===`), including origin and axes. Equal names do
not connect frames or convert coordinates. `FrameTree.convert` selects the mapping
from the state's frame identity. Low-level `FrameTransform.applying` has no stored
source-frame identity, so the caller must associate it with the correct coordinates.
`pseudoInertial` is a classification flag, not validation of a physical inertial model.

## Validation and error handling

The shared `CelMek` constructors and queries use throwing APIs. Provider errors propagate
without being replaced by a generic error or a fallback state. Catch errors at the
application boundary; select an alternate provider or approximation explicitly.

| Error family | Meaning and caller action |
| --- | --- |
| `TimeArithmeticError` | Nonfinite input or unsupported arithmetic range; validate numeric inputs and simulation intervals. |
| `TDBModelError` | Outside selected coverage, invalid correction, wrong sample epoch, or inversion did not converge; select covered data and correct the model contract. |
| `UTCError` | Invalid table/data, invalid UTC second label, or outside table coverage; correct the label or supply an appropriate table. |
| `CalendarConversionError` | Invalid day/clock, outside calendar range, or missing UTC table; correct the reading or provide the required data. |
| `CartesianStateError` | Nonfinite position, velocity or acceleration, including unit-conversion overflow; correct the state inputs. |
| `FrameTransformError` | Invalid quaternion, nonfinite transform/provider configuration, epoch mismatch, or disconnected/invalid registration; correct data or connect the intended identities. |
| `FrameQueryError` | Provider result uses the wrong frame, or an adapter's weakly referenced tree was released; fix the provider contract or retain the tree. |
| `BodyFixedFrameError` | Origin/orientation frames differ, origin acceleration is missing, or an estimation step is invalid; align frames and choose an explicit acceleration policy. |
| `IAUOrientationError` | Unknown body, legacy model requires opt-in, or query outside selected coverage; choose a supported model and coverage. |
| `EphemerisError` | Invalid metadata/derivative methods, mismatched frame or acceleration claims, coverage failure, duplicate or unsupported body; align provider metadata and use covered queries. |
| `BodyModelError` | Mismatched model identities, invalid surface directions/positions, singular point-mass gravity or nonfinite results; select compatible models and valid local coordinates. |
| `CelestialBodyError` | Invalid identity, physical values or provenance; duplicate catalogue identity or absent body; provide validated metadata and an explicit catalogue. |
| `LocalOrbitalFrameError` | Undefined relative position/velocity/plane, invalid step or plane threshold, or incompatible reference frame; supply nondegenerate motion in one pseudo-inertial frame. |
| `EarthOrientationError` | Invalid EOP/configuration, missing realization label, or outside data coverage; supply valid parameters and covered samples. |

A provider must return the requested epoch and its declared native frame. Frame-tree
edges must return the requested epoch. A query may require neighbouring samples:
Earth-frame derivatives need ±2 derivative steps, and estimated origin acceleration
needs ±1 estimation step. Coverage for the central sample alone is insufficient.
Custom data providers define their own coverage/interpolation errors and provenance.

Finite values, successful conversion, round trips and wide storage ranges do not
establish physical accuracy. Model approximation limits, EOP provenance and coverage
remain part of the caller's selection. Numeric equality is exact for shared state and
epoch values; use appropriate tolerances for floating-point conversions and equivalent
rotations. There is no universal library-wide tolerance.

## Existing verification

[RotationAndScalingTests](Tests/MathTests/RotationAndScalingTests.swift) checks the +z
quarter-turn, axis normalization and conjugate inverse. [FrameTransformTests](Tests/CelMekTests/FrameTransformTests.swift)
checks rotating signs, derivative terms, composition/inversion, identity and provider
errors. Time, calendar, body-fixed and Earth-frame suites check their own boundary and
reference cases. These checks support the implemented conventions rather than a claim
that every legacy routine follows the shared contract.

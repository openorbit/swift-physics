# Cartesian states and providers

`CartesianState` carries an `Epoch`, a `Math.Frame` identity, Cartesian position in metres, velocity in metres per elapsed atomic second, and optional acceleration in metres per second squared. Missing acceleration means unavailable; it does not mean zero.

```swift
import CelMek
import Math

let state = try earth.state(at: .j2000)
let localFrame = Frame(name: "Local inertial", pseudoInertial: true)
let initial = try CartesianState(
  epoch: .j2000, frame: localFrame,
  position: SIMD3<Double>(7000, 0, 0),
  velocity: SIMD3<Double>(0, 7.5, 0),
  units: .kilometresAndSeconds)
// initial.position.x is 7,000,000 m; initial.velocity.y is 7,500 m/s.
```

The initializer normalizes units once at the boundary:

| Input convention | Position multiplier | Velocity multiplier | Acceleration multiplier |
| --- | --- | --- | --- |
| Metres and seconds (default) | 1 | 1 | 1 |
| Kilometres and seconds | 1000 | 1000 | 1000 |
| Astronomical units and days | 149,597,870,700 | AU / 86,400 | AU / 86,400² |

Here AU is exactly 149,597,870,700 metres and a day is 86,400 seconds. Nonfinite components and conversion overflow throw `CartesianStateError`. Input vectors must already be Cartesian coordinates in the supplied frame; units do not convert spherical coordinates or time scales.

Frame equality uses object identity, so two frames with the same name are distinct. A frame defines both origin and axes. [FrameTree](FRAMES.md) supplies transformations between connected frames; changing a frame label cannot convert coordinates. These types do not promise concurrent mutation or Sendable conformance.

## VSOP87 provider

`StateProvider` exposes a native `frame` and `state(at:)`. Existing VSOP87 body instances conform, allowing queries through `any StateProvider`. Results carry the requested epoch and the shared `CelestialFrames.vsop87E` identity.

The native convention is barycentric Cartesian coordinates with dynamical ecliptic and equinox J2000 axes. The original [VSOP87 evaluator](https://github.com/ctdk/vsop87/blob/master/vsop87.f) documents the E variant and its TDB input, AU positions, and AU/day velocities. This frame is distinct from an equatorial ICRF frame.

The adapter evaluates the complete implemented series using split TDB date components, converts output to SI, and applies the TDB-to-atomic time rate so velocity follows `Epoch` intervals. It returns no acceleration. The provisional TDB approximation and its limits are described in [TIME.md](TIME.md); independent characterization and scientific coverage bounds remain outstanding. Precise epoch storage does not establish model accuracy or preserve arbitrarily small changes through every floating-point evaluation.

Legacy `position(jd:)` retains its TDB Julian-date input and AU/AU-day outputs for compatibility. Other legacy ephemerides retain their documented units: the lunar model and spherical outputs require their own coordinate and derivative adapters, and DE441 integration remains separate work. SI normalization currently applies to the new shared state boundary, not every legacy API.

Public-import Swift Testing checks cover input normalization, missing acceleration, frame identity, invalid components and overflow, a planetary reference sample, and all nine VSOP87 bodies' velocities against position differences over elapsed seconds.

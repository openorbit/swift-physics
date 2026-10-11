# Epochs and durations

`Epoch` represents an instant. `CelMek.Duration` represents an elapsed interval in SI seconds. `TimeScale` selects a representation of an instant; changing that representation does not change the instant's physical meaning.

The types live in `CelMek` and are immutable, `Sendable`, `Hashable`, and comparable. `JulianDateComponents` is an immutable, scale-tagged split Julian date. The existing `JulianDate` type remains a Julian **calendar date**, not a time-scale-aware instant.

## Reference instant and storage

`Epoch.j2000` is JD 2451545.0 **TT**, or 2000-01-01 at 12:00 TT. `offsetFromJ2000` is elapsed continuous atomic time in SI seconds from that instant. At this instant the TAI clock representation is 32.184 seconds earlier than the TT representation.

Durations store an `Int64` whole-second part and a `Double` fractional-second part. Fractions are normalized to `0..<1`; for example, −0.25 seconds is represented as −1 whole second plus 0.75 fractional seconds. An epoch stores one such duration, keeping a large epoch offset separate from fine fractional steps.

The supported interval and epoch-offset range is ±9,007,199,254,740,991 seconds (2^53 − 1), approximately 285 million years. This is a storage/arithmetic range, not a scientific validity guarantee for an ephemeris or time-scale model. Fractional values retain ordinary `Double` rounding limits; operations are not exact decimal arithmetic.

## Using the public API

Qualify `CelMek.Duration` when necessary to distinguish it from Swift's standard-library `Duration`.

```swift
import CelMek

let step = try CelMek.Duration(seconds: 0.02)
let later = try Epoch.j2000.advanced(by: step)
let elapsed = try later.duration(since: .j2000)

let tt = try later.julianDate(in: .tt)
let tai = try later.julianDate(in: .tai)
let restored = try Epoch(julianDate: tt)

// Bridge to an existing astronomical algorithm when its conventions match.
let legacyJulianDate: Double = tt.value
```

`Duration.adding`, `subtracting`, and `negated` provide interval arithmetic. `Epoch.advanced(by:)` supports positive and negative intervals. Comparisons and equality use canonical stored values; equality is exact, so conversions through floating-point date representations may require a duration tolerance rather than `==`.

Construction and arithmetic throw `TimeArithmeticError.nonFinite` for NaN/infinite numeric inputs and `.outOfRange` for unsupported storage or result ranges. Differences between two individually valid epochs can exceed the supported duration range and throw. A scale conversion near a range endpoint can also throw when its offset crosses that endpoint.

## Julian dates and precision

Use `JulianDateComponents(wholeDay:fractionalDay:scale:)` to preserve fine date information. The fraction is normalized to `0..<1`; whole Julian day boundaries occur at noon. The scalar initializer `Epoch(julianDate:scale:)` accepts existing `Double` Julian dates, but cannot recover precision already lost by the caller.

At J2000 a single `Double` Julian date has a spacing of about 40 microseconds. Both `JulianDateComponents.value` and `Duration.seconds` are explicitly lossy conversions when the whole part is large. Keep simulation arithmetic in `Epoch`/`Duration`, and convert at the algorithm boundary.

Regression tests demonstrate nanosecond increments on offsets of ±500 billion seconds (roughly ±16,000 years). Split Julian date round trips are checked within 30 picoseconds at the tested epochs, including these large offsets. Those results are tests of the chosen inputs, not universal error bounds; converting fractional days and seconds still rounds.

## Supported scales

| Scale | Adapter |
| --- | --- |
| TAI | Continuous atomic representation. |
| TT | Fixed 32.184-second offset from TAI. |
| TDB | Existing low-order periodic TT/TDB approximation, evaluated from split epoch offsets and inverted iteratively. |

TDB is provisional: its physical error and valid astronomical interval have not been independently characterized. It must not be treated as a precision relativistic conversion merely because the epoch storage is precise. `Epoch` storage and TT/TAI arithmetic are independent of that optional adapter.

UTC is intentionally deferred until a versioned leap-second provider and explicit coverage policy are implemented. UT1 requires separate Earth-orientation data. No local timezone, wall clock, leap-second table, or implicit UTC conversion is used by these types. Calendar parsing and scale-aware civil date presentation remain follow-up work.

## Validation

[EpochTests](Tests/CelMekTests/EpochTests.swift) uses a normal `import CelMek`, exercising public construction and result access. It covers normalization, arithmetic, comparisons/hashing, nonfinite inputs, range errors, J2000 scale offsets, fine steps at distant epochs, split-date precision, conversion round trips, and the millisecond scale of the TDB approximation.

```sh
swift test --filter EpochTests
```

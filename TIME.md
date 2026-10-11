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

The TDB adapter is the two-term approximation
`0.001658*sin(g) + 0.000014*sin(2*g)` seconds, with
`g = (357.53 + 0.9856003*daysSinceJ2000TT)` degrees. It represents a geocentric
conversion, without observer-dependent topocentric terms.

We compared it with the full Fairhead–Bretagnon model in
[ERFA 2.0.1 `eraDtdb`](https://github.com/liberfa/erfa/blob/v2.0.1/src/dtdb.c),
using zero observer coordinates. The reference permits TT as its date argument
and reports accuracy better than ±3 ns relative to DE405 time ephemerides over
1950–2050. That reference accuracy does **not** apply to our two-term adapter.

| Sampled TT interval (inclusive January 1 endpoints) | Cadence / samples | Maximum absolute offset error | RMS error |
| --- | --- | --- | --- |
| 1950–2050 | 6 hours / 146,101 | 36.483 µs | 16.831 µs |
| 1900–2100 | 6 hours / 292,197 | 38.763 µs | 17.011 µs |

These are empirical sampled errors, not continuous error bounds or certification
against a modern DE time ephemeris. The 1950–2050 interval is the characterized
use interval with the reference's stated accuracy; the wider survey is a numerical
stress comparison. No accuracy claim is made outside these surveys, even though
`Epoch` permits a much larger storage range. Existing calls remain available without
an implicit coverage restriction; callers must choose whether this model suits them.

For scale, a 36.483 µs epoch error corresponds to about 1.1 m of position error for
motion at 30 km/s, before other model errors. This can suit game ephemerides, but
metre-level work should select a more accurate time model. A configurable precision
TDB model remains follow-up work. TT/TAI arithmetic and `Epoch` storage are independent
of the TDB approximation.

Regression tests use 4,874 independently generated ERFA samples, including all four
six-hour phases every 30 days, the final endpoint and J2000. The sampled maximum rate
error is 9.919e-12 s/s and acceleration error is 1.121e-17 1/s. Reference derivatives
use centred differences with a 300-second step; these derivative figures include
that numerical approximation and are not analytic error bounds. The analytic
Swift derivatives remain consistent with the implemented two-term model. TDB round
trips check numerical inversion, separately from physical accuracy.

[Survey results](TDB_CHARACTERIZATION.json) record sample counts and worst TT Julian
dates. `Tools/characterize_tdb_adapter.py` regenerates the survey and fixtures using
an external ERFA shared library. This developer reference adds no C target or runtime
dependency to the package.

UTC uses `UTCDate` with an integer modified Julian day (midnight) and split
`secondOfDay`, plus an explicitly supplied `LeapSecondTable`. It is separate from
`TimeScale`: an ordinary 86400-second Julian day cannot unambiguously label the
extra UTC second. No implicit leap-second smearing or wall-clock conversion occurs.

```swift
let table = LeapSecondTable.iersBulletinC72
let leap = try UTCDate(modifiedJulianDay: 57753,
                       secondOfDay: CelMek.Duration(seconds: 86400.25))
let instant = try Epoch(utc: leap, using: table) // 2016-12-31 23:59:60.25 UTC
let restored = try instant.utc(using: table)
let tai = try instant.julianDate(in: .tai)
```

The frozen table follows [IERS Leap_Second.dat](https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat),
updated through Bulletin C 72 (July 2026). Its UTC coverage starts at 1972-01-01
(MJD 41317) and ends exclusively at 2027-06-28 midnight (MJD 61584), the stated
expiry date. Pre-1972 UTC frequency adjustments are unsupported. The table never
extrapolates beyond coverage and never refreshes itself over the network.

Supply updated entries or parse IERS text with `LeapSecondTable(iersData:...)`;
version, source and exclusive coverage are required caller metadata. Data loading
is separate from evaluation. Row MJD values define effective midnights; the redundant
calendar columns are checked for integer syntax. Tables must be ordered and changes
must be ±1 second. Negative leap days are supported, although none occur in the
bundled history. Invalid leap labels and out-of-coverage dates throw `UTCError`.

UTC days contain 86399, 86400 or 86401 labelled seconds. `secondOfDay` is validated
against the chosen table, so 86400 is accepted only on positive leap days. Advance
`Epoch` by SI durations for simulation arithmetic, rather than adding nominal UTC
calendar days. Calendar parsing/presentation and UTC-to-UT1 via DUT1 remain follow-up
work; UT1 still requires separate Earth-orientation data.

## Validation

[EpochTests](Tests/CelMekTests/EpochTests.swift) uses a normal `import CelMek`, exercising public construction and result access. It covers normalization, arithmetic, comparisons/hashing, nonfinite inputs, range errors, J2000 scale offsets, fine steps at distant epochs, split-date precision, conversion round trips, and the millisecond scale of the TDB approximation.

```sh
swift test --filter EpochTests
```

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
metre-level work should select a more accurate time model. The optional full model
below supports callers needing a more accurate conversion. TT/TAI arithmetic and `Epoch` storage are independent
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

### Selecting a TDB model

`TDBTimeModel` supplies a correction at an atomic/TT `Epoch`, a coverage interval and
provenance. `TDBCorrection` carries TDB−TT as a split duration, `dTDB/dTT` as a positive
finite rate, and `d²TDB/dTT²` in 1/s. Returned epochs are validated. Models should be
deterministic for fixed inputs; custom providers define their own data errors and ownership.

`FairheadTDBModel` evaluates the full 787-term geocentric Fairhead–Bretagnon series
and the planetary-mass adjustments used by ERFA 2.0.1, in native Swift. The supported
coverage is 1950-01-01 TT through 2050-01-01 TT inclusive; callers may narrow that
interval. Wider intervals are rejected rather than borrowing the reference's accuracy
claim beyond its stated range. Topocentric observer terms are absent. Analytic first
and second derivatives include the series amplitudes' time dependence.

```swift
let model = try FairheadTDBModel()
let date = try Epoch.j2000.julianDate(in: .tdb, using: model)
let epoch = try Epoch(julianDate: date, using: model)
let reading = try epoch.calendar(in: .tdb, tdbModel: model)
let recovered = try reading.epoch(tdbModel: model)

let earthProvider = VSOP87StateProvider(body: earth, tdbModel: model)
let tree = try CelestialFrames.makeTree()
let icrf = try earthProvider.state(at: epoch, in: CelestialFrames.icrf, using: tree)
let marsAxes = try IAUOrientationProvider(
  body: "mars", coverage: model.coverage, tdbModel: model)
```

Use the same selected model for forward/inverse dates and for providers that share
those epochs. `earth.state(at:using:)` also accepts a model directly. The configured
VSOP provider and IAU orientation provider use one correction sample for both the date
and derivative chain rule. Orientation queries must satisfy both orientation and model
coverage. `TwoTermTDBModel(coverage:)` offers the inexpensive approximation with an
explicit caller-selected interval; that choice does not establish an accuracy guarantee.

Existing calls without a model retain their previous two-term behavior. The optional
model is consulted only for TDB representations; TT/TAI/UTC arithmetic does not need it.
Reverse conversion uses Newton iteration, a 12-iteration limit and explicit errors for
nonconvergence, wrong-epoch corrections and out-of-coverage queries. At covered boundaries,
TDB labels within 30 picoseconds of the boundary's predicted label return that boundary
to accommodate split-JD rounding. No uncovered model sample is evaluated.

The Swift full-series offsets agree with 4,874 independent ERFA fixtures within the
2-femtosecond test tolerance. This is numerical agreement with a reference model, **not**
a femtosecond physical accuracy claim. Derivative checks use independently differenced
reference values; round trips also cover both coverage boundaries and a synthetic model
with a deliberately non-unit rate. The reference's nanosecond-level physical statement
and its model limitations described above still apply. Coefficient extraction is reproducible
with `Tools/generate_tdb_series.py`; BSD attribution is retained in
`Documentation/licenses/ERFA.txt`. No C target or runtime dependency is added.

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
calendar days. UT1 requires separate Earth-orientation data.

`EarthOrientationRecord.sample(using:)` converts UTC/DUT1 to an atomic EOP sample:
UT1−TT = DUT1 − (TAI−UTC) − 32.184 seconds. `sample.dut1(using:)` reverses that
conversion, including inside a labelled leap second. The tabulated EOP provider
interpolates UT1−TT rather than the discontinuous DUT1. `epoch.ut1JulianDate(using:)`
returns distinct two-part UT1 components; UT1 is not a model-free `TimeScale` value.
See [FRAMES.md](FRAMES.md) for measured file ingestion and interpolation policies.

## Gregorian calendar readings

`GregorianDay` validates a proleptic Gregorian date with astronomical year numbering:
year 0 is 1 BCE. Its range is −1,000,000 through +1,000,000 inclusive, with no historical
calendar cutover. It converts to/from integer midnight MJD using integer 400-year eras
([algorithm derivation](https://howardhinnant.github.io/date_algorithms.html)).

`CalendarDateTime` carries a `CalendarTimeScale` (`tai`, `tt`, `tdb`, or `utc`) and
validated hour, minute and split `Duration` seconds. Seconds of 60 are allowed only
in UTC at 23:59; the chosen leap-second table validates the actual date on conversion.
A missing UTC table throws `CalendarConversionError.leapSecondTableRequired`.

```swift
let reading = try CalendarDateTime(
  date: GregorianDay(year: 2016, month: 12, day: 31),
  hour: 23, minute: 59, second: CelMek.Duration(seconds: 60.25), scale: .utc)
let epoch = try Epoch(calendar: reading, using: .iersBulletinC72)
let ttClock = try epoch.calendar(in: .tt)
let utcClock = try epoch.calendar(in: .utc, using: .iersBulletinC72)
let ttJD = try epoch.julianDate(in: .tt)
```

There is no timezone interpretation, implicit UTC default, or wall-clock dependency.
Continuous-scale dates use split Julian dates; UTC retains its civil day and leap label.
Day conversions are integer-exact within the calendar range; fractional clock/Julian-date
conversion rounds in `Double` (tests use a 30-picosecond tolerance). The calendar range
is an API arithmetic range, not a model accuracy guarantee, and UTC table coverage is
usually much narrower. The characterized TDB limits above still apply.

The legacy `Date`, `Time`, `DateTime`, and astronomical calendar algorithms retain their
existing conventions. The new types provide an explicit alternative; locale formatting,
ISO string parsing, timezones and historical-calendar presentation remain separate work.

## Validation

[EpochTests](Tests/CelMekTests/EpochTests.swift) uses a normal `import CelMek`, exercising public construction and result access. It covers normalization, arithmetic, comparisons/hashing, nonfinite inputs, range errors, J2000 scale offsets, fine steps at distant epochs, split-date precision, conversion round trips, and the millisecond scale of the TDB approximation.

```sh
swift test --filter EpochTests
```

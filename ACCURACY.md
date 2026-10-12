# Accuracy hints

`AccuracyLevel` expresses the approximation an end user tolerates. It is a hint:
a model may evaluate at higher accuracy than requested, including full evaluation.
A lower hint does not require degraded results or promise a particular runtime.

| Hint | Intent |
| --- | --- |
| `low` | Allow the model's largest supported approximation to reduce cost |
| `balanced` | Permit a model-specific compromise between cost and approximation error |
| `high` | Prefer smaller approximation error while allowing documented simplifications |
| `full` | Evaluate the complete implemented model |

Levels are not universal physical error bounds. Full evaluation can still mean a
truncated theory, old coefficients, or approximate physics. Each algorithm must
document its selection rule and limitations. Independent physical accuracy and
runtime characterization remain separate from coefficient truncation bounds.

`AccuracyAwareEphemerisProvider.accuracy` is the requested hint;
`effectiveAccuracy` reports the evaluation mode used. `availableAccuracyLevels`
lists the implemented evaluation modes, not the accepted hints. Every hint is
accepted by the analytic adapters. Ephemeris metadata records both requested and
effective levels; nil means unspecified. Composed providers leave those fields
unspecified rather than inventing a single accuracy level for their inputs.

## VSOP87 subsets

`VSOP87Ephemeris` supports all four evaluation modes. It prepares its coefficient
subset once during construction and reuses it for position and analytic velocity.
The original coefficient tables remain intact. Full evaluation retains every term
and remains the default.

For a term A cos(B+Ct) tⁿ, selection bounds its individual positional contribution
by |A| max(|t|)ⁿ over the caller's coverage. Here t is TDB Julian millennia from
J2000. Coverage endpoints use the selected TDB model. The adapter adds a small
margin to max(|t|), preserves the original order of retained terms, and retains
terms at or above these cutoffs:

| Hint | Individual contribution cutoff, AU |
| --- | --- |
| `low` | 1e-6 |
| `balanced` | 1e-8 |
| `high` | 1e-10 |
| `full` | No terms removed |

These are initial algorithm-specific settings, not certified position tolerances.
They depend on coverage: larger intervals retain higher polynomial-order terms.
`coefficientSelection` exposes original/retained term counts, the cutoff, and a
conservative per-axis `positionErrorBoundMetres`. The bound sums omitted maximum
amplitudes relative to the complete implemented series; cancellation generally
makes observed error smaller. Ordinary floating evaluation error remains separate.
It does not bound error against observations or a modern DE solution.

If selection removes no terms, effective accuracy is full even for a lower hint.
Velocity analytically differentiates exactly the selected terms and uses the same
TDB-to-atomic rate conversion as full evaluation. The position bound is not a
velocity error bound: small, fast oscillating terms can have larger derivatives.
[Measured results](VSOP87_ACCURACY_CHARACTERIZATION.md) cover all nine bodies at
4,567 epochs from 1950 to 2050 (eight-day spacing plus the endpoint), using the
same Fairhead TDB model in every mode. These are sampled maxima relative to full
VSOP87E and local release-build timings, not independent physical error bounds.

| Mode | Worst sampled position difference, km | Worst sampled velocity difference, m/s | Query speedup across bodies |
| --- | ---: | ---: | ---: |
| low | 2.9e+03 (Saturn system) | 2.69 | 3.07–10.03× |
| balanced | 55.8 (Mars system) | 0.114 | 1.29–6.88× |
| high | 0.513 (Mercury) | 0.000905 | 1.03–2.86× |

The complete report includes per-body RMS errors, retained counts, query timings,
preparation costs and conservative coefficient bounds. Run the characterization
executable to measure your actual workload and machine. The eight-day sampling
can miss short-period extrema; velocity accuracy deserves its own assessment.

```swift
let tdb = try FairheadTDBModel()
let model = try VSOP87Ephemeris(
  body: earth, coverage: coverage, tdbModel: tdb, accuracy: .balanced)
let effective = model.effectiveAccuracy
let retained = model.coefficientSelection.retainedTermCount
let positionBound = model.coefficientSelection.positionErrorBoundMetres
```

## Other models and independent settings

The shorter Meeus lunar and Goffin Pluto adapters accept every hint and currently
use full evaluation. Their requested hint remains visible while effective accuracy
is full. No artificial reduction or unsupported-level error occurs.

TDB conversion, frame differentiation, coefficient reduction, body orientation,
force selection and integrator tolerances are independent sources of error.
Accuracy hints do not change explicit TDB models, derivative intervals, coverage,
or missing-body handling. Legacy low-level functions keep complete-series behavior.
Selecting full lunar evaluation still evaluates the existing truncated Meeus series,
not the full ELP theory.

Tests verify that lower hints can return full results, full VSOP87 output remains
unchanged, retained term counts increase with tighter hints, sampled position
errors respect the coefficient bound, and selected analytic velocities match
independently differenced positions. Reduced-versus-full checks measure truncation
within one theory, not independent physical validation.

## Numerical spacecraft integration

`NumericalPropagationOptions.preset(for:)` maps shared hints to version 1
Dormand–Prince 5(4) local-error settings:

| Hint | Absolute position, m | Absolute velocity, m/s | Absolute propellant, kg | Relative tolerance |
| --- | ---: | ---: | ---: | ---: |
| `low` | 10 | 1e-2 | 1e-4 | 1e-5 |
| `balanced` | 1e-1 | 1e-4 | 1e-6 | 1e-7 |
| `high` | 1e-3 | 1e-6 | 1e-9 | 1e-9 |
| `full` | 1e-5 | 1e-8 | 1e-11 | 1e-11 |

Every preset uses initial/minimum/maximum steps of 10 / 1e-6 / 300 seconds,
100,000 maximum attempts and 700,001 maximum derivative evaluations. These work
controls remain explicit and do not vary with the hint. The existing default
`NumericalPropagationOptions()` equals the `high` preset, preserving its behavior.

```swift
let provider = try NumericalSpacecraftPropagator(
    initialState: seed, dynamics: dynamics, accuracy: .balanced)
let requested = provider.accuracy
let preset = provider.effectiveAccuracy
let actualSettings = provider.options
```

`accuracy` reports requested intent; `effectiveAccuracy` reports the selected
preset, not measured trajectory error. All levels evaluate every supplied force
model with the same numerical method. `full` selects the tightest preset; it does
not imply exact integration, new force models, complete physics or unlimited work.
Low hints may obtain much better results than their tolerances, particularly for
unforced or polynomial motion; there is no artificial degradation.

For direct control, construct validated options and pass them through either
initializer. Explicit options are authoritative even when a hint is present:

```swift
let settings = try NumericalPropagationOptions(
    absolutePositionToleranceMetres: 1e-4,
    absoluteVelocityToleranceMetresPerSecond: 1e-7,
    relativeTolerance: 1e-11,
    maximumStepSeconds: 60)
let custom = try NumericalSpacecraftPropagator(
    initialState: seed, dynamics: dynamics, accuracy: .low, options: settings)
// custom.accuracy == .low; custom.options == settings.
// custom.effectiveAccuracy == nil: explicit settings are not classified as a preset.
```

The direct-options initializer reports nil for both hint fields. Nil means
unspecified, not full. To adjust a preset's work controls, obtain it with
`preset(for:)` and copy its tolerance fields into the explicit options initializer
with the desired step limits; the options are immutable.

Presets are initial numerical policies, not universal physical or final-position
bounds. Acceptance scales each component by absolute tolerance plus relative
tolerance times its state magnitude. Coordinate origin, trip duration, model
smoothness and supplied forces all matter. Tightening tolerances can increase
cost or fail at work/step/precision limits; there is no automatic downgrade to a
looser preset. See [spacecraft dynamics](SPACECRAFT_DYNAMICS.md) for the error norm,
propellant handling and numerical limits.

Tests verify every preset's values tighten in order, unchanged defaults/work
controls, requested/effective reporting and identical results for explicit
settings regardless of a lower hint. An analytic exponential-velocity problem
shows smaller error and greater evaluation cost for full versus low. This is a
representative convergence check, not a benchmark or an error envelope covering
all spacecraft scenarios; use propagation diagnostics for the actual workload.

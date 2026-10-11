# VSOP87 reduced-mode characterization

These results compare reduced coefficient sets against the complete built-in VSOP87E series. They measure truncation within one model, not physical error against observations or DE ephemerides.

Coverage: TT JD 2433282.5 through 2469807.5 inclusive (1950-01-01 through 2050-01-01). Sample spacing: 8 days, plus the final endpoint; 4567 epochs per body/mode. Maxima are sampled maxima, not guarantees between samples. Norms are Euclidean in the native VSOP87E Cartesian frame.

Time model: Fairhead–Bretagnon 1990, full geocentric ERFA 2.0.1 adaptation. State velocities use atomic seconds. All modes use the same coverage and time model.

## Query cost and approximation differences

| Body | Mode | Retained / full terms | Max position, km | RMS position, km | Max velocity, m/s | RMS velocity, m/s | Query, µs | Speedup |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Sun | low | 65 / 6586 | 755.775 | 329.381 | 0.0300737 | 0.0134068 | 10.25 | 9.16× |
| Sun | balanced | 314 / 6586 | 16.4459 | 7.985 | 0.00314338 | 0.00186409 | 13.64 | 6.88× |
| Sun | high | 1819 / 6586 | 0.368516 | 0.183113 | 4.08916e-05 | 2.17339e-05 | 32.83 | 2.86× |
| Sun | full | 6586 / 6586 | 0 | 0 | 0 | 0 | 93.85 | 1× |
| Mercury | low | 126 / 7871 | 1574.03 | 534.884 | 2.69275 | 0.829168 | 11 | 10× |
| Mercury | balanced | 734 / 7871 | 45.3228 | 12.4263 | 0.114374 | 0.0278137 | 18.74 | 5.89× |
| Mercury | high | 3978 / 7871 | 0.512669 | 0.1869 | 0.000904898 | 0.000342936 | 60.86 | 1.81× |
| Mercury | full | 7871 / 7871 | 0 | 0 | 0 | 0 | 110.4 | 1× |
| Venus | low | 124 / 4657 | 1385.77 | 529.643 | 0.531114 | 0.199825 | 10.98 | 6.41× |
| Venus | balanced | 696 / 4657 | 29.9123 | 11.3224 | 0.0247763 | 0.00549464 | 18.27 | 3.85× |
| Venus | high | 2778 / 4657 | 0.251128 | 0.0730238 | 0.000174514 | 2.92237e-05 | 45.05 | 1.56× |
| Venus | full | 4657 / 4657 | 0 | 0 | 0 | 0 | 70.42 | 1× |
| Earth | low | 143 / 5556 | 1926 | 624.167 | 1.24327 | 0.722269 | 11.19 | 7.23× |
| Earth | balanced | 890 / 5556 | 40.5881 | 13.6216 | 0.0478871 | 0.0174986 | 20.69 | 3.91× |
| Earth | high | 3425 / 5556 | 0.313781 | 0.0779132 | 0.000193568 | 4.93775e-05 | 53.45 | 1.51× |
| Earth | full | 5556 / 5556 | 0 | 0 | 0 | 0 | 80.9 | 1× |
| Mars system | low | 270 / 7575 | 1846.2 | 721.821 | 0.379441 | 0.141431 | 12.82 | 8.27× |
| Mars system | balanced | 1605 / 7575 | 55.8276 | 18.2887 | 0.0155703 | 0.00492511 | 29.91 | 3.54× |
| Mars system | high | 4931 / 7575 | 0.450582 | 0.107177 | 0.000233988 | 3.44021e-05 | 73.15 | 1.45× |
| Mars system | full | 7575 / 7575 | 0 | 0 | 0 | 0 | 106 | 1× |
| Jupiter system | low | 387 / 4542 | 1738.22 | 828.737 | 0.0693785 | 0.0279537 | 14.3 | 4.77× |
| Jupiter system | balanced | 2024 / 4542 | 50.743 | 19.7451 | 0.00965156 | 0.00361728 | 35.34 | 1.93× |
| Jupiter system | high | 3785 / 4542 | 0.174369 | 0.0422236 | 1.38185e-05 | 2.4206e-06 | 57.83 | 1.18× |
| Jupiter system | full | 4542 / 4542 | 0 | 0 | 0 | 0 | 68.13 | 1× |
| Saturn system | low | 687 / 7522 | 2900.42 | 1558.36 | 0.108128 | 0.0322942 | 18.34 | 5.77× |
| Saturn system | balanced | 3907 / 7522 | 52.2638 | 20.4943 | 0.0163451 | 0.00630954 | 59.57 | 1.78× |
| Saturn system | high | 6309 / 7522 | 0.167389 | 0.0520039 | 1.16847e-05 | 2.46345e-06 | 90.17 | 1.17× |
| Saturn system | full | 7522 / 7522 | 0 | 0 | 0 | 0 | 105.8 | 1× |
| Uranus system | low | 704 / 5125 | 2103.97 | 1209 | 0.0826466 | 0.0317825 | 18.12 | 4.13× |
| Uranus system | balanced | 3597 / 5125 | 19.181 | 7.23625 | 0.000666761 | 0.000167514 | 56.13 | 1.33× |
| Uranus system | high | 4801 / 5125 | 0.12759 | 0.0380069 | 3.22576e-06 | 9.20925e-07 | 71.1 | 1.05× |
| Uranus system | full | 5125 / 5125 | 0 | 0 | 0 | 0 | 74.78 | 1× |
| Neptune system | low | 299 / 2389 | 1821.71 | 994.074 | 0.0907455 | 0.0384916 | 12.89 | 3.07× |
| Neptune system | balanced | 1695 / 2389 | 12.3751 | 4.61282 | 0.00155808 | 0.000595579 | 30.58 | 1.29× |
| Neptune system | high | 2257 / 2389 | 0.1013 | 0.0265321 | 2.28581e-06 | 4.98261e-07 | 38.47 | 1.03× |
| Neptune system | full | 2389 / 2389 | 0 | 0 | 0 | 0 | 39.57 | 1× |

## Method and limits

Timing: Release public state query; includes Fairhead TDB conversion, allocations and SI conversion; warm providers; median of rotated-order batches; excludes provider preparation. Each timing batch contains 64 queries; 5 repetitions per mode. Query epochs span the sampled coverage. A checksum consumes returned positions and velocities to prevent unused-result elimination.

Host: Version 27.0.1 (Build 26A434), 8 logical processors; architecture: arm64. Compiler: Apple Swift 6.4 (swiftlang-6.4.0.34.1); release configuration. Query timings are one local measurement and depend on compiler, hardware, system load and selected TDB model. They are not a CI pass/fail threshold or performance guarantee.

Provider preparation is excluded from query timing and recorded separately in the JSON. Preparation selects subsets and allocates their tables; results are reused for every query. Full evaluation also performs metadata/selection preparation, so preparation costs should not be confused with steady-state query cost.

Position differences are checked against each model’s conservative per-axis omitted-term bound, with 0.1 m allowance for floating evaluation differences. The full mode has zero truncation bound. No velocity bound is inferred from positional amplitudes. A small fast term may have a significant derivative. The JSON records per-axis bounds, preparation time and TT epochs of sampled position/velocity maxima (zero for full mode, whose differences are identically zero).

The eight-day grid does not resolve every short-period coefficient. More densely sampled intervals and independent physical reference validation may be appropriate for a particular game. Changing the selected coverage changes the subsets; these numbers do not apply unchanged to every interval.

No cutoff was tuned against these samples. Initial cutoffs remain low 1e-6 AU, balanced 1e-8 AU, high 1e-10 AU per maximum term contribution. The shared levels remain tolerance hints and permit a higher effective evaluation mode. Lunar and Pluto reduced modes are not implemented and still evaluate their full existing approximation.

## Reproduce

```sh
swift run --disable-sandbox -c release VSOP87Characterization VSOP87_CHARACTERIZATION.json 8
python3 Tools/summarize_vsop_characterization.py VSOP87_CHARACTERIZATION.json
```

The optional final argument controls positive integer sample spacing in days. Keep the release configuration and supply --compiler-description and --architecture to the summarizer when publishing a new report. The executable writes a report only at the supplied path, performs no network access and does not modify coefficient tables.

Machine-readable results: [VSOP87_CHARACTERIZATION.json](VSOP87_CHARACTERIZATION.json). Policy and selection rules: [ACCURACY.md](ACCURACY.md).

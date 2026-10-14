# Fast nutation accuracy

`fastNutation(jd:)` is the existing four-term approximation with all trigonometric arguments converted to radians. `nutation(jd:)` evaluates the existing 63-term coefficient table. Both return longitude and obliquity corrections in radians.

## Measured differences

On 2026-10-11, the corrected fast functions were compared with the full functions at 73,414 epochs, one Julian day apart, from JD 2415020.5 through JD 2488433.5 inclusive (calendar dates 1900-01-01 through 2100-12-31). Both models receive exactly the same JD; this comparison does not test time-scale conversion.

| Quantity | Maximum absolute difference (arcseconds) | RMS difference (arcseconds) | JD of sampled maximum |
| --- | ---: | ---: | ---: |
| Nutation in longitude | 0.330237 | 0.114803 | 2481630.5 |
| Nutation in obliquity | 0.088571 | 0.026546 | 2465910.5 |
| True obliquity, fast versus full | 0.092367 | 0.026771 | 2418793.5 |

True obliquity compares `fastTrueObliquityOfTheEcliptic` with `trueObliquityOfTheEcliptic`, so its difference includes both the mean-obliquity polynomial choice and the nutation approximation.

These are empirical differences on a daily grid. They do not bound intermediate epochs, epochs outside the sampled interval, or error relative to observed orientation or a modern reference model. The full coefficient table remains subject to the independent model validation listed in the roadmap.

## Regression budgets and reproduction

[FastNutationAccuracyTests](https://github.com/lorrden/swift-celmek/blob/master/Tests/CelMekTests/FastNutationAccuracyTests.swift) repeats the daily sweep, measures maxima and RMS, and checks these budgets:

| Quantity | Maximum budget (arcseconds) | RMS budget (arcseconds) |
| --- | ---: | ---: |
| Nutation in longitude | 0.35 | 0.12 |
| Nutation in obliquity | 0.10 | 0.03 |
| True obliquity | 0.10 | 0.03 |

Separate fractional-epoch cases check points away from the daily grid. The existing Meeus example at JD 2446895.5 now uses approximation tolerances stated in arcseconds; the original tolerance was approximately two arcseconds and could hide the degree/radian defect.

Run the characterization tests with a Swift 6+ toolchain:

```sh
swift test --filter FastNutationAccuracyTests
```

Use the fast function when these measured differences fit the application's needs, and the full function when a tighter approximation is needed. The shared accuracy policy and independent astronomical reference validation remain future work.

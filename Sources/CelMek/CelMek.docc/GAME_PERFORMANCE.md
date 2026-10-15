# Game integration performance

Run the public-API benchmark with an optimized build:

```sh
swift run -c release GameIntegrationBenchmark > game-benchmark.json
```

The executable requires no DE441 or SPICE assets. It emits JSON with operating
system, architecture, logical processor count, setup costs, fidelity comparisons,
seven timing samples summarized as minimum/median/maximum, and a consumed checksum.
Run on an otherwise quiet target device and retain the Swift version alongside the
report. Rerun after compiler, hardware, model or configuration changes. Timings are
observations, not API promises or pass/fail thresholds.

Ephemeris workloads evaluate all nine built-in VSOP87E bodies at 64 equally spaced
epochs from J2000 through 365.25 TT days, using the same full Fairhead TDB conversion
for every accuracy hint. Modes are prepared separately for this coverage and queried
in low/balanced/high/full order. Each workload warms once, then runs seven sequential
batches using monotonic uptime. Query timings exclude model preparation. Accuracy
rows report total retained terms and the largest sampled Euclidean position/velocity
differences across bodies against full VSOP87E. Zero difference in full mode means
self-comparison, not zero physical error. Differences and retained terms depend on
coverage; this one-year sample does not bound arbitrary dates.

Cache-miss batches invalidate before timing; cache-hit batches reuse all 64 exact
Earth epochs in a 128-entry shared memory cache. Both include batch-array allocation.
The two-body workload evaluates a 7000 km circular Newtonian Earth orbit at 256 epochs
one minute apart. The orbit-sampling workload adds the sampler's validation and array
construction. Frame conversion uses precomputed low-mode Earth states and the
standard VSOP-to-ICRF frame tree. Renderer timing adapts 1024 prebuilt positions around
a common moving origin at 1e6 metres per render unit. Setup, provider queries and GPU
work are excluded from those conversion workloads. Reported cost is per state/sample,
not per entire batch.

For broader model truncation comparisons, see the existing
<doc:VSOP87_ACCURACY_CHARACTERIZATION>, including
per-body sampled differences and coefficient bounds. Those timings used a different
run and configuration and should not be compared directly with this report.

## Practical tradeoffs

Choose fidelity from the model's documented errors and the game's positional needs,
then measure the actual target platform. Full evaluation preserves all implemented
terms but still carries the model's intrinsic physical limitations. Reduced modes
save coefficient work; the shared time conversion and state construction cost remains.
The aggregate nine-body maximum may be dominated by a different body in each mode.

Cache repeated exact epochs when useful, especially when several subsystems query
one playback instant. Continuously changing epochs generally miss; caching does not
provide interpolation. A shared cache also has capacity and synchronization costs.
Orbit samples are independent queries, so path resolution directly changes workload;
reuse sampled buffers until their inputs change under the
<doc:GAME_INTEGRATION#Cache-invalidation>.

Keep native states in Double and convert at the render boundary. Frame and renderer
benchmarks measure CPU conversion only, not draw calls, uploads or application
scheduling. They cannot establish a frame budget for a complete game.

These workloads do not measure numerical force integration, attitude dynamics,
Lambert searches, gravity-assist planning, file-backed kernels, cold filesystem IO,
concurrent contention or end-to-end scene loading. Benchmark those separately with
representative asset and mission configurations before extrapolating these results.
No wall-clock assertions are added to the correctness test suite.

## Measured baseline

Recorded 2026-10-13 on this host: arm64, eight logical processors,
macOS version 27.0.1 (build 26A434), Apple Swift 6.4
(swiftlang-6.4.0.34.1). The optimized binary was linked successfully; the sandbox
blocked the subsequent dSYM generation, so this run invoked the linked release
binary directly. Raw results are in `GAME_PERFORMANCE.json`.
The hardware CPU model was unavailable in the sandbox. No claim is made that this
host represents a supported target platform.

| Workload | Median µs per state | Min–max µs per state |
| --- | ---: | ---: |
| nine-body VSOP87E low | 10.338 | 10.205–10.385 |
| nine-body VSOP87E balanced | 27.361 | 26.608–29.612 |
| nine-body VSOP87E high | 46.300 | 45.840–46.624 |
| nine-body VSOP87E full | 83.241 | 82.897–89.386 |
| Earth low cache miss | 11.254 | 11.189–11.366 |
| Earth low cache hit | 0.550 | 0.546–0.553 |
| two-body state | 0.271 | 0.271–0.272 |
| 256-point orbit sampling | 0.299 | 0.298–0.300 |
| VSOP axes to ICRF conversion | 1.500 | 1.458–1.558 |
| moving-origin Float adapter | 0.010 | 0.010–0.010 |

| VSOP hint | Total retained/full terms | Sampled max position difference, km | Sampled max velocity difference, m/s | Nine-provider setup, ms |
| --- | ---: | ---: | ---: | ---: |
| low | 2558/51823 | 1524.49 | 1.90078 | 3.217 |
| balanced | 13758/51823 | 41.4018 | 0.0787766 | 0.996 |
| high | 26773/51823 | 0.355955 | 0.000643356 | 0.906 |
| full | 51823/51823 | 0 | 0 | 0.714 |

In this run, low-mode native queries were about eight times faster than full-mode
queries, while their largest sampled position difference was about 1524 km. High
mode reduced that difference to about 356 metres while costing roughly half of full
evaluation. These are sampled cross-body maxima, not universal accuracy guarantees.
The exact-epoch cache substantially reduced repeated-query cost, while misses added
storage overhead. Setup measurements are single observations in a fixed order and
include first-use effects; do not interpret their ordering as a mode speed ranking.
The renderer workload lasts only microseconds per batch and is especially sensitive
to timer overhead and optimizer effects; use it as a small-operation baseline.

For repeatability assumptions, platform limits and replay/comparison policy, see
<doc:REPRODUCIBILITY>.

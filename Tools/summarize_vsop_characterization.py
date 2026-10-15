"""Render a VSOP87 characterization JSON report as Markdown.

Usage: python3 Tools/summarize_vsop_characterization.py VSOP87_CHARACTERIZATION.json
"""
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('report')
parser.add_argument('--output',default='VSOP87_ACCURACY_CHARACTERIZATION.md')
parser.add_argument('--compiler-description',default='not recorded')
parser.add_argument('--architecture',default='not recorded')
args = parser.parse_args()
r = json.loads(Path(args.report).read_text())
rows = r['measurements']
assert r['schemaVersion'] == 1 and len(rows) == 36
assert all(row['sampleCount'] == rows[0]['sampleCount'] for row in rows)
lines = [
    '# VSOP87 reduced-mode characterization', '',
    'These results compare reduced coefficient sets against the complete built-in VSOP87E series. They measure truncation within one model, not physical error against observations or DE ephemerides.', '',
    f"Coverage: TT JD {r['coverageStartTTJulianDate']} through {r['coverageEndTTJulianDate']} inclusive (1950-01-01 through 2050-01-01). Sample spacing: {r['samplingStepDays']} days, plus the final endpoint; {rows[0]['sampleCount']} epochs per body/mode. Maxima are sampled maxima, not guarantees between samples. Norms are Euclidean in the native VSOP87E Cartesian frame.", '',
    f"Time model: {r['tdbModel']}. State velocities use atomic seconds. All modes use the same coverage and time model.", '',
    '## Query cost and approximation differences', '',
    '| Body | Mode | Retained / full terms | Max position, km | RMS position, km | Max velocity, m/s | RMS velocity, m/s | Query, µs | Speedup |',
    '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |',
]
for row in rows:
    lines.append(f"| {row['body']} | {row['mode']} | {row['retainedTerms']} / {row['originalTerms']} | {row['maximumPositionErrorMetres']/1000:.6g} | {row['rmsPositionErrorMetres']/1000:.6g} | {row['maximumVelocityErrorMetresPerSecond']:.6g} | {row['rmsVelocityErrorMetresPerSecond']:.6g} | {row['medianQueryNanoseconds']/1000:.4g} | {row['speedupAgainstFull']:.3g}× |")
lines += ['', '## Method and limits', '',
    f"Timing: {r['timingDescription']}. Each timing batch contains {r['queriesPerRepetition']} queries; {r['timingRepetitions']} repetitions per mode. Query epochs span the sampled coverage. A checksum consumes returned positions and velocities to prevent unused-result elimination.", '',
    f"Host: {r['operatingSystem']}, {r['processors']} logical processors; architecture: {args.architecture}. Compiler: {args.compiler_description}; release configuration. Query timings are one local measurement and depend on compiler, hardware, system load and selected TDB model. They are not a CI pass/fail threshold or performance guarantee.", '',
    'Provider preparation is excluded from query timing and recorded separately in the JSON. Preparation selects subsets and allocates their tables; results are reused for every query. Full evaluation also performs metadata/selection preparation, so preparation costs should not be confused with steady-state query cost.', '',
    'Position differences are checked against each model’s conservative per-axis omitted-term bound, with 0.1 m allowance for floating evaluation differences. The full mode has zero truncation bound. No velocity bound is inferred from positional amplitudes. A small fast term may have a significant derivative. The JSON records per-axis bounds, preparation time and TT epochs of sampled position/velocity maxima (zero for full mode, whose differences are identically zero).', '',
    'The eight-day grid does not resolve every short-period coefficient. More densely sampled intervals and independent physical reference validation may be appropriate for a particular game. Changing the selected coverage changes the subsets; these numbers do not apply unchanged to every interval.', '',
    'No cutoff was tuned against these samples. Initial cutoffs remain low 1e-6 AU, balanced 1e-8 AU, high 1e-10 AU per maximum term contribution. The shared levels remain tolerance hints and permit a higher effective evaluation mode. Lunar and Pluto reduced modes are not implemented and still evaluate their full existing approximation.', '',
    '## Reproduce', '', '```sh',
    'swift run --disable-sandbox -c release VSOP87Characterization VSOP87_CHARACTERIZATION.json 8',
    'python3 Tools/summarize_vsop_characterization.py VSOP87_CHARACTERIZATION.json', '```', '',
    'The optional final argument controls positive integer sample spacing in days. Keep the release configuration and supply --compiler-description and --architecture to the summarizer when publishing a new report. The executable writes a report only at the supplied path, performs no network access and does not modify coefficient tables.', '',
    'Machine-readable results: [VSOP87_CHARACTERIZATION.json](VSOP87_CHARACTERIZATION.json). Policy and selection rules: [ACCURACY.md](ACCURACY.md).', '',
]
Path(args.output).write_text('\n'.join(lines))

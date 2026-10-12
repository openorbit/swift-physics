#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Select independent JPL testpo values; never evaluate the Swift reader."""
import argparse
import hashlib
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('testpo', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
raw = args.testpo.read_bytes()
text = raw.decode('ascii')
start, end, record_days = -3027215.5, 7930192.5, 32
rows = []
for line in text.split('EOT', 1)[1].splitlines():
    if not line.strip():
        continue
    fields = line[15:].split()
    if len(fields) != 5:
        raise ValueError(f'Malformed reference row: {line}')
    jd, target, center, component, value = fields
    row = dict(jd=float(jd), target=int(target), center=int(center),
               component=int(component), expected=float(value))
    if start <= row['jd'] <= end and 1 <= row['target'] <= 14:
        rows.append(row)
selected = {}
def add(row):
    selected[(row['jd'], row['target'], row['center'], row['component'])] = row
# Every vector/nutation component for every supported target, near J2000.
for target in range(1, 15):
    for component in range(1, 5 if target == 14 else 7):
        add(min((r for r in rows if r['target'] == target and r['component'] == component),
                key=lambda r: abs(r['jd']-2451545)))
# Independently supplied rows exactly on record boundaries and in interiors,
# including the Earth/Moon/EMB reconstruction and both position and velocity.
for target in [3, 10, 13, 14]:
    for boundary in [False, True]:
        for derivative in [False, True]:
            position_components = 2 if target == 14 else 3
            candidates = [r for r in rows if r['target'] == target
                          and ((r['jd']-start) % record_days == 0) == boundary
                          and (r['component'] > position_components) == derivative]
            add(min(candidates, key=lambda r: abs(r['jd']-2451545)))
# Modern relative Earth/Moon/EMB values also verify the SI adapter.
for target in [3, 10, 13]:
    for component in range(1, 7):
        candidates = [r for r in rows if r['target'] == target
                      and r['component'] == component and 2433282.5 <= r['jd'] <= 2469807.5]
        if not candidates:
            raise ValueError(f'Missing modern reference for {target}/{component}')
        add(min(candidates, key=lambda r: abs(r['jd']-2451545)))
# Widely separated epochs from the supplied subset, not only modern dates.
for jd in [start, start+(end-start)/4, start+(end-start)/2, start+3*(end-start)/4, end]:
    add(min(rows, key=lambda r: abs(r['jd']-jd)))
fixture = dict(source='JPL DE441 testpo.441 independently supplied reference values',
    sourceURL='https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de441/testpo.441',
    sourceSHA256=hashlib.sha256(raw).hexdigest(),
    binaryCoverage=[start,end], recordDays=record_days, tolerance=1e-12,
    units='vectors AU and AU/TDB-day; nutations radians and radians/TDB-day',
    rows=sorted(selected.values(), key=lambda r: (r['jd'], r['target'], r['component'])))
args.output.write_text(json.dumps(fixture, indent=2)+'\n')
print(f"Selected {len(fixture['rows'])} independent rows")

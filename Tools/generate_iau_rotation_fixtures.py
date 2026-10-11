#!/usr/bin/env python3
"""Generate independent rotation/derivative fixtures with CSPICE N0067+.

Usage: python generate_iau_rotation_fixtures.py pck00011.tpc output.json
Requires spiceypy and numpy only when regenerating fixtures, never in Swift CI.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import spiceypy as spice

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('kernel', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
generator = ast.parse(Path(__file__).with_name('generate_iau_rotation_data.py').read_text())
constants = {node.targets[0].id: ast.literal_eval(node.value) for node in generator.body
             if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)
             and node.targets[0].id in ('BODIES', 'EXPECTED_SHA256')}
sha = hashlib.sha256(args.kernel.read_bytes()).hexdigest()
if sha != constants['EXPECTED_SHA256']:
    raise SystemExit('Kernel hash differs from reviewed pck00011.tpc')
spice.kclear()
spice.furnsh(str(args.kernel.resolve()))
rows = []
for name, body_id in constants['BODIES'].items():
    for days in [0, 1, -1, 1234.5, -1234.5, 36525, -36525]:
        seconds = days * 86400
        transform = spice.tisbod('J2000', body_id, seconds)
        # Reference acceleration from a centred derivative of CSPICE Rdot.
        # A five-point derivative with a 10 s step limits both cancellation
        # at distant epochs and rotational-curvature truncation error.
        second = (8*(spice.tisbod('J2000', body_id, seconds + 10)[3:, :3]
                     - spice.tisbod('J2000', body_id, seconds - 10)[3:, :3])
                  - (spice.tisbod('J2000', body_id, seconds + 20)[3:, :3]
                     - spice.tisbod('J2000', body_id, seconds - 20)[3:, :3])) / 120
        rows.append(dict(body=name, days=days, rotation=transform[:3, :3].flatten().tolist(),
                         derivative=transform[3:, :3].flatten().tolist(),
                         secondDerivative=second.flatten().tolist()))
args.output.write_text(json.dumps(dict(toolkit=spice.tkvrsn('TOOLKIT'),
    kernelSHA256=sha, frame='J2000', secondDerivativeStepSeconds=10, secondDerivativeStencil="five-point", samples=rows), indent=2)+'\n')
spice.kclear()
print(f'Generated {len(rows)} independent samples')

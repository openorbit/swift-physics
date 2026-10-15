"""Extract the bundled body's GM reference subset from NAIF gm_de440.tpc.

Usage: python3 Tools/generate_body_gm_fixtures.py /path/to/gm_de440.tpc
Source: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc
Expected snapshot: 2022-12-14. All fixture values retain source km^3/s^2 units.
"""
import argparse
import json
import re
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('kernel')
args = parser.parse_args()
text = Path(args.kernel).read_text()
selected = set(range(1,11)) | {199,299,301,399,499,599,699,799,899,999}
values = {}
for code, value in re.findall(r'BODY(\d+)_GM\s*=\s*\(\s*([\d.+\-EDed]+)\s*\)',text):
    code = int(code)
    if code in selected:
        if code in values:
            raise ValueError(f'duplicate GM for {code}')
        values[code] = float(value.upper().replace('D','E'))
if set(values) != selected:
    raise ValueError('Missing requested GM entries')
rows = [dict(naif=code,km3PerSecond2=values[code]) for code in sorted(values)]
Path('Tests/CelMekTests/Fixtures/jpl-gm-2022.json').write_text(json.dumps(rows,indent=2)+'\n')

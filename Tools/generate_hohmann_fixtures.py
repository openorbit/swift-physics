#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Independent 80-digit energy/vis-viva transfer references, standard library only."""
import json
from decimal import Decimal, localcontext
from pathlib import Path
PI = Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628')
rows = []
with localcontext() as ctx:
    ctx.prec = 80
    for r1, r2, mu in [(7e6,42164e3,3.986004418e14),(42164e3,7e6,3.986004418e14),
                       (1,2,1),(2,1,1),(7e6,7e6,3.986004418e14),(7e6,7e6+1e-6,3.986004418e14),
                       (0.001,1000,1)]:
        x,y,m = map(Decimal.from_float,[float(r1),float(r2),float(mu)])
        a = (x+y)/2
        circular1,circular2 = (m/x).sqrt(),(m/y).sqrt()
        transfer1,transfer2 = (m*(2/x-1/a)).sqrt(),(m*(2/y-1/a)).sqrt()
        dv1,dv2 = transfer1-circular1,circular2-transfer2
        rows.append(dict(r1=r1,r2=r2,gm=mu,a=float(a),e=float(abs(x-y)/(x+y)),
                         dv1=float(dv1),dv2=float(dv2),total=float(abs(dv1)+abs(dv2)),
                         seconds=float(PI*(a**3/m).sqrt())))
Path('Tests/CelMekTests/Fixtures/hohmann-reference.json').write_text(json.dumps(rows,indent=2)+'\n')

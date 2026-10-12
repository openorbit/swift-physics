#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""80-digit independent piecewise constant-acceleration integration references."""
import json
from decimal import Decimal, localcontext
from pathlib import Path
rows = []
with localcontext() as ctx:
    ctx.prec = 80
    for start,end,a in [([0,0,0],[100,0,0],4),([12,-30,7],[-288,370,-1193],2),
                         ([1e11,-2e11,3e11],[1.1e11,-1.8e11,2.7e11],9.80665),
                         ([1,2,3],[1,2,3],1)]:
        p,q = [list(map(lambda x: Decimal.from_float(float(x)),v)) for v in [start,end]]
        acc = Decimal.from_float(float(a))
        delta = [y-x for x,y in zip(p,q)]
        distance = sum(x*x for x in delta).sqrt()
        half = (distance/acc).sqrt()
        direction = [x/distance if distance else Decimal(0) for x in delta]
        samples = []
        for fraction in [0,.125,.25,.5,.75,.875,1]:
            # Use exact binary64 elapsed inputs, not fractions inside Swift.
            seconds = float(2*half*Decimal.from_float(fraction))
            t = Decimal.from_float(seconds)
            if fraction == 1 or not distance:
                position,velocity,acceleration = q,[Decimal(0)]*3,[Decimal(0)]*3
            elif t < half:
                position = [x+u*acc*t*t/2 for x,u in zip(p,direction)]
                velocity = [u*acc*t for u in direction]
                acceleration = [u*acc for u in direction]
            else:
                remaining = 2*half-t
                position = [y-u*acc*remaining*remaining/2 for y,u in zip(q,direction)]
                velocity = [u*acc*remaining for u in direction]
                acceleration = [-u*acc for u in direction]
            # At midpoint, explicitly use the braking-side convention.
            if fraction == .5 and distance:
                acceleration = [-u*acc for u in direction]
            samples.append(dict(seconds=seconds,position=list(map(float,position)),
                                velocity=list(map(float,velocity)),acceleration=list(map(float,acceleration))))
        rows.append(dict(start=start,end=end,acceleration=a,distance=float(distance),
                         seconds=float(2*half),peak=float(acc*half),samples=samples))
Path('Tests/CelMekTests/Fixtures/brachistochrone-reference.json').write_text(json.dumps(rows,indent=2)+'\n')

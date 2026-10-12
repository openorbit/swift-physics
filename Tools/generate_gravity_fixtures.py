#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Independent 80-digit direct Newtonian acceleration differences."""
import json
from decimal import Decimal,localcontext
from pathlib import Path
rows=[]
with localcontext() as ctx:
    ctx.prec=80
    for p,R,mu in [([7e6,0,0],None,3.986004418e14),([1,2,3],None,1),
                    ([7e6,0,0],[384400e3,0,0],4.9048695e12),
                    ([1e-3,2e-3,-3e-3],[1e12,0,0],1e20),
                    ([2,3,4],[-4,8,1],10),([0,0,0],[3,4,5],1)]:
        pD=[Decimal.from_float(float(x)) for x in p]; m=Decimal.from_float(float(mu))
        def attraction(v):
            length=sum(x*x for x in v).sqrt()
            return [m*x/length**3 for x in v]
        if R is None:
            result=attraction([-x for x in pD])
        else:
            RD=[Decimal.from_float(float(x)) for x in R]
            direct=attraction([y-x for x,y in zip(pD,RD)])
            indirect=attraction(RD)
            result=[x-y for x,y in zip(direct,indirect)]
        rows.append(dict(position=p,bodyPosition=R,gm=mu,acceleration=list(map(float,result))))
Path('Tests/CelMekTests/Fixtures/gravity-reference.json').write_text(json.dumps(rows,indent=2)+'\n')

#!/usr/bin/env python3
"""Optional black-box fixture generator. Requires binary SpiceyPy/CSPICE and NumPy.
Uses public writer/query APIs only; do not consult toolkit implementation source.
Run with PYTHONPATH pointing to an optional external binary-wheel installation.
"""
from pathlib import Path
import hashlib, json, math, platform
import numpy as np
import spiceypy as s

out = Path(__file__).resolve().parents[1] / 'Tests/CelMekTests/Fixtures/spice-oracle'
out.mkdir(parents=True, exist_ok=True)
s.kclear()
for name in ['states.bsp','attitudes.bc','rotation.bpc']:
    (out/name).unlink(missing_ok=True)
refs = dict(schema=1, toolkit=s.tkvrsn('TOOLKIT'), wrapper=s.__version__,
            generator='Scripts/generate_spice_oracle.py', platform=platform.platform(),
            basis='Public NAIF format/behavior/API documentation; binary oracle only',
            spk=[], ck=[], pck=[], utc=[], sclk=[])
# Position coefficients deliberately differ from supplied type-3 velocity coefficients.
p = np.array([[2,.3,-.2],[4,-.1,.05],[-1,.2,.01]])
v = np.array([[.4,.02,0],[.1,0,.003],[-.2,.01,0]])
h = s.spkopn(str(out/'states.bsp'),'CELMEK ORACLE SPK',0)
for body,kind in [(399,2),(499,3)]:
    coeff = p if kind == 2 else np.vstack([p,v])
    s.__getattribute__('spkw%02d'%kind)(h,body,3,'J2000',-20,20,'TYPE %d'%kind,20,2,2,
                                     np.vstack([coeff,coeff+0.01]).ravel(),-20)
s.spkw02(h,3,0,'J2000',-20,20,'CENTRE CHAIN',40,1,1,[10,.5,-3,.2,2,-.1],-20)
s.spkcls(h); s.furnsh(str(out/'states.bsp'))
for body in [399,499]:
    for t in [-20,-17.125,-0.000001,0,6.25,20]:
        state, _ = s.spkez(body,t,'J2000','NONE',0)
        refs['spk'].append(dict(body=body,centre=0,t=t,state=state.tolist()))

h = s.ckopn(str(out/'attitudes.bc'),'CELMEK ORACLE CK',0)
times = np.array([0.,20.,60.,80.,100.]); angles=[.3,.7,1.1,1.4,1.6]
# Use oracle matrix-to-quaternion conversion to avoid sharing Swift quaternion assumptions.
quats=np.array([s.m2q(s.rotate(a,3) @ s.rotate(.4,1)) for a in angles])
av=np.array([[i*.01,.02,-.03] for i in range(5)])
s.ckw01(h,0,100,-77001,'J2000',True,'TYPE 1',5,times,quats,av)
s.ckw02(h,0,100,-77002,'J2000','TYPE 2',2,np.array([0.,60.]),np.array([40.,100.]),
        quats[[0,2]],av[[0,2]],np.array([2.,3.]))
s.ckw03(h,0,100,-77003,'J2000',True,'TYPE 3',5,times,quats,av,3,np.array([0.,60.,100.]))
s.ckw03(h,0,100,-77004,'J2000',False,'TYPE 3 NO RATE',5,times,quats,av,3,np.array([0.,60.,100.]))
s.ckcls(h); s.furnsh(str(out/'attitudes.bc'))
for inst,queries in [(-77001,[(0,0),(20,0),(19,2),(100,0)]),
                     (-77002,[(0,0),(15,0),(40,0),(60,0),(77,0),(100,0)]),
                     (-77003,[(0,0),(7,0),(20,0),(35,16),(60,0),(73,0),(100,0)]),
                     (-77004,[(7,0),(73,0)])]:
    for tick,tol in queries:
        if inst == -77004:
            matrix,selected=s.ckgp(inst,tick,tol,'J2000'); rate=None
        else:
            matrix,rate,selected=s.ckgpav(inst,tick,tol,'J2000'); rate=rate.tolist()
        refs['ck'].append(dict(instrument=inst,tick=tick,tolerance=tol,selected=selected,matrix=matrix.tolist(),rate=rate))

h=s.pckopn(str(out/'rotation.bpc'),'CELMEK ORACLE PCK',0)
s.pckw02(h,3000,'J2000',-20,20,'PCK TYPE 2',20,2,2,np.vstack([p*.1,p*.1+.01]).ravel(),-20)
s.pckcls(h)
(out/'frames.tf').write_text("""KPL/FK
\\begindata
FRAME_ORACLE_PCK = 13000
FRAME_13000_NAME = 'ORACLE_PCK'
FRAME_13000_CLASS = 2
FRAME_13000_CLASS_ID = 3000
FRAME_13000_CENTER = 301
""")
s.furnsh(str(out/'frames.tf')); s.furnsh(str(out/'rotation.bpc'))
for t in [-20,-17.125,-.000001,0,6.25,20]:
    refs['pck'].append(dict(t=t,matrix=s.pxform('J2000','ORACLE_PCK',t).tolist()))

(out/'time.tls').write_text("""KPL/LSK
\\begindata
DELTET/DELTA_T_A=32.184
DELTET/K=1.657D-3
DELTET/EB=1.671D-2
DELTET/M=(6.239996 1.99096871D-7)
DELTET/DELTA_AT=(10 @1972-JAN-1 11 @1972-JUL-1 12 @1973-JAN-1)
""")
s.furnsh(str(out/'time.tls'))
for label,mjd,seconds in [('1972 JAN 01 00:00:00 UTC',41317,0),('1972 JUN 30 23:59:59 UTC',41498,86399),
                          ('1972 JUN 30 23:59:60 UTC',41498,86400),('1972 JUL 01 00:00:00 UTC',41499,0),
                          ('1973 JAN 01 00:00:00 UTC',41683,0)]:
    refs['utc'].append(dict(label=label,mjd=mjd,seconds=seconds,tt=s.unitim(s.str2et(label),'ET','TT')))
(out/'clock.tsc').write_text("""KPL/SCLK
\\begindata
SCLK_DATA_TYPE_77=1
SCLK01_N_FIELDS_77=2
SCLK01_MODULI_77=(1000 10)
SCLK01_OFFSETS_77=(0 0)
SCLK01_OUTPUT_DELIM_77=2
SCLK01_TIME_SYSTEM_77=1
SCLK_PARTITION_START_77=(0 0)
SCLK_PARTITION_END_77=(1000 1000)
SCLK01_COEFFICIENTS_77=(0 0 20 1000 2000 10)
""")
s.furnsh(str(out/'clock.tsc'))
for tick in [0,17.125,999.5,1000,1000.5,1750,2000]:
    refs['sclk'].append(dict(tick=tick,tdb=s.sct2e(-77,tick)))
refs['assets']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.suffix in ['.bsp','.bc','.bpc','.tf','.tls','.tsc']}
(out/'reference.json').write_text(json.dumps(refs,indent=2)+'\n')
s.kclear()
print('Wrote pinned oracle fixtures:',out)

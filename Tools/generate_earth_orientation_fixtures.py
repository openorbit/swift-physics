import ctypes,json
from pathlib import Path
import argparse
parser = argparse.ArgumentParser(description="Generate reference matrices using an external ERFA 2.0.1 library.")
parser.add_argument("library", help="Path to an external ERFA shared library exporting eraPnm06a, eraC2ixys, etc.")
args = parser.parse_args()
lib=ctypes.CDLL(args.library)
Matrix = (ctypes.c_double * 3) * 3
for name in ["eraPnm06a", "eraPom00"]:
 getattr(lib,name).restype = None
lib.eraPnm06a.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.POINTER(ctypes.c_double*3)]
lib.eraPom00.argtypes = [ctypes.c_double]*3+[ctypes.POINTER(ctypes.c_double*3)]
lib.eraBpn2xy.argtypes = [ctypes.POINTER(ctypes.c_double*3)]+[ctypes.POINTER(ctypes.c_double)]*2
for name in ["eraS06","eraEra00","eraSp00","eraEors"]:
 getattr(lib,name).restype = ctypes.c_double
lib.eraS06.argtypes = [ctypes.c_double]*4
lib.eraEra00.argtypes = lib.eraSp00.argtypes = [ctypes.c_double]*2
lib.eraEors.argtypes = [ctypes.POINTER(ctypes.c_double*3),ctypes.c_double]
lib.eraC2ixys.argtypes = [ctypes.c_double]*3+[ctypes.POINTER(ctypes.c_double*3)]
lib.eraRz.argtypes = [ctypes.c_double,ctypes.POINTER(ctypes.c_double*3)]
def reference(edge, tt1, tt2, ut11, ut12, xp, yp):
 r,b = Matrix(),Matrix()
 if edge == 2:
  lib.eraPom00(xp,yp,lib.eraSp00(tt1,tt2),r)
 else:
  lib.eraPnm06a(tt1,tt2,b)
  x,y = ctypes.c_double(),ctypes.c_double()
  lib.eraBpn2xy(b,ctypes.byref(x),ctypes.byref(y))
  s = lib.eraS06(tt1,tt2,x.value,y.value)
  if edge == 0: lib.eraC2ixys(x.value,y.value,s,r)
  elif edge == 3: r = b
  else:
   angle = lib.eraEra00(ut11,ut12)
   if edge == 4: angle -= lib.eraEors(b,s)
   lib.eraIr(r);lib.eraRz(angle,r)
 return [r[i][j] for i in range(3) for j in range(3)]
rows=[]
for day in [-365250,-73050,-36525,-1,0,1,1234,36525,73050,365250]:
 for fraction in [0.000001,0.499999,0.999999]:
  for edge in range(5):
   matrix=reference(edge,2451545+day,fraction,2451545+day,fraction-64.184/86400,2.55060238e-7,1.860359247e-6)
   rows.append(dict(day=2451545+day,fraction=fraction,edge=edge,matrix=list(matrix)))
Path('Tests/CelMekTests/Fixtures/earth-orientation-erfa.json').write_text(json.dumps(rows,indent=2)+'\n')

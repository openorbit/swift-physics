"""Generate fixtures and an empirical survey against an external ERFA 2.0.1.
The shared library is a developer-only reference, never a package dependency.
Run: python3 Tools/characterize_tdb_adapter.py /path/to/erfa/shared/library
"""
import argparse,ctypes,datetime,hashlib,json,math
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('library')
args=parser.parse_args()
lib=ctypes.CDLL(args.library)
reference=lib.eraDtdb
reference.argtypes=[ctypes.c_double]*6
reference.restype=ctypes.c_double
def ref(day):return reference(2451545.,day,0.,0.,0.,0.)
def approximate(day):
 g=math.fmod(357.53+0.9856003*day,360)*math.pi/180
 return 0.001658*math.sin(g)+0.000014*math.sin(2*g)
def days(year):return (datetime.datetime(year,1,1)-datetime.datetime(2000,1,1,12)).total_seconds()/86400
surveys=[]
fixtures=[]
for start,end in [(1950,2050),(1900,2100)]:
 count=0;squares=0.;maximum=-1;worst=None
 first,last=days(start),days(end)
 for tick in range(int((last-first)*4)+1):
  day=first+tick/4
  full=ref(day);error=abs(approximate(day)-full)
  count+=1;squares+=error*error
  if error>maximum:maximum=error;worst=day+2451545
  # Every 30 days at all four phases, plus endpoints and the maximum below.
  if start==1950 and (tick%120<4 or day==last):
   h=300/86400
   fixtures.append(dict(daysSinceJ2000=day,tdbMinusTT=full,
    rate=(ref(day+h)-ref(day-h))/600,
    acceleration=(ref(day+h)-2*full+ref(day-h))/(300*300)))
 surveys.append(dict(start=f'{start}-01-01 TT',end=f'{end}-01-01 TT',stepSeconds=21600,
  samples=count,maximumSeconds=maximum,rmsSeconds=math.sqrt(squares/count),worstJulianDateTT=worst))
fixtures.append(dict(daysSinceJ2000=0.,tdbMinusTT=ref(0.),rate=(ref(300/86400)-ref(-300/86400))/600,
 acceleration=(ref(300/86400)-2*ref(0.)+ref(-300/86400))/(300*300)))
rateErrors=[];accelerationErrors=[]
for row in fixtures:
 g=math.fmod(357.53+0.9856003*row['daysSinceJ2000'],360)*math.pi/180
 k=0.9856003*math.pi/(180*86400)
 rateErrors.append(abs((0.001658*math.cos(g)+0.000028*math.cos(2*g))*k-row['rate']))
 accelerationErrors.append(abs(-(0.001658*math.sin(g)+0.000056*math.sin(2*g))*k*k-row['acceleration']))
Path('Tests/CelMekTests/Fixtures/tdb-erfa.json').write_text(json.dumps(fixtures,separators=(',',':'))+'\n')
Path('TDB_CHARACTERIZATION.json').write_text(json.dumps(dict(reference='ERFA 2.0.1 eraDtdb, geocentric (u=v=0), TT supplied as permitted by reference documentation',
 generator='Tools/characterize_tdb_adapter.py',referenceLibrarySHA256=hashlib.sha256(Path(args.library).read_bytes()).hexdigest(),
 fixtureCount=len(fixtures),derivativeReferenceStepSeconds=300,
 maximumFixtureRateError=max(rateErrors),maximumFixtureAccelerationError=max(accelerationErrors),surveys=surveys),indent=2)+'\n')
print(json.dumps(surveys,indent=2));print('Fixtures:',len(fixtures))

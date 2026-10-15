"""Regenerate fixtures using NOAA's independent reference Fortran (developer-only).
Usage: python3 generate_igrf_fixtures.py /path/to/igrf14.f
Source: https://www.ngdc.noaa.gov/IAGA/vmod/igrf14.f
Requires gfortran. The reference implementation is never linked into the package.
"""
from pathlib import Path
import subprocess
import sys
import tempfile

source = Path(sys.argv[1]).read_text()
source = source[source.lower().index("      subroutine igrf14syn"):]
driver = """program oracle
implicit none
real(8) :: year, height, latitude, longitude, x, y, z, f
integer :: status
 do
 read(*,*,iostat=status) year,height,latitude,longitude
 if (status /= 0) exit
 call igrf14syn(0,year,1,height,90.0d0-latitude,longitude,x,y,z,f)
 write(*,'(8(f20.8,1x))') year,height,latitude,longitude,x,y,z,f
 end do
end program
"""
rows = [(year, height, latitude, longitude)
        for year in [1900, 1945, 1997.5, 2000, 2020, 2025, 2027.5, 2030]
        for height, latitude, longitude in
        [(0, 0, 0), (0, 80, 120), (100, -80, 240), (0, 90, 0), (0, -90, 90), (850, 42, -71)]]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    (root / "synthesis.f").write_text(source)
    (root / "driver.f90").write_text(driver)
    executable = root / "oracle"
    subprocess.run(["gfortran", "-std=legacy", str(root / "synthesis.f"),
                    str(root / "driver.f90"), "-o", str(executable)], check=True)
    output = subprocess.run([str(executable)], input="".join(
        " ".join(map(str, row)) + "\n" for row in rows), text=True, capture_output=True, check=True).stdout
fixture = Path(__file__).resolve().parents[2] / "Tests/MagneticFieldsTests/Fixtures/IGRF14_REFERENCE.txt"
fixture.write_text("# year WGS84-height-km latitude-deg longitude-deg north-nT east-nT down-nT total-nT\n" + output)

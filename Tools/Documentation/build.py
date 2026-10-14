"""Build five module DocC archives and merge a package archive.
Usage: python3 Tools/Documentation/build.py [--disable-sandbox]
Requires a Swift toolchain with DocC merge support (Xcode 16 or newer).
"""
import argparse
from pathlib import Path
import shutil
import shlex
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--disable-sandbox", action="store_true", help="Disable SwiftPM's nested plugin sandbox when running inside an existing sandbox.")
parser.add_argument("--output", type=Path, help="Output directory (default: .build/documentation).")
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
output = (args.output or root / ".build/documentation").resolve()
output.mkdir(parents=True, exist_ok=True)
docc = shutil.which("docc")
if not docc:
    docc = subprocess.check_output(["xcrun", "--find", "docc"], text=True).strip()
help_text = subprocess.check_output([docc, "convert", "--help"], text=True)
# Older plugin releases inject --index, which newer DocC toolchains removed.
index_flags = ["--disable-indexing", "--emit-lmdb-index"] if "--emit-lmdb-index" in help_text else []
modules = ["CelMek", "Math", "Aerodynamics", "Atmospheres", "MagneticFields"]
archives = []
for module in modules:
    archive = output / f"{module}.doccarchive"
    command = ["swift", "package"]
    if args.disable_sandbox:
        command.append("--disable-sandbox")
    # The command plugin needs explicit permission for a custom output directory.
    command += ["--allow-writing-to-directory", str(output), "generate-documentation",
                "--target", module, "--output-path", str(archive), "--warnings-as-errors", *index_flags]
    print(f"Building {module} documentation", flush=True)
    with (output / f"{module}.log").open("w") as log:
        result = subprocess.run(command, cwd=root, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise SystemExit(f"{module} documentation failed; see {output / (module + '.log')}")
    archives.append(str(archive))
combined = output / "Physics.doccarchive"
subprocess.run([docc, "merge", *archives, "--output-path", str(combined),
                "--synthesized-landing-page-name", "Physics", "--synthesized-landing-page-kind", "Package"], check=True, cwd=root)
print(f"Combined documentation: {combined}", flush=True)

print(f"Browser preview: python3 Tools/Documentation/preview.py --archive {shlex.quote(str(combined))}", flush=True)

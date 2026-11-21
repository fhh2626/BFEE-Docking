"""Install and exercise the Linux wheel without Python package dependencies."""

from email.parser import BytesParser
from hashlib import sha256
import os
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile


# The upstream bin/vina is a QuickVina binary that reports the same version
# string as the real classic Vina. Verify the selected artifacts by content.
ENGINES = {
    "qvina2": "bc1f908869181d17b48c9a9b4c2a28cce5264ccf3053e15bce32bb5baf83b43d",
    "qvinaw": "bc0ad6d2b91e61d4bb878ef13fe524e55642bf264d48da75605dd31ff56f6cd4",
    "vina": "5986ea839e306a0130a5816a1e8359d2b6c05e0654fa32aad5e963c5a21166b7",
}
SPLITTER_SHA256 = "4d43eb3a70ae41412e036d784ac3e77c42c8a9dbe0010e495293908475ad227e"
WHEEL_TAG = "py3-none-manylinux2014_x86_64"


def check_wheel(wheel):
    if not wheel.name.endswith(f"-{WHEEL_TAG}.whl"):
        raise RuntimeError(f"Wheel filename must use the {WHEEL_TAG} tag: {wheel.name}")
    with ZipFile(wheel) as archive:
        names = archive.namelist()
        wheel_files = [name for name in names if name.endswith(".dist-info/WHEEL") and name.count("/") == 1]
        if len(wheel_files) != 1:
            raise RuntimeError("Expected exactly one .dist-info/WHEEL metadata file")
        metadata = BytesParser().parsebytes(archive.read(wheel_files[0]))
        if metadata.get_all("Tag") != [WHEEL_TAG]:
            raise RuntimeError(f"WHEEL metadata must use the {WHEEL_TAG} tag")
        if [value.strip().lower() for value in metadata.get_all("Root-Is-Purelib", [])] != ["false"]:
            raise RuntimeError("WHEEL metadata must mark this wheel as non-pure")
        if any(".data/purelib/" in name for name in names):
            raise RuntimeError("Platform wheel contains files under .data/purelib")

    env = os.environ.copy()
    for key in ("LD_LIBRARY_PATH", "LD_PRELOAD", "LD_AUDIT"):
        env.pop(key, None)
    env["LC_ALL"] = "C"

    with TemporaryDirectory(prefix="bfee-wheel-") as tmp:
        work = Path(tmp)
        site = work / "site"
        env["HOME"] = str(work)
        env["XDG_CACHE_HOME"] = str(work / "cache")
        env["PIP_CONFIG_FILE"] = os.devnull
        subprocess.run(
            [sys.executable, "-m", "pip", "--isolated", "install", "--no-index", "--no-deps", "--target", site, wheel],
            check=True, env=env, timeout=120,
        )
        tools = site / "bfee_docking" / "third_party"
        qvina = tools / "qvina"
        if any(qvina.glob("*.so*")) or any((tools / "lib").glob("*.so*")):
            raise RuntimeError("Unexpected QVina shared libraries in the installed wheel")
        for name, digest in ENGINES.items():
            for filename, expected in ((name, digest), (f"{name}_split", SPLITTER_SHA256)):
                artifact = qvina / filename
                if not artifact.is_file() or sha256(artifact.read_bytes()).hexdigest() != expected:
                    raise RuntimeError(f"Missing or unexpected QVina binary: qvina/{filename}")

        def run(*args):
            result = subprocess.run(args, cwd=work, env=env, timeout=30, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if result.returncode:
                raise RuntimeError(f"Command failed ({result.returncode}): {args}\n{result.stdout}")
            return result.stdout

        executables = [tools / "smina" / "smina", tools / "vina" / "vina"]
        executables.extend(qvina / f"{name}{suffix}" for name in ENGINES for suffix in ("", "_split"))
        for executable in executables:
            elf = run("readelf", "--program-headers", "--dynamic", "--wide", executable)
            if "(NEEDED)" in elf or re.search(r"^\s*INTERP\s", elf, re.MULTILINE):
                raise RuntimeError(f"{executable} is not statically linked")
            run(executable, "--help")
            print(f"{executable.relative_to(tools)}: static ELF and --help OK")

        def atom(i, x, y, z):
            return f"ATOM  {i:5d}  C   UNK A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00     0.000 C\n"

        (work / "receptor.pdbqt").write_text(atom(1, 2, 0, 0) + atom(2, -2, 0, 0) + atom(3, 0, 3, 0), encoding="ascii")
        (work / "ligand.pdbqt").write_text("ROOT\n" + atom(1, 0, 0, 0) + atom(2, 1.4, 0, 0) + "ENDROOT\nTORSDOF 0\n", encoding="ascii")
        for name in ENGINES:
            result = work / f"{name}-out.pdbqt"
            run(qvina / name, "--receptor", "receptor.pdbqt", "--ligand", "ligand.pdbqt",
                "--center_x", "0", "--center_y", "0", "--center_z", "0",
                "--size_x", "10", "--size_y", "10", "--size_z", "10",
                "--cpu", "1", "--seed", "42", "--num_modes", "1", "--exhaustiveness", "1",
                "--energy_range", "3", "--out", result)
            if "REMARK VINA RESULT:" not in result.read_text(encoding="ascii"):
                raise RuntimeError(f"{name} did not produce a scored docking pose")
            run(qvina / f"{name}_split", "--input", result)
            pose = work / f"{result.stem}_ligand_1.pdbqt"
            if not pose.is_file() or not any(
                line.startswith(("ATOM  ", "HETATM")) for line in pose.read_text(encoding="ascii").splitlines()
            ):
                raise RuntimeError(f"{name}_split did not produce a ligand pose")
            print(f"{name}: docking and splitting OK")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python scripts/check_linux_wheel.py <linux-wheel>")
    check_wheel(Path(sys.argv[1]).resolve(strict=True))

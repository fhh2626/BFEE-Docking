# PyPI packaging

This project publishes platform wheels only. Do not upload an sdist or a
`py3-none-any` wheel if macOS users should be prevented from installing the
package.

## Windows wheel

Run from the repository root in the build environment:

```powershell
micromamba activate auto_ff
python -m pip install --upgrade build wheel twine
.\scripts\build_windows_wheel.ps1
```

The output is written to:

```text
dist/windows/bfee_docking-<version>-py3-none-win_amd64.whl
```

The Windows build script removes the Linux `vina` and `smina` binaries from its
temporary build copy before creating the wheel.

## Linux wheel

Run from the repository root on Linux:

```bash
micromamba activate auto_ff
python -m pip install --upgrade build wheel twine
bash scripts/build_linux_wheel.sh
```

The output is written to:

```text
dist/linux/bfee_docking-<version>-py3-none-manylinux2014_x86_64.whl
```

The Linux build script removes Windows executables, DLLs, batch files, jars,
and the bundled Windows Open Babel directory from its temporary build copy
before creating the wheel.

## Check and upload

Check both wheels:

```bash
python -m twine check dist/windows/*.whl dist/linux/*.whl
```

Upload to TestPyPI:

```bash
python -m twine upload --repository testpypi dist/windows/*.whl dist/linux/*.whl
```

Upload to PyPI:

```bash
python -m twine upload dist/windows/*.whl dist/linux/*.whl
```

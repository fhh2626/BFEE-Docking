# conda-forge packaging

This repository includes a conda-forge recipe in `recipe/`.

## Strategy

- `qvina`, `pymol-open-source`, `openbabel`, `rdkit`, and the Python scientific
  stack are declared as conda-forge runtime dependencies.
- Bundled `vina` and `smina` are kept inside `bfee_docking/third_party`.
- macOS is skipped.
- Windows and Linux builds prune the other platform's bundled binaries before
  installation.

## Files

```text
recipe/meta.yaml
recipe/build.sh
recipe/bld.bat
```

## Local test build

For local testing before submitting to `conda-forge/staged-recipes`, either
temporarily change `source:` in `recipe/meta.yaml` to:

```yaml
source:
  path: ..
```

or create a GitHub tag and fill in the real release URL and `sha256`.

Then build:

```bash
conda build recipe -c conda-forge
```

On Windows:

```powershell
conda build recipe -c conda-forge
```

## staged-recipes submission

For the conda-forge PR, keep `source.url` and replace `PUT_SHA256_HERE` with the
SHA256 of the GitHub tag archive:

```bash
curl -L -o bfee-docking.tar.gz https://github.com/fhh2626/BFEE-Docking/archive/refs/tags/v1.0.0rc7.tar.gz
sha256sum bfee-docking.tar.gz
```

Copy `recipe/` into:

```text
staged-recipes/recipes/bfee-docking/
```

and open the PR.

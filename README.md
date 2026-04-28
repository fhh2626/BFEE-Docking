# BFEE-Docking

Author: Haohao Fu (fhh2626@nankai.edu.cn, fhh2626@gmail.com)

## Introduction

BFEE-Docking is an open-source companion platform to BFEE3 for automated structure-based drug discovery workflows. It integrates virtual screening, flexible redocking, MD preparation, MM-GBSA calculations, and BFEE3-based absolute binding free-energy calculations into a single graphical workflow.

Starting from a standard protein PDB file and a ligand database in SMI/SDF/MOL2 format, BFEE-Docking can automate many routine preparation steps, including protein cleaning, pH-dependent protonation-state assignment, ligand preprocessing, binding-site definition, docking/redocking, file-format conversion, and preparation of downstream MD or free-energy calculation files.

BFEE-Docking also interfaces with PyMOL, allowing users to visualize binding pockets, docking boxes, docked ligand conformations, and flexible residues during the workflow.

Main features:
- Automated protein and ligand preprocessing
- pH-dependent protonation-state handling
- Binding-site definition by blind docking, pocket prediction, co-crystallized ligand, or manual input
- Rigid docking and flexible redocking
- Support for DSDP, smina, AutoDock Vina, QVina2, and QVinaW
- PyMOL-based visualization of intermediate results
- Preparation of MD-ready files for OpenMM simulations
- Support for MM-GBSA calculations
- Connection to BFEE3 for absolute binding free-energy calculations

## License

- Main code: GPL-3.0 or later

This distribution bundles third-party components. Each component remains under its own license. They are distributed as separate works by mere aggregation and are not relicensed by this project.

The bundled third-party programs, except the MIT-licensed pocketeer, are used as independent executables and are invoked via the command line. They are not linked into the main program, and no third-party source code is incorporated into the main codebase.

All third-party license texts are included in the `licenses/` folder.

Bundled third-party components:
- OpenBabel 3.1.1 (binary): GPL-2.0
- smina 2020.12.10 (binary): GPL-2.0
- AutoDock Vina 1.2.7 (binary): Apache-2.0
- DSDP-redocking (source code): Apache-2.0
- pocketeer (source code): MIT

Their source code corresponding to the bundled version:
- OpenBabel 3.1.1: https://github.com/openbabel/openbabel/releases/tag/openbabel-3-1-1
- smina 2020.12.10: https://github.com/jaimergp/smina/releases/tag/2020.12.10
- AutoDock Vina 1.2.7: https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7
- DSDP-redocking: https://github.com/PKUGaoGroup/DSDP
- pocketeer: https://github.com/cch1999/pocketeer

## Installation and Usage

```bash
# optional
conda create -c conda-forge -n bfee_dock
conda activate bfee_dock

conda install -c conda-forge biotite biopython dimorphite-dl mdanalysis numba pdb2pqr pdb-tools pyside6 pymol-open-source qvina rdkit scipy

python main_gui.py
```

On Linux platforms, BFEE-Docking will prefer `obabel` from the current environment if available, and fall back to the bundled `third_party/obabel/obabel` otherwise.

Hence, one usually needs to install OpenBabel:

```bash
sudo apt install openbabel
```

One also needs to give execution permission to the smina and vina binaries:

```bash
chmod +x ./third_party/vina/vina
chmod +x ./third_party/smina/smina
```

For servers running GPU-based docking, unzip `third_party/DSDP_redocking.zip` and enter the DSDP folder. Then type:

```bash
make
```

For servers running MD simulations:

```bash
# optional
conda create -c conda-forge -n openmm python=3.13
conda activate openmm

conda install -c conda-forge numpy=2.2 openmm parmed mdanalysis cuda-version=12
```

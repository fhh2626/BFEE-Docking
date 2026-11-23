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
- Support for DSDP, smina, AutoDock Vina, QVina2, QVinaW and GNINA
- PyMOL-based visualization of intermediate results
- Preparation of MD-ready files for OpenMM simulations
- Support for MM-GBSA calculations
- Connection to BFEE3 for absolute binding free-energy calculations

## Citation

If you use BFEE-Docking in your research, please cite the following article.

Sun, X., Li, Z., Shao, X., Cai, W., Li, Z., and Fu, H. “BFEE-Docking: A User-Friendly and Customizable End-to-End Tool from High-Throughput Virtual Screening to Binding Free-Energy Calculations.” *Journal of Chemical Theory and Computation* **2026**, *22*, 6976–6985. [doi:10.1021/acs.jctc.6c00830](https://doi.org/10.1021/acs.jctc.6c00830).

## License

- Main code: GPL-3.0 or later

This distribution bundles third-party components. Each component remains under its own license. They are distributed as separate works by mere aggregation and are not relicensed by this project.

The bundled third-party programs, except the MIT-licensed pocketeer, are used as independent executables and are invoked via the command line. They are not linked into the main program, and no third-party source code is incorporated into the main codebase.

All third-party license texts are included in the `licenses/` folder.

Bundled third-party components:
- OpenBabel 3.1.1 (binary): GPL-2.0
- smina 2020.12.10 (binary): Apache-2.0
- AutoDock Vina 1.2.7 (binary): Apache-2.0
- QuickVina-W/2 (binary): Apache-2.0
- DSDP-redocking (source code): Apache-2.0
- pocketeer (source code): MIT

Their source code corresponding to the bundled version:
- OpenBabel 3.1.1: https://github.com/openbabel/openbabel/releases/tag/openbabel-3-1-1
- smina 2020.12.10: https://github.com/jaimergp/smina/releases/tag/2020.12.10
- AutoDock Vina 1.2.7: https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7
- QuickVina-W/2 and classic AutoDock Vina 1.1.2: https://github.com/QVina/qvina/tree/f4bb3b1073a0d50bb2f1fdd14d38594f937602ee
- DSDP-redocking: https://github.com/PKUGaoGroup/DSDP
- pocketeer: https://github.com/cch1999/pocketeer

## Installation and Usage

BFEE-Docking requires Python 3.11 or later.

```bash
pip install bfee-docking
```

### Ligand preparation for MD

BFEE-Docking uses the CGenFF 5.0 data bundled with [AutoCGENFF](https://github.com/fhh2626/AutoCGENFF) 20261123.1 or later. No force-field directory setup is needed.

After exporting a docking pose, click **Process Ligand**. BFEE-Docking prepares the ligand PDB and MOL2, runs AutoCGENFF locally, and fills **Ligand str File** with the absolute path to `md_prepper/fixed_files/ligand.str`. Review the assigned parameters before running simulations. You can edit the path or browse to a custom CHARMM `.str` or `.top` file. Generate PSF uses the file currently selected in the input box. A GROMACS `.top` file is not a CHARMM topology and cannot be used there.

### Optional tools for docking and simulation

For servers running GPU-based docking, unzip `third_party/DSDP_redocking.zip` and enter the DSDP folder. Then type:

```bash
make
```

For servers running MD simulations:

```bash
pip install parmed mdanalysis
pip install "openmm[cuda13]>=8.5"
```

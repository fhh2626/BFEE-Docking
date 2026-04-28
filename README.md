BFEE-Docking
Author: Haohao Fu (fhh2626@nankai.edu.cn, fhh2626@gmail.com)

License
- Main code: GPL-3.0 or later
This distribution bundles third-party components. Each component remains under its own license. They are distributed as separate works (mere aggregation) and are not relicensed by this project.
The bundled third-party programs, except MIT-licensed pocketeer, are used as independent executables and are invoked via the command line. They are not linked into the main program and no third-party source code is incorporated into the main codebase.
All third-party license texts are included in the licenses/ folder.

Bundled third-party components:
- OpenBabel 3.1.1 (binary): GPL-2.0
- smina 2020.12.10 (binary): GPL-2.0
- AutoDock Vina 1.2.7 (binary): Apache-2.0
- DSDP-redocking (source code): Apache-2.0
- pocketeer (source code): MIT

Their source code corresponding to the boundled version:
- OpenBabel 3.1.1: https://github.com/openbabel/openbabel/releases/tag/openbabel-3-1-1 
- smina 2020.12.10: https://github.com/jaimergp/smina/releases/tag/2020.12.10
- AutoDock Vina 1.2.7: https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7
- DSDP-redocking: https://github.com/PKUGaoGroup/DSDP
- pocketeer: https://github.com/cch1999/pocketeer

Installation and Usage:
(optional) conda create -c conda-forge -n bfee_dock
(optional) conda activate bfee_dock
conda install -c conda-forge biotite biopython dimorphite-dl mdanalysis numba pdb2pqr pdb-tools pyside6 pymol-open-source qvina rdkit scipy
python main_gui.py

On Linux platforms, BFEE-Docking will prefer `obabel` from the current environment if available, and fall back to the bundled `third_party/obabel/obabel` otherwise.
Hence, one usually need to:
sudo apt install openbabel

One also need to give the permission to smina and vina binaries:
chmod +x ./third_party/vina/vina
chmod +x ./third_party/smina/smina

For servers running GPU-based docking:
just unzip third_party/DSDP_redocking.zip and cd to the DSDP folder. Then type
make

For servers running MD simulations:
(optional) conda create -c conda-forge -n openmm python=3.13
(optional) conda activate openmm
conda install -c conda-forge numpy=2.2 openmm parmed mdanalysis cuda-version=12

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Haohao Fu <fhh2626@nankai.edu.cn, fhh2626@gmail.com>

# Standard library imports
import collections
import collections.abc
import os
import pathlib
import shutil
import string
import subprocess
import sys
import textwrap

from . import third_party_tools

# Add third_party directory to Python path for bundled packages (e.g., pocketeer)
third_party_tools.add_third_party_to_path()

# Module-level reference for backward compatibility
_SCRIPT_DIR = third_party_tools.get_third_party_dir().parent

# Third-party library imports
import Bio.PDB
import MDAnalysis
import numpy as np
import pocketeer

class PDBParser:
    """Manage paired input/output PDB paths and convenience helpers."""

    # Common metal ion residue names in PDB files
    _METAL_ION_RESNAMES = {
        "LI", "NA", "K", "RB", "CS",           # Alkali metals
        "MG", "CA", "CAL", "SR", "BA",         # Alkaline earth metals (CAL = CHARMM calcium)
        "MN", "FE", "CO", "NI", "CU", "ZN",    # Transition metals (row 1)
        "CD", "HG",                            # Transition metals (row 2-3)
        "AL", "GA", "IN", "SN", "PB", "BI",    # Post-transition metals
        "CR", "MO", "W", "V",                  # Other transition metals
        "ZN2", "FE2", "FE3", "MN2", "CU1", "CU2",  # Common ion states
        "CA2", "MG2", "NA1", "K1", "POT", "SOD", "CLA", "CL",  # Force field variants
    }

    _WATER_RESNAMES = {"HOH", "H2O", "WAT", "WT", "SOL", "TIP3", "TP3"}

    def __init__(
        self,
        input_pdb: str,
        output_pdb: str,
        pdb_id: str | None = None,
        unit_structure: bool = False,
    ) -> None:
        """Store both the source (input) and generated (output) paths."""
        self._input_pdb = pathlib.Path(input_pdb).resolve()
        self._output_pdb = pathlib.Path(output_pdb).resolve()
        self._docking_region: list[np.ndarray] | None = None
        self._hetatm_labels: list[str] = []
        self._hetatm_selections: list[str] = []
        self._generated_pdbqt_file: pathlib.Path | None = None
        self._detected_pockets: list | None = None
        
        # Get absolute path to obabel executable
        self._obabel_executable = third_party_tools.get_obabel_executable()

        if pdb_id is not None:
            self._input_pdb.parent.mkdir(parents=True, exist_ok=True)
            self._fetch_pdb(pdb_id, unit_structure)

        self._collect_hetatm_molecules()

    def _fetch_pdb(self, pdb_id: str, unit_structure: bool) -> None:
        """Run pdb_fetch and write its stdout to the input PDB file."""
        cmd = ["pdb_fetch"]
        if unit_structure:
            cmd.append("-biounit")
        cmd.append(pdb_id)
        
        try:
            completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
            self._input_pdb.write_text(completed.stdout, encoding="utf-8")
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to fetch PDB '{pdb_id}'"
            if e.stderr:
                error_msg += f"\nError details: {e.stderr.strip()}"
            raise RuntimeError(error_msg) from e
        except FileNotFoundError:
            raise FileNotFoundError(
                "pdb_fetch command not found. Please ensure pdb-tools is installed."
            )
            

    def renew_output(self) -> None:
        """Overwrite the output PDB by copying the input file in place."""
        if not self._input_pdb.exists():
            raise FileNotFoundError(f"Input PDB not found: {self._input_pdb}")
        self._output_pdb.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(self._input_pdb, self._output_pdb)

    def _run_pdb_tool(
        self,
        tool_name: str,
        args: list[str],
        error_description: str
    ) -> None:
        """
        Run a pdb-tools command and write output to the output PDB file.
        
        This is a helper method to reduce code duplication for all pdb-tools
        commands that follow the pattern: tool [args] file.pdb > file.pdb
        
        Args:
            tool_name: Name of the pdb-tools command (e.g., 'pdb_selchain').
            args: Additional arguments for the command (not including the file path).
            error_description: Human-readable description of the operation for error messages.
        
        Raises:
            FileNotFoundError: If output PDB doesn't exist or the tool is not installed.
            RuntimeError: If the command fails.
        """
        if not self._output_pdb.exists():
            raise FileNotFoundError(f"Output PDB not found: {self._output_pdb}")
        
        # Build command: tool_name [args...] output_pdb
        cmd = [tool_name] + args + [str(self._output_pdb)]
        
        try:
            completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
            # Overwrite output file with command stdout
            self._output_pdb.write_text(completed.stdout, encoding="utf-8")
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to {error_description}"
            if e.stderr:
                error_msg += f"\nError details: {e.stderr.strip()}"
            raise RuntimeError(error_msg) from e
        except FileNotFoundError:
            raise FileNotFoundError(
                f"{tool_name} command not found. Please ensure pdb-tools is installed."
            )

    def _get_available_chains(self) -> set[str]:
        """
        Get the set of available chain IDs in the output PDB file.
        
        Returns:
            set[str]: Set of chain IDs found in the PDB. Empty set if no chains are found.
        
        Raises:
            FileNotFoundError: If output PDB file does not exist.
        """
        if not self._output_pdb.exists():
            raise FileNotFoundError(f"Output PDB not found: {self._output_pdb}")
        
        chains = set()
        with open(self._output_pdb, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')) and len(line) > 21:
                    chain_id = line[21].strip()
                    if chain_id:
                        chains.add(chain_id)
        return chains

    def _select_chains(self, chains: list[str] | str) -> None:
        """
        Select specific chains in the output PDB file.

        Uses `pdb_selchain` command: `pdb_selchain -A,C output.pdb > output.pdb`
        
        If the PDB file has no chain IDs, this step is skipped.
        If the requested chains don't exist but other chains do, all available chains are selected.
        """
        # Get available chains in the PDB file
        available_chains = self._get_available_chains()
        
        # If PDB has no chain IDs, skip chain selection
        if not available_chains:
            print("Warning: PDB file has no chain IDs. Skipping chain selection.")
            return
        
        # Parse requested chains
        if isinstance(chains, list):
            requested_chains = set(chains)
            chain_arg = ",".join(chains)
        else:
            requested_chains = set(c.strip() for c in chains.split(',') if c.strip())
            chain_arg = chains
        
        # Check if any requested chain exists in the PDB
        matching_chains = requested_chains & available_chains
        
        if not matching_chains:
            # No matching chains found - warn user and keep all chains (skip selection)
            print(f"Warning: Requested chain(s) '{chain_arg}' not found in PDB. "
                  f"Available chains: {', '.join(sorted(available_chains))}. "
                  f"Keeping all chains.")
            return
        
        # If only some chains match, use only the matching ones
        if matching_chains != requested_chains:
            missing = requested_chains - available_chains
            print(f"Warning: Chain(s) '{', '.join(sorted(missing))}' not found in PDB. "
                  f"Selecting only available chains: {', '.join(sorted(matching_chains))}")
            chain_arg = ",".join(sorted(matching_chains))

        self._run_pdb_tool(
            "pdb_selchain",
            [f"-{chain_arg}"],
            f"select chains '{chain_arg}'"
        )

    def _remove_useless_lines(self) -> None:
        """
        Remove useless lines from the output PDB file (keep only coordinates).

        Uses `pdb_keepcoord` command: `pdb_keepcoord output.pdb > output.pdb`
        """
        self._run_pdb_tool(
            "pdb_keepcoord",
            [],
            "remove useless lines from PDB"
        )

    def _remove_hetatm(self) -> None:
        """
        Remove HETATM records from the output PDB file.

        Uses `pdb_delhetatm` command: `pdb_delhetatm output.pdb > output.pdb`
        """
        self._run_pdb_tool(
            "pdb_delhetatm",
            [],
            "remove HETATM records from PDB"
        )

    def _save_ligands_to_mol2(self, output_path: pathlib.Path) -> None:
        """
        Extract non-metal and non-water HETATM from output PDB and save as mol2.
        
        This method selects atoms with HETATM record type that are not in the 
        predefined list of waters and metal ions, then uses Open Babel to
        convert them to a MOL2 file with guessed bonds and atom types.
        
        Atom names are ensured to be unique by renaming duplicates with
        element symbol + incrementing number format.
        """
        # Load the universe from current output PDB
        u = MDAnalysis.Universe(str(self._output_pdb))
        
        # Build exclusion list
        excluded = self._WATER_RESNAMES | self._METAL_ION_RESNAMES
        resname_sel = " ".join(excluded)
        
        # Select ligands: HETATM and not in excluded list
        ligands = u.select_atoms(f"record_type HETATM and not resname {resname_sel}")
        
        if ligands.n_atoms == 0:
            print("No non-metal/non-water HETATM found to save to mol2.")
            return

        # Ensure unique atom names to avoid issues in MOL2 format
        # Track used names and element counters
        used_names: set[str] = set()
        element_counters: dict[str, int] = {}
        
        for atom in ligands:
            original_name = atom.name.strip()
            
            if original_name not in used_names:
                # Name is unique, just add it to the set
                used_names.add(original_name)
            else:
                # Name is duplicated, generate a new unique name
                # Use element symbol + incrementing number
                element = atom.element.strip() if hasattr(atom, 'element') and atom.element else original_name[0]
                
                # Initialize counter for this element if not exists
                if element not in element_counters:
                    # Start from a number higher than typical atom names
                    element_counters[element] = 100
                
                # Find a unique name
                while True:
                    new_name = f"{element}{element_counters[element]}"
                    element_counters[element] += 1
                    # Ensure new name is max 4 characters (PDB atom name limit)
                    if len(new_name) > 4:
                        new_name = new_name[:4]
                    if new_name not in used_names:
                        break
                
                atom.name = new_name
                used_names.add(new_name)
        
        # Save to temporary PDB first to let obabel handle bond guessing
        temp_pdb = output_path.with_suffix(".temp_ligands.pdb")
        ligands.write(str(temp_pdb))
        
        try:
            # obabel -ipdb temp.pdb -omol2 -O output_path
            cmd = [
                str(self._obabel_executable),
                str(temp_pdb),
                "-O", str(output_path)
            ]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            # Post-process MOL2 to fix unsupported bond types for pdb2pqr
            # pdb2pqr doesn't support: am (amide), ar (aromatic), etc.
            # Replace them with standard bond type "1" (single bond)
            self._fix_mol2_bond_types(output_path)
            
            print(f"✓ Saved {ligands.n_atoms} ligand atoms to {output_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error converting ligands to mol2: {e.stderr.strip()}")
        finally:
            if temp_pdb.exists():
                temp_pdb.unlink()

    def _fix_mol2_bond_types(self, mol2_path: pathlib.Path) -> None:
        """
        Fix unsupported bond types in MOL2 file for pdb2pqr compatibility.
        
        pdb2pqr only supports these bond types: 1 (single), 2 (double), 3 (triple)
        This method replaces unsupported types:
          - am (amide) -> 1 (single)
          - ar (aromatic) -> ar (keep, usually supported as 1.5)
          - du (dummy) -> 1 (single)
          - un (unknown) -> 1 (single)
          - nc (not connected) -> 1 (single)
        """
        if not mol2_path.exists():
            return
        
        # Map of unsupported bond types to their replacements
        # pdb2pqr specifically doesn't support 'am'
        bond_type_map = {
            'am': '1',   # amide -> single
            'du': '1',   # dummy -> single
            'un': '1',   # unknown -> single
            'nc': '1',   # not connected -> single
        }
        
        content = mol2_path.read_text(encoding='utf-8')
        lines = content.splitlines()
        modified_lines = []
        in_bond_section = False
        
        for line in lines:
            # Check for section markers
            if line.startswith('@<TRIPOS>BOND'):
                in_bond_section = True
                modified_lines.append(line)
                continue
            elif line.startswith('@<TRIPOS>'):
                in_bond_section = False
                modified_lines.append(line)
                continue
            
            if in_bond_section and line.strip():
                # Bond line format: bond_id atom1_id atom2_id bond_type [status_bits]
                parts = line.split()
                if len(parts) >= 4:
                    bond_type = parts[3]
                    if bond_type in bond_type_map:
                        parts[3] = bond_type_map[bond_type]
                        # Reconstruct line preserving spacing approximately
                        line = '     '.join(parts[:4])
                        if len(parts) > 4:
                            line += '     ' + '     '.join(parts[4:])
            
            modified_lines.append(line)
        
        # Write back
        mol2_path.write_text('\n'.join(modified_lines) + '\n', encoding='utf-8')

    def _remove_hydrogens(self) -> None:
        """
        Remove all Hydrogen atoms from the output PDB file.

        Uses `pdb_delelem` command: `pdb_delelem -H output.pdb > output.pdb`
        """
        self._run_pdb_tool(
            "pdb_delelem",
            ["-H"],
            "remove hydrogen atoms from PDB"
        )

    def _select_altloc(self) -> None:
        """
        Select the highest occupancy alternative location for each atom.

        Uses `pdb_selaltloc` command: `pdb_selaltloc output.pdb > output.pdb`
        """
        self._run_pdb_tool(
            "pdb_selaltloc",
            [],
            "select altloc in PDB"
        )

    def add_chain_identifiers(self) -> None:
        """
        Add chain identifiers to the output PDB file.

        Uses `pdb_chainbows` command: `pdb_chainbows output.pdb > output.pdb`
        This assigns chain IDs (A, B, C, ...) to segments without chain identifiers.
        """
        self._run_pdb_tool(
            "pdb_chainbows",
            [],
            "add chain identifiers to PDB"
        )

    def has_chain_identifiers(self) -> bool:
        """
        Check if all protein residues in the output PDB have chain identifiers.
        
        Returns:
            bool: True if all protein residues have chain IDs, False otherwise.
        
        Raises:
            FileNotFoundError: If output PDB file does not exist.
        """
        if not self._output_pdb.exists():
            raise FileNotFoundError(f"Output PDB not found: {self._output_pdb}")
        
        u = MDAnalysis.Universe(str(self._output_pdb))
        all_protein_atoms = u.select_atoms("protein")
        
        for residue in all_protein_atoms.residues:
            first_atom = residue.atoms[0]
            chain = getattr(first_atom, "chainID", "").strip()
            if not chain:
                return False
        
        return True


    def _add_hydrogens_with_ph(
        self, 
        ph: float = 7.0, 
        ff: str = 'AMBER',
        ligand_mol2: pathlib.Path | None = None
    ) -> None:
        """
        Add hydrogens according to pH using pdb2pqr and convert back to PDB.

        1. pdb2pqr --with-ph=<ph> --ff=AMBER <output_pdb> <temp_pqr>
        2. Filter <temp_pqr> to <output_pdb> removing lines with "CHARGE" or "RADIUS".
        
        Args:
            ph: pH value for protonation. Default: 7.0
            ff: Force field to use. Default: 'AMBER'
            ligand_mol2: Path to ligand MOL2 file. If provided, pdb2pqr will use
                        --ligand to generate parameters for the small molecule.
        """
        if not self._output_pdb.exists():
            raise FileNotFoundError(f"Output PDB not found: {self._output_pdb}")

        temp_pqr = self._output_pdb.with_suffix(".pqr")

        # pdb2pqr --with-ph=7.0 --ff=AMBER out.pdb out.pqr
        # Use third_party_tools.run_pdb2pqr which has the propka monkey patch applied
        
        
        try:
            third_party_tools.run_pdb2pqr(
                input_pdb=str(self._output_pdb),
                output_pqr=str(temp_pqr),
                ph=ph,
                ff="CHARMM_WITH_ION",
                titration_state_method="propka",
                keep_chain=True,
                ligand_mol2=str(ligand_mol2) if ligand_mol2 else None
            )

            # Use BioPython to convert PQR to PDB
            if temp_pqr.exists():
                parser = Bio.PDB.PDBParser(QUIET=True)
                structure = parser.get_structure("temp", str(temp_pqr))
                io = Bio.PDB.PDBIO()
                io.set_structure(structure)
                io.save(str(self._output_pdb))
        except RuntimeError as e:
            print(f"Error running pdb2pqr: {e}")
            raise


    def _convert_to_pdbqt(self) -> None:
        """
        Convert the output PDB to PDBQT format using Open Babel.
        
        Uses obabel command to generate Vina-compatible PDBQT files for
        rigid receptors (proteins), preserving existing hydrogen atoms
        (added by pdb2pqr) without re-protonation.
        
        Note:
            This method only performs format conversion and does NOT add or 
            modify hydrogen atoms. Ensure hydrogens have been added beforehand
            using _add_hydrogens_with_ph().
            
            The -xr flag ensures the protein is treated as a rigid receptor,
            preventing the addition of ROOT/BRANCH tags that would cause
            Vina parsing errors.
        """
        if not self._output_pdb.exists():
            raise FileNotFoundError(f"Output PDB not found: {self._output_pdb}")

        output_pdbqt = self._output_pdb.with_suffix(".pdbqt")

        # Use Open Babel to convert PDB to PDBQT format (with absolute path)
        # -xr: Treat as rigid receptor (no ROOT/BRANCH/TORSDOF tags)
        # This preserves existing hydrogens from pdb2pqr without modification
        cmd = [
            str(self._obabel_executable),
            str(self._output_pdb),
            "-O", str(output_pdbqt),
            "-xr"  # Rigid receptor flag - critical for Vina compatibility
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running obabel: {e}")
            print(f"STDERR:\n{e.stderr}")
            raise
        
        # Store the generated PDBQT file path
        self._generated_pdbqt_file = output_pdbqt

    def process_protein_pipeline(
        self,
        ph: float = 7.0,
        chains: str | None = None,
        preserve_metal: bool = False,
        preserve_coord_water: bool | int = False,
        preserve_metal_coord_water: bool = False,
        ligand_mol2: pathlib.Path | None = None,
    ) -> None:
        """
        Execute the full protein processing pipeline.

        Sequence:
        1. Renew output from input
        2. Select chains (if specified)
        3. Select altloc
        4. Remove HETATM/Save non-metal and non-water HETATM to a mol2
        5. Remove existing hydrogens
        6. (Optional) Re-add metal ions if preserve_metal=True
        7. (Optional) Re-add metal-coordinating waters if preserve_metal_coord_water=True
        8. (Optional) Re-add hydrogen-bonded waters if preserve_coord_water=True/int
        9. Add hydrogens at specified pH
        10. Remove useless lines
        11. Convert to PDBQT
        
        Args:
            ph: pH value for protonation. Default: 7.0
            chains: Comma-separated chain IDs to select. Default: None (all chains)
            preserve_metal: If True, re-add metal ions from input PDB that are within
                          2.5Å of the processed protein. Default: False
            preserve_coord_water: If True or int, re-add water molecules that form
                                hydrogen bonds with the protein. If True, uses 3 as
                                the minimum number of hydrogen bonds. If int (2 or 3),
                                uses that value. Default: False
            preserve_metal_coord_water: If True, re-add water molecules that coordinate
                                       with metal ions (within 2.5Å). Requires
                                       preserve_metal=True. Default: False
            ligand_mol2: Path to the ligand mol2 file. Default: None
        
        Raises:
            ValueError: If preserve_metal_coord_water=True but preserve_metal=False.
        """
        # Validate parameters
        if preserve_metal_coord_water and not preserve_metal:
            raise ValueError(
                "preserve_metal_coord_water=True requires preserve_metal=True. "
                "Metal ions must be preserved first before their coordinating waters can be added."
            )
        
        print("Starting protein processing...")
        
        self.renew_output()
        print("✓ Output file renewed")
        
        if chains:
            self._select_chains(chains)
            print(f"✓ Selected chains: {chains}")
        
        self._select_altloc()
        print("✓ Selected highest occupancy altloc")
        
        # Save non-metal and non-water HETATM to a mol2 if requested
        if ligand_mol2:
            self._save_ligands_to_mol2(pathlib.Path(ligand_mol2))
        
        self._remove_hetatm()
        print("✓ Removed HETATM records")
        
        self._remove_hydrogens()
        print("✓ Removed old hydrogens")
        
        # Re-add metal ions if requested (before adding new hydrogens)
        if preserve_metal:
            n_ions = self.readd_ions(distance_cutoff=2.5)
            print(f"✓ Re-added {n_ions} metal ion(s)")
        
        # Re-add metal-coordinating waters if requested (requires preserve_metal)
        if preserve_metal_coord_water:
            n_coord_waters = self.readd_metal_coor_waters(distance_cutoff=2.5)
            print(f"✓ Re-added {n_coord_waters} metal-coordinating water(s)")
        
        # Re-add hydrogen-bonded waters if requested
        if preserve_coord_water:
            # If True, use default of 3 hydrogen bonds; if int, use that value
            if preserve_coord_water is True:
                min_hbonds = 3
            else:
                min_hbonds = int(preserve_coord_water)
                if min_hbonds not in (2, 3):
                    raise ValueError("preserve_coord_water must be True, False, 2, or 3")
            n_hbond_waters = self.readd_waters(hydrogens=min_hbonds)
            print(f"✓ Re-added {n_hbond_waters} hydrogen-bonded water(s) (>={min_hbonds} H-bonds)")
        
        self._add_hydrogens_with_ph(
            ph=ph, 
            ligand_mol2=pathlib.Path(ligand_mol2) if ligand_mol2 else None
        )

        print(f"✓ Added hydrogens at pH {ph}")
        
        self._remove_useless_lines()
        print("✓ Removed useless lines")
        
        self._convert_to_pdbqt()
        print("✓ Converted to PDBQT format")

    def set_blind_docking_range(self, margin: float | collections.abc.Sequence[float] = 0.0) -> None:
        """Calculate center and side lengths then store docking region with margin."""
        if not self._output_pdb.exists():
            raise FileNotFoundError(f"Output PDB not found: {self._output_pdb}")

        u = MDAnalysis.Universe(str(self._output_pdb))
        atoms = u.atoms
        if atoms is None:
            raise ValueError("No atoms found in output PDB.")
        positions = atoms.positions
        if positions.size == 0:
            raise ValueError("No atoms found in output PDB.")

        min_corner = positions.min(axis=0)
        max_corner = positions.max(axis=0)
        center = (min_corner + max_corner) / 2.0

        margin_arr = np.array(margin, dtype=float)
        if margin_arr.ndim == 0:
            margin_arr = np.full(3, margin_arr)
        if margin_arr.shape != (3,):
            raise ValueError("margin must be a scalar or a length-3 sequence.")

        side_length = (max_corner - min_corner) + margin_arr
        self._docking_region = [center, side_length]

    def set_known_ligand_range(
        self,
        hetatm_label: str,
        margin: float | collections.abc.Sequence[float] = 0.0
    ) -> None:
        """
        Calculate docking range based on a known ligand from the input PDB.

        1. Find selection string for the given hetatm_label.
        2. Select atoms using MDAnalysis on the input PDB.
        3. Calculate center and dimensions + margin.
        4. Store in self._docking_region.
        """
        try:
            index = self._hetatm_labels.index(hetatm_label)
        except ValueError:
            raise ValueError(
                f"Label '{hetatm_label}' not found. Available: {self._hetatm_labels}"
            )

        selection_str = self._hetatm_selections[index]

        if not self._input_pdb.exists():
            raise FileNotFoundError(f"Input PDB not found: {self._input_pdb}")

        u = MDAnalysis.Universe(str(self._input_pdb))
        ligand_atoms = u.select_atoms(selection_str)

        if ligand_atoms.n_atoms == 0:
            raise ValueError(f"No atoms found for selection: {selection_str}")

        positions = ligand_atoms.positions
        min_corner = positions.min(axis=0)
        max_corner = positions.max(axis=0)
        center = (min_corner + max_corner) / 2.0

        margin_arr = np.array(margin, dtype=float)
        if margin_arr.ndim == 0:
            margin_arr = np.full(3, margin_arr)
        if margin_arr.shape != (3,):
            raise ValueError("margin must be a scalar or a length-3 sequence.")

        side_length = (max_corner - min_corner) + margin_arr
        self._docking_region = [center, side_length]

    def set_custom_docking_range(
        self,
        center: collections.abc.Sequence[float] | np.ndarray,
        side_length: collections.abc.Sequence[float] | np.ndarray,
    ) -> None:
        """Directly set docking box center and side lengths (each length-3)."""
        center_arr = np.asarray(center, dtype=float)
        side_arr = np.asarray(side_length, dtype=float)

        if center_arr.shape != (3,):
            raise ValueError("center must be a length-3 sequence of coordinates.")
        if side_arr.shape != (3,):
            raise ValueError("side_length must be a length-3 sequence.")
        if np.any(side_arr <= 0):
            raise ValueError("side_length values must be positive.")

        self._docking_region = [center_arr, side_arr]

    def get_generated_pdbqt_path(self) -> pathlib.Path | None:
        """
        Get the full path of the generated PDBQT file.
        
        Returns:
            pathlib.Path | None: The absolute path of the generated PDBQT file as a Path object,
                                or None if convert_to_pdbqt() has not been called yet.
        """
        return self._generated_pdbqt_file

    def get_generated_docking_range(self) -> list[np.ndarray] | None:
        """
        Get the docking region that was previously set.
        
        Returns:
            list[np.ndarray] | None: A list containing [center, side_length] where:
                - center: numpy array of shape (3,) with x, y, z coordinates
                - side_length: numpy array of shape (3,) with dimensions in x, y, z
                Returns None if no docking range has been set yet.
        """
        return self._docking_region

    def _open_pdb_only(
        self,
        pdb_path: pathlib.Path,
    ) -> None:
        """Open the PDB in PyMOL showing only the protein."""
        if not pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")

        script_content = "set orthoscopic, 1\n"
        pml_path = pdb_path.with_name(f"{pdb_path.stem}_view.pml")
        pml_path.write_text(script_content, encoding="utf-8")

        cmd = ["pymol", str(pdb_path), str(pml_path)]
        try:
            subprocess.run(cmd, check=True)
        finally:
            if pml_path.exists():
                pml_path.unlink()

    def _open_pdb_with_box(
        self,
        pdb_path: pathlib.Path,
        line_width: float = 10.0,
    ) -> None:
        """
        Open the PDB in PyMOL and draw the docking box if available.
        Falls back to protein-only view when `_docking_region` is None.
        """
        if self._docking_region is None:
            self._open_pdb_only(pdb_path)
            return

        center, side_length = self._docking_region
        cx, cy, cz = (float(c) for c in center)
        lx, ly, lz = (float(s) for s in side_length)
        hx, hy, hz = lx / 2, ly / 2, lz / 2
        v = [
            (cx - hx, cy - hy, cz - hz), (cx + hx, cy - hy, cz - hz),
            (cx + hx, cy + hy, cz - hz), (cx - hx, cy + hy, cz - hz),
            (cx - hx, cy - hy, cz + hz), (cx + hx, cy - hy, cz + hz),
            (cx + hx, cy + hy, cz + hz), (cx - hx, cy + hy, cz + hz),
        ]
        script_content = "set orthoscopic, 1\n"
        script_content += textwrap.dedent(f"""
            python
            from pymol.cgo import *
            from pymol import cmd

            verts = {v}
            edges = [
                (0,1), (1,2), (2,3), (3,0),
                (4,5), (5,6), (6,7), (7,4),
                (0,4), (1,5), (2,6), (3,7)
            ]

            obj = [BEGIN, LINES, COLOR, 1.0, 0.0, 0.0]
            for (i, j) in edges:
                v1 = verts[i]
                v2 = verts[j]
                obj.extend([VERTEX, v1[0], v1[1], v1[2], VERTEX, v2[0], v2[1], v2[2]])
            obj.append(END)

            cmd.load_cgo(obj, 'box')
            cmd.set('cgo_line_width', {line_width}, 'box')
            python end
        """)
        pml_path = pdb_path.with_name(f"{pdb_path.stem}_view.pml")
        pml_path.write_text(script_content, encoding="utf-8")

        cmd = ["pymol", str(pdb_path), str(pml_path)]
        try:
            subprocess.run(cmd, check=True)
        finally:
            if pml_path.exists():
                pml_path.unlink()

    def _collect_hetatm_molecules(self) -> None:
        """
        Populate `_hetatm_labels` and `_hetatm_selections` with non-water HETATM groups.
        Groups are defined by (resname, chainid) or (resname, resid) when chainid is missing.
        Labels look like "(RESNAME:CHAIN)" or "(RESNAME:@RESID)" when no chain is available.
        """
        if not self._input_pdb.exists():
            raise FileNotFoundError(f"Input PDB not found: {self._input_pdb}")

        u = MDAnalysis.Universe(str(self._input_pdb))
        het_atoms = u.select_atoms("record_type HETATM")
        excluded = {"HOH", "H2O", "WAT", "WT"}
        
        # Use dict to track unique groups
        # Key: (resname, chain, resid_or_none)
        # - If chain exists: (resname, chain, None) - grouped by chain
        # - If no chain: (resname, "", resid) - grouped by resid
        groups: dict[tuple[str, str, int | None], None] = {}

        for atom in het_atoms:
            resname = atom.resname.strip()
            chain = (getattr(atom, "chainID", "") or getattr(atom, "segid", "") or "").strip()
            resid = getattr(atom, "resid", None)
            
            if resname.upper() in excluded:
                continue
            
            if chain:
                # Chain ID exists, group by (resname, chain)
                groups[(resname, chain, None)] = None
            else:
                # No chain ID, use resid to distinguish same-named ligands
                groups[(resname, "", resid)] = None

        self._hetatm_labels = []
        self._hetatm_selections = []
        for resname, chain, resid in sorted(groups):
            if chain:
                # Has chain ID
                label = f"({resname}:{chain})"
                selection = f"chainid {chain} and resname {resname}"
            elif resid is not None:
                # No chain, but has resid
                label = f"({resname}:@{resid})"
                selection = f"resname {resname} and resid {resid}"
            else:
                # No chain, no resid (edge case)
                label = f"({resname}:)"
                selection = f"resname {resname}"
            
            self._hetatm_labels.append(label)
            self._hetatm_selections.append(selection)

    def get_hetatm_labels(self) -> list[str]:
        """Get the list of available HETATM labels (e.g., '(LIG:A)')."""
        return self._hetatm_labels.copy()

    def open_input_pdb(self) -> None:
        """Open the input PDB file in PyMOL (protein only)."""
        self._open_pdb_only(self._input_pdb)

    def open_output_pdb_only(self) -> None:
        """Open the output PDB file in PyMOL (protein only)."""
        self._open_pdb_only(self._output_pdb)

    def open_output_pdb_with_box(
        self,
        line_width: float = 3.0,
    ) -> None:
        """Open the output PDB file in PyMOL, with docking box if available."""
        self._open_pdb_with_box(self._output_pdb, line_width=line_width)

    def detect_flexible_residues(self, margin: float | collections.abc.Sequence[float] = 0.0) -> list[str]:
        """
        Detect residues within the docking region + margin.
        
        Args:
            margin: Additional margin to expand the docking region. Can be a scalar
                   or a length-3 sequence for different margins in x, y, z.
        
        Returns:
            list[str]: List of residue identifiers in "chain:resid" format,
                      e.g., ["A:45", "A:82", "B:120"]
        
        Raises:
            ValueError: If docking region has not been set yet.
            FileNotFoundError: If output PDB file does not exist.
        """
        if self._docking_region is None:
            raise ValueError("Docking region has not been set. Call set_blind_docking_range, "
                           "set_known_ligand_range, or set_custom_docking_range first.")
        
        if not self._output_pdb.exists():
            raise FileNotFoundError(f"Output PDB not found: {self._output_pdb}")
        
        # Get docking region center and side length
        center, side_length = self._docking_region
        
        # Convert margin to numpy array
        margin_arr = np.array(margin, dtype=float)
        if margin_arr.ndim == 0:
            margin_arr = np.full(3, margin_arr)
        if margin_arr.shape != (3,):
            raise ValueError("margin must be a scalar or a length-3 sequence.")
        
        # Calculate expanded box boundaries
        expanded_side_length = side_length + margin_arr
        half_side = expanded_side_length / 2.0
        
        min_corner = center - half_side
        max_corner = center + half_side
        
        # Load the universe from output PDB
        u = MDAnalysis.Universe(str(self._output_pdb))
        
        # Select all protein atoms
        all_protein_atoms = u.select_atoms("protein")
        
        # Helper function to check if an atom is within the box
        def is_atom_in_box(atom):
            pos = atom.position
            return (min_corner[0] <= pos[0] <= max_corner[0] and
                    min_corner[1] <= pos[1] <= max_corner[1] and
                    min_corner[2] <= pos[2] <= max_corner[2])
        
        # Use MDAnalysis residues property for efficient grouping
        fully_enclosed_residues = []
        for residue in all_protein_atoms.residues:
            # Get chain ID from first atom of the residue
            first_atom = residue.atoms[0]
            chain = getattr(first_atom, "chainID", "").strip()
            if not chain:
                raise ValueError(f"Residue {residue.resid} does not have a chainID. "
                               "The PDB file must contain chain identifiers for flexible residue detection.")
            
            # Check if ALL atoms in this residue are within the box
            if all(is_atom_in_box(atom) for atom in residue.atoms):
                fully_enclosed_residues.append((chain, residue.resid))
        
        # Format as "chain:resid" and sort
        result = [f"{chain}:{resid}" for chain, resid in sorted(fully_enclosed_residues)]
        
        return result

    def detect_pockets(
        self,
        *,
        r_min: float = 3.0,
        r_max: float = 6.0,
        polar_probe_radius: float = 1.4,
        sasa_threshold: float = 20.0,
        merge_distance: float = 1.75,
        min_spheres: int = 35,
        ignore_hydrogens: bool = True,
        ignore_water: bool = True,
        ignore_hetero: bool = True,
    ) -> list:
        """
        Detect binding pockets in the output PDB structure using pocketeer.
        
        This method uses the alpha-sphere algorithm to detect potential binding
        pockets in the protein structure. Detected pockets are stored internally
        and can be retrieved using get_detected_pockets().
        
        Args:
            r_min: Minimum alpha-sphere radius (Å). Default: 3.0
            r_max: Maximum alpha-sphere radius (Å). Default: 6.0
            polar_probe_radius: Probe radius for SASA calculation (Å). Default: 1.4
            sasa_threshold: Threshold for mean SASA value to determine if a sphere
                is buried (Å²). Default: 20.0
            merge_distance: Distance threshold for merging nearby sphere clusters (Å).
                Default: 1.75
            min_spheres: Minimum number of spheres per pocket cluster. Default: 35
            ignore_hydrogens: Ignore hydrogen atoms. Default: True
            ignore_water: Ignore water molecules. Default: True
            ignore_hetero: Ignore hetero atoms (ligands). Default: True
        
        Returns:
            list: List of pocket regions, each as [center, side_length] where:
                  - center: numpy array of shape (3,) with x, y, z coordinates
                  - side_length: numpy array of shape (3,) with bounding box dimensions
                  Pockets are sorted by score (highest first).
        
        Raises:
            FileNotFoundError: If the output PDB file does not exist.
        """
        if not self._output_pdb.exists():
            raise FileNotFoundError(f"Output PDB not found: {self._output_pdb}")
        
        # Load structure using pocketeer's loader
        atomarray = pocketeer.load_structure(str(self._output_pdb))
        
        # Normalize non-standard histidine residue names to standard HIS
        # Biotite's SASA calculation doesn't recognize CHARMM/AMBER-specific
        # histidine names (HSD, HSE, HSP from CHARMM; HID, HIE, HIP from AMBER)
        # which causes "Residue 'HSD' does not contain an atom named N" errors
        his_variants = {"HSD", "HSE", "HSP", "HID", "HIE", "HIP"}
        for i, res_name in enumerate(atomarray.res_name):
            if res_name in his_variants:
                atomarray.res_name[i] = "HIS"
        
        # Detect pockets
        pockets = pocketeer.find_pockets(
            atomarray,
            r_min=r_min,
            r_max=r_max,
            polar_probe_radius=polar_probe_radius,
            sasa_threshold=sasa_threshold,
            merge_distance=merge_distance,
            min_spheres=min_spheres,
            ignore_hydrogens=ignore_hydrogens,
            ignore_water=ignore_water,
            ignore_hetero=ignore_hetero,
        )
        
        # Convert Pocket objects to [center, side_length] format
        pocket_regions = []
        for pocket in pockets:
            # Collect all sphere centers and radii
            sphere_centers = np.array([sphere.center for sphere in pocket.spheres])
            sphere_radii = np.array([sphere.radius for sphere in pocket.spheres])
            
            # Calculate bounding box considering sphere radii
            # Each sphere extends from (center - radius) to (center + radius)
            min_corner = (sphere_centers - sphere_radii[:, np.newaxis]).min(axis=0)
            max_corner = (sphere_centers + sphere_radii[:, np.newaxis]).max(axis=0)
            
            center = (min_corner + max_corner) / 2.0
            side_length = max_corner - min_corner
            
            pocket_regions.append([center, side_length])
        
        # Print diagnostic information
        if pocket_regions:
            print(f"Pocketeer detected {len(pocket_regions)} pocket(s)")
        else:
            print("Pocketeer detected no pockets in the protein structure")
        
        # Store detected pockets
        self._detected_pockets = pocket_regions
        
        return pocket_regions

    def get_detected_pockets(self) -> list[list[np.ndarray]] | None:
        """
        Get the previously detected binding pockets.
        
        Returns:
            list | None: List of pocket regions if detect_pockets() was called,
                        None otherwise. Each pocket region is [center, side_length]:
                        - center: numpy array (3,) with x, y, z coordinates
                        - side_length: numpy array (3,) with bounding box dimensions
        """
        return self._detected_pockets

    def readd_ions(self, distance_cutoff: float = 2.5) -> int:
        """
        Re-add metal ions from the input PDB that are within a distance cutoff of the output PDB.
        
        This method finds metal ions in the input PDB that are within the specified
        distance of any atom in the output PDB, and appends them to the output PDB file.
        
        Common metal ions detected: Li, Na, K, Rb, Cs, Mg, Ca, Sr, Ba, Mn, Fe, Co, Ni, Cu, Zn,
        Cd, Hg, Al, Ga, In, Sn, Pb, Bi, Cr, Mo, W, V.
        
        Args:
            distance_cutoff: Maximum distance (Å) for an ion to be considered close to
                           the output structure. Default: 2.5
        
        Returns:
            int: Number of ions added to the output PDB.
        
        Raises:
            FileNotFoundError: If input or output PDB files do not exist.
        """
        if not self._input_pdb.exists():
            raise FileNotFoundError(f"Input PDB not found: {self._input_pdb}")
        if not self._output_pdb.exists():
            raise FileNotFoundError(f"Output PDB not found: {self._output_pdb}")
        
        # Load both universes
        u_input = MDAnalysis.Universe(str(self._input_pdb))
        u_output = MDAnalysis.Universe(str(self._output_pdb))

        # Build selection string for metal ions from class constant
        resname_list = " ".join(self._METAL_ION_RESNAMES)
        metal_selection = f"resname {resname_list}"
        
        try:
            input_metals = u_input.select_atoms(metal_selection)
        except Exception:
            # If selection fails (e.g., no matching atoms), return 0
            print("No metal ions found in input PDB.")
            return 0
        
        if input_metals.n_atoms == 0:
            print("No metal ions found in input PDB.")
            return 0
        
        output_atoms = u_output.atoms
        if output_atoms.n_atoms == 0:
            print("Output PDB has no atoms.")
            return 0
        
        # Find metal ions within distance cutoff of any output atom
        ions_to_add = []
        for metal_atom in input_metals:
            metal_pos = metal_atom.position
            # Calculate distances to all output atoms
            distances = np.linalg.norm(output_atoms.positions - metal_pos, axis=1)
            min_distance = distances.min()
            
            if min_distance <= distance_cutoff:
                ions_to_add.append(metal_atom)
        
        if not ions_to_add:
            print(f"No metal ions found within {distance_cutoff} Å of the output structure.")
            return 0
        
        # Read existing output PDB content
        with open(self._output_pdb, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Remove trailing END/ENDMDL lines for appending
        content_lines = []
        end_lines = []
        for line in lines:
            if line.strip().startswith(('END', 'ENDMDL')):
                end_lines.append(line)
            else:
                content_lines.append(line)
        
        # Generate HETATM lines for ions
        # Get the highest atom serial number from existing content
        max_serial = 0
        for line in content_lines:
            if line.startswith(('ATOM', 'HETATM')) and len(line) > 11:
                try:
                    serial = int(line[6:11].strip())
                    max_serial = max(max_serial, serial)
                except ValueError:
                    pass
        
        ion_lines = []
        for i, atom in enumerate(ions_to_add):
            serial = max_serial + i + 1
            atom_name = getattr(atom, 'name', atom.resname)[:4].ljust(4)
            resname = atom.resname[:3].ljust(3)
            chain = getattr(atom, 'chainID', '') or getattr(atom, 'segid', '') or 'X'
            chain = chain[:1]
            resid = getattr(atom, 'resid', 1)
            x, y, z = atom.position
            occupancy = getattr(atom, 'occupancy', 1.0)
            tempfactor = getattr(atom, 'tempfactor', 0.0)
            element = atom.element if hasattr(atom, 'element') and atom.element else atom_name.strip()[:2]
            
            # Format HETATM line according to PDB specification
            hetatm_line = (
                f"HETATM{serial:5d} {atom_name} {resname} {chain}{resid:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{tempfactor:6.2f}          {element:>2s}\n"
            )
            ion_lines.append(hetatm_line)
        
        # Write updated content
        with open(self._output_pdb, 'w', encoding='utf-8') as f:
            f.writelines(content_lines)
            f.writelines(ion_lines)
            if end_lines:
                f.writelines(end_lines)
            else:
                f.write("END\n")
        
        print(f"Added {len(ions_to_add)} metal ion(s) to the output PDB.")
        return len(ions_to_add)

    def readd_waters(self, hydrogens: int = 2, distance_cutoff: float = 3.5) -> int:
        """
        Re-add water molecules from the input PDB that form potential hydrogen bonds 
        with the protein in the output PDB.
        
        This method uses a distance-based approach to identify water molecules that
        could form hydrogen bonds with the processed protein structure. It checks
        distances between water oxygen atoms (from input PDB) and polar atoms 
        (N, O from output PDB protein).
        
        A water molecule is considered to form a potential hydrogen bond if its
        oxygen atom is within the distance_cutoff of a polar protein atom.
        
        Args:
            hydrogens: Minimum number of potential hydrogen bonds (polar contacts)
                      required for a water molecule to be added. Must be 2 or 3. 
                      Default: 2
            distance_cutoff: Maximum distance (Å) between water oxygen and protein
                           polar atoms to be considered a potential hydrogen bond.
                           Default: 3.5
        
        Returns:
            int: Number of water molecules added to the output PDB.
        
        Raises:
            ValueError: If hydrogens is not 2 or 3.
            FileNotFoundError: If input or output PDB files do not exist.
        """
        if hydrogens not in (2, 3):
            raise ValueError("hydrogens must be 2 or 3")
        
        if not self._input_pdb.exists():
            raise FileNotFoundError(f"Input PDB not found: {self._input_pdb}")
        if not self._output_pdb.exists():
            raise FileNotFoundError(f"Output PDB not found: {self._output_pdb}")
        
        # Load both universes
        u_input = MDAnalysis.Universe(str(self._input_pdb))
        u_output = MDAnalysis.Universe(str(self._output_pdb))
        
        # Select water molecules from input PDB
        water_resnames = "HOH H2O WAT WT TIP3 TIP4 SPC"
        input_waters = u_input.select_atoms(f"resname {water_resnames}")
        
        if input_waters.n_atoms == 0:
            print("No water molecules found in input PDB.")
            return 0
        
        # Select polar atoms (N, O) from protein in output PDB
        # These are the potential hydrogen bond donors/acceptors
        try:
            output_protein_polar = u_output.select_atoms("protein and (name N* or name O*)")
        except Exception:
            print("No protein atoms found in output PDB.")
            return 0
        
        if output_protein_polar.n_atoms == 0:
            print("No polar atoms (N, O) found in output PDB protein.")
            return 0
        
        # Get positions of protein polar atoms
        protein_polar_positions = output_protein_polar.positions
        
        # Select only oxygen atoms from water molecules for distance calculation
        input_water_oxygens = u_input.select_atoms(f"resname {water_resnames} and name O*")
        
        if input_water_oxygens.n_atoms == 0:
            print("No water oxygen atoms found in input PDB.")
            return 0
        
        # Count potential hydrogen bonds (polar contacts) per water residue
        water_contact_counts = collections.Counter()
        
        for water_o in input_water_oxygens:
            water_pos = water_o.position
            
            # Calculate distances to all protein polar atoms
            distances = np.linalg.norm(protein_polar_positions - water_pos, axis=1)
            
            # Count contacts within cutoff distance
            n_contacts = np.sum(distances <= distance_cutoff)
            
            if n_contacts > 0:
                # Store unique water residue identifier
                chain = getattr(water_o, 'chainID', '') or getattr(water_o, 'segid', '')
                water_resid = (water_o.resname, water_o.resid, chain)
                water_contact_counts[water_resid] = n_contacts
        
        # Filter waters with enough potential hydrogen bonds
        qualified_waters = {resid for resid, count in water_contact_counts.items() 
                          if count >= hydrogens}
        
        if not qualified_waters:
            print(f"No water molecules found with >= {hydrogens} polar contacts to protein.")
            return 0
        
        # Collect atoms from qualified water residues
        waters_to_add = []
        for resname, resid, chain in qualified_waters:
            if chain:
                selection = f"resname {resname} and resid {resid} and (chainid {chain} or segid {chain})"
            else:
                selection = f"resname {resname} and resid {resid}"
            
            try:
                water_atoms = u_input.select_atoms(selection)
                if water_atoms.n_atoms > 0:
                    waters_to_add.extend(water_atoms)
            except Exception:
                pass
        
        if not waters_to_add:
            print("No water molecules to add.")
            return 0
        
        # Read existing output PDB content
        with open(self._output_pdb, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Remove trailing END/ENDMDL lines for appending
        content_lines = []
        end_lines = []
        for line in lines:
            if line.strip().startswith(('END', 'ENDMDL')):
                end_lines.append(line)
            else:
                content_lines.append(line)
        
        # Get the highest atom serial number from existing content
        max_serial = 0
        for line in content_lines:
            if line.startswith(('ATOM', 'HETATM')) and len(line) > 11:
                try:
                    serial = int(line[6:11].strip())
                    max_serial = max(max_serial, serial)
                except ValueError:
                    pass
        
        # Generate HETATM lines for water atoms
        water_lines = []
        for i, atom in enumerate(waters_to_add):
            serial = max_serial + i + 1
            atom_name = atom.name[:4].ljust(4)
            resname = atom.resname[:3].ljust(3)
            chain = getattr(atom, 'chainID', '') or getattr(atom, 'segid', '') or 'W'
            chain = chain[:1]
            resid = getattr(atom, 'resid', 1)
            x, y, z = atom.position
            occupancy = getattr(atom, 'occupancy', 1.0)
            tempfactor = getattr(atom, 'tempfactor', 0.0)
            element = atom.element if hasattr(atom, 'element') and atom.element else atom_name.strip()[:1]
            
            # Format HETATM line according to PDB specification
            hetatm_line = (
                f"HETATM{serial:5d} {atom_name} {resname} {chain}{resid:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{tempfactor:6.2f}          {element:>2s}\n"
            )
            water_lines.append(hetatm_line)
        
        # Write updated content
        with open(self._output_pdb, 'w', encoding='utf-8') as f:
            f.writelines(content_lines)
            f.writelines(water_lines)
            if end_lines:
                f.writelines(end_lines)
            else:
                f.write("END\n")
        
        n_waters = len(qualified_waters)
        print(f"Added {n_waters} water molecule(s) ({len(waters_to_add)} atoms) to the output PDB.")
        return n_waters

    def readd_metal_coor_waters(self, distance_cutoff: float = 2.5) -> int:
        """
        Re-add water molecules from the input PDB that coordinate with metal ions in the output PDB.
        
        This method finds water molecules in the input PDB that are within the specified
        distance of any metal ion in the output PDB, and appends them to the output PDB file.
        These waters typically serve as coordination ligands for metal ions.
        
        Args:
            distance_cutoff: Maximum distance (Å) for a water molecule to be considered
                           coordinating with a metal ion. Default: 2.5
        
        Returns:
            int: Number of water molecules added to the output PDB.
        
        Raises:
            FileNotFoundError: If input or output PDB files do not exist.
        
        Note:
            This function should typically be called after readd_ions() to ensure
            metal ions are present in the output PDB.
        """
        if not self._input_pdb.exists():
            raise FileNotFoundError(f"Input PDB not found: {self._input_pdb}")
        if not self._output_pdb.exists():
            raise FileNotFoundError(f"Output PDB not found: {self._output_pdb}")
        
        # Common metal ion residue names in PDB files (same as readd_ions)
        metal_ion_resnames = {
            "LI", "NA", "K", "RB", "CS",           # Alkali metals
            "MG", "CA", "CAL", "SR", "BA",         # Alkaline earth metals
            "MN", "FE", "CO", "NI", "CU", "ZN",    # Transition metals (row 1)
            "CD", "HG",                            # Transition metals (row 2-3)
            "AL", "GA", "IN", "SN", "PB", "BI",    # Post-transition metals
            "CR", "MO", "W", "V",                  # Other transition metals
            "ZN2", "FE2", "FE3", "MN2", "CU1", "CU2",  # Common ion states
            "CA2", "MG2", "NA1", "K1", "POT", "SOD",   # Force field variants
        }
        
        # Water residue names
        water_resnames = {"HOH", "H2O", "WAT", "WT", "TIP3", "TIP4", "SPC"}
        
        # Load both universes
        u_input = MDAnalysis.Universe(str(self._input_pdb))
        u_output = MDAnalysis.Universe(str(self._output_pdb))
        
        # Find metal ions in the output PDB
        resname_list = " ".join(metal_ion_resnames)
        metal_selection = f"resname {resname_list}"
        
        try:
            output_metals = u_output.select_atoms(metal_selection)
        except Exception:
            print("No metal ions found in output PDB.")
            return 0
        
        if output_metals.n_atoms == 0:
            print("No metal ions found in output PDB. Call readd_ions() first.")
            return 0
        
        # Find water molecules in the input PDB
        water_resname_list = " ".join(water_resnames)
        water_selection = f"resname {water_resname_list}"
        
        try:
            input_waters = u_input.select_atoms(water_selection)
        except Exception:
            print("No water molecules found in input PDB.")
            return 0
        
        if input_waters.n_atoms == 0:
            print("No water molecules found in input PDB.")
            return 0
        
        # Get metal ion positions from output PDB
        metal_positions = output_metals.positions
        
        # Find water residues within distance cutoff of any metal ion
        coordinating_water_resids = set()
        
        for water_atom in input_waters:
            water_pos = water_atom.position
            # Calculate distances to all metal ions
            distances = np.linalg.norm(metal_positions - water_pos, axis=1)
            min_distance = distances.min()
            
            if min_distance <= distance_cutoff:
                # Store unique water residue identifier
                chain = getattr(water_atom, 'chainID', '') or getattr(water_atom, 'segid', '')
                water_resid = (water_atom.resname, water_atom.resid, chain)
                coordinating_water_resids.add(water_resid)
        
        if not coordinating_water_resids:
            print(f"No water molecules found within {distance_cutoff} Å of metal ions.")
            return 0
        
        # Collect all atoms from coordinating water residues
        waters_to_add = []
        for resname, resid, chain in coordinating_water_resids:
            if chain:
                selection = f"resname {resname} and resid {resid} and (chainid {chain} or segid {chain})"
            else:
                selection = f"resname {resname} and resid {resid}"
            
            try:
                water_atoms = u_input.select_atoms(selection)
                if water_atoms.n_atoms > 0:
                    waters_to_add.extend(water_atoms)
            except Exception:
                pass
        
        if not waters_to_add:
            print("No water molecules to add.")
            return 0
        
        # Read existing output PDB content
        with open(self._output_pdb, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Remove trailing END/ENDMDL lines for appending
        content_lines = []
        end_lines = []
        for line in lines:
            if line.strip().startswith(('END', 'ENDMDL')):
                end_lines.append(line)
            else:
                content_lines.append(line)
        
        # Get the highest atom serial number from existing content
        max_serial = 0
        for line in content_lines:
            if line.startswith(('ATOM', 'HETATM')) and len(line) > 11:
                try:
                    serial = int(line[6:11].strip())
                    max_serial = max(max_serial, serial)
                except ValueError:
                    pass
        
        # Generate HETATM lines for water atoms
        water_lines = []
        for i, atom in enumerate(waters_to_add):
            serial = max_serial + i + 1
            atom_name = atom.name[:4].ljust(4)
            resname = atom.resname[:3].ljust(3)
            chain = getattr(atom, 'chainID', '') or getattr(atom, 'segid', '') or 'W'
            chain = chain[:1]
            resid = getattr(atom, 'resid', 1)
            x, y, z = atom.position
            occupancy = getattr(atom, 'occupancy', 1.0)
            tempfactor = getattr(atom, 'tempfactor', 0.0)
            element = atom.element if hasattr(atom, 'element') and atom.element else atom_name.strip()[:1]
            
            # Format HETATM line according to PDB specification
            hetatm_line = (
                f"HETATM{serial:5d} {atom_name} {resname} {chain}{resid:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{tempfactor:6.2f}          {element:>2s}\n"
            )
            water_lines.append(hetatm_line)
        
        # Write updated content
        with open(self._output_pdb, 'w', encoding='utf-8') as f:
            f.writelines(content_lines)
            f.writelines(water_lines)
            if end_lines:
                f.writelines(end_lines)
            else:
                f.write("END\n")
        
        n_waters = len(coordinating_water_resids)
        print(f"Added {n_waters} metal-coordinating water molecule(s) ({len(waters_to_add)} atoms) to the output PDB.")
        return n_waters

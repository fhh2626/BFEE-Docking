# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Haohao Fu <fhh2626@nankai.edu.cn, fhh2626@gmail.com>

import collections
import os
import pathlib
import shutil
import subprocess
import sys
import traceback
import typing

import Bio.PDB
import MDAnalysis
import numpy
import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.rdFMCS
import scipy.optimize

import third_party_tools

# Add third_party directory to Python path for bundled packages (e.g., pocketeer)
third_party_tools.add_third_party_to_path()

# Module-level reference for backward compatibility
_SCRIPT_DIR = third_party_tools.get_third_party_dir().parent

class MDPrepper:
    """
    MD Prepper class for preparing molecular dynamics simulation input files.
    
    This class handles:
    1. Fixing flexible residue atom names by matching to reference receptor
    2. Merging rigid protein with corrected flexible residues
    3. Outputting complete protein PDB (without ligand)
    """
    
    def __init__(
        self,
        receptor_pdb_path: str,
        ligand_pdb_path: str,
        flex_res_pdb_path: str,
        ligand_smiles: str,
        ph: float,
        output_dir: str,
        vmd_path: str
    ):
        """
        Initialize MDPrepper.
        
        Args:
            receptor_pdb_path: Path to the reference receptor PDB file (complete protein)
            ligand_pdb_path: Path to the docked ligand PDB file
            flex_res_pdb_path: Path to the flexible residues PDB file from docking
            ligand_smiles: SMILES string of the ligand
            ph: pH value for protonation
            output_dir: Directory for output files
            vmd_path: Path to the VMD executable
        """
        self._receptor_pdb_path = str(pathlib.Path(receptor_pdb_path).resolve())
        self._ligand_pdb_path = str(pathlib.Path(ligand_pdb_path).resolve())
        self._flex_res_pdb_path = str(pathlib.Path(flex_res_pdb_path).resolve())
        self._ligand_smiles = ligand_smiles
        self._ph = ph
        self._output_dir = str(pathlib.Path(output_dir).resolve())
        
        # Cache for parsed PDB data
        self._receptor_lines: typing.Optional[list] = None
        self._flex_res_lines: typing.Optional[list] = None
        
        # Cache for disulfide bonds detection
        # Each bond is stored as ((chain1, resid1, resname1), (chain2, resid2, resname2), distance)
        self._disulfide_bonds: typing.Optional[typing.List[typing.Tuple[tuple, tuple, float]]] = None
        
        # Cache for protonation state patches
        # Dictionary mapping patch names to list of residues: {patch_name: [(chain, resid, resname, charge), ...]}
        self._protonation_patches: typing.Optional[typing.Dict[str, typing.List[typing.Tuple[str, int, str, float]]]] = None
        
        # Get absolute path to obabel executable (for mol2 generation)
        self._obabel_executable = third_party_tools.get_obabel_executable()
        
        # Get VMD executable path and verify it exists
        self._vmd_executable = third_party_tools.get_vmd_executable(vmd_path)
        if not self._vmd_executable.exists():
            raise FileNotFoundError(f"VMD executable not found at: {self._vmd_executable}")
    
    # ========== Private Helper Functions for PDB Parsing ==========
    
    @staticmethod
    def _parse_pdb_line(line: str) -> typing.Optional[dict]:
        """
        Parse a single PDB ATOM/HETATM line into a dictionary.
        
        Args:
            line: A single line from a PDB file
            
        Returns:
            Dictionary with atom properties, or None if not an ATOM/HETATM line
        """
        if not line.startswith(("ATOM", "HETATM")):
            return None
        
        try:
            record = {
                "record_type": line[0:6].strip(),
                "atom_serial": int(line[6:11].strip()) if line[6:11].strip() else 0,
                "atom_name": line[12:16].strip(),
                "alt_loc": line[16:17].strip(),
                "res_name": line[17:20].strip(),
                "chain_id": line[21:22].strip(),
                "res_seq": int(line[22:26].strip()) if line[22:26].strip() else 0,
                "icode": line[26:27].strip(),
                "x": float(line[30:38].strip()) if line[30:38].strip() else 0.0,
                "y": float(line[38:46].strip()) if line[38:46].strip() else 0.0,
                "z": float(line[46:54].strip()) if line[46:54].strip() else 0.0,
                "occupancy": float(line[54:60].strip()) if line[54:60].strip() else 1.0,
                "temp_factor": float(line[60:66].strip()) if line[60:66].strip() else 0.0,
                "segment_id": line[72:76].strip() if len(line) > 72 else "",
                "element": line[76:78].strip() if len(line) > 76 else "",
                "charge": line[78:80].strip() if len(line) > 78 else "",
                "original_line": line,
            }
            return record
        except (ValueError, IndexError):
            return None
    
    @staticmethod
    def _format_pdb_line(atom: dict) -> str:
        """
        Format a PDB atom record back to a PDB line.
        
        PDB atom name field is columns 13-16 (4 characters).
        Standard alignment rules:
        - 1-char names (C, N, O): position 2, e.g., " C  "
        - 2-char names (CA, CB): position 2, e.g., " CA "
        - 3-char names (CG1): position 2, e.g., " CG1"
        - 4-char names (1HG1, HD21): position 1, e.g., "1HG1"
        - Names starting with digit: position 1, e.g., "1HG1"
        """
        atom_name = atom["atom_name"]
        name_len = len(atom_name)
        
        if name_len == 4:
            pass  # 4-char names stay as-is
        elif atom_name[0].isdigit():
            atom_name = f"{atom_name:<4}"
        else:
            atom_name = f" {atom_name:<3}"
        
        # Ensure single-char fields output correct width even when empty
        alt_loc = atom['alt_loc'] if atom['alt_loc'] else ' '
        chain_id = atom['chain_id'] if atom['chain_id'] else ' '
        icode = atom['icode'] if atom['icode'] else ' '
        
        # Handle residue name: standard PDB uses 3 chars (cols 18-20), CHARMM allows 4 chars
        # 4-char residue names (e.g., TIP3) occupy cols 18-21 (including the space at col 21)
        res_name = atom['res_name']
        if len(res_name) >= 4:
            # 4-char residue name: alt_loc(col17) + 4-char resname(cols18-21) + chain_id(col22)
            # No space at col 21
            alt_loc_resname_space = f"{alt_loc}{res_name[:4]}"
        else:
            # Standard 3-char residue name: alt_loc(col17) + 3-char resname(cols18-20) + space(col21)
            alt_loc_resname_space = f"{alt_loc}{res_name:>3} "
        
        line = (
            f"{atom['record_type']:<6}"      # cols 1-6:   record type
            f"{atom['atom_serial']:>5} "     # cols 7-11:  atom serial + col 12 space
            f"{atom_name}"                   # cols 13-16: atom name (4 chars)
            f"{alt_loc_resname_space}"       # cols 17-21: alt_loc + residue name (+ space)
            f"{chain_id}"                    # col 22:     chain ID
            f"{atom['res_seq']:>4}"          # cols 23-26: residue seq
            f"{icode}   "                    # col 27:     icode + cols 28-30 spaces
            f"{atom['x']:>8.3f}"             # cols 31-38: X
            f"{atom['y']:>8.3f}"             # cols 39-46: Y
            f"{atom['z']:>8.3f}"             # cols 47-54: Z
            f"{atom['occupancy']:>6.2f}"     # cols 55-60: occupancy
            f"{atom['temp_factor']:>6.2f}"   # cols 61-66: temp factor
        )
        
        line = f"{line:<76}"
        
        if atom.get("segment_id"):
            line = line[:72] + f"{atom['segment_id']:<4}" + line[76:]
        if atom.get("element"):
            line = line[:76] + f"{atom['element']:>2}"
        else:
            line = line[:76] + "  "
        if atom.get("charge"):
            line = line[:78] + f"{atom['charge']:>2}"
        
        return line.rstrip()
    
    @staticmethod
    def _calculate_distance(atom1: dict, atom2: dict) -> float:
        """Calculate Euclidean distance between two atoms."""
        dx = atom1["x"] - atom2["x"]
        dy = atom1["y"] - atom2["y"]
        dz = atom1["z"] - atom2["z"]
        return (dx * dx + dy * dy + dz * dz) ** 0.5
    
    @staticmethod
    def _get_element_from_atom_name(atom_name: str) -> str:
        """
        Extract element from atom name.
        
        PDB atom names follow specific conventions:
        - Standard atoms: element is the first 1-2 characters
        - Hydrogen atoms often start with a digit (e.g., 1HG1 -> H)
        """
        name = atom_name.strip()
        if not name:
            return ""
        
        if len(name) == 1 and name.isalpha():
            return name.upper()
        
        if name[0].isdigit():
            for c in name:
                if c.isalpha():
                    return c.upper()
            return ""
        
        two_char_elements = ["CL", "BR", "FE", "ZN", "MG", "NA", "CU", "MN", "CO", "NI", "SE"]
        if len(name) >= 2:
            first_two = name[0:2].upper()
            if first_two in two_char_elements:
                return first_two
        
        return name[0].upper()
    
    @staticmethod
    def _greedy_matching(cost_matrix: list) -> list:
        """
        Greedy matching algorithm for bipartite assignment.
        
        Args:
            cost_matrix: 2D list of distances between atoms
            
        Returns:
            List of (row, col) pairs representing the matching
        """
        n = len(cost_matrix)
        m = len(cost_matrix[0]) if n > 0 else 0
        
        if n == 0 or m == 0:
            return []
        
        edges = []
        for i in range(n):
            for j in range(m):
                edges.append((cost_matrix[i][j], i, j))
        
        edges.sort(key=lambda x: x[0])
        
        used_rows = set()
        used_cols = set()
        assignment = []
        
        for dist, row, col in edges:
            if row not in used_rows and col not in used_cols:
                assignment.append((row, col))
                used_rows.add(row)
                used_cols.add(col)
                
                if len(used_rows) >= n:
                    break
        
        return assignment
    
    # ========== Private Functions for Flex Residue Processing ==========
    
    def _parse_flex_residue_atoms(self) -> dict:
        """
        Parse flexible residue atoms from the flex_res PDB file.
        
        Returns:
            Dictionary of {(chain_id, res_seq, res_name): [atoms]}
        """
        if self._flex_res_lines is None:
            path = pathlib.Path(self._flex_res_pdb_path)
            if not path.exists() or not path.is_file():
                return {}
                
            with open(self._flex_res_pdb_path, 'r', encoding='utf-8') as f:
                self._flex_res_lines = f.readlines()
        
        flex_atoms = collections.defaultdict(list)
        
        for line in self._flex_res_lines:
            atom = self._parse_pdb_line(line)
            if atom is None:
                continue
            
            key = (atom["chain_id"], atom["res_seq"], atom["res_name"])
            flex_atoms[key].append(atom)
        
        return dict(flex_atoms)
    
    def _get_receptor_residue_atoms(
        self, 
        chain: str, 
        res_seq: int, 
        res_name: str,
        sidechain_only: bool = True,
        backbone_only: bool = False
    ) -> list:
        """
        Extract atoms for a specific residue from the receptor PDB.
        
        Args:
            chain: Chain ID
            res_seq: Residue sequence number
            res_name: Residue name
            sidechain_only: If True, only return side chain atoms (excludes backbone)
            backbone_only: If True, only return backbone atoms (N, CA, C, O)
            
        Returns:
            List of atom dictionaries
        """
        if self._receptor_lines is None:
            with open(self._receptor_pdb_path, 'r', encoding='utf-8') as f:
                self._receptor_lines = f.readlines()
        
        backbone_atoms = {"N", "CA", "C", "O"}
        atoms = []
        
        for line in self._receptor_lines:
            atom = self._parse_pdb_line(line)
            if atom is None:
                continue
            
            if (atom["chain_id"] == chain and 
                atom["res_seq"] == res_seq and 
                atom["res_name"] == res_name):
                
                is_backbone = atom["atom_name"] in backbone_atoms
                
                if backbone_only:
                    if is_backbone:
                        atoms.append(atom)
                elif sidechain_only:
                    if not is_backbone:
                        atoms.append(atom)
                else:
                    atoms.append(atom)
        
        return atoms
    
    def _filter_backbone_duplicates(
        self, 
        flex_atoms: list,
        receptor_backbone_atoms: list,
        overlap_threshold: float = 0.1
    ) -> list:
        """
        Filter out backbone atoms from flexible residue atoms by coordinate matching.
        
        Since flexible residue atom names may be incorrect (e.g., all named 'C'),
        we use coordinate comparison with receptor backbone atoms to identify
        which flex atoms are actually backbone atoms.
        
        Args:
            flex_atoms: List of flexible residue atoms
            receptor_backbone_atoms: List of backbone atoms (N, CA, C, O) from receptor
            overlap_threshold: Distance threshold for backbone detection
            
        Returns:
            List of atoms without backbone duplicates
        """
        filtered = []
        
        for flex_atom in flex_atoms:
            is_backbone = False
            
            # Check if this flex atom overlaps with any receptor backbone atom
            for bb_atom in receptor_backbone_atoms:
                dist = self._calculate_distance(flex_atom, bb_atom)
                if dist < overlap_threshold:
                    # This flex atom is a backbone atom, skip it
                    is_backbone = True
                    break
            
            if not is_backbone:
                filtered.append(flex_atom)
        
        return filtered
    
    def _match_atoms_and_assign_names(
        self, 
        flex_atoms: list, 
        receptor_atoms: list
    ) -> list:
        """
        Match flexible residue atoms to receptor atoms and assign correct names.
        
        This method uses topology-based matching that does NOT rely on coordinates,
        since flexible docking changes atom positions significantly.
        
        Strategy:
        1. MCS-based topological matching using RDKit
        2. Element-order matching as fallback (atoms of same element matched in order)
        
        Note: Coordinate-based matching (distance comparison, Hungarian algorithm) 
        is NOT used because flexible docking changes sidechain coordinates.
        
        Args:
            flex_atoms: List of flexible residue atoms (with incorrect names)
            receptor_atoms: List of reference receptor atoms (with correct names)
            
        Returns:
            List of atoms with corrected atom names
        """
        if not flex_atoms or not receptor_atoms:
            return [atom.copy() for atom in flex_atoms]
        
        # Initialize result with None
        result = [None] * len(flex_atoms)
        used_receptor_indices = set()
        mcs_mapped_flex_indices = set()
        
        # Helper function to get element from atom (robust extraction)
        def get_element_robust(atom: dict) -> str:
            """Extract element from atom using multiple fallbacks."""
            # Try element field first
            if atom.get("element"):
                return atom["element"].strip().upper()
            # Fall back to atom name parsing
            return self._get_element_from_atom_name(atom["atom_name"]).upper()
        
        # ========== Strategy 1: MCS-based Topological Matching ==========
        try:
            
            def atoms_to_mol(atom_list: list):
                """Convert atom list to RDKit Mol object."""
                lines = []
                for i, at in enumerate(atom_list):
                    at_copy = at.copy()
                    at_copy['atom_serial'] = i + 1
                    lines.append(self._format_pdb_line(at_copy))
                block = "\n".join(lines)
                return rdkit.Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False)
            
            def get_mol_element(mol, idx: int) -> str:
                """Get element symbol from RDKit Mol atom."""
                if mol is None:
                    return ""
                try:
                    atom = mol.GetAtomWithIdx(idx)
                    return atom.GetSymbol().upper()
                except Exception:
                    return ""

            mol_flex = atoms_to_mol(flex_atoms)
            mol_rec = atoms_to_mol(receptor_atoms)
            
            if mol_flex and mol_rec:
                # Try MCS with strict element matching
                mcs = rdkit.Chem.rdFMCS.FindMCS(
                    [mol_flex, mol_rec], 
                    atomCompare=rdkit.Chem.rdFMCS.AtomCompare.CompareElements,
                    bondCompare=rdkit.Chem.rdFMCS.BondCompare.CompareAny,
                    ringMatchesRingOnly=False,
                    completeRingsOnly=False,
                    timeout=2.0
                )
                
                if mcs.smartsString and mcs.numAtoms > 0:
                    patt = rdkit.Chem.MolFromSmarts(mcs.smartsString)
                    if patt:
                        match_flex = mol_flex.GetSubstructMatch(patt)
                        match_rec = mol_rec.GetSubstructMatch(patt)
                        
                        # Apply MCS matches with element validation
                        valid_matches = 0
                        for k in range(len(match_flex)):
                            flex_idx = match_flex[k]
                            rec_idx = match_rec[k]
                            
                            # Validate element consistency (using both mol and dict)
                            mol_flex_elem = get_mol_element(mol_flex, flex_idx)
                            mol_rec_elem = get_mol_element(mol_rec, rec_idx)
                            dict_flex_elem = get_element_robust(flex_atoms[flex_idx])
                            dict_rec_elem = get_element_robust(receptor_atoms[rec_idx])
                            
                            # Element must match from at least one source
                            if (mol_flex_elem and mol_rec_elem and mol_flex_elem != mol_rec_elem):
                                continue  # Skip mismatched elements
                            if (not mol_flex_elem or not mol_rec_elem):
                                # Use dict-based elements as fallback
                                if dict_flex_elem and dict_rec_elem and dict_flex_elem != dict_rec_elem:
                                    continue
                            
                            # Valid match, assign name
                            rec_atom = receptor_atoms[rec_idx]
                            new_atom = flex_atoms[flex_idx].copy()
                            new_atom["atom_name"] = rec_atom["atom_name"]
                            new_atom["element"] = rec_atom.get("element", dict_rec_elem)
                            
                            # Debug output: show the mapping
                            flex_coord = f"({flex_atoms[flex_idx]['x']:.3f}, {flex_atoms[flex_idx]['y']:.3f}, {flex_atoms[flex_idx]['z']:.3f})"
                            rec_coord = f"({rec_atom['x']:.3f}, {rec_atom['y']:.3f}, {rec_atom['z']:.3f})"
                            print(f"      MCS Map: flex[{flex_idx}] '{flex_atoms[flex_idx].get('atom_name', '?')}' {flex_coord} -> rec[{rec_idx}] '{rec_atom['atom_name']}' {rec_coord}")
                            
                            result[flex_idx] = new_atom
                            mcs_mapped_flex_indices.add(flex_idx)
                            used_receptor_indices.add(rec_idx)
                            valid_matches += 1
                        
                        if valid_matches > 0:
                            print(f"    MCS matched {valid_matches}/{len(flex_atoms)} atoms (validated).")
                            
        except ImportError:
            print("    MCS matching skipped: RDKit not available.")
        except Exception as e:
            print(f"    MCS matching skipped: {e}")
        
        # ========== Strategy 2: Element-Order Matching for Remaining Atoms ==========
        # This strategy matches atoms of the same element type in order of appearance.
        # It does NOT use coordinates, which is important for flexible docking results.
        unmatched_flex_indices = [i for i in range(len(flex_atoms)) if i not in mcs_mapped_flex_indices]
        remaining_rec_indices = [j for j in range(len(receptor_atoms)) if j not in used_receptor_indices]
        
        if unmatched_flex_indices and remaining_rec_indices:
            # Group by element
            flex_by_element = collections.defaultdict(list)
            for i in unmatched_flex_indices:
                elem = get_element_robust(flex_atoms[i])
                flex_by_element[elem].append(i)
            
            rec_by_element = collections.defaultdict(list)
            for j in remaining_rec_indices:
                elem = get_element_robust(receptor_atoms[j])
                rec_by_element[elem].append(j)
            
            order_matches = 0
            
            for elem, flex_indices_for_elem in flex_by_element.items():
                rec_indices_for_elem = rec_by_element.get(elem, [])
                
                if not rec_indices_for_elem:
                    # No receptor atoms of this element, keep original names
                    for idx in flex_indices_for_elem:
                        if result[idx] is None:
                            result[idx] = flex_atoms[idx].copy()
                    continue
                
                # Match in order: first flex atom of element X matches first rec atom of element X
                # This is a simple heuristic that works when atom order is consistent
                for i, flex_idx in enumerate(flex_indices_for_elem):
                    if i < len(rec_indices_for_elem):
                        rec_idx = rec_indices_for_elem[i]
                        rec_atom = receptor_atoms[rec_idx]
                        new_atom = flex_atoms[flex_idx].copy()
                        new_atom["atom_name"] = rec_atom["atom_name"]
                        new_atom["element"] = rec_atom.get("element", "")
                        result[flex_idx] = new_atom
                        order_matches += 1
                    else:
                        # More flex atoms than receptor atoms for this element
                        if result[flex_idx] is None:
                            result[flex_idx] = flex_atoms[flex_idx].copy()
            
            if order_matches > 0:
                print(f"    Element-order matched {order_matches} remaining atoms.")
        
        # ========== Final Pass: Ensure All Atoms Have Values ==========
        for i in range(len(flex_atoms)):
            if result[i] is None:
                result[i] = flex_atoms[i].copy()
                print(f"    Warning: Atom {i} ({flex_atoms[i].get('atom_name', '?')}) could not be matched.")
        
        # Calculate and report match quality
        total_matched = sum(1 for i, r in enumerate(result) 
                          if r and r.get("atom_name") != flex_atoms[i].get("atom_name"))
        print(f"    Total: {total_matched}/{len(flex_atoms)} atoms renamed.")
        
        return result
    
    def _sort_atoms_by_receptor_order(
        self, 
        corrected_atoms: list, 
        receptor_atoms: list
    ) -> list:
        """
        Sort corrected atoms according to the order they appear in receptor.
        
        Args:
            corrected_atoms: List of atoms with corrected names
            receptor_atoms: Reference atoms from receptor
            
        Returns:
            Sorted list of atoms
        """
        receptor_order = {}
        for i, atom in enumerate(receptor_atoms):
            receptor_order[atom["atom_name"]] = i
        
        def sort_key(atom):
            name = atom["atom_name"]
            if name in receptor_order:
                return receptor_order[name]
            return 9999
        
        return sorted(corrected_atoms, key=sort_key)
    
    # ========== Public Interface ==========
    
    def fix_pdb(self, output_path: typing.Optional[str] = None) -> pathlib.Path:
        """
        Fix the PDB by merging rigid protein with corrected flexible residues.
        
        This method:
        1. Parses flexible residue atoms from flex_res_pdb_path
        2. Matches them to receptor to get correct atom names
        3. Replaces corresponding residues in rigid protein with corrected atoms
        4. Outputs a complete protein PDB (without ligand)
        
        Args:
            output_path: Optional custom output path. If None, uses default.
            
        Returns:
            Path to the output fixed PDB file
        """
        output_dir = pathlib.Path(self._output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if output_path is None:
            output_path = output_dir / "protein_fixed.pdb"
        else:
            output_path = pathlib.Path(output_path).resolve()
        
        # Parse flexible residue atoms
        flex_atoms_dict = self._parse_flex_residue_atoms()
        
        # If no flexible residues, just copy the receptor
        if not flex_atoms_dict:
            with open(self._receptor_pdb_path, 'r', encoding='utf-8') as f:
                receptor_content = f.read()
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(receptor_content)
            return output_path
        
        # Load receptor lines
        if self._receptor_lines is None:
            with open(self._receptor_pdb_path, 'r', encoding='utf-8') as f:
                self._receptor_lines = f.readlines()
        
        # Build corrected residues
        corrected_residues = {}  # {(chain, res_seq): [corrected_atoms]}
        
        for key, atoms in flex_atoms_dict.items():
            chain, res_seq, res_name = key
            print(f"\n  Processing flexible residue: {chain}:{res_seq} ({res_name})")
            print(f"    Input atoms from PDBQT ({len(atoms)} total):")
            for i, at in enumerate(atoms):
                print(f"      [{i}] '{at.get('atom_name', '?')}' elem={at.get('element', '?')} coord=({at['x']:.3f}, {at['y']:.3f}, {at['z']:.3f})")
            
            # Get backbone atoms from receptor for coordinate-based filtering
            receptor_backbone = self._get_receptor_residue_atoms(
                chain, res_seq, res_name, sidechain_only=False, backbone_only=True
            )
            
            # Filter backbone atoms using coordinate comparison
            sidechain_atoms = self._filter_backbone_duplicates(atoms, receptor_backbone)
            
            print(f"    After backbone filtering ({len(sidechain_atoms)} sidechain atoms):")
            for i, at in enumerate(sidechain_atoms):
                print(f"      [{i}] '{at.get('atom_name', '?')}' elem={at.get('element', '?')} coord=({at['x']:.3f}, {at['y']:.3f}, {at['z']:.3f})")
            
            if not sidechain_atoms:
                continue
            
            # Get reference atoms from receptor
            receptor_atoms = self._get_receptor_residue_atoms(
                chain, res_seq, res_name, sidechain_only=True
            )
            
            print(f"    Receptor reference sidechain atoms ({len(receptor_atoms)} atoms):")
            for i, at in enumerate(receptor_atoms):
                print(f"      [{i}] '{at.get('atom_name', '?')}' elem={at.get('element', '?')} coord=({at['x']:.3f}, {at['y']:.3f}, {at['z']:.3f})")
            
            if not receptor_atoms:
                continue
            
            # Match and correct atom names
            corrected = self._match_atoms_and_assign_names(
                sidechain_atoms, receptor_atoms
            )
            
            # Sort by receptor order
            corrected = self._sort_atoms_by_receptor_order(corrected, receptor_atoms)
            
            corrected_residues[(chain, res_seq)] = corrected
        
        # Build output: replace residues in receptor with corrected ones
        backbone_names = {"N", "CA", "C", "O"}
        new_lines = []
        current_residue_key = None
        residue_backbone_written = False
        
        for line in self._receptor_lines:
            atom = self._parse_pdb_line(line)
            
            if atom is None:
                # Non-ATOM line, keep as-is (TER, END, etc.)
                new_lines.append(line)
                continue
            
            key = (atom["chain_id"], atom["res_seq"])
            
            if key in corrected_residues:
                # This residue needs replacement
                if key != current_residue_key:
                    # New flexible residue, reset state
                    current_residue_key = key
                    residue_backbone_written = False
                
                # Keep backbone atoms from receptor
                if atom["atom_name"] in backbone_names:
                    new_lines.append(line)
                    
                    # After writing backbone, insert corrected sidechain
                    if atom["atom_name"] == "O" and not residue_backbone_written:
                        residue_backbone_written = True
                        # Get last atom serial
                        last_serial = atom["atom_serial"]
                        
                        for corrected_atom in corrected_residues[key]:
                            last_serial += 1
                            corrected_atom["atom_serial"] = last_serial
                            new_line = self._format_pdb_line(corrected_atom)
                            new_lines.append(new_line + "\n")
                # Skip original sidechain atoms (they will be replaced)
            else:
                # Not a flexible residue, keep as-is
                new_lines.append(line)
                current_residue_key = None
        
        # Write output
        with open(str(output_path), 'w', encoding='utf-8') as f:
            f.writelines(new_lines)
        
        return output_path
    
    def add_hydrogens(
        self, 
        input_pdb_path: str,
        output_pdb_path: typing.Optional[str] = None,
        cleanup_temp_files: bool = True
    ) -> str:
        """
        Add hydrogens to a PDB file using PDB2PQR with CHARMM force field.
        
        This method uses PDB2PQR to:
        1. Clean the PDB file using BioPython (remove existing H, fix format)
        2. Determine protonation states based on self.ph using PROPKA
        3. Add hydrogen atoms according to CHARMM force field
        4. Convert PQR back to PDB using BioPython
        
        Args:
            input_pdb_path: Path to input PDB file (e.g., from fix_pdb)
            output_pdb_path: Optional output path. If None, uses default naming.
            cleanup_temp_files: If True, remove temporary PQR file after processing.
                               If False, keep the PQR file for analysis. Default: True.
            
        Returns:
            Path to the output protonated PDB file
        """
        input_path = pathlib.Path(input_pdb_path).resolve()
        output_dir = pathlib.Path(self._output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if output_pdb_path is None:
            output_path = output_dir / "protein_protonated.pdb"
        else:
            output_path = pathlib.Path(output_pdb_path).resolve()
        
        # Temporary files
        temp_clean_pdb = output_path.with_name("temp_clean.pdb")
        temp_pqr = output_path.with_suffix(".pqr")
        
        # Step 1: Clean the PDB file using BioPython
        # - This fixes the AutoDock-style element types
        # - Also remove existing hydrogen atoms for re-protonation
        try:
            parser = Bio.PDB.PDBParser(QUIET=True)
            structure = parser.get_structure("protein", str(input_path))
            
            # Create a Select class to filter out hydrogens
            class NoHydrogen(Bio.PDB.Select):
                def accept_atom(self, atom):
                    # Check element first if available (most robust)
                    element = atom.element.strip().upper() if atom.element else ""
                    if element == "H":
                        return False
                    
                    # Fallback to name check
                    # Note: We don't use startswith("H") directly to avoid filtering out 
                    # atoms like Hg (Mercury), Hf (Hafnium), He (Helium), etc.
                    name = atom.get_name().strip().upper()
                    if name.startswith("H") and (len(name) == 1 or not name[1].isalpha()):
                        return False
                        
                    return True
            
            io = Bio.PDB.PDBIO()
            io.set_structure(structure)
            io.save(str(temp_clean_pdb), NoHydrogen())
            
        except Exception as e:
            print(f"Error cleaning PDB file: {e}")
            raise
        
        # Step 2: Run pdb2pqr on the cleaned PDB
        # Note: --titration-state-method propka is required to use PROPKA for pKa prediction
        # Without this, PDB2PQR will use standard protonation states regardless of pH
        # Use third_party_tools.run_pdb2pqr which has the propka monkey patch applied
        try:
            third_party_tools.run_pdb2pqr(
                input_pdb=str(temp_clean_pdb),
                output_pqr=str(temp_pqr),
                ph=self._ph,
                ff="CHARMM_WITH_ION",
                titration_state_method="propka",
                keep_chain=True
            )
            
            # Step 3: Convert PQR back to PDB using BioPython
            if temp_pqr.exists():
                parser = Bio.PDB.PDBParser(QUIET=True)
                structure = parser.get_structure("protein", str(temp_pqr))
                io = Bio.PDB.PDBIO()
                io.set_structure(structure)
                io.save(str(output_path))
                
                # Step 4: Analyze protonation states from PQR file
                print("\nAnalyzing protonation states from PDB2PQR output...")
                self.analyze_protonation_states(str(temp_pqr))
                
        except RuntimeError as e:
            print(f"Error running pdb2pqr: {e}")
            raise
        finally:
            # Clean up temporary files
            if temp_clean_pdb.exists():
                temp_clean_pdb.unlink()
            if cleanup_temp_files and temp_pqr.exists():
                temp_pqr.unlink()
        
        return str(output_path)

    
    def fix_histidine_names(
        self,
        input_pdb_path: str,
        output_pdb_path: typing.Optional[str] = None
    ) -> str:
        """
        Fix histidine residue names based on protonation state (CHARMM naming).
        
        Detects HIS protonation state by checking hydrogen atom names:
        - HD1 present, HE2 absent → HSD (Nδ protonated, neutral)
        - HE2 present, HD1 absent → HSE (Nε protonated, neutral)
        - Both HD1 and HE2 present → HSP (doubly protonated, +1)
        
        Args:
            input_pdb_path: Path to input PDB file (e.g., from add_hydrogens)
            output_pdb_path: Optional output path. If None, overwrites input file.
            
        Returns:
            Path to the output PDB file with corrected HIS names
        """
        input_path = pathlib.Path(input_pdb_path).resolve()
        
        if output_pdb_path is None:
            output_path = input_path
        else:
            output_path = pathlib.Path(output_pdb_path).resolve()
        
        # Read input PDB file
        with open(input_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # First pass: Identify HIS residues and their hydrogen atoms
        # Key: (chain_id, res_seq) -> {"has_HD1": bool, "has_HE2": bool}
        his_residues: dict = {}
        
        for line in lines:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            
            # Parse residue name and check if it's HIS
            res_name = line[17:20].strip()
            if res_name != "HIS":
                continue
            
            chain_id = line[21].strip()
            try:
                res_seq = int(line[22:26].strip())
            except ValueError:
                continue
            
            atom_name = line[12:16].strip()
            key = (chain_id, res_seq)
            
            if key not in his_residues:
                his_residues[key] = {"has_HD1": False, "has_HE2": False}
            
            if atom_name == "HD1":
                his_residues[key]["has_HD1"] = True
            elif atom_name == "HE2":
                his_residues[key]["has_HE2"] = True
        
        # Determine protonation state for each HIS
        # CHARMM naming: HSD (Nδ), HSE (Nε), HSP (both)
        his_new_names: dict = {}
        
        for key, state in his_residues.items():
            has_HD1 = state["has_HD1"]
            has_HE2 = state["has_HE2"]
            
            if has_HD1 and has_HE2:
                # Doubly protonated
                his_new_names[key] = "HSP"
            elif has_HD1:
                # Nδ protonated only
                his_new_names[key] = "HSD"
            elif has_HE2:
                # Nε protonated only
                his_new_names[key] = "HSE"
            else:
                # No protons on ring nitrogens (unusual), keep as HIS
                his_new_names[key] = "HIS"
        
        # Second pass: Rewrite lines with corrected residue names
        new_lines = []
        
        for line in lines:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                new_lines.append(line)
                continue
            
            res_name = line[17:20].strip()
            if res_name != "HIS":
                new_lines.append(line)
                continue
            
            chain_id = line[21].strip()
            try:
                res_seq = int(line[22:26].strip())
            except ValueError:
                new_lines.append(line)
                continue
            
            key = (chain_id, res_seq)
            new_name = his_new_names.get(key, "HIS")
            
            # Replace residue name in the line (columns 18-20, 1-indexed)
            # PDB format: columns 18-20 are residue name (right-justified in 3 chars)
            new_line = line[:17] + f"{new_name:>3}" + line[20:]
            new_lines.append(new_line)
        
        # Write output
        with open(str(output_path), 'w', encoding='utf-8') as f:
            f.writelines(new_lines)
        
        # Print summary
        for key, new_name in his_new_names.items():
            chain, resid = key
            print(f"  HIS {chain}:{resid} -> {new_name}")
        
        return str(output_path)
    
    def analyze_protonation_states(
        self,
        pqr_path: str,
        charge_threshold: float = 0.3
    ) -> typing.Dict[str, typing.List[typing.Tuple[str, int, str, float]]]:
        """
        Analyze protonation states of ionizable residues from PDB2PQR output PQR file.
        
        This method reads the PQR file (which contains partial charges) and calculates
        the total charge for each residue to determine its protonation state according
        to CHARMM force field conventions.
        
        CHARMM patches for non-standard protonation states:
        - CYM: Deprotonated CYS (negatively charged, ~-1.0) - for ionized cysteine
        - GLUP: Protonated GLU (neutral, ~0.0) - applies patch to add hydrogen to carboxyl
        - ASPP: Protonated ASP (neutral, ~0.0) - applies patch to add hydrogen to carboxyl
        - LSN: Deprotonated LYS (neutral, ~0.0) - removes charge from lysine
        
        Note: Standard states don't need patches:
        - CYS (protonated, neutral, ~0.0) - default
        - GLU (deprotonated, ~-1.0) - default
        - ASP (deprotonated, ~-1.0) - default  
        - LYS (protonated, ~+1.0) - default
        
        Args:
            pqr_path: Path to the PQR file from PDB2PQR
            charge_threshold: Threshold for determining charge state (default 0.3)
            
        Returns:
            Dictionary with CHARMM patch names as keys and residue lists as values
        """
        pqr_file = pathlib.Path(pqr_path).resolve()
        
        if not pqr_file.exists():
            print(f"Warning: PQR file not found: {pqr_file}")
            return {}
        
        # Read PQR file
        with open(pqr_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Parse PQR format: similar to PDB but with charge and radius in place of occupancy/B-factor
        # PQR format: ATOM serial name resName chainID resSeq x y z charge radius
        residue_charges = collections.defaultdict(lambda: {"atoms": [], "total_charge": 0.0})
        
        for line in lines:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            
            try:
                # Parse PQR line
                record_type = line[0:6].strip()
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                res_seq = int(line[22:26].strip())
                
                # PQR format: charge is after coordinates
                # Typically: x(8.3f) y(8.3f) z(8.3f) charge(8.4f) radius(7.4f)
                parts = line[30:].split()
                if len(parts) >= 5:
                    charge = float(parts[3])  # 4th element after coordinates
                else:
                    continue
                
                key = (res_name, chain_id, res_seq)
                residue_charges[key]["atoms"].append((atom_name, charge))
                residue_charges[key]["total_charge"] += charge
                
            except (ValueError, IndexError):
                continue
        
        # Analyze protonation states (CHARMM patch names)
        results = {
            'CYM': [],    # Deprotonated CYS (charge ~ -1, anionic)
            'GLUP': [],   # Protonated GLU (charge ~ 0, neutral)
            'ASPP': [],   # Protonated ASP (charge ~ 0, neutral)
            'LSN': [],    # Deprotonated LYS (charge ~ 0, neutral)
        }
        
        for key, data in residue_charges.items():
            res_name, chain_id, res_seq = key
            total_charge = data["total_charge"]
            
            # Check each residue type
            if res_name == "CYS":
                # CYS: neutral (~0) vs CYM (~-1)
                if total_charge < -charge_threshold:
                    results['CYM'].append((chain_id, res_seq, res_name, total_charge))
                    
            elif res_name == "GLU":
                # GLU: deprotonated (~-1, default) vs GLUP patch for protonated (~0, neutral)
                if abs(total_charge) < charge_threshold:
                    results['GLUP'].append((chain_id, res_seq, res_name, total_charge))
                    
            elif res_name == "ASP":
                # ASP: deprotonated (~-1, default) vs ASPP patch for protonated (~0, neutral)
                if abs(total_charge) < charge_threshold:
                    results['ASPP'].append((chain_id, res_seq, res_name, total_charge))
                    
            elif res_name == "LYS":
                # LYS: protonated (~+1, default) vs LSN patch for deprotonated (~0, neutral)
                if abs(total_charge) < charge_threshold:
                    results['LSN'].append((chain_id, res_seq, res_name, total_charge))
        
        # Store results in private variable
        self._protonation_patches = results
        
        # Print summary
        print("\n=== Protonation State Analysis (CHARMM Patches) ===")
        for patch_name, residues in results.items():
            if residues:
                print(f"\n{patch_name} patch needed for ({len(residues)} residue(s)):")
                for chain, resid, resname, charge in residues:
                    print(f"  {resname} {chain}:{resid} (total charge: {charge:+.3f})")
        
        if not any(results.values()):
            print("No unusual protonation states detected.")
            print("All CYS, GLU, ASP, LYS residues are in their standard states.")
        
        return results
    
    def detect_disulfide_bonds(
        self,
        pdb_path: typing.Optional[str] = None,
        distance_threshold: float = 2.5
    ) -> typing.List[typing.Tuple[tuple, tuple, float]]:
        """
        Detect disulfide bonds in the protein structure using MDAnalysis.
        
        A disulfide bond is identified when two cysteine (CYS) residues have their
        SG (sulfur gamma) atoms within the specified distance threshold.
        
        Args:
            pdb_path: Path to the PDB file to analyze. If None, uses receptor_pdb_path.
            distance_threshold: Maximum distance (in Angstroms) between SG atoms to 
                                consider a disulfide bond. Default is 2.5 Å.
        
        Returns:
            List of tuples, each containing:
                ((chain1, resid1, resname1), (chain2, resid2, resname2), distance)
            where chain is the chain ID, resid is the residue number, 
            resname is the residue name (typically 'CYS'), and distance is 
            the distance between the SG atoms in Angstroms.
        """

        
        
        # Use receptor PDB if no path is provided
        if pdb_path is None:
            pdb_path = self._receptor_pdb_path
        else:
            pdb_path = str(pathlib.Path(pdb_path).resolve())
        
        # Load the structure using MDAnalysis
        u = MDAnalysis.Universe(pdb_path)
        
        # Select all cysteine SG atoms
        # CYS residues might also be named CYX in some force fields when involved in disulfide bonds
        cys_sg = u.select_atoms("(resname CYS or resname CYX) and name SG")
        
        if len(cys_sg) == 0:
            print("No cysteine residues with SG atoms found in the structure.")
            self._disulfide_bonds = []
            return self._disulfide_bonds
        
        # Find disulfide bonds
        disulfide_bonds = []
        
        # Check all pairs of SG atoms
        for i in range(len(cys_sg)):
            for j in range(i + 1, len(cys_sg)):
                atom1 = cys_sg[i]
                atom2 = cys_sg[j]
                
                # Calculate distance between the two SG atoms
                distance = numpy.linalg.norm(atom1.position - atom2.position)
                
                # If distance is below threshold, it's a disulfide bond
                if distance <= distance_threshold:
                    # Extract residue information
                    res1_info = (
                        atom1.segid if atom1.segid else atom1.chainID,
                        atom1.resid,
                        atom1.resname
                    )
                    res2_info = (
                        atom2.segid if atom2.segid else atom2.chainID,
                        atom2.resid,
                        atom2.resname
                    )
                    
                    disulfide_bonds.append((res1_info, res2_info, distance))
        
        # Store in private variable
        self._disulfide_bonds = disulfide_bonds
        
        # Print summary
        if disulfide_bonds:
            print(f"Found {len(disulfide_bonds)} disulfide bond(s):")
            for bond in disulfide_bonds:
                res1, res2, dist = bond
                chain1, resid1, resname1 = res1
                chain2, resid2, resname2 = res2
                print(f"  {resname1} {chain1}:{resid1} - {resname2} {chain2}:{resid2} (distance: {dist:.2f} Å)")
        else:
            print("No disulfide bonds found.")
        
        return self._disulfide_bonds
    
    def generate_psf(
        self,
        input_pdb_path: str,
        output_psf_path: typing.Optional[str] = None,
        output_pdb_path: typing.Optional[str] = None,
        topology_file: typing.Optional[str] = None
    ) -> typing.Tuple[str, str]:
        """
        Generate PSF file using VMD psfgen plugin from a PDB file.
        
        This method:
        1. Analyzes the PDB structure to identify segments and residues
        2. Detects disulfide bonds if not already detected
        3. Identifies C-terminal residues and determines if CTFS patch is needed
        4. Generates a TCL script for psfgen
        5. Executes VMD to run psfgen and generate the PSF and final PDB files
        
        Args:
            input_pdb_path: Path to input PDB file (should be protonated with correct residue names)
            output_psf_path: Optional output PSF path. If None, uses default naming.
            output_pdb_path: Optional output PDB path. If None, uses default naming.
            topology_file: Path to topology RTF file. Default is third_party/force_field/top_all36_prot.rtf
            
        Returns:
            Tuple of (psf_path, pdb_path) for the generated files
        """
        # Convert all paths to absolute pathlib.Path objects
        input_path = pathlib.Path(input_pdb_path).resolve()
        output_dir = pathlib.Path(self._output_dir).resolve()
        
        # Ensure output directory exists
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create temp directory for intermediate files
        temp_dir = output_dir / "md_prepper" / "temp_files"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Set output paths
        if output_psf_path is None:
            output_psf_path = output_dir / "protein_final.psf"
        else:
            output_psf_path = pathlib.Path(output_psf_path).resolve()
        
        if output_pdb_path is None:
            output_pdb_path = output_dir / "protein_final.pdb"
        else:
            output_pdb_path = pathlib.Path(output_pdb_path).resolve()
        
        # Convert topology file to absolute path (use default if not provided)
        if topology_file is None:
            topology_path = third_party_tools.get_protein_topology_file()
        else:
            topology_path = pathlib.Path(topology_file).resolve()
        
        # Detect disulfide bonds if not already done
        if self._disulfide_bonds is None:
            print("Detecting disulfide bonds...")
            self.detect_disulfide_bonds(pdb_path=str(input_path))
        
        # Parse PDB to get segment information
        segments_info = self._analyze_pdb_structure(str(input_path))
        
        # Split PDB by chain/segment for multi-chain proteins
        print("Splitting PDB by chain/segment...")
        split_pdb_paths = self._split_pdb_by_chain(
            pdb_path=str(input_path),
            output_dir=temp_dir,
            segments_info=segments_info
        )
        
        # Generate TCL script with split PDB files
        tcl_script_path = temp_dir / "generate_psf.tcl"
        self._generate_psfgen_tcl(
            tcl_script_path=str(tcl_script_path),
            topology_file=str(topology_path),
            split_pdb_paths=split_pdb_paths,
            output_psf_path=str(output_psf_path),
            output_pdb_path=str(output_pdb_path),
            segments_info=segments_info
        )
        
        # Execute VMD with psfgen TCL script
        print(f"Executing VMD psfgen...")
        print(f"  VMD Executable: {self._vmd_executable}")
        print(f"  TCL Script: {tcl_script_path}")
        
        stdout = third_party_tools.run_vmd(
            vmd_executable=self._vmd_executable,
            tcl_script=tcl_script_path,
            cwd=output_dir
        )
        print("VMD psfgen execution successful!")
        if stdout:
            print(f"STDOUT:\n{stdout}")
        
        return (str(output_psf_path), str(output_pdb_path))
    
    def _analyze_pdb_structure(self, pdb_path: str) -> dict:
        """
        Analyze PDB structure to extract segment information.
        
        Args:
            pdb_path: Path to PDB file
            
        Returns:
            Dictionary with segment information:
            {
                segment_id: {
                    'chain': chain_id,
                    'residues': [(resid, resname), ...],
                    'first_residue': (resid, resname),
                    'last_residue': (resid, resname)
                }
            }
        """
        analyze_pdb_path = pathlib.Path(pdb_path).resolve()
        with open(analyze_pdb_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Group residues by chain/segment
        segments = collections.defaultdict(lambda: {'chain': '', 'residues': []})
        residue_set = set()  # To track unique residues per segment
        
        for line in lines:
            atom = self._parse_pdb_line(line)
            if atom is None:
                continue
            
            # Use chain_id as segment_id if segment_id is empty
            segment_id = atom.get('segment_id', '').strip()
            chain_id = atom['chain_id'].strip() if atom['chain_id'] else ''
            
            # Handle empty segment_id and chain_id
            if not segment_id and not chain_id:
                print(f"  Warning: Atom at res {res_seq} has no chain/segment ID, defaulting to 'PROT'")
            if not segment_id:
                segment_id = chain_id if chain_id else 'PROT'
            if not chain_id:
                chain_id = segment_id
            
            res_seq = atom['res_seq']
            res_name = atom['res_name']
            
            # Track unique residues
            res_key = (segment_id, res_seq, res_name)
            if res_key not in residue_set:
                residue_set.add(res_key)
                segments[segment_id]['residues'].append((res_seq, res_name))
                segments[segment_id]['chain'] = chain_id
        
        # Sort residues and identify first/last
        for seg_id in segments:
            segments[seg_id]['residues'].sort(key=lambda x: x[0])
            if segments[seg_id]['residues']:
                segments[seg_id]['first_residue'] = segments[seg_id]['residues'][0]
                segments[seg_id]['last_residue'] = segments[seg_id]['residues'][-1]
        
        return dict(segments)
    
    # Common water molecule residue names
    _WATER_RESNAMES = {'HOH', 'WAT', 'TIP', 'TIP3', 'TIP4', 'TIP5', 'SOL', 'H2O', 'DOD', 'OHH'}
    
    # Common metal ion residue names (single-atom residues)
    _METAL_RESNAMES = {
        'ZN', 'ZN2', 'MG', 'CA', 'CAL', 'FE', 'FE2', 'FE3', 'CU', 'CU1', 'CU2',
        'MN', 'MN2', 'MN3', 'CO', 'CO2', 'NI', 'NI2', 'NA', 'SOD', 'K', 'POT',
        'CD', 'CD2', 'HG', 'HG2', 'PB', 'PB2', 'SR', 'BA', 'BAR', 'AG', 'AU',
        'LI', 'LIT', 'RB', 'RUB', 'CS', 'CES', 'AL', 'AL3', 'CR', 'CR3',
        'MO', 'W', 'V', 'SE', 'AS', 'CL', 'CLA'
    }
    
    # Mapping from common PDB metal/ion residue names to CHARMM (resname, atomname) tuples
    # Based on toppar_water_ions.str definitions
    _METAL_CHARMM_MAPPING = {
        # Lithium
        'LI': ('LIT', 'LIT'),
        'LIT': ('LIT', 'LIT'),
        # Sodium
        'NA': ('SOD', 'SOD'),
        'SOD': ('SOD', 'SOD'),
        # Magnesium
        'MG': ('MG', 'MG'),
        # Potassium
        'K': ('POT', 'POT'),
        'POT': ('POT', 'POT'),
        # Calcium
        'CA': ('CAL', 'CAL'),
        'CAL': ('CAL', 'CAL'),
        # Rubidium
        'RB': ('RUB', 'RUB'),
        'RUB': ('RUB', 'RUB'),
        # Cesium
        'CS': ('CES', 'CES'),
        'CES': ('CES', 'CES'),
        # Barium
        'BA': ('BAR', 'BAR'),
        'BAR': ('BAR', 'BAR'),
        # Zinc (II) - note: resname is ZN2, atomname is ZN
        'ZN': ('ZN2', 'ZN'),
        'ZN2': ('ZN2', 'ZN'),
        # Cadmium (II) - note: resname is CD2, atomname is CD
        'CD': ('CD2', 'CD'),
        'CD2': ('CD2', 'CD'),
        # Chloride
        'CL': ('CLA', 'CLA'),
        'CLA': ('CLA', 'CLA'),
    }
    
    def _is_water_residue(self, resname: str) -> bool:
        """Check if a residue name corresponds to a water molecule."""
        return resname.upper() in self._WATER_RESNAMES
    
    def _is_metal_residue(self, resname: str) -> bool:
        """Check if a residue name corresponds to a metal ion."""
        return resname.upper() in self._METAL_RESNAMES
    
    def _get_charmm_metal_names(self, resname: str) -> typing.Tuple[str, str]:
        """
        Get CHARMM-compatible residue name and atom name for a metal ion.
        
        Args:
            resname: Original residue name from PDB
            
        Returns:
            Tuple of (charmm_resname, charmm_atomname)
            If not found in mapping, returns the original name for both
        """
        resname_upper = resname.upper()
        if resname_upper in self._METAL_CHARMM_MAPPING:
            return self._METAL_CHARMM_MAPPING[resname_upper]
        else:
            # For metals not in mapping (like FE, CU, etc.), keep original name
            # User may need to add custom topology for these
            print(f"  Warning: Metal ion '{resname}' not in CHARMM mapping, keeping original name")
            return (resname_upper, resname_upper)
    
    def _rename_water_atom(self, atom_name: str, element: str) -> str:
        """
        Rename water atom to TIP3 format.
        
        Original names can be: O/OW/OH2 -> OH2, H1/HW1/H -> H1, H2/HW2/H -> H2
        """
        atom_name_upper = atom_name.strip().upper()
        element_upper = element.strip().upper() if element else ''
        
        # Oxygen atom
        if element_upper == 'O' or atom_name_upper in ('O', 'OW', 'OH2', 'OT', 'O1'):
            return 'OH2'
        # Hydrogen atoms
        elif element_upper == 'H' or atom_name_upper.startswith('H'):
            # Try to determine if it's H1 or H2 based on original name
            if '1' in atom_name_upper or atom_name_upper == 'H':
                return 'H1'
            elif '2' in atom_name_upper:
                return 'H2'
            else:
                # Default to H1 if we can't determine
                return 'H1'
        return atom_name
    
    def _process_water_molecules(
        self, 
        water_atoms: typing.List[dict], 
        output_dir: pathlib.Path
    ) -> typing.Optional[str]:
        """
        Process water molecules: rename to TIP3 format and output to a single PDB file.
        
        Each water molecule gets a unique resid, and atoms are renamed to OH2, H1, H2.
        
        Args:
            water_atoms: List of atom dictionaries for water molecules
            output_dir: Directory to write the output PDB file
            
        Returns:
            Path to the water PDB file, or None if no water molecules
        """
        if not water_atoms:
            return None
        
        # Group atoms by original residue
        water_residues: typing.Dict[typing.Tuple[str, int], typing.List[dict]] = collections.defaultdict(list)
        for atom in water_atoms:
            key = (atom['chain_id'], atom['res_seq'])
            water_residues[key].append(atom)
        
        # Process each water molecule
        processed_atoms = []
        new_resid = 0
        
        for res_key in sorted(water_residues.keys()):
            atoms = water_residues[res_key]
            new_resid += 1
            
            # Track which hydrogen positions are used
            h_count = 0
            
            for atom in atoms:
                new_atom = atom.copy()
                new_atom['res_name'] = 'TIP3'
                new_atom['res_seq'] = new_resid
                
                # Rename atom based on element
                element = atom.get('element', '')
                old_name = atom['atom_name'].strip().upper()
                
                if element.upper() == 'O' or old_name in ('O', 'OW', 'OH2', 'OT', 'O1'):
                    new_atom['atom_name'] = 'OH2'
                    new_atom['element'] = 'O'
                elif element.upper() == 'H' or old_name.startswith('H'):
                    h_count += 1
                    new_atom['atom_name'] = f'H{h_count}'
                    new_atom['element'] = 'H'
                
                processed_atoms.append(new_atom)
        
        if not processed_atoms:
            return None
        
        # Write water PDB file
        water_pdb_path = output_dir / "segment_SOLV.pdb"
        
        with open(water_pdb_path, 'w', encoding='utf-8') as f:
            atom_serial = 0
            for atom in processed_atoms:
                atom_serial += 1
                atom['atom_serial'] = atom_serial
                atom['segment_id'] = 'SOLV'
                line = self._format_pdb_line(atom)
                f.write(line + '\n')
            f.write("END\n")
        
        print(f"  Split water molecules: {len(water_residues)} waters ({len(processed_atoms)} atoms) -> segment_SOLV.pdb")
        
        return str(water_pdb_path)
    
    def _process_metal_ions(
        self, 
        metal_atoms: typing.List[dict], 
        output_dir: pathlib.Path
    ) -> typing.Dict[str, str]:
        """
        Process metal ions: each metal ion gets its own PDB file with unique segment.
        
        Metal ion residue names and atom names are converted to CHARMM-compatible format.
        
        Args:
            metal_atoms: List of atom dictionaries for metal ions
            output_dir: Directory to write the output PDB files
            
        Returns:
            Dictionary mapping segment_id to the path of its PDB file
        """
        metal_pdb_paths = {}
        
        if not metal_atoms:
            return metal_pdb_paths
        
        # Group metal atoms by their original residue (chain, resid)
        metal_residues: typing.Dict[typing.Tuple[str, int, str], typing.List[dict]] = collections.defaultdict(list)
        for atom in metal_atoms:
            key = (atom['chain_id'], atom['res_seq'], atom['res_name'])
            metal_residues[key].append(atom)
        
        # Process each metal ion
        metal_count = 0
        for res_key in sorted(metal_residues.keys()):
            atoms = metal_residues[res_key]
            chain_id, res_seq, res_name = res_key
            metal_count += 1
            
            # Get CHARMM-compatible names
            charmm_resname, charmm_atomname = self._get_charmm_metal_names(res_name)
            
            # Create segment ID for this metal ion (e.g., ION1, ION2, ...)
            seg_id = f"ION{metal_count}"
            
            # Write metal PDB file
            metal_pdb_path = output_dir / f"segment_{seg_id}.pdb"
            
            with open(metal_pdb_path, 'w', encoding='utf-8') as f:
                atom_serial = 0
                for atom in atoms:
                    atom_serial += 1
                    new_atom = atom.copy()
                    new_atom['atom_serial'] = atom_serial
                    new_atom['segment_id'] = seg_id
                    new_atom['res_seq'] = 1  # Single residue in segment
                    new_atom['res_name'] = charmm_resname  # CHARMM residue name
                    new_atom['atom_name'] = charmm_atomname  # CHARMM atom name
                    line = self._format_pdb_line(new_atom)
                    f.write(line + '\n')
                f.write("END\n")
            
            metal_pdb_paths[seg_id] = str(metal_pdb_path)
            if res_name.upper() != charmm_resname:
                print(f"  Split metal ion {res_name} -> {charmm_resname} (Chain {chain_id}, Res {res_seq}): -> segment_{seg_id}.pdb")
            else:
                print(f"  Split metal ion {charmm_resname} (Chain {chain_id}, Res {res_seq}): -> segment_{seg_id}.pdb")
        
        return metal_pdb_paths
    
    def _split_pdb_by_chain(
        self, 
        pdb_path: str, 
        output_dir: pathlib.Path,
        segments_info: dict
    ) -> typing.Dict[str, str]:
        """
        Split a multi-chain PDB file into separate PDB files, one per segment/chain.
        Also handles water molecules and metal ions specially.
        
        This is necessary because psfgen's `pdb` command reads all atoms from a file
        into a segment. For multi-chain proteins, we need separate files for each chain.
        
        Water molecules are:
        - Converted to TIP3 residue name
        - Each water molecule gets a unique resid
        - Atoms are renamed to OH2, H1, H2
        - All waters are output to a single "SOLV" segment PDB
        
        Metal ions are:
        - Each metal ion gets its own segment (ION1, ION2, ...)
        - Each metal is output to a separate PDB file
        
        Args:
            pdb_path: Path to the input PDB file
            output_dir: Directory to write the split PDB files
            segments_info: Segment information from _analyze_pdb_structure
            
        Returns:
            Dictionary mapping segment_id to the path of its split PDB file
        """
        split_pdb_input_path = pathlib.Path(pdb_path).resolve()
        
        with open(split_pdb_input_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Group lines by segment/chain, and collect water/metal atoms separately
        segment_lines: typing.Dict[str, typing.List[str]] = {seg_id: [] for seg_id in segments_info}
        water_atoms: typing.List[dict] = []
        metal_atoms: typing.List[dict] = []
        
        for line in lines:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                # Keep non-ATOM lines like TER, END for all segments
                continue
            
            atom = self._parse_pdb_line(line)
            if atom is None:
                continue
            
            res_name = atom['res_name'].strip().upper()
            
            # Check if this is a water molecule
            if self._is_water_residue(res_name):
                water_atoms.append(atom)
                continue
            
            # Check if this is a metal ion
            if self._is_metal_residue(res_name):
                metal_atoms.append(atom)
                continue
            
            # Regular protein/other atoms - determine which segment
            segment_id = atom.get('segment_id', '').strip()
            chain_id = atom['chain_id'].strip() if atom['chain_id'] else ''
            
            if not segment_id:
                segment_id = chain_id if chain_id else 'PROT'
            
            if segment_id in segment_lines:
                segment_lines[segment_id].append(line)
        
        # Write separate PDB files for each protein segment
        split_pdb_paths = {}
        for seg_id, atom_lines in segment_lines.items():
            if not atom_lines:
                print(f"  Warning: No atoms found for segment {seg_id}")
                continue
            
            # Create output file path
            split_pdb_path = output_dir / f"segment_{seg_id}.pdb"
            
            with open(split_pdb_path, 'w', encoding='utf-8') as f:
                for line in atom_lines:
                    f.write(line)
                f.write("END\n")
            
            split_pdb_paths[seg_id] = str(split_pdb_path)
            print(f"  Split segment {seg_id}: {len(atom_lines)} atoms -> {split_pdb_path.name}")
        
        # Process water molecules
        water_pdb_path = self._process_water_molecules(water_atoms, output_dir)
        if water_pdb_path:
            split_pdb_paths['SOLV'] = water_pdb_path
        
        # Process metal ions
        metal_pdb_paths = self._process_metal_ions(metal_atoms, output_dir)
        split_pdb_paths.update(metal_pdb_paths)
        
        return split_pdb_paths
    
    def _generate_psfgen_tcl(
        self,
        tcl_script_path: str,
        topology_file: str,
        split_pdb_paths: typing.Dict[str, str],
        output_psf_path: str,
        output_pdb_path: str,
        segments_info: dict
    ):
        """
        Generate TCL script for psfgen.
        
        Args:
            tcl_script_path: Path to output TCL script
            topology_file: Path to topology RTF file
            split_pdb_paths: Dictionary mapping segment_id to path of its PDB file
            output_psf_path: Path to output PSF file
            output_pdb_path: Path to output PDB file
            segments_info: Segment information from _analyze_pdb_structure
        """
        # Convert all paths to POSIX format (forward slashes) for TCL
        topology_posix = pathlib.Path(topology_file).as_posix()
        output_psf_posix = pathlib.Path(output_psf_path).as_posix()
        output_pdb_posix = pathlib.Path(output_pdb_path).as_posix()
        
        # Get the path to toppar_water_ions.str for water and ion topology
        water_ions_topology = third_party_tools.get_third_party_dir() / "force_field" / "toppar_water_ions.str"
        water_ions_topology_posix = water_ions_topology.as_posix()
        
        # Check if we have water or metal segments
        has_water = 'SOLV' in split_pdb_paths
        metal_segments = [seg for seg in split_pdb_paths.keys() if seg.startswith('ION')]
        has_metals = len(metal_segments) > 0
        
        lines = []
        
        # Header
        lines.append("# PSF Generation Script")
        lines.append("# Generated by MDPrepper")
        lines.append("")
        
        # Load psfgen plugin in VMD
        lines.append("package require psfgen")
        lines.append("")
        
        # Load topology - always load protein topology first
        lines.append("# Load topology files")
        lines.append(f"topology {topology_posix}")
        
        # Load water and ions topology if we have water or metal segments
        if has_water or has_metals:
            lines.append(f"topology {water_ions_topology_posix}")
            print(f"  Loading water/ions topology: {water_ions_topology.name}")
        lines.append("")
        
        # Create a mapping from Chain ID to Segment ID for looking up correct segments for patches
        chain_to_seg_map = {}
        
        # Set up pdbalias for common residue name mappings
        lines.append("# Set up residue aliases")
        lines.append("pdbalias residue HIS HSE  ; # Alias HIS to default tautomer if not already renamed")
        lines.append("")
        
        # Build chain_to_seg_map
        for seg_id, info in segments_info.items():
            chain_id = info['chain']
            chain_to_seg_map[chain_id] = seg_id

        # Process each protein segment - each segment uses its own split PDB file
        lines.append("# Protein segments")
        for seg_id, info in segments_info.items():
            chain_id = info['chain']
            
            # Get the split PDB file for this segment
            if seg_id not in split_pdb_paths:
                print(f"  Warning: No split PDB file for segment {seg_id}, skipping")
                continue
            
            seg_pdb_posix = pathlib.Path(split_pdb_paths[seg_id]).as_posix()
            
            lines.append(f"# Segment {seg_id} (Chain {chain_id})")
            lines.append(f"segment {seg_id} {{")
            lines.append(f"  first NTER")
            lines.append(f"  last CTER")
            lines.append(f"  pdb {seg_pdb_posix}")
            lines.append("}")
            lines.append("")
        
        # Process water segment if present
        if has_water:
            water_pdb_posix = pathlib.Path(split_pdb_paths['SOLV']).as_posix()
            lines.append("# Water segment")
            lines.append("segment SOLV {")
            lines.append("  first NONE")
            lines.append("  last NONE")
            lines.append(f"  pdb {water_pdb_posix}")
            lines.append("}")
            lines.append("")
            print(f"  Creating water segment: SOLV")
        
        # Process metal ion segments if present
        if has_metals:
            lines.append("# Metal ion segments")
            for seg_id in sorted(metal_segments):
                metal_pdb_posix = pathlib.Path(split_pdb_paths[seg_id]).as_posix()
                lines.append(f"# Metal ion segment {seg_id}")
                lines.append(f"segment {seg_id} {{")
                lines.append("  first NONE")
                lines.append("  last NONE")
                lines.append(f"  pdb {metal_pdb_posix}")
                lines.append("}")
                lines.append("")
                print(f"  Creating metal ion segment: {seg_id}")
        
        # Coordinate reading - each segment reads from its own split PDB file
        lines.append("# Read coordinates from split PDB files")
        
        # Protein segments
        for seg_id in segments_info.keys():
            if seg_id not in split_pdb_paths:
                continue
            seg_pdb_posix = pathlib.Path(split_pdb_paths[seg_id]).as_posix()
            lines.append(f"coordpdb {seg_pdb_posix} {seg_id}")
        
        # Water segment
        if has_water:
            water_pdb_posix = pathlib.Path(split_pdb_paths['SOLV']).as_posix()
            lines.append(f"coordpdb {water_pdb_posix} SOLV")
        
        # Metal ion segments
        for seg_id in sorted(metal_segments):
            metal_pdb_posix = pathlib.Path(split_pdb_paths[seg_id]).as_posix()
            lines.append(f"coordpdb {metal_pdb_posix} {seg_id}")
        
        lines.append("")
        
        # Apply protonation state patches
        if self._protonation_patches:
            lines.append("# Apply protonation state patches")
            for patch_name, residues in self._protonation_patches.items():
                if residues:
                    lines.append(f"# {patch_name} patch")
                    for chain, resid, resname, charge in residues:
                        # Find the correct segment for this chain
                        # If chain not found (shouldn't happen), default to first segment or use chain as fallback
                        seg = chain_to_seg_map.get(chain, list(segments_info.keys())[0])
                        
                        lines.append(f"patch {patch_name} {seg}:{resid}")
                        print(f"  Applying {patch_name} patch to {resname} {seg}:{resid} (charge: {charge:+.3f})")
            lines.append("")
        
        # Apply disulfide bond patches
        if self._disulfide_bonds:
            lines.append("# Apply disulfide bond patches")
            for bond in self._disulfide_bonds:
                res1, res2, dist = bond
                chain1, resid1, resname1 = res1
                chain2, resid2, resname2 = res2
                
                # Get correct segment IDs
                seg1 = chain_to_seg_map.get(chain1, list(segments_info.keys())[0])
                seg2 = chain_to_seg_map.get(chain2, list(segments_info.keys())[0])
                
                lines.append(f"patch DISU {seg1}:{resid1} {seg2}:{resid2}")
                print(f"  Disulfide bond: {resname1} {seg1}:{resid1} - {resname2} {seg2}:{resid2} (distance: {dist:.2f} Å)")
            lines.append("")
        
        # Regenerate angles and dihedrals
        lines.append("# Regenerate angles and dihedrals")
        lines.append("regenerate angles dihedrals")
        lines.append("")
        
        # Guess missing coordinates (if any)
        lines.append("# Guess missing coordinates")
        lines.append("guesscoord")
        lines.append("")
        
        # Write output files
        lines.append("# Write PSF and PDB files")
        lines.append(f"writepsf {output_psf_posix}")
        lines.append(f"writepdb {output_pdb_posix}")
        lines.append("")
        
        lines.append("# Done")
        lines.append("exit")
        
        # Write TCL script
        with open(tcl_script_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        
        print(f"TCL script generated: {tcl_script_path}")
    
    def _generate_ligand_psfgen_tcl(
        self,
        tcl_script_path: str,
        topology_files: typing.List[str],
        input_pdb_path: str,
        output_psf_path: str,
        output_pdb_path: str,
        segment_id: str = "LIG",
        residue_name: str = "UNL"
    ):
        """
        Generate TCL script for psfgen to create ligand PSF.
        
        This is a simplified version of _generate_psfgen_tcl for small molecules.
        It supports loading multiple topology files including CGenFF .str files.
        
        Args:
            tcl_script_path: Path to output TCL script
            topology_files: List of paths to topology files (RTF and/or STR files)
            input_pdb_path: Path to input ligand PDB file
            output_psf_path: Path to output PSF file
            output_pdb_path: Path to output PDB file
            segment_id: Segment ID for the ligand (default: "LIG")
            residue_name: Residue name in the PDB/topology (default: "UNL")
        """
        lines = []
        
        # Header
        lines.append("# Ligand PSF Generation Script")
        lines.append("# Generated by MDPrepper")
        lines.append("")
        
        # Load psfgen plugin in VMD
        lines.append("package require psfgen")
        lines.append("")
        
        # Load topology files (convert to POSIX format for TCL)
        lines.append("# Load topology files")
        for topo_file in topology_files:
            topo_posix = pathlib.Path(topo_file).as_posix()
            lines.append(f"topology {topo_posix}")
        lines.append("")
        
        # Convert paths to POSIX format
        input_pdb_posix = pathlib.Path(input_pdb_path).as_posix()
        output_psf_posix = pathlib.Path(output_psf_path).as_posix()
        output_pdb_posix = pathlib.Path(output_pdb_path).as_posix()
        
        # Create segment for ligand
        lines.append(f"# Create ligand segment")
        lines.append(f"segment {segment_id} {{")
        lines.append(f"  pdb {input_pdb_posix}")
        lines.append("}")
        lines.append("")
        
        # Read coordinates
        lines.append("# Read coordinates from PDB")
        lines.append(f"coordpdb {input_pdb_posix} {segment_id}")
        lines.append("")
        
        # Regenerate angles and dihedrals
        lines.append("# Regenerate angles and dihedrals")
        lines.append("regenerate angles dihedrals")
        lines.append("")
        
        # Guess missing coordinates (if any, e.g., hydrogens)
        lines.append("# Guess missing coordinates")
        lines.append("guesscoord")
        lines.append("")
        
        # Write output files
        lines.append("# Write PSF and PDB files")
        lines.append(f"writepsf {output_psf_posix}")
        lines.append(f"writepdb {output_pdb_posix}")
        lines.append("")
        
        lines.append("# Done")
        lines.append("exit")
        
        # Write TCL script
        with open(tcl_script_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        
        print(f"Ligand TCL script generated: {tcl_script_path}")
    
    def generate_ligand_psf(
        self,
        input_pdb_path: str,
        topology_files: typing.Optional[typing.List[str]] = None,
        output_psf_path: typing.Optional[str] = None,
        output_pdb_path: typing.Optional[str] = None,
        segment_id: str = "LIG"
    ) -> typing.Tuple[str, str]:
        """
        Generate PSF file for a ligand using VMD psfgen plugin.
        
        This method generates a PSF file for small molecules using CGenFF topology.
        It reads topology from RTF and STR files and creates the PSF/PDB pair.
        
        Args:
            input_pdb_path: Path to input ligand PDB file
            topology_files: List of topology files to load. If None, uses default:
                            [top_all36_cgenff.rtf, ligand.str]
            output_psf_path: Optional output PSF path. If None, uses default naming.
            output_pdb_path: Optional output PDB path. If None, uses default naming.
            segment_id: Segment ID for the ligand (default: "LIG")
            
        Returns:
            Tuple of (psf_path, pdb_path) for the generated files
        """
        # Convert paths to absolute pathlib.Path objects
        input_path = pathlib.Path(input_pdb_path).resolve()
        output_dir = pathlib.Path(self._output_dir).resolve()
        
        # Ensure output directory exists
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create temp directory for intermediate files
        temp_dir = output_dir / "md_prepper" / "temp_files"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Set default topology files if not provided
        if topology_files is None:
            topology_files = [
                str(third_party_tools.get_cgenff_topology_file()),
                str(output_dir / "ligand.str")
            ]
        
        # Convert all topology files to absolute paths
        topology_paths = [pathlib.Path(f).resolve() for f in topology_files]
        
        # Verify topology files exist
        for topo_path in topology_paths:
            if not topo_path.exists():
                raise FileNotFoundError(f"Topology file not found: {topo_path}")
        
        # Set output paths
        if output_psf_path is None:
            output_psf_path = output_dir / "ligand_final.psf"
        else:
            output_psf_path = pathlib.Path(output_psf_path).resolve()
        
        if output_pdb_path is None:
            output_pdb_path = output_dir / "ligand_final.pdb"
        else:
            output_pdb_path = pathlib.Path(output_pdb_path).resolve()
        
        print(f"Generating ligand PSF file...")
        print(f"  Input PDB: {input_path}")
        print(f"  Topology files:")
        for topo_path in topology_paths:
            print(f"    - {topo_path}")
        
        # Generate TCL script
        tcl_script_path = temp_dir / "generate_ligand_psf.tcl"
        self._generate_ligand_psfgen_tcl(
            tcl_script_path=str(tcl_script_path),
            topology_files=[str(p) for p in topology_paths],
            input_pdb_path=str(input_path),
            output_psf_path=str(output_psf_path),
            output_pdb_path=str(output_pdb_path),
            segment_id=segment_id
        )
        
        # Execute VMD with psfgen TCL script
        print(f"Executing VMD psfgen...")
        print(f"  VMD Executable: {self._vmd_executable}")
        print(f"  TCL Script: {tcl_script_path}")
        
        stdout = third_party_tools.run_vmd(
            vmd_executable=self._vmd_executable,
            tcl_script=tcl_script_path,
            cwd=output_dir
        )
        print("VMD psfgen execution successful!")
        if stdout:
            print(f"STDOUT:\n{stdout}")
        
        print(f"  Output PSF: {output_psf_path}")
        print(f"  Output PDB: {output_pdb_path}")
        
        return (str(output_psf_path), str(output_pdb_path))
    
    def _generate_merge_tcl(
        self,
        tcl_script_path: str,
        psf_pdb_pairs: typing.List[typing.Tuple[str, str]],
        output_psf_path: str,
        output_pdb_path: str
    ):
        """
        Generate TCL script for merging multiple PSF/PDB pairs using VMD TopoTools.
        
        Args:
            tcl_script_path: Path to output TCL script
            psf_pdb_pairs: List of (psf_path, pdb_path) tuples to merge
            output_psf_path: Path for merged output PSF file
            output_pdb_path: Path for merged output PDB file
        """
        lines = []
        
        # Header
        lines.append("# PSF/PDB Merge Script")
        lines.append("# Generated by MDPrepper")
        lines.append("# Using TopoTools to merge multiple PSF/PDB pairs")
        lines.append("")
        
        # Load required packages
        lines.append("# Load required packages")
        lines.append("package require topotools 1.6")
        lines.append("")
        
        # Convert output paths to POSIX format for TCL
        output_psf_posix = pathlib.Path(output_psf_path).as_posix()
        output_pdb_posix = pathlib.Path(output_pdb_path).as_posix()
        
        # Initialize molecule list
        lines.append("# Load to be merged molecules into VMD")
        lines.append("set midlist {}")
        lines.append("")
        
        # Load each PSF/PDB pair and append to list
        for i, (psf_path, pdb_path) in enumerate(psf_pdb_pairs):
            psf_posix = pathlib.Path(psf_path).as_posix()
            pdb_posix = pathlib.Path(pdb_path).as_posix()
            
            lines.append(f"# Load pair {i + 1}: {pathlib.Path(psf_path).name}")
            lines.append(f'set mol [mol new "{psf_posix}" waitfor all]')
            lines.append(f'mol addfile "{pdb_posix}" $mol')
            lines.append("lappend midlist $mol")
            lines.append("")
        
        # Merge all molecules using TopoTools
        lines.append("# Do the magic - merge all molecules")
        lines.append("set mol [::TopoTools::mergemols $midlist]")
        lines.append("")
        
        # Write merged output files
        lines.append("# Write merged PSF and PDB files")
        lines.append(f'animate write psf "{output_psf_posix}" $mol')
        lines.append(f'animate write pdb "{output_pdb_posix}" $mol')
        lines.append("")
        
        lines.append("# Done")
        lines.append("exit")
        
        # Write TCL script
        with open(tcl_script_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        
        print(f"Merge TCL script generated: {tcl_script_path}")
    
    
    def merge_psf_pdb(
        self,
        psf_pdb_pairs: typing.List[typing.Tuple[str, str]],
        output_psf_path: typing.Optional[str] = None,
        output_pdb_path: typing.Optional[str] = None
    ) -> typing.Tuple[str, str]:
        """
        Merge multiple PSF/PDB file pairs into a single PSF/PDB pair using VMD TopoTools.
        
        This method is useful for combining protein and ligand structures, or
        merging multiple molecules into a single system for MD simulation.
        
        Args:
            psf_pdb_pairs: List of (psf_path, pdb_path) tuples to merge.
                          Example: [("protein.psf", "protein.pdb"), ("ligand.psf", "ligand.pdb")]
            output_psf_path: Optional output PSF path. If None, uses "complex.psf" in output_dir.
            output_pdb_path: Optional output PDB path. If None, uses "complex.pdb" in output_dir.
            
        Returns:
            Tuple of (merged_psf_path, merged_pdb_path) for the generated files
            
        Raises:
            ValueError: If less than 2 PSF/PDB pairs are provided
            FileNotFoundError: If any input PSF or PDB file doesn't exist
        """
        # Validate input
        if len(psf_pdb_pairs) < 2:
            raise ValueError("At least 2 PSF/PDB pairs are required for merging")
        
        # Convert paths to absolute pathlib.Path objects and verify existence
        validated_pairs = []
        for i, (psf_path, pdb_path) in enumerate(psf_pdb_pairs):
            psf_abs = pathlib.Path(psf_path).resolve()
            pdb_abs = pathlib.Path(pdb_path).resolve()
            
            if not psf_abs.exists():
                raise FileNotFoundError(f"PSF file not found for pair {i + 1}: {psf_abs}")
            if not pdb_abs.exists():
                raise FileNotFoundError(f"PDB file not found for pair {i + 1}: {pdb_abs}")
            
            validated_pairs.append((str(psf_abs), str(pdb_abs)))
        
        output_dir = pathlib.Path(self._output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set output paths
        if output_psf_path is None:
            output_psf_path = output_dir / "complex.psf"
        else:
            output_psf_path = pathlib.Path(output_psf_path).resolve()
        
        if output_pdb_path is None:
            output_pdb_path = output_dir / "complex.pdb"
        else:
            output_pdb_path = pathlib.Path(output_pdb_path).resolve()
        
        print(f"Merging {len(validated_pairs)} PSF/PDB pairs...")
        for i, (psf_path, pdb_path) in enumerate(validated_pairs):
            print(f"  Pair {i + 1}: {pathlib.Path(psf_path).name}, {pathlib.Path(pdb_path).name}")
        
        # Create temp directory for intermediate files
        temp_dir = output_dir / "md_prepper" / "temp_files"
        temp_dir.mkdir(parents=True, exist_ok=True)

        # Generate TCL script
        tcl_script_path = temp_dir / "merge_psf_pdb.tcl"
        self._generate_merge_tcl(
            tcl_script_path=str(tcl_script_path),
            psf_pdb_pairs=validated_pairs,
            output_psf_path=str(output_psf_path),
            output_pdb_path=str(output_pdb_path)
        )
        
        # Execute VMD with merge TCL script
        print(f"Executing VMD merge...")
        print(f"  VMD Executable: {self._vmd_executable}")
        print(f"  TCL Script: {tcl_script_path}")
        
        stdout = third_party_tools.run_vmd(
            vmd_executable=self._vmd_executable,
            tcl_script=tcl_script_path,
            cwd=output_dir
        )
        print("VMD merge execution successful!")
        if stdout:
            print(f"STDOUT:\n{stdout}")
        
        print(f"  Merged PSF: {output_psf_path}")
        print(f"  Merged PDB: {output_pdb_path}")
        
        return (str(output_psf_path), str(output_pdb_path))
    
    def _generate_solvate_tcl(
        self,
        tcl_script_path: str,
        input_psf_path: str,
        input_pdb_path: str,
        output_prefix: str,
        padding: float = 15.0,
        use_cubic_box: bool = True
    ):
        """
        Generate TCL script for solvating a system using VMD solvate package.
        
        Args:
            tcl_script_path: Path to output TCL script
            input_psf_path: Path to input PSF file
            input_pdb_path: Path to input PDB file
            output_prefix: Prefix for output files (without extension)
            padding: Padding distance in Angstroms to add around the molecule
            use_cubic_box: If True, adjust box to be cubic (equal sides)
        """
        lines = []
        
        # Header
        lines.append("# Solvation Script")
        lines.append("# Generated by MDPrepper")
        lines.append("")
        
        # Load required packages
        lines.append("# Load required packages")
        lines.append("package require solvate")
        lines.append("")
        
        # Convert paths to POSIX format for TCL
        input_psf_posix = pathlib.Path(input_psf_path).as_posix()
        input_pdb_posix = pathlib.Path(input_pdb_path).as_posix()
        output_prefix_posix = pathlib.Path(output_prefix).as_posix()
        
        if use_cubic_box:
            # Calculate cubic box dimensions
            lines.append("# Load molecule to calculate bounding box")
            lines.append(f'set mol [mol new "{input_psf_posix}" type psf waitfor all]')
            lines.append(f'mol addfile "{input_pdb_posix}" type pdb molid $mol waitfor all')
            lines.append("")
            
            lines.append("# Calculate molecule bounding box")
            lines.append("set sel [atomselect top all]")
            lines.append("set minmax [measure minmax $sel]")
            lines.append("set min [lindex $minmax 0]")
            lines.append("set max [lindex $minmax 1]")
            lines.append("")
            
            lines.append(f"# Add padding of {padding} Angstroms")
            lines.append(f"set padding {padding}")
            lines.append("")
            
            lines.append("# Calculate box dimensions with padding")
            lines.append("set minX [expr [lindex $min 0] - $padding]")
            lines.append("set minY [expr [lindex $min 1] - $padding]")
            lines.append("set minZ [expr [lindex $min 2] - $padding]")
            lines.append("set maxX [expr [lindex $max 0] + $padding]")
            lines.append("set maxY [expr [lindex $max 1] + $padding]")
            lines.append("set maxZ [expr [lindex $max 2] + $padding]")
            lines.append("")
            
            lines.append("# Calculate box sizes")
            lines.append("set sizeX [expr $maxX - $minX]")
            lines.append("set sizeY [expr $maxY - $minY]")
            lines.append("set sizeZ [expr $maxZ - $minZ]")
            lines.append("")
            
            lines.append("# Find maximum dimension for cubic box")
            lines.append("set maxSize $sizeX")
            lines.append("if {$sizeY > $maxSize} {set maxSize $sizeY}")
            lines.append("if {$sizeZ > $maxSize} {set maxSize $sizeZ}")
            lines.append("")
            
            lines.append("# Calculate molecule center")
            lines.append("set centerX [expr ([lindex $min 0] + [lindex $max 0]) / 2.0]")
            lines.append("set centerY [expr ([lindex $min 1] + [lindex $max 1]) / 2.0]")
            lines.append("set centerZ [expr ([lindex $min 2] + [lindex $max 2]) / 2.0]")
            lines.append("")
            
            lines.append("# Calculate cubic box boundaries (centered on molecule)")
            lines.append("set halfSize [expr $maxSize / 2.0]")
            lines.append("set cubeMinX [expr $centerX - $halfSize]")
            lines.append("set cubeMinY [expr $centerY - $halfSize]")
            lines.append("set cubeMinZ [expr $centerZ - $halfSize]")
            lines.append("set cubeMaxX [expr $centerX + $halfSize]")
            lines.append("set cubeMaxY [expr $centerY + $halfSize]")
            lines.append("set cubeMaxZ [expr $centerZ + $halfSize]")
            lines.append("")
            
            lines.append("# Print box information")
            lines.append('puts "Molecule bounding box: $min to $max"')
            lines.append('puts "Cubic box size: $maxSize Angstroms"')
            lines.append('puts "Cubic box: ($cubeMinX, $cubeMinY, $cubeMinZ) to ($cubeMaxX, $cubeMaxY, $cubeMaxZ)"')
            lines.append("")
            
            lines.append("# Delete the loaded molecule (solvate will reload it)")
            lines.append("mol delete top")
            lines.append("")
            
            lines.append("# Solvate with cubic water box")
            lines.append("set boxMin [list $cubeMinX $cubeMinY $cubeMinZ]")
            lines.append("set boxMax [list $cubeMaxX $cubeMaxY $cubeMaxZ]")
            lines.append("set minmaxBox [list $boxMin $boxMax]")
            lines.append(f'solvate "{input_psf_posix}" "{input_pdb_posix}" -minmax $minmaxBox -o "{output_prefix_posix}"')
        else:
            # Simple solvation with padding on each side
            lines.append("# Solvate with rectangular water box")
            lines.append(f'solvate "{input_psf_posix}" "{input_pdb_posix}" -t {padding} -o "{output_prefix_posix}"')
        
        lines.append("")
        lines.append("# Done")
        lines.append("exit")
        
        # Write TCL script
        with open(tcl_script_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        
        print(f"Solvate TCL script generated: {tcl_script_path}")
    
    def solvate_system(
        self,
        input_psf_path: str,
        input_pdb_path: str,
        output_psf_path: typing.Optional[str] = None,
        output_pdb_path: typing.Optional[str] = None,
        padding: float = 15.0,
        use_cubic_box: bool = True
    ) -> typing.Tuple[str, str]:
        """
        Add a water box around the system using VMD solvate package.
        
        This method adds TIP3P water molecules around the solute. The water box
        extends 'padding' Angstroms beyond the molecule boundaries. If use_cubic_box
        is True, the box is adjusted to be cubic (equal sides) centered on the molecule.
        
        Args:
            input_psf_path: Path to input PSF file
            input_pdb_path: Path to input PDB file
            output_psf_path: Optional output PSF path. If None, uses "solvated.psf" in output_dir.
            output_pdb_path: Optional output PDB path. If None, uses "solvated.pdb" in output_dir.
            padding: Padding distance in Angstroms to add around the molecule (default: 15.0)
            use_cubic_box: If True, adjust box to be cubic with equal sides (default: True)
            
        Returns:
            Tuple of (solvated_psf_path, solvated_pdb_path) for the generated files
            
        Note:
            The VMD solvate plugin uses TIP3P water model by default.
            For cubic boxes, the box size is determined by the largest dimension
            of the molecule plus 2*padding.
        """
        # Convert paths to absolute pathlib.Path objects and verify existence
        input_psf = pathlib.Path(input_psf_path).resolve()
        input_pdb = pathlib.Path(input_pdb_path).resolve()
        
        if not input_psf.exists():
            raise FileNotFoundError(f"Input PSF file not found: {input_psf}")
        if not input_pdb.exists():
            raise FileNotFoundError(f"Input PDB file not found: {input_pdb}")
        
        output_dir = pathlib.Path(self._output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create temp directory for intermediate files
        temp_dir = output_dir / "md_prepper" / "temp_files"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Set output paths
        if output_psf_path is None:
            output_psf_path = output_dir / "solvated.psf"
        else:
            output_psf_path = pathlib.Path(output_psf_path).resolve()
        
        if output_pdb_path is None:
            output_pdb_path = output_dir / "solvated.pdb"
        else:
            output_pdb_path = pathlib.Path(output_pdb_path).resolve()
        
        # Determine output prefix (solvate command adds .psf and .pdb extensions)
        output_prefix = str(output_psf_path.with_suffix(""))
        
        box_type = "cubic" if use_cubic_box else "rectangular"
        print(f"Solvating system with {box_type} water box...")
        print(f"  Input PSF: {input_psf}")
        print(f"  Input PDB: {input_pdb}")
        print(f"  Padding: {padding} Angstroms")
        
        # Generate TCL script
        tcl_script_path = temp_dir / "solvate_system.tcl"
        self._generate_solvate_tcl(
            tcl_script_path=str(tcl_script_path),
            input_psf_path=str(input_psf),
            input_pdb_path=str(input_pdb),
            output_prefix=output_prefix,
            padding=padding,
            use_cubic_box=use_cubic_box
        )
        
        # Execute VMD with solvate TCL script
        print(f"Executing VMD solvate...")
        print(f"  VMD Executable: {self._vmd_executable}")
        print(f"  TCL Script: {tcl_script_path}")
        
        stdout = third_party_tools.run_vmd(
            vmd_executable=self._vmd_executable,
            tcl_script=tcl_script_path,
            cwd=output_dir
        )
        print("VMD solvate execution successful!")
        if stdout:
            print(f"STDOUT:\n{stdout}")
        
        print(f"  Solvated PSF: {output_psf_path}")
        print(f"  Solvated PDB: {output_pdb_path}")
        
        return (str(output_psf_path), str(output_pdb_path))
    
    def _generate_ionize_tcl(
        self,
        tcl_script_path: str,
        input_psf_path: str,
        input_pdb_path: str,
        output_prefix: str,
        salt_type: str = "NaCl",
        salt_concentration: float = 0.10,
        neutralize: bool = True
    ):
        """
        Generate TCL script for ionizing a system using VMD autoionize package.
        
        Args:
            tcl_script_path: Path to output TCL script
            input_psf_path: Path to input PSF file
            input_pdb_path: Path to input PDB file
            output_prefix: Prefix for output files (without extension)
            salt_type: Type of salt to add ("NaCl", "KCl", "CaCl2", "MgCl2")
            salt_concentration: Salt concentration in mol/L (M)
            neutralize: If True, neutralize the system charge first
        """
        # Define ion types for different salts (VMD naming convention)
        salt_ions = {
            "NaCl": {"cation": "SOD", "anion": "CLA", "cation_charge": 1},
            "KCl": {"cation": "POT", "anion": "CLA", "cation_charge": 1},
            "CaCl2": {"cation": "CAL", "anion": "CLA", "cation_charge": 2},
            "MgCl2": {"cation": "MG", "anion": "CLA", "cation_charge": 2}
        }
        
        if salt_type not in salt_ions:
            raise ValueError(f"Unsupported salt type: {salt_type}. "
                           f"Supported types: {list(salt_ions.keys())}")
        
        ion_info = salt_ions[salt_type]
        cation = ion_info["cation"]
        anion = ion_info["anion"]
        
        lines = []
        
        # Header
        lines.append("# Ionization Script")
        lines.append("# Generated by MDPrepper")
        lines.append(f"# Salt type: {salt_type}, Concentration: {salt_concentration} M")
        lines.append("")
        
        # Load required packages
        lines.append("# Load required packages")
        lines.append("package require autoionize")
        lines.append("")
        
        # Convert paths to POSIX format for TCL
        input_psf_posix = pathlib.Path(input_psf_path).as_posix()
        input_pdb_posix = pathlib.Path(input_pdb_path).as_posix()
        output_prefix_posix = pathlib.Path(output_prefix).as_posix()
        
        # Build autoionize command
        lines.append("# Add ions to neutralize and set salt concentration")
        
        if neutralize:
            # autoionize with -neutralize flag and salt concentration
            lines.append(f'autoionize -psf "{input_psf_posix}" -pdb "{input_pdb_posix}" '
                        f'-sc {salt_concentration} -cation {cation} -anion {anion} '
                        f'-o "{output_prefix_posix}"')
        else:
            # autoionize without neutralization (just add salt)
            lines.append(f'autoionize -psf "{input_psf_posix}" -pdb "{input_pdb_posix}" '
                        f'-sc {salt_concentration} -cation {cation} -anion {anion} '
                        f'-nna 0 -ncl 0 -o "{output_prefix_posix}"')
        
        lines.append("")
        lines.append("# Done")
        lines.append("exit")
        
        # Write TCL script
        with open(tcl_script_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        
        print(f"Ionize TCL script generated: {tcl_script_path}")
    
    def ionize_system(
        self,
        input_psf_path: str,
        input_pdb_path: str,
        output_psf_path: typing.Optional[str] = None,
        output_pdb_path: typing.Optional[str] = None,
        salt_type: str = "NaCl",
        salt_concentration: float = 0.10,
        neutralize: bool = True
    ) -> typing.Tuple[str, str]:
        """
        Add ions to the system using VMD autoionize package.
        
        This method adds ions to neutralize the system charge and achieve a
        specified salt concentration. Supports common salt types used in MD.
        
        Args:
            input_psf_path: Path to input PSF file (typically solvated system)
            input_pdb_path: Path to input PDB file (typically solvated system)
            output_psf_path: Optional output PSF path. If None, uses "ionized.psf" in output_dir.
            output_pdb_path: Optional output PDB path. If None, uses "ionized.pdb" in output_dir.
            salt_type: Type of salt to add. Options:
                      - "NaCl" (default): Sodium chloride (Na+, Cl-)
                      - "KCl": Potassium chloride (K+, Cl-)
                      - "CaCl2": Calcium chloride (Ca2+, 2Cl-)
                      - "MgCl2": Magnesium chloride (Mg2+, 2Cl-)
            salt_concentration: Salt concentration in mol/L (M). Default: 0.10 M
            neutralize: If True, neutralize the system charge first. Default: True
            
        Returns:
            Tuple of (ionized_psf_path, ionized_pdb_path) for the generated files
            
        Note:
            The autoionize plugin will:
            1. Calculate system net charge
            2. Add counter-ions to neutralize the charge
            3. Add additional salt ions to reach the specified concentration
            
        Example:
            # Add 0.15 M NaCl and neutralize
            ionized_psf, ionized_pdb = prepper.ionize_system(
                input_psf_path="solvated.psf",
                input_pdb_path="solvated.pdb",
                salt_type="NaCl",
                salt_concentration=0.15
            )
            
            # Add 0.1 M KCl
            ionized_psf, ionized_pdb = prepper.ionize_system(
                input_psf_path="solvated.psf",
                input_pdb_path="solvated.pdb",
                salt_type="KCl",
                salt_concentration=0.1
            )
        """
        # Validate salt type
        valid_salt_types = ["NaCl", "KCl", "CaCl2", "MgCl2"]
        if salt_type not in valid_salt_types:
            raise ValueError(f"Invalid salt type: {salt_type}. "
                           f"Valid options: {valid_salt_types}")
        
        # Convert paths to absolute pathlib.Path objects and verify existence
        input_psf = pathlib.Path(input_psf_path).resolve()
        input_pdb = pathlib.Path(input_pdb_path).resolve()
        
        if not input_psf.exists():
            raise FileNotFoundError(f"Input PSF file not found: {input_psf}")
        if not input_pdb.exists():
            raise FileNotFoundError(f"Input PDB file not found: {input_pdb}")
        
        output_dir = pathlib.Path(self._output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create temp directory for intermediate files
        temp_dir = output_dir / "md_prepper" / "temp_files"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Set output paths
        if output_psf_path is None:
            output_psf_path = output_dir / "ionized.psf"
        else:
            output_psf_path = pathlib.Path(output_psf_path).resolve()
        
        if output_pdb_path is None:
            output_pdb_path = output_dir / "ionized.pdb"
        else:
            output_pdb_path = pathlib.Path(output_pdb_path).resolve()
        
        # Determine output prefix (autoionize command adds .psf and .pdb extensions)
        output_prefix = str(output_psf_path.with_suffix(""))
        
        print(f"Ionizing system with {salt_type}...")
        print(f"  Input PSF: {input_psf}")
        print(f"  Input PDB: {input_pdb}")
        print(f"  Salt type: {salt_type}")
        print(f"  Salt concentration: {salt_concentration} M")
        print(f"  Neutralize: {neutralize}")
        
        # Generate TCL script
        tcl_script_path = temp_dir / "ionize_system.tcl"
        self._generate_ionize_tcl(
            tcl_script_path=str(tcl_script_path),
            input_psf_path=str(input_psf),
            input_pdb_path=str(input_pdb),
            output_prefix=output_prefix,
            salt_type=salt_type,
            salt_concentration=salt_concentration,
            neutralize=neutralize
        )
        
        # Execute VMD with ionize TCL script
        print(f"Executing VMD autoionize...")
        print(f"  VMD Executable: {self._vmd_executable}")
        print(f"  TCL Script: {tcl_script_path}")
        
        stdout = third_party_tools.run_vmd(
            vmd_executable=self._vmd_executable,
            tcl_script=tcl_script_path,
            cwd=output_dir
        )
        print("VMD autoionize execution successful!")
        if stdout:
            print(f"STDOUT:\n{stdout}")
        
        print(f"  Ionized PSF: {output_psf_path}")
        print(f"  Ionized PDB: {output_pdb_path}")
        
        return (str(output_psf_path), str(output_pdb_path))
    
    def _create_template_from_smiles(self, smiles: str) -> typing.Tuple[typing.Any, int]:
        """
        Create RDKit template molecule from SMILES string.
        
        Args:
            smiles: SMILES string
            
        Returns:
            Tuple of (template molecule with hydrogens, formal charge)
            
        Raises:
            ValueError: If SMILES parsing fails
        """
        template = rdkit.Chem.MolFromSmiles(smiles)
        if template is None:
            raise ValueError(f"Failed to parse SMILES: {smiles}")
        
        template = rdkit.Chem.AddHs(template)
        formal_charge = rdkit.Chem.GetFormalCharge(template)
        
        print(f"  Template molecule: {template.GetNumAtoms()} atoms (with hydrogens)")
        print(f"  Formal charge from SMILES: {formal_charge}")
        
        return template, formal_charge
    
    def _load_pdb_with_mdanalysis(self, pdb_path: pathlib.Path, expected_atom_count: int):
        """
        Load PDB file using MDAnalysis to preserve atom names.
        
        Args:
            pdb_path: Path to PDB file
            expected_atom_count: Expected number of atoms from SMILES template
            
        Returns:
            MDAnalysis AtomGroup
            
        Raises:
            ValueError: If atom count doesn't match expected count
        """
        u = MDAnalysis.Universe(str(pdb_path))
        atoms = u.select_atoms("all")
        
        print(f"  PDB molecule: {len(atoms)} atoms")
        
        if expected_atom_count != len(atoms):
            print(
                f"  Warning: Atom count mismatch: SMILES has {expected_atom_count} atoms, "
                f"PDB has {len(atoms)} atoms.\n"
                f"  Continuing, but this may require fixing the structure (e.g. re-adding hydrogens)."
            )
        
        return atoms
    
    def _clean_pdb_for_rdkit(self, pdb_path: pathlib.Path) -> pathlib.Path:
        """
        Clean PDB file for RDKit parsing (remove PDBQT-specific lines).
        
        Args:
            pdb_path: Path to original PDB file
            
        Returns:
            Path to cleaned temporary PDB file
        """
        temp_pdb = pdb_path.parent / f"temp_rdkit_{pdb_path.name}"
        
        with open(pdb_path, 'r') as f:
            lines = f.readlines()
        
        cleaned_lines = []
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                # Clear element column (columns 77-78) to let RDKit infer from atom name
                if len(line) > 76:
                    line = line[:76] + '  ' + (line[78:] if len(line) > 78 else '')
                cleaned_lines.append(line)
            elif not line.startswith(('ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF')):
                cleaned_lines.append(line)
        
        with open(temp_pdb, 'w') as f:
            f.writelines(cleaned_lines)
        
        return temp_pdb
    
    def _assign_bond_orders_from_template(self, template, pdb_path: pathlib.Path):
        """
        Assign bond orders from SMILES template to PDB structure.
        
        Args:
            template: RDKit molecule from SMILES with correct bond orders
            pdb_path: Path to cleaned PDB file
            
        Returns:
            RDKit molecule with assigned bond orders
            
        Raises:
            ValueError: If PDB parsing fails
            RuntimeError: If bond order assignment fails
        """
        pdb_mol = rdkit.Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
        
        if pdb_mol is None:
            raise ValueError(f"Failed to read PDB file with RDKit: {pdb_path}")
        
        print("  Assigning bond orders from SMILES template...")
        
        try:
            mol_with_bonds = rdkit.Chem.AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)
        except Exception as e:
            print(f"  Error during bond order assignment: {e}")
            print("  Attempting alternative method with sanitization...")
            try:
                rdkit.Chem.SanitizeMol(pdb_mol)
                mol_with_bonds = rdkit.Chem.AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)
            except Exception as e2:
                print(f"  Standard assignment failed: {e2}")
                print("  FALLBACK: Forcing coordinates from PDB to SMILES template (Trusting SMILES Topology)...")
                
                # Clone template to avoid modifying original
                mol_with_bonds = rdkit.Chem.Mol(template)
                mol_with_bonds.RemoveAllConformers()
                
                pdb_conf = pdb_mol.GetConformer()
                new_conf = rdkit.Chem.Conformer(mol_with_bonds.GetNumAtoms())
                
                # Helper to transfer PDB info
                def transfer_pdb_info(src_atom, dst_atom):
                    info = src_atom.GetPDBResidueInfo()
                    if info:
                        dst_atom.SetProp("AtomName", info.GetName().strip())
                        dst_atom.SetProp("ResName", info.GetResidueName().strip())
                        dst_atom.SetProp("ResId", str(info.GetResidueNumber()))
                        dst_atom.SetProp("ChainId", info.GetChainId().strip())

                # Check if we can do 1-to-1 mapping
                if template.GetNumAtoms() == pdb_mol.GetNumAtoms():
                    # 1-to-1 Mapping
                    for i in range(mol_with_bonds.GetNumAtoms()):
                        pos = pdb_conf.GetAtomPosition(i)
                        new_conf.SetAtomPosition(i, pos)
                        transfer_pdb_info(pdb_mol.GetAtomWithIdx(i), mol_with_bonds.GetAtomWithIdx(i))
                    mol_with_bonds.AddConformer(new_conf, assignId=True)
                    print("  Successfully forced PDB coordinates onto SMILES structure (1-to-1).")
                    
                else:
                    # MCS-based Smart Alignment (handles atom order differences)
                    print(f"  Atom count mismatch or alignment check needed. Using MCS to map PDB structure to SMILES...")
                    from rdkit.Chem import rdFMCS
                    
                    # 1. Find Maximum Common Substructure (MCS) strict on elements
                    mcs_result = rdFMCS.FindMCS(
                        [mol_with_bonds, pdb_mol],
                        atomCompare=rdFMCS.AtomCompare.CompareElements,
                        bondCompare=rdFMCS.BondCompare.CompareAny,
                        matchValences=False,
                        ringMatchesRingOnly=False,
                        completeRingsOnly=False,
                        timeout=5
                    )
                    
                    if not mcs_result.smartsString:
                         raise RuntimeError("Could not find any common substructure between SMILES and PDB!")

                    common_mol = rdkit.Chem.MolFromSmarts(mcs_result.smartsString)
                    print(f"  MCS found (SMARTS len: {len(mcs_result.smartsString)}), mapping atoms...")
                    
                    match_smi = mol_with_bonds.GetSubstructMatch(common_mol)
                    match_pdb = pdb_mol.GetSubstructMatch(common_mol)
                    
                    if len(match_smi) != len(match_pdb):
                        # Should rarely happen if derived from MCS
                        print("  Warning: MCS match length mismatch.")
                    
                    # 2. Build Coordinate Map: SMILES_Index -> PDB_Coordinate
                    coord_map = {}
                    fixed_indices = [] # Keep track for naming/info transfer
                    
                    for i in range(len(match_smi)):
                        smi_idx = match_smi[i]
                        pdb_idx = match_pdb[i]
                        
                        smi_atom = mol_with_bonds.GetAtomWithIdx(smi_idx)
                        pdb_atom = pdb_mol.GetAtomWithIdx(pdb_idx)
                        
                        # Double check element match (MCS should guarantee this but safe is good)
                        if smi_atom.GetSymbol() != pdb_atom.GetSymbol():
                             # E.g. MCS general match on aromatic might match c to n?
                             continue

                        pos = pdb_conf.GetAtomPosition(pdb_idx)
                        coord_map[smi_idx] = pos
                        fixed_indices.append(smi_idx)
                        
                        # Transfer Props
                        transfer_pdb_info(pdb_atom, smi_atom)
                    
                    print(f"  Mapped {len(coord_map)} atoms via MCS.")
                    
                    if len(coord_map) < mol_with_bonds.GetNumAtoms() * 0.5:
                        print("  Warning: Less than 50% of atoms mapped! Structure might be very wrong.")

                    # Add empty conformer to hold result
                    mol_with_bonds.AddConformer(new_conf, assignId=True)

                    # Now handle unmapped atoms via Constrained Embedding (below)
                    print("  Regenerating coordinates for unmapped atoms (including Hydrogens)...")
                    
                    try:
                        
                        # 2. Embed molecule using constraints
                        # This generates a 3D structure that respects the fixed coordinates in coord_map
                        # forceTol=0.01 makes it strictly respect the map
                        res = rdkit.Chem.AllChem.EmbedMolecule(
                            mol_with_bonds,
                            coordMap=coord_map,
                            useExpTorsionAnglePrefs=True,
                            useBasicKnowledge=True,
                            enforceChirality=True
                        )
                        
                        if res == -1:
                             print("  Warning: Constrained Embedding failed. Retrying with random coordinates...")
                             res = rdkit.Chem.AllChem.EmbedMolecule(mol_with_bonds, coordMap=coord_map, useRandomCoords=True)
                        
                        if res != -1:
                             # 3. STRICTLY restore heavy atom coordinates from PDB
                             # Embedding might slightly drift or relax heavy atoms, but valid docking requires exact PDB coords
                             conf = mol_with_bonds.GetConformer()
                             for smi_idx, pos in coord_map.items():
                                 conf.SetAtomPosition(smi_idx, pos)
                             
                             # 4. Minimize ONLY Hydrogens to fix bond lengths/clashes
                             # Fix heavy atoms
                             ff = rdkit.Chem.AllChem.MMFFGetMoleculeForceField(mol_with_bonds, confId=0)
                             if ff is None: # Fallback to UFF
                                 ff = rdkit.Chem.AllChem.UFFGetMoleculeForceField(mol_with_bonds)
                             
                             if ff:
                                 for smi_idx in coord_map.keys():
                                     ff.AddFixedPoint(smi_idx)
                                 ff.Initialize()
                                 ff.Minimize(maxIts=500)
                                 print("  Regenerated H coordinates and refined geometry.")
                             else:
                                 print("  Warning: Could not set up ForceField for hydrogen refinement.")
                        else:
                             print("  Error: Failed to embed molecule dimensions.")
                             
                    except Exception as e_gen:
                         print(f"  Warning: Hydrogen generation failed: {e_gen}")
                         traceback.print_exc()

        print("  Successfully assigned bond orders (and coordinates)")
        return mol_with_bonds
    
    def _transfer_formal_charges(self, template, mol):
        """
        Transfer formal charges from template to molecule.
        
        Args:
            template: RDKit molecule (source of charges)
            mol: RDKit molecule (destination)
        """
        print("  Transferring formal charges from template...")
        try:
            # Find match between mol (pdb derived) and template (smiles)
            # Both should now have consistent bond orders/connectivity
            match = mol.GetSubstructMatch(template)
            if not match or len(match) != mol.GetNumAtoms():
                # Try fallback: if sizes match, assume 1-to-1 if AssignBondOrders worked? 
                # No, AssignBondOrders matched based on structure. 
                # If GetSubstructMatch fails here, something is weird.
                print(f"  Warning: Could not find perfect substructure match for charge transfer (Match len: {len(match) if match else 0})")
                return

            for template_idx, mol_idx in enumerate(match):
                template_atom = template.GetAtomWithIdx(template_idx)
                mol_atom = mol.GetAtomWithIdx(mol_idx)
                if template_atom.GetFormalCharge() != 0:
                    mol_atom.SetFormalCharge(template_atom.GetFormalCharge())
                
            print(f"  Transferred formal charges based on substructure match")
            
        except Exception as e:
            print(f"  Warning: Failed to transfer formal charges: {e}")
    
    def _extract_bond_information(self, mol):
        """
        Extract bond information from RDKit molecule.
        
        Args:
            mol: RDKit molecule with bond information
            
        Returns:
            Tuple of (bonds list, bond_dict) where:
            - bonds: List of (begin_idx, end_idx, bond_order_str)
            - bond_dict: Dict mapping atom index to list of (neighbor_idx, bond_type)
        """
        bonds = []
        bond_dict = {}  # Map atom index to list of (neighbor_idx, bond_type)
        
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType()
            
            # Convert bond type to MOL2 format
            if bond_type == rdkit.Chem.BondType.SINGLE:
                bond_order = '1'
            elif bond_type == rdkit.Chem.BondType.DOUBLE:
                bond_order = '2'
            elif bond_type == rdkit.Chem.BondType.TRIPLE:
                bond_order = '3'
            elif bond_type == rdkit.Chem.BondType.AROMATIC:
                bond_order = 'ar'
            else:
                bond_order = '1'  # Default to single
            
            bonds.append((begin_idx, end_idx, bond_order))
            
            # Build bond dictionary for atom type inference
            if begin_idx not in bond_dict:
                bond_dict[begin_idx] = []
            if end_idx not in bond_dict:
                bond_dict[end_idx] = []
            bond_dict[begin_idx].append((end_idx, bond_type))
            bond_dict[end_idx].append((begin_idx, bond_type))
        
        print(f"  Extracted {len(bonds)} bonds from molecule")
        return bonds, bond_dict
    
    def _determine_mol2_atom_types(self, mol, bond_dict: dict) -> typing.List[str]:
        """
        Determine MOL2 atom types from RDKit molecule.
        
        Args:
            mol: RDKit molecule
            bond_dict: Bond dictionary from _extract_bond_information
            
        Returns:
            List of MOL2 atom type strings
        """
        print("  Determining MOL2 atom types from RDKit...")
        mol2_atom_types = []
        
        for i in range(mol.GetNumAtoms()):
            rdkit_atom = mol.GetAtomWithIdx(i)
            element = rdkit_atom.GetSymbol()
            is_aromatic = rdkit_atom.GetIsAromatic()
            hybridization = rdkit_atom.GetHybridization()
            
            # Determine MOL2 atom type based on element, hybridization, and aromaticity
            if element == 'C':
                mol2_type = self._determine_carbon_type(is_aromatic, hybridization)
            elif element == 'N':
                mol2_type = self._determine_nitrogen_type(
                    i, is_aromatic, hybridization, rdkit_atom, mol, bond_dict
                )
            elif element == 'O':
                mol2_type = self._determine_oxygen_type(i, hybridization, bond_dict)
            elif element == 'S':
                mol2_type = self._determine_sulfur_type(i, hybridization, bond_dict)
            elif element == 'P':
                mol2_type = 'P.3'
            elif element == 'H':
                mol2_type = 'H'
            elif element in ['F', 'Cl', 'Br', 'I']:
                mol2_type = element
            else:
                mol2_type = element
                print(f"  Warning: Unknown element '{element}' at atom {i}, using element symbol as type")
            
            mol2_atom_types.append(mol2_type)
        
        print(f"  Assigned {len(mol2_atom_types)} atom types")
        return mol2_atom_types
    
    def _determine_carbon_type(self, is_aromatic: bool, hybridization) -> str:
        """Determine MOL2 type for carbon atom."""
        if is_aromatic:
            return 'C.ar'
        elif hybridization == rdkit.Chem.HybridizationType.SP3:
            return 'C.3'
        elif hybridization == rdkit.Chem.HybridizationType.SP2:
            return 'C.2'
        elif hybridization == rdkit.Chem.HybridizationType.SP:
            return 'C.1'
        else:
            return 'C.3'  # Default
    
    def _determine_nitrogen_type(
        self, 
        atom_idx: int,
        is_aromatic: bool, 
        hybridization, 
        rdkit_atom,
        mol,
        bond_dict: dict
    ) -> str:
        """Determine MOL2 type for nitrogen atom."""
        if is_aromatic:
            return 'N.ar'
        elif hybridization == rdkit.Chem.HybridizationType.SP3:
            # Check if it's ammonium (N+)
            formal_charge = rdkit_atom.GetFormalCharge()
            return 'N.4' if formal_charge > 0 else 'N.3'
        elif hybridization == rdkit.Chem.HybridizationType.SP2:
            # Check if it's amide nitrogen
            is_amide = self._is_amide_nitrogen(atom_idx, mol, bond_dict)
            return 'N.am' if is_amide else 'N.2'
        elif hybridization == rdkit.Chem.HybridizationType.SP:
            return 'N.1'
        else:
            return 'N.3'  # Default
    
    def _is_amide_nitrogen(self, atom_idx: int, mol, bond_dict: dict) -> bool:
        """Check if nitrogen is part of an amide group."""
        if atom_idx not in bond_dict:
            return False
        
        for neighbor_idx, bond_type in bond_dict[atom_idx]:
            neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
            if neighbor_atom.GetSymbol() == 'C':
                # Check if carbon is bonded to O with double bond (C=O)
                for c_neighbor_idx, c_bond_type in bond_dict.get(neighbor_idx, []):
                    c_neighbor = mol.GetAtomWithIdx(c_neighbor_idx)
                    if (c_neighbor.GetSymbol() == 'O' and 
                        c_bond_type == rdkit.Chem.BondType.DOUBLE):
                        return True
        return False
    
    def _determine_oxygen_type(self, atom_idx: int, hybridization, bond_dict: dict) -> str:
        """Determine MOL2 type for oxygen atom."""
        if hybridization == rdkit.Chem.HybridizationType.SP3:
            return 'O.3'
        elif hybridization == rdkit.Chem.HybridizationType.SP2:
            return 'O.2'
        else:
            # Check valence to distinguish O.2 vs O.3
            num_bonds = len(bond_dict.get(atom_idx, []))
            return 'O.2' if num_bonds == 1 else 'O.3'
    
    def _determine_sulfur_type(self, atom_idx: int, hybridization, bond_dict: dict) -> str:
        """Determine MOL2 type for sulfur atom."""
        if hybridization == rdkit.Chem.HybridizationType.SP3:
            # Check oxidation state
            num_bonds = len(bond_dict.get(atom_idx, []))
            return 'S.O2' if num_bonds >= 4 else 'S.3'  # Sulfone/sulfate
        elif hybridization == rdkit.Chem.HybridizationType.SP2:
            return 'S.2'
        else:
            return 'S.3'  # Default
    
    def _compute_gasteiger_charges(self, mol, num_atoms: int) -> typing.List[float]:
        """
        Compute Gasteiger partial charges for molecule.
        
        Args:
            mol: RDKit molecule
            num_atoms: Number of atoms (for fallback zero charges)
            
        Returns:
            List of partial charges
        """
        print("  Computing Gasteiger charges...")
        try:
            rdkit.Chem.AllChem.ComputeGasteigerCharges(mol)
            charges = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge')
                      for i in range(mol.GetNumAtoms())]
            total_charge = sum(charges)
            print(f"  Total molecular charge (Gasteiger): {total_charge:.4f}")
            return charges
        except Exception as e:
            print(f"  Warning: Failed to compute Gasteiger charges: {e}")
            print("  Using zero charges...")
            return [0.0] * num_atoms
    
    
    def _write_mol2_file(
        self,
        output_path: pathlib.Path,
        pdb_path: pathlib.Path,
        atoms,
        bonds: list,
        mol2_atom_types: list,
        charges: list,
        mol=None
    ):
        """
        Write MOL2 file with preserved atom names.
        
        Args:
            output_path: Path to output MOL2 file
            pdb_path: Original PDB path (for molecule name)
            atoms: MDAnalysis AtomGroup (legacy source of names)
            bonds: List of bond tuples (begin_idx, end_idx, bond_order)
            mol2_atom_types: List of MOL2 atom types
            charges: List of partial charges
            mol: RDKit molecule (optional, preferred source of coords/names)
        """
        print(f"  Writing MOL2 file: {output_path}")
        
        # Determine number of atoms
        num_atoms = 0
        if mol:
            num_atoms = mol.GetNumAtoms()
        elif atoms:
            num_atoms = len(atoms)
            
        with open(output_path, 'w') as f:
            # Header - MOLECULE record
            f.write("@<TRIPOS>MOLECULE\n")
            molecule_name = pdb_path.stem
            f.write(f"{molecule_name}\n")
            f.write(f"{num_atoms} {len(bonds)} 1 0 0\n")
            f.write("SMALL\n")
            f.write("USER_CHARGES\n")
            f.write("\n")
            
            # ATOM record
            f.write("@<TRIPOS>ATOM\n")
            
            # Track H count for naming if needed
            h_count = 0
            
            for i in range(num_atoms):
                atom_id = i + 1
                
                # Default values
                atom_name = f"A{atom_id}"
                res_name = "UNL"
                res_id = 1
                x, y, z = 0.0, 0.0, 0.0
                symbol = "C"
                
                # Try to get data from RDKit Mol first
                if mol:
                    m_atom = mol.GetAtomWithIdx(i)
                    symbol = m_atom.GetSymbol()
                    pos = mol.GetConformer().GetAtomPosition(i)
                    x, y, z = pos.x, pos.y, pos.z
                    
                    # Try PDBResidueInfo (standard)
                    info = m_atom.GetPDBResidueInfo()
                    if info:
                        atom_name = info.GetName().strip()
                        res_name = info.GetResidueName().strip()
                        res_id = info.GetResidueNumber()
                    
                    # Try Properties (FALLBACK set these)
                    if m_atom.HasProp("AtomName"):
                        atom_name = m_atom.GetProp("AtomName")
                    if m_atom.HasProp("ResName"):
                        res_name = m_atom.GetProp("ResName")
                    if m_atom.HasProp("ResId"):
                        try:
                            res_id = int(m_atom.GetProp("ResId"))
                        except:
                            pass
                            
                # Fallback to MDAnalysis atoms if alignment seems valid and no mol/info
                # Use strict check: if mol provided, only use MDA if counts match exactly
                use_mda = False
                if atoms and i < len(atoms):
                    if not mol:
                         use_mda = True
                    elif mol and len(atoms) == mol.GetNumAtoms():
                         # Only use MDA if names are missing in RDKit
                         if atom_name.startswith("A") and atom_name[1:].isdigit():
                             use_mda = True
                
                if use_mda:
                    atom = atoms[i]
                    atom_name = atom.name
                    if not mol:
                        x, y, z = atom.position
                    res_name = atom.resname if hasattr(atom, 'resname') else 'UNL'
                    res_id = atom.resid if hasattr(atom, 'resid') else 1
                
                # Final cleanup for generated Hydrogens
                if atom_name.startswith("A") and atom_name[1:].isdigit():
                     # Likely default name
                     if symbol == 'H':
                         h_count += 1
                         atom_name = f"H{h_count}"
                     else:
                         atom_name = f"{symbol}{atom_id}"

                atom_type = mol2_atom_types[i]
                charge = charges[i] if i < len(charges) else 0.0
                
                f.write(f"{atom_id:7d} {atom_name:<8s} {x:10.4f} {y:10.4f} {z:10.4f} "
                       f"{atom_type:<8s} {res_id:>4d} {res_name:<8s} {charge:>10.6f}\n")
            
            # BOND record
            f.write("@<TRIPOS>BOND\n")
            for bond_id, (begin_idx, end_idx, bond_order) in enumerate(bonds, 1):
                atom1_id = begin_idx + 1
                atom2_id = end_idx + 1
                f.write(f"{bond_id:6d} {atom1_id:>5d} {atom2_id:>5d} {bond_order:<4s}\n")
            
            # SUBSTRUCTURE record
            f.write("@<TRIPOS>SUBSTRUCTURE\n")
            # Use data from first atom
            first_res_id = 1
            first_res_name = "UNL"
            
            # Try to extract from first written atom (simplification)
            # Ideally we wrote these, but let's just grab reasonable defaults or from atoms[0]
            if mol:
                 m_atom = mol.GetAtomWithIdx(0)
                 info = m_atom.GetPDBResidueInfo()
                 if info:
                     first_res_name = info.GetResidueName().strip()
                     first_res_id = info.GetResidueNumber()
                 elif m_atom.HasProp("ResName"):
                     first_res_name = m_atom.GetProp("ResName")
                     try: 
                         first_res_id = int(m_atom.GetProp("ResId"))
                     except: pass
            elif atoms and len(atoms) > 0:
                 atom = atoms[0]
                 first_res_name = atom.resname if hasattr(atom, 'resname') else 'UNL'
                 first_res_id = atom.resid if hasattr(atom, 'resid') else 1

            f.write(f"{first_res_id:>5d} {first_res_name:<8s} {1:>5d} ****               0 ****  **** \n")
    
    def generate_mol2(
        self,
        pdb_path: str,
        smiles: str,
        output_mol2_path: typing.Optional[str] = None
    ) -> str:
        """
        Generate mol2 file from PDB coordinates using MDAnalysis to preserve atom names.
        
        This method uses MDAnalysis to read the PDB file, which preserves the original
        atom names. It then uses RDKit to assign bond orders based on SMILES, and finally
        writes a MOL2 file with the original atom names preserved.
        
        Args:
            pdb_path: Path to input PDB file with 3D coordinates
            smiles: SMILES string (not a file path)
            output_mol2_path: Optional output path. If None, uses default naming.
            
        Returns:
            Path to the output mol2 file
            
        Raises:
            ValueError: If SMILES is invalid or structure doesn't match PDB
            RuntimeError: If mol2 file generation fails
        """
        # Prepare paths
        pdb_path = pathlib.Path(pdb_path).resolve()
        output_dir = pathlib.Path(self._output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if output_mol2_path is None:
            output_mol2_path = output_dir / "ligand.mol2"
        else:
            output_mol2_path = pathlib.Path(output_mol2_path).resolve()
        
        print(f"Generating mol2 file from PDB using MDAnalysis (preserving atom names)...")
        print(f"  SMILES: {smiles}")
        print(f"  PDB: {pdb_path}")
        
        # Step 1: Create template from SMILES
        template, formal_charge = self._create_template_from_smiles(smiles)
        
        # Step 2: Load PDB with MDAnalysis
        atoms = self._load_pdb_with_mdanalysis(pdb_path, template.GetNumAtoms())
        
        # Step 3: Clean PDB and assign bond orders
        temp_pdb_for_rdkit = self._clean_pdb_for_rdkit(pdb_path)
        
        try:
            mol_with_bonds = self._assign_bond_orders_from_template(template, temp_pdb_for_rdkit)
            
            # Transfer formal charges from template to ensure correct atom tying
            self._transfer_formal_charges(template, mol_with_bonds)
            
            # Step 4: Extract bond information
            bonds, bond_dict = self._extract_bond_information(mol_with_bonds)
            
            # Step 5: Determine atom types
            mol2_atom_types = self._determine_mol2_atom_types(mol_with_bonds, bond_dict)
            
            # Step 6: Compute charges
            charges = self._compute_gasteiger_charges(mol_with_bonds, mol_with_bonds.GetNumAtoms())
            
            # Step 7: Write MOL2 file
            self._write_mol2_file(
                output_mol2_path, pdb_path, atoms, bonds, mol2_atom_types, charges, mol=mol_with_bonds
            )
            
            print(f"  Successfully generated MOL2 file with preserved atom names")
            print(f"  Formal charge: {formal_charge}")
            
            return str(output_mol2_path)
        
        finally:
            # Clean up temporary file
            if temp_pdb_for_rdkit.exists():
                temp_pdb_for_rdkit.unlink()
    
    @staticmethod
    def fix_str_residue_name(
        input_str_path: str, 
        output_str_path: typing.Optional[str] = None,
        new_residue_name: str = "UNL"
    ) -> str:
        """
        Modify the residue name in RESI lines of a CHARMM str file.
        
        This function reads a CHARMM topology stream file (.str), finds RESI lines,
        and replaces the residue name (can be any name like "ligand") with the 
        specified name (default "UNL").
        
        Example RESI line format:
            RESI ligand        -1.000 ! param penalty=   0.000 ; charge penalty=   0.000
            RESI MOL           0.000 ! ...
        
        Will be replaced with:
            RESI UNL           -1.000 ! param penalty=   0.000 ; charge penalty=   0.000
            RESI UNL           0.000 ! ...
        
        Args:
            input_str_path: Path to input str file
            output_str_path: Path to output str file. If None, overwrites the input file
            new_residue_name: New residue name, default is "UNL"
            
        Returns:
            Path to the output file
            
        Raises:
            FileNotFoundError: If input file does not exist
            ValueError: If no RESI line is found in the file
        """
        # Convert to pathlib.Path object
        input_path = pathlib.Path(input_str_path).resolve()
        
        if not input_path.exists():
            raise FileNotFoundError(f"Input file does not exist: {input_path}")
        
        # If output path is not specified, overwrite the input file
        if output_str_path is None:
            output_path = input_path
        else:
            output_path = pathlib.Path(output_str_path).resolve()
        
        print(f"Processing str file: {input_path}")
        print(f"New residue name: {new_residue_name}")
        
        # Read file contents
        with open(input_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        modified_lines = []
        resi_found = False
        
        for line in lines:
            # Check if this is a RESI line (case insensitive)
            if line.strip().upper().startswith('RESI'):
                # Split the line to get its parts
                parts = line.split()
                if len(parts) >= 3:  # Must have at least RESI, residue name, and charge
                    # Modify residue name (index 1)
                    parts[1] = new_residue_name
                    # Rejoin the line, preserving the line ending (\n)
                    new_line = ' '.join(parts) + '\n'
                    modified_lines.append(new_line)
                    resi_found = True
                    print(f"  Found and modified: {line.strip()}")
                    print(f"  Changed to:         {new_line.strip()}")
                else:
                    # RESI line format is incorrect, keep as-is
                    modified_lines.append(line)
            else:
                modified_lines.append(line)
        
        if not resi_found:
            raise ValueError(f"No RESI line found in file {input_path}")
        
        # Write to output file
        with open(output_path, 'w', encoding='utf-8') as f:
            f.writelines(modified_lines)
        
        print(f"Modified file saved to: {output_path}")
        
        return str(output_path)


    def process_protein(self) -> str:
        """
        Process protein: fix flexible residues, add hydrogens, fix histidine names.
        
        Returns:
            Path to the final processed protein PDB file.
        """
        # Define output directory for fixed files
        md_prepper_dir = pathlib.Path(self._output_dir) / "md_prepper"
        fixed_files_dir = md_prepper_dir / "fixed_files"
        temp_files_dir = md_prepper_dir / "temp_files"
        
        fixed_files_dir.mkdir(parents=True, exist_ok=True)
        temp_files_dir.mkdir(parents=True, exist_ok=True)
        
        print("\nProcessing protein pipeline...")
        print(f"Output directory (final): {fixed_files_dir}")
        print(f"Output directory (temp): {temp_files_dir}")
        
        # 1. Fix PDB (merge flexible residues if any)
        print("\nStep 1: Fixing PDB structure...")
        fix_pdb_out = temp_files_dir / "protein_fixed.pdb"
        fixed_pdb_path = self.fix_pdb(output_path=str(fix_pdb_out))
        print(f"Fixed PDB: {fixed_pdb_path}")
        
        # 2. Add hydrogens using PDB2PQR
        print("\nStep 2: Adding hydrogens (protonation at pH {self._ph})...")
        protonated_out = temp_files_dir / "protein_protonated.pdb"
        protonated_pdb_path = self.add_hydrogens(
            input_pdb_path=fixed_pdb_path,
            output_pdb_path=str(protonated_out)
        )
        print(f"Protonated PDB: {protonated_pdb_path}")
        
        # 3. Fix histidine names (CHARMM compatibility)
        print("\nStep 3: Fixing histidine residue names...")
        final_out = fixed_files_dir / "protein_final.pdb"
        final_pdb_path = self.fix_histidine_names(
            input_pdb_path=protonated_pdb_path,
            output_pdb_path=str(final_out)
        )
        print(f"Final processed protein PDB: {final_pdb_path}")
        
        return final_pdb_path
    
    def process_ligand(self) -> str:
        """
        Process ligand: generate MOL2 and copy PDB to fixed_files.
        
        Returns:
            Path to the generated MOL2 file.
        """
        # Define output directory for fixed files
        md_prepper_dir = pathlib.Path(self._output_dir) / "md_prepper"
        fixed_files_dir = md_prepper_dir / "fixed_files"
        fixed_files_dir.mkdir(parents=True, exist_ok=True)
        
        output_mol2 = fixed_files_dir / "ligand.mol2"
        output_pdb = fixed_files_dir / "ligand.pdb"
        
        print("\nProcessing ligand pipeline...")
        print(f"Output directory: {fixed_files_dir}")
        
        # 1. Generate MOL2
        print("\nStep 1: Generating MOL2 file...")
        self.generate_mol2(
            pdb_path=self._ligand_pdb_path,
            smiles=self._ligand_smiles,
            output_mol2_path=str(output_mol2)
        )
        
        # 2. Copy PDB
        print("\nStep 2: Copying PDB file...")
        shutil.copy(self._ligand_pdb_path, output_pdb)
        print(f"Copied ligand PDB to: {output_pdb}")
        

        return str(output_mol2)

    def generate_system_psf(
        self,
        protein_pdb_path: str,
        ligand_pdb_path: str,
        ligand_str_path: typing.Optional[str],
        padding: float,
        salt_type: str,
        salt_concentration: float,
        temp_output_path: str,
        output_path: str
    ) -> typing.Tuple[str, str]:
        """
        Generate final system PSF and PDB by running the full pipeline.
        
        Steps:
        1. Generate protein PSF
        2. Generate ligand PSF (if ligand provided)
        3. Merge protein and ligand (if ligand provided)
        4. Solvate system
        5. Ionize system
        
        Args:
            protein_pdb_path: Path to prepared protein PDB
            ligand_pdb_path: Path to prepared ligand PDB (optional/empty if none)
            ligand_str_path: Path to ligand topology STR file (optional/empty if none)
            padding: padding for water box in Angstroms
            salt_type: Salt type (e.g., "NaCl")
            salt_concentration: Salt concentration in M
            temp_output_path: Directory for intermediate files
            output_path: Directory for final output files
            
        Returns:
            Tuple of (final_psf_path, final_pdb_path)
        """
        temp_dir = pathlib.Path(temp_output_path).resolve()
        temp_dir.mkdir(parents=True, exist_ok=True)
        final_dir = pathlib.Path(output_path).resolve()
        final_dir.mkdir(parents=True, exist_ok=True)
        
        # 1. Generate Protein PSF
        print(f"Generating protein PSF from {protein_pdb_path}")
        prot_psf, prot_pdb = self.generate_psf(
            input_pdb_path=protein_pdb_path,
            output_psf_path=str(temp_dir / "protein.psf"),
            output_pdb_path=str(temp_dir / "protein.pdb")
        )
        
        current_psf = prot_psf
        current_pdb = prot_pdb
        
        # 2. Generate Ligand PSF (if applicable)
        if ligand_pdb_path and pathlib.Path(ligand_pdb_path).exists():
            print(f"Generating ligand PSF from {ligand_pdb_path}")
            
            topo_files = [str(third_party_tools.get_cgenff_topology_file())]
            if ligand_str_path and pathlib.Path(ligand_str_path).exists():
                # Copy STR file to temp dir to avoid modifying the original
                temp_str_path = temp_dir / "ligand.str"
                shutil.copy(ligand_str_path, temp_str_path)
                
                # Fix residue name in the copy to match "UNL"
                print(f"Fixing residue name in {temp_str_path}...")
                self.fix_str_residue_name(str(temp_str_path), new_residue_name="UNL")
                
                topo_files.append(str(temp_str_path))
            
            lig_psf, lig_pdb = self.generate_ligand_psf(
                input_pdb_path=ligand_pdb_path,
                topology_files=topo_files,
                output_psf_path=str(temp_dir / "ligand.psf"),
                output_pdb_path=str(temp_dir / "ligand.pdb")
            )
            
            # 3. Merge PSF/PDB
            print("Merging protein and ligand...")
            psf_pdb_pairs = [(prot_psf, prot_pdb), (lig_psf, lig_pdb)]
            merged_psf, merged_pdb = self.merge_psf_pdb(
                psf_pdb_pairs=psf_pdb_pairs,
                output_psf_path=str(temp_dir / "complex.psf"),
                output_pdb_path=str(temp_dir / "complex.pdb")
            )
            current_psf = merged_psf
            current_pdb = merged_pdb
        
        # 4. Solvate System
        print(f"Solvating system with padding {padding}...")
        solvated_psf, solvated_pdb = self.solvate_system(
            input_psf_path=current_psf,
            input_pdb_path=current_pdb,
            output_psf_path=str(temp_dir / "solvated.psf"),
            output_pdb_path=str(temp_dir / "solvated.pdb"),
            padding=padding,
            use_cubic_box=True
        )
        
        # 5. Ionize System
        print(f"Ionizing system with {salt_concentration} M {salt_type}...")
        ionized_psf, ionized_pdb = self.ionize_system(
            input_psf_path=solvated_psf,
            input_pdb_path=solvated_pdb,
            output_psf_path=str(final_dir / "ionized.psf"),
            output_pdb_path=str(final_dir / "ionized.pdb"),
            salt_type=salt_type,
            salt_concentration=salt_concentration
        )
        
        return ionized_psf, ionized_pdb
    
    def copy_openmm_scripts(self, destination_dir: str, enhanced_sampling: bool = False) -> None:
        """
        Copy OpenMM scripts (equilibration.py, MMGBSA.py) to the specified destination.
        
        Args:
            destination_dir: Directory where scripts should be copied
            enhanced_sampling: If False, copy equilibration.py and MMGBSA.py.
                               If True, copy WTMtABF.py and run_WTM-TABF.py.
        """
        dest_path = pathlib.Path(destination_dir).resolve()
        dest_path.mkdir(parents=True, exist_ok=True)
        
        scripts_dir = third_party_tools.get_third_party_dir() / "openmm_scripts"
        
        if not scripts_dir.exists():
            print(f"Warning: OpenMM scripts directory not found at {scripts_dir}")
            return
            
        if enhanced_sampling:
            scripts_to_copy = ["WTMtABF.py", "run_WTM-TABF.py"]
        else:
            scripts_to_copy = ["equilibration.py", "MMGBSA.py"]
            
        for script_name in scripts_to_copy:
            src_file = scripts_dir / script_name
            if src_file.exists():
                dest_file = dest_path / script_name
                shutil.copy(src_file, dest_file)
                print(f"Copied {script_name} to {dest_file}")
            else:
                print(f"Warning: Script {script_name} not found in {scripts_dir}")



# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Haohao Fu <fhh2626@nankai.edu.cn, fhh2626@gmail.com>

# Standard library imports
import collections.abc
import csv
import pathlib
import subprocess

# Third-party imports
import MDAnalysis
import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.rdmolfiles
import rdkit.DataStructs
import rdkit.Chem.rdFingerprintGenerator
import rdkit.ML.Cluster.Butina

from . import third_party_tools

class Docking:
    """
    A class for performing molecular docking using AutoDock Vina-compatible engines.
    
    Supports multiple docking engines including smina, gnina, vina, qvina, and DSDP.
    
    Attributes:
        _protein_pdbqt: Path to the protein PDBQT file.
        _ligand_pdbqts: List of paths to ligand PDBQT files.
        _ligand_smiles: List of SMILES strings corresponding to ligand PDBQT files.
        _docking_range: Docking range containing center and size information.
        _exhaustiveness: Exhaustiveness of the global search.
        _num_modes: Maximum number of binding modes to generate.
        _scoring: Scoring function to use ('vina' or 'vinardo').
        _output_dir: Directory to save docking results.
        _energy_range: Maximum energy difference between best and worst binding mode.
        _docking_executable: Path to the docking engine executable.
        _docking_results: List of dictionaries containing docking results (populated after run_docking).
        _affinity_results: List of lists containing parsed affinity values for each ligand (populated after run_docking).
        _engine: Docking engine being used.
    """
    
    def __init__(
        self,
        protein_pdbqt: str | pathlib.Path,
        ligand_pdbqts: collections.abc.Sequence[str | pathlib.Path],
        ligand_smiles: collections.abc.Sequence[str],
        docking_range: collections.abc.Sequence,
        exhaustiveness: int = 32,
        num_modes: int = 10,
        scoring: str = "vina",
        output_dir: str | pathlib.Path = "docking_results",
        energy_range: float = 3.0,
        engine: str = "smina",
        flexible_residues: str = "",
        cnn_scoring: str = "rescore"
    ):
        """
        Initialize the Docking class.
        
        Args:
            protein_pdbqt: Path to the protein PDBQT file.
            ligand_pdbqts: List of paths to ligand PDBQT files.
            ligand_smiles: List of SMILES strings corresponding to ligand_pdbqts.
            docking_range: Docking range as [center, size] where:
                - center is a sequence of 3 floats [x, y, z]
                - size is a sequence of 3 floats [size_x, size_y, size_z]
            exhaustiveness: Exhaustiveness of the global search (default: 32).
            num_modes: Maximum number of binding modes to generate (default: 10).
            scoring: Scoring function to use, either 'vina' or 'vinardo' (default: 'vina').
            output_dir: Directory to save docking results (default: 'docking_results').
            energy_range: Maximum energy difference between best and worst binding mode (default: 3.0 kcal/mol).
            engine: Docking engine to use ('qvina2', 'qvinaw', 'vina-classic', 'vina-new', 'smina', 'DSDP', 'gnina').
                    - 'smina' (default): Uses bundled smina executable.
                    - 'vina-new': Uses bundled Vina with custom scoring support.
                    - 'gnina': Generates job package for server submission (supports CNN scoring).
                    - 'DSDP': Generates job package for server submission.
                    - others: Use system installed executables, forced to 'vina' scoring.
            flexible_residues: Flexible residues specification (e.g., "A:45,A:82,B:120").
                             Only used when engine is 'smina' or 'gnina'. Empty string means no flexible residues.
            cnn_scoring: CNN scoring mode for gnina ('none', 'rescore', 'refinement', 'all').
                        Default is 'rescore'. Only used when engine is 'gnina'.
        
        Raises:
            FileNotFoundError: If protein PDBQT file doesn't exist.
            ValueError: If ligand_pdbqts is empty, ligand_smiles length doesn't match, or docking_range format is invalid.
        """
        self._protein_pdbqt = pathlib.Path(protein_pdbqt)
        self._ligand_pdbqts = [pathlib.Path(ligand) for ligand in ligand_pdbqts]
        self._ligand_smiles = list(ligand_smiles)
        
        # Validate and store docking range
        if len(docking_range) != 2:
            raise ValueError(f"docking_range must contain [center, size], got {len(docking_range)} elements")
        
        center, size = docking_range
        if len(center) != 3:
            raise ValueError(f"Center must have 3 coordinates (x, y, z), got {len(center)}")
        if len(size) != 3:
            raise ValueError(f"Size must have 3 values (size_x, size_y, size_z), got {len(size)}")
        
        self._docking_range = docking_range
        
        # Store docking parameters
        self._exhaustiveness = exhaustiveness
        self._num_modes = num_modes
        self._energy_range = energy_range
        
        # Validate engine
        valid_engines = ["qvina2", "qvinaw", "vina-classic", "vina-new", "smina", "DSDP", "gnina"]
        if engine not in valid_engines:
            raise ValueError(f"engine must be one of {valid_engines}, got '{engine}'")
        self._engine = engine

        # Validate and store scoring function based on engine
        if self._engine in ["qvina2", "qvinaw", "vina-classic", "DSDP"]:
            # For these engines, force scoring to 'vina' (as per requirements)
            # The actual --scoring argument is omitted during execution for these
            self._scoring = "vina"
        else:
            # vina-new, smina, and gnina support custom scoring

            if scoring not in ["vina", "vinardo"]:
                raise ValueError(f"Scoring must be 'vina' or 'vinardo', got '{scoring}'")
            self._scoring = scoring
        
        # Store flexible residues parameter (only for smina and gnina)
        self._flexible_residues = flexible_residues if engine in ["smina", "gnina"] else ""
        
        # Validate and store CNN scoring mode (only for gnina)
        valid_cnn_modes = ["none", "rescore", "refinement", "all"]
        if cnn_scoring not in valid_cnn_modes:
            raise ValueError(f"cnn_scoring must be one of {valid_cnn_modes}, got '{cnn_scoring}'")
        self._cnn_scoring = cnn_scoring if engine == "gnina" else "rescore"
        
        # Store output directory as absolute path
        self._output_dir = pathlib.Path(output_dir).resolve()
        
        # Determine docking executable path or command using third_party_tools
        # For DSDP and gnina, we don't need a local executable
        if self._engine not in ["DSDP", "gnina"]:
            self._docking_executable = third_party_tools.get_docking_executable(self._engine)
        elif self._engine == "DSDP":
            self._docking_executable = "DSDP"
        elif self._engine == "gnina":
            self._docking_executable = "gnina"
        
        # Validate inputs
        if not self._protein_pdbqt.exists():
            raise FileNotFoundError(f"Protein PDBQT file not found: {self._protein_pdbqt}")
        
        # Determine and store docking type
        self._last_docking_type = "flexible" if (self._engine in ["smina", "gnina"] and self._flexible_residues) else "rigid"
        
        if not self._ligand_pdbqts:
            raise ValueError("Ligand PDBQT list cannot be empty")
        
        if len(self._ligand_smiles) != len(self._ligand_pdbqts):
            raise ValueError(f"Length of ligand_smiles ({len(self._ligand_smiles)}) must match ligand_pdbqts ({len(self._ligand_pdbqts)})")
        
        for ligand in self._ligand_pdbqts:
            if not ligand.exists():
                raise FileNotFoundError(f"Ligand PDBQT file not found: {ligand}")
        
        # Initialize docking results storage
        self._docking_results = []
        
        # Initialize affinity results storage (list of lists, one per ligand)
        self._affinity_results = []
        
        # Cache for combined results
        self._combined_results = None
        self._clustered_combined_results = None
        
        # Track the last docking type (rigid or flexible)
        self._last_docking_type = None
    
    def cluster_results(self, similarity_threshold: float = 0.7) -> None:
        """
        Cluster results based on chemical structure clustering.
        
        Uses RDKit to generate Morgan fingerprints and performs Butina clustering
        based on Tanimoto similarity.
        
        Args:
            similarity_threshold: Tanimoto similarity threshold for clustering (default: 0.7).
                                 Must be in range (0, 1].
        
        Raises:
            ValueError: If similarity_threshold is not in range (0, 1].
        """
        # Validate similarity_threshold parameter
        if not (0 < similarity_threshold <= 1):
            raise ValueError(
                f"similarity_threshold must be in range (0, 1], got {similarity_threshold}"
            )
        
        if self._combined_results is None:
            return

        # Generate fingerprints for valid molecules
        mols = []
        valid_indices = []
        for idx, item in enumerate(self._combined_results):
            if not item.get('smiles'):
                continue
            mol = rdkit.Chem.MolFromSmiles(item['smiles'])
            if mol:
                mols.append(mol)
                valid_indices.append(idx)
        
        if not mols:
            self._clustered_combined_results = [item.copy() for item in self._combined_results]
            for item in self._clustered_combined_results:
                item['cluster'] = 0
            return

        try:
            # New RDKit API (fixes deprecation warning)
            gen = rdkit.Chem.rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
            fps = [gen.GetFingerprint(x) for x in mols]
        except AttributeError:
            # Fallback for older RDKit versions
            fps = [rdkit.Chem.AllChem.GetMorganFingerprintAsBitVect(x, 2, 1024) for x in mols]
        
        # Calculate distance matrix (1 - similarity) for Butina clustering
        dists = []
        nfps = len(fps)
        
        for i in range(1, nfps):
            sims = rdkit.DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
            dists.extend([1 - x for x in sims])
            
        # Perform clustering
        cutoff = 1.0 - similarity_threshold
        clusters = rdkit.ML.Cluster.Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
        
        # Initialize clustered results with original data
        clustered_results = [item.copy() for item in self._combined_results]
        
        # Default cluster is 0 (unclustered/invalid)
        for item in clustered_results:
            item['cluster'] = 0
        
        # Assign cluster IDs
        for cluster_rank, cluster_indices in enumerate(clusters):
            cluster_id = cluster_rank + 1
            for idx_in_fps in cluster_indices:
                original_idx = valid_indices[idx_in_fps]
                clustered_results[original_idx]['cluster'] = cluster_id
                
        self._clustered_combined_results = clustered_results

    def save_results_to_csv(self, file_path: str) -> None:
        """
        Save the clustered combined results to a CSV file.
        
        Args:
            file_path: Full path to the output CSV file.
        """
        # Use clustered results if available, otherwise fallback to combined results
        data = getattr(self, '_clustered_combined_results', None)
        if data is None:
            data = self._combined_results
            
        if not data:
            return
            
        # Collect all fields
        # Ensure 'cluster' is included if present
        fieldnames = list(data[0].keys())
        
        # Ensure directory exists
        pathlib.Path(file_path).parent.mkdir(parents=True, exist_ok=True)
        
        try:
            with open(file_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(data)
        except OSError as e:
            print(f"Error saving CSV: {e}")

    
    def run_docking(
        self, 
        progress_callback: collections.abc.Callable[[int, int], None] | None = None
    ) -> list[dict]:

        """
        Perform docking for all ligands using AutoDock Vina.
        
        Args:
            progress_callback (callable, optional): A callback function to report progress.
                The callback will be called with (current, total) where:
                - current (int): Current ligand index being processed
                - total (int): Total number of ligands
        
        Returns:
            A list of dictionaries containing docking results for each ligand.
            Each dictionary contains:
                - 'ligand': Path to the input ligand
                - 'output': Path to the output file
                - 'log': Path to the log file
                - 'return_code': Return code from Vina
                - 'success': Boolean indicating if docking succeeded
        
        Raises:
            FileNotFoundError: If Vina executable is not found (checked in __init__).
        """
        # Reset combined results cache
        self._combined_results = None
        
        # Create output directory
        self._output_dir.mkdir(parents=True, exist_ok=True)
        
        # Extract center and size from docking_range
        center, size = self._docking_range
        
        results = []
        affinity_results = []
        total_ligands = len(self._ligand_pdbqts)
        
        # Run docking for each ligand
        for idx, ligand_pdbqt in enumerate(self._ligand_pdbqts):
            # Report progress
            if progress_callback:
                progress_callback(idx, total_ligands)
            
            result = self._dock_single_ligand(
                ligand_pdbqt=ligand_pdbqt,
                center=center,
                size=size,
                output_dir=self._output_dir,
                energy_range=self._energy_range
            )
            results.append(result)
            
            # Parse affinity values from the log file
            if result['success']:
                affinity_values = self._parse_affinity_from_log(result['log'])
                affinity_results.append(affinity_values)
            else:
                affinity_results.append([])
        
        # Report completion
        if progress_callback:
            progress_callback(total_ligands, total_ligands)
        
        # Store results for later retrieval
        self._docking_results = results
        self._affinity_results = affinity_results
        
        # Update combined results cache
        self._update_combined_results()
        
        return results
    
    def _dock_single_ligand(
        self,
        ligand_pdbqt: pathlib.Path,
        center: collections.abc.Sequence[float],
        size: collections.abc.Sequence[float],
        output_dir: pathlib.Path,
        energy_range: float
    ) -> dict:
        """
        Perform docking for a single ligand.
        
        Args:
            ligand_pdbqt: Path to the ligand PDBQT file.
            center: Center of the search space (x, y, z).
            size: Size of the search space (size_x, size_y, size_z).
            output_dir: Directory to save results.
            energy_range: Energy range parameter for the docking engine.
        
        Returns:
            Dictionary with docking results for this ligand.
        """
        # Prepare output filenames
        ligand_name = ligand_pdbqt.stem
        output_file = output_dir / f"{ligand_name}_out.pdbqt"
        log_file = output_dir / f"{ligand_name}_log.txt"
        
        # Construct docking command
        docking_cmd = [
            str(self._docking_executable),
            "--receptor", str(self._protein_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "--out", str(output_file),
            "--exhaustiveness", str(self._exhaustiveness),
            "--num_modes", str(self._num_modes),
            "--energy_range", str(energy_range)
        ]
        
        # Add scoring parameter for vina-new, smina, and gnina engines
        if self._engine in ["vina-new", "smina", "gnina"]:
            docking_cmd.extend(["--scoring", self._scoring])
        
        # Determine docking type and flexible residue output file
        is_flexible = self._engine in ["smina", "gnina"] and self._flexible_residues
        docking_type = "flexible" if is_flexible else "rigid"
        out_flex_file = None
        
        # Add flexible residues parameter for smina/gnina
        if is_flexible:
            docking_cmd.extend(["--flexres", self._flexible_residues])
            out_flex_file = output_dir / f"{ligand_name}_out_flex_res.pdbqt"
            docking_cmd.extend(["--out_flex", str(out_flex_file)])
        
        # Run docking
        try:
            result = subprocess.run(
                docking_cmd,
                capture_output=True,
                text=True,
                check=False,
                env=third_party_tools.get_docking_subprocess_env(self._docking_executable)
            )
            
            success = result.returncode == 0
            
            # Manually write log file from stdout (some docking engines don't support --log)
            if success or result.stdout:
                with open(log_file, 'w') as f:
                    f.write(result.stdout)
            
            return {
                'ligand': ligand_pdbqt,
                'output': output_file,
                'log': log_file,
                'return_code': result.returncode,
                'success': success,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'type': docking_type,
                'flexible_residue_output': out_flex_file
            }
        
        except Exception as e:
            return {
                'ligand': ligand_pdbqt,
                'output': output_file,
                'log': log_file,
                'return_code': -1,
                'success': False,
                'error': str(e),
                'type': docking_type,
                'flexible_residue_output': out_flex_file
            }
    
    def _parse_affinity_from_log(self, log_file_path: str | pathlib.Path) -> list[float]:
        """
        Parse affinity values from a docking log file.
        
        Args:
            log_file_path: Path to the log file to parse.
        
        Returns:
            list[float]: List of affinity values (in kcal/mol).
        """
        try:
            with open(log_file_path, 'r') as f:
                lines = f.readlines()
            
            # Check for Vina/Smina format
            if self._is_vina_format(lines):
                return self._parse_vina_log(lines)
            
            # Check for DSDP format
            if self._is_dsdp_format(lines):
                return self._parse_dsdp_log(lines)
                
            return []
            
        except FileNotFoundError:
            return []
        except Exception as e:
            print(f"Warning: Failed to parse affinity from {log_file_path}: {e}")
            return []

    def _is_vina_format(self, lines: list[str]) -> bool:
        """Check if log lines match Vina/Smina format."""
        for line in lines:
            if 'mode' in line and 'affinity' in line:
                return True
        return False
        
    def _is_dsdp_format(self, lines: list[str]) -> bool:
        """Check if log lines match DSDP format."""
        # DSDP format: check first few non-empty lines
        count = 0
        for line in lines:
            parts = line.strip().split()
            if not parts:
                continue
            
            # Check if format is: filename.pdbqt float_value
            if len(parts) >= 2 and parts[0].endswith('.pdbqt'):
                try:
                    float(parts[1])
                    return True
                except ValueError:
                    pass
            
            count += 1
            if count >= 5: # Only check first 5 non-empty lines
                break
        return False

    def _parse_vina_log(self, lines: list[str]) -> list[float]:
        """Parse Vina/Smina style log."""
        affinity_values = []
        table_start_idx = -1
        
        for i, line in enumerate(lines):
             if 'mode' in line and 'affinity' in line:
                 table_start_idx = i
                 break
        
        if table_start_idx == -1:
            return []
            
        for line in lines[table_start_idx + 2:]:
            line = line.strip()
            if not line or line.startswith('===') or line.startswith('---'):
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                try:
                    affinity_values.append(float(parts[1]))
                except ValueError:
                    continue
            else:
                break
        
        return affinity_values

    def _parse_dsdp_log(self, lines: list[str]) -> list[float]:
        """Parse DSDP style log."""
        affinity_values = []
        
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0].endswith('.pdbqt'):
                try:
                    # Second column is affinity
                    affinity_values.append(float(parts[1]))
                except ValueError:
                    continue
        
        return affinity_values

    
    def get_output_logs(self) -> list[pathlib.Path]:
        """
        Get the list of output log file paths from docking results.
        
        Returns:
            list[pathlib.Path]: List of paths to log files. Returns empty list if run_docking() hasn't been called.
        """
        if not self._docking_results:
            return []
        return [result['log'] for result in self._docking_results]
    
    def get_output_pdbqts(self) -> list[pathlib.Path]:
        """
        Get the list of output PDBQT file paths from docking results.
        
        Returns:
            list[pathlib.Path]: List of paths to output PDBQT files. Returns empty list if run_docking() hasn't been called.
        """
        if not self._docking_results:
            return []
        return [result['output'] for result in self._docking_results]
    
    def get_ligand_smiles(self) -> list[str]:
        """
        Get the list of input ligand SMILES strings.
        
        Returns:
            list[str]: List of SMILES strings corresponding to the input ligands.
        """
        return self._ligand_smiles.copy()
    
    def get_affinity_values(self) -> list[list[float]]:
        """
        Get the parsed affinity values from docking results.
        
        Returns:
            list[list[float]]: List of affinity value lists, one for each ligand.
                              Each inner list contains affinity values (in kcal/mol) 
                              for all binding modes of that ligand.
                              Returns empty list if run_docking() hasn't been called.
        
        Example:
            [[-1.913, -1.725, -1.721], [-2.5, -2.3]]  # 2 ligands with multiple modes each
        """
        return [affinities.copy() for affinities in self._affinity_results]
    
    def _update_combined_results(self) -> None:
        """Update the internal combined results cache."""
        smiles_list = self.get_ligand_smiles()
        affinity_values = self.get_affinity_values()
        # Get filenames from ligand pdbqt paths
        filenames = [p.name for p in self._ligand_pdbqts]
        
        combined_results = []
        for idx, (smiles, affinities, filename) in enumerate(zip(smiles_list, affinity_values, filenames)):
            best_affinity = min(affinities) if affinities else "N/A"
            
            combined_results.append({
                "number": idx + 1,
                "smiles": smiles,
                "filename": filename,
                "best_affinity": best_affinity,
                "affinity": affinities,
                "molecule_number": idx
            })
            
        self._combined_results = combined_results

    def get_combined_results(self) -> list[dict]:
        """
        Get combined results for the Results Groupbox table.
        
        Returns:
            list[dict]: A list of dictionaries, where each dictionary represents a row:
                - 'number': 1-based index (int)
                - 'smiles': Ligand SMILES (str)
                - 'filename': Ligand PDBQT filename (str)
                - 'best_affinity': Best affinity value or "N/A" (float | str)
                - 'affinity': List of affinity values (list[float])
                - 'molecule_number': 0-based index (int)
        """
        if self._combined_results is None:
            return []
        return self._combined_results
    
    def get_sorted_combined_results(self) -> list[dict]:
        """
        Get combined results sorted by best affinity (ascending / most negative first).
        
        Returns:
            list[dict]: Sorted list of combined result dictionaries.
        """
        results = self._combined_results if self._combined_results is not None else []
        
        def get_affinity_key(item):
            val = item['best_affinity']
            if val == "N/A":
                return float('inf')
            return float(val)
            
        return sorted(results, key=get_affinity_key)
    
    def get_last_docking_type(self) -> str | None:
        """
        Get the type of the last docking run.
        
        Returns:
            str | None: "rigid" if no flexible residues were used,
                       "flexible" if flexible residues were specified,
                       None if run_docking() hasn't been called yet.
        """
        return self._last_docking_type

    def get_redocking_range(self, ligand_index: int, pose_index: int) -> tuple[list[float], list[float]]:
        """
        Get the redocking range (center and length) for a specific pose.
        
        This method parses the output PDBQT file to extract atomic coordinates
        for the specified binding mode, then calculates a bounding box suitable
        for redocking.
        
        Args:
            ligand_index: Index of the ligand (0-based, corresponding to ligand_pdbqts).
            pose_index: Index of the pose/binding mode for the ligand (0-based).
        
        Returns:
            tuple[list[float], list[float]]: A tuple containing:
                - center: List of 3 floats [x, y, z] representing the center of the bounding box
                - length: List of 3 floats [size_x, size_y, size_z] representing the dimensions
        
        Raises:
            IndexError: If ligand_index or pose_index is out of range.
            FileNotFoundError: If output PDBQT file doesn't exist.
            RuntimeError: If run_docking() hasn't been called yet.
            ValueError: If no valid coordinates are found in the specified pose.
        
        Example:
            center, length = docking.get_redocking_range(0, 2)
            # Returns: ([10.5, 20.3, 15.7], [5.2, 8.3, 6.1])
        """
        # Check if docking has been run
        if not self._docking_results:
            raise RuntimeError("run_docking() must be called before getting redocking range")
        
        # Validate ligand index
        if ligand_index < 0 or ligand_index >= len(self._docking_results):
            raise IndexError(
                f"ligand_index {ligand_index} out of range [0, {len(self._docking_results) - 1}]"
            )
        
        # Get the output PDBQT file for the specified ligand
        output_pdbqt = self._docking_results[ligand_index]['output']
        
        # Validate pose index
        if pose_index < 0:
            raise IndexError(f"pose_index must be non-negative, got {pose_index}")
        
        # Use shared parsing function to get ATOM lines
        atom_lines = self._parse_pdbqt_model_atoms(output_pdbqt, pose_index)
        
        # Parse coordinates from ATOM/HETATM lines
        # PDB/PDBQT format: columns 31-38 (x), 39-46 (y), 47-54 (z)
        coordinates = []
        for line in atom_lines:
            if line.startswith(('ATOM', 'HETATM')) and len(line) >= 54:
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coordinates.append([x, y, z])
                except ValueError:
                    continue
        
        # Check if we have valid positions
        if len(coordinates) == 0:
            raise ValueError(
                f"No atoms found for ligand {ligand_index}, pose {pose_index}"
            )
        
        # Calculate bounding box
        min_x = min(c[0] for c in coordinates)
        min_y = min(c[1] for c in coordinates)
        min_z = min(c[2] for c in coordinates)
        max_x = max(c[0] for c in coordinates)
        max_y = max(c[1] for c in coordinates)
        max_z = max(c[2] for c in coordinates)
        
        # Calculate center and size
        center = [
            (min_x + max_x) / 2.0,
            (min_y + max_y) / 2.0,
            (min_z + max_z) / 2.0
        ]
        length = [
            max_x - min_x,
            max_y - min_y,
            max_z - min_z
        ]
        
        return center, length

    def _create_protein_without_flex_residues(self, output_path: pathlib.Path) -> pathlib.Path:
        """
        Create a temporary protein PDB file with flexible residue sidechains removed.
        
        Uses MDAnalysis to read the protein PDBQT, remove only the sidechain atoms 
        (non-backbone atoms) belonging to flexible residues, and write the result 
        to a new PDB file for PyMOL visualization. Backbone atoms (N, CA, C, O) 
        are preserved to maintain protein chain connectivity.
        
        Args:
            output_path: Directory to save the temporary protein file.
        
        Returns:
            pathlib.Path: Path to the temporary protein PDB file.
        
        Note:
            The flexible residues are parsed from self._flexible_residues format:
            "A:45,A:82,B:120" -> [(A, 45), (A, 82), (B, 120)]
            
            Backbone atoms preserved: N, CA, C, O
            Sidechain atoms removed: All other atoms in the flexible residues
        """
        # Parse flexible residues from format "A:45,A:82,B:120"
        flex_res_list = []
        for res_spec in self._flexible_residues.split(","):
            res_spec = res_spec.strip()
            if ":" in res_spec:
                chain_id, res_num = res_spec.split(":")
                flex_res_list.append((chain_id.strip(), int(res_num.strip())))
        
        if not flex_res_list:
            # No flexible residues to remove, return original protein
            return self._protein_pdbqt
        
        # Load protein with MDAnalysis
        u = MDAnalysis.Universe(str(self._protein_pdbqt))
        
        # Build selection string to exclude sidechain atoms of flexible residues
        # Keep backbone atoms (N, CA, C, O) but remove sidechain atoms
        # Format: "not ((resid 45 and (segid A or chainID A)) and not backbone)"
        # Note: In PDBQT, chain ID is often in 'segid' or 'chainID' field
        exclusion_parts = []
        for chain_id, res_num in flex_res_list:
            # Select sidechain atoms only (not backbone) for each flexible residue
            # Try both segid and chainID as MDAnalysis may use either depending on file format
            exclusion_parts.append(
                f"((resid {res_num} and (segid {chain_id} or chainID {chain_id})) and not backbone)"
            )
        
        if exclusion_parts:
            selection_str = "not (" + " or ".join(exclusion_parts) + ")"
        else:
            selection_str = "all"
        
        # Select atoms excluding sidechain atoms of flexible residues
        try:
            atoms_to_keep = u.select_atoms(selection_str)
        except Exception:
            # Fallback: try simpler selection without chainID
            exclusion_parts = []
            for chain_id, res_num in flex_res_list:
                # Exclude sidechain atoms only
                exclusion_parts.append(f"(resid {res_num} and not backbone)")
            selection_str = "not (" + " or ".join(exclusion_parts) + ")"
            atoms_to_keep = u.select_atoms(selection_str)
        
        # Write to temporary PDB file (PDB is more compatible with PyMOL)
        temp_protein_path = output_path / f"protein_no_flex_sidechain_{self._protein_pdbqt.stem}.pdb"
        atoms_to_keep.write(str(temp_protein_path))
        
        return temp_protein_path

    def _parse_pdbqt_model_atoms(
        self, 
        pdbqt_file: pathlib.Path, 
        pose_index: int
    ) -> list[str]:
        """
        Parse a PDBQT file and extract ATOM/HETATM lines for a specific pose.
        
        This method handles PDBQT files where:
        1. MODEL numbering may not start at 1 (e.g., Smina output starts at 1 for ligands but
           may start at different numbers for flexible residues depending on the pose)
        2. Multiple MODEL blocks may share the same ID (for flexible residues, each residue
           gets its own MODEL block with the same pose number)
        
        Args:
            pdbqt_file: Path to the PDBQT file.
            pose_index: Index of the pose to extract (0-based).
        
        Returns:
            list[str]: List of ATOM/HETATM lines for the specified pose.
        """
        if not pdbqt_file.exists():
            raise FileNotFoundError(f"PDBQT file not found: {pdbqt_file}")
        
        # First pass: find the minimum MODEL ID in the file
        min_model_id = None
        found_model_tag = False
        
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('MODEL'):
                    found_model_tag = True
                    parts = line.split()
                    if len(parts) >= 2 and parts[1].isdigit():
                        model_id = int(parts[1])
                        if min_model_id is None or model_id < min_model_id:
                            min_model_id = model_id
        
        # Handle file with no MODEL tags (single model)
        if not found_model_tag and pose_index == 0:
             with open(pdbqt_file, 'r') as f:
                 return [l for l in f if l.startswith(('ATOM', 'HETATM'))]
        
        if min_model_id is None:
            return []
        
        # Calculate target MODEL ID: min_id + pose_index
        # e.g., if file starts with MODEL 3 and pose_index=0, target is 3
        # if pose_index=1, target is 4, etc.
        target_model_id = min_model_id + pose_index
        
        # Second pass: collect all ATOM/HETATM lines from blocks with target MODEL ID
        lines = []
        in_desired_model = False
        
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('MODEL'):
                    parts = line.split()
                    if len(parts) >= 2 and parts[1].isdigit():
                        model_id = int(parts[1])
                        if model_id == target_model_id:
                            in_desired_model = True
                        elif model_id > target_model_id:
                            # Models are sorted, we can stop
                            break
                        else:
                            in_desired_model = False
                    else:
                        in_desired_model = False
                        
                elif line.startswith('ENDMDL'):
                    # Do not break or reset here - continue to find other blocks
                    # with the same MODEL ID (for flexible residues)
                    pass
                
                elif in_desired_model:
                    if line.startswith(('ATOM', 'HETATM')):
                        lines.append(line)
                        
        return lines

    def generate_dsdp_job_package(self) -> str:
        """
        Generate a job package for DSDP to be run on Linux.
        
        Creates a directory containing:
        - The receptor PDBQT.
        - All ligand PDBQT files.
        - A shell script (run_dsdp.sh) to run the docking.
        - An index file (index.csv) mapping filenames to SMILES.
        
        Returns:
            str: Path to the generated job directory.
        """
        import shutil
        import csv
        
        # Create job directory
        job_dir = self._output_dir / "dsdp_job"
        job_dir.mkdir(parents=True, exist_ok=True)
        
        # Create ligands subdirectory
        ligands_dir = job_dir / "ligands"
        ligands_dir.mkdir(exist_ok=True)
        
        # Create results subdirectory (to keep root clean)
        results_dir = job_dir / "results"
        results_dir.mkdir(exist_ok=True)
        
        # 1. Copy protein file
        # Use simple name for the script
        new_protein_name = "receptor.pdbqt"
        new_protein_path = job_dir / new_protein_name
        shutil.copy2(self._protein_pdbqt, new_protein_path)
        
        # 2. Prepare script and index
        script_lines = [
            "#!/bin/bash",
            "# Generated by BFEE-Docking for DSDP",
            "# Make sure 'DSDP' is in your PATH",
            "",
            "mkdir -p results",
            ""
        ]
        
        index_rows = [["Filename", "OriginalPath", "SMILES"]]
        
        center, size = self._docking_range
        
        # Calculate box min/max
        box_min = [c - s/2 for c, s in zip(center, size)]
        box_max = [c + s/2 for c, s in zip(center, size)]
        
        # Calculate padding width based on total ligands to ensure correct sorting
        # e.g., if 10000 ligands, padding will be 5 digits (ligand_00001)
        total_ligands = len(self._ligand_pdbqts)
        pad_width = len(str(total_ligands))
        # Ensure at least 4 digits for consistency with small sets
        pad_width = max(pad_width, 4)
        
        for idx, (ligand_path, smiles) in enumerate(zip(self._ligand_pdbqts, self._ligand_smiles)):
            # Use safe filename with dynamic padding
            # f"{val:0{width}d}" is the syntax for dynamic width
            safe_name = f"ligand_{idx + 1:0{pad_width}d}_{ligand_path.name}"
            dest_path = ligands_dir / safe_name
            shutil.copy2(ligand_path, dest_path)
            
            # Record in index
            index_rows.append([safe_name, str(ligand_path), smiles])
            
            # Construct command
            # Using relative paths from the script location
            # Use .stem to safely remove extension
            stem_name = pathlib.Path(safe_name).stem
            out_pdbqt = f"results/{stem_name}_out.pdbqt"
            out_log = f"results/{stem_name}_log.txt"
            
            # DSDP command construction
            # DSDP uses --protein instead of --receptor
            # DSDP uses --box_min/max instead of center/size
            # DSDP uses --top_n instead of --num_modes
            # DSDP has explicit --log argument
            # DSDP does not appear to support --energy_range in the provided help
            # Quote paths to handle spaces
            cmd = (
                f"DSDP "
                f"--protein \"{new_protein_name}\" "
                f"--ligand \"ligands/{safe_name}\" "
                f"--box_min {box_min[0]:.3f} {box_min[1]:.3f} {box_min[2]:.3f} "
                f"--box_max {box_max[0]:.3f} {box_max[1]:.3f} {box_max[2]:.3f} "
                f"--out \"{out_pdbqt}\" "
                f"--log \"{out_log}\" "
                f"--exhaustiveness {self._exhaustiveness} "
                f"--top_n {self._num_modes}"
            )
            
            script_lines.append(f"echo 'Processing {safe_name}...'")
            script_lines.append(cmd)
        
        # 3. Write Shell Script
        script_path = job_dir / "run_dsdp.sh"
        with open(script_path, "w", encoding="utf-8", newline='\n') as f:
            f.write("\n".join(script_lines))
            
        # 4. Write Index CSV
        index_path = job_dir / "index.csv"
        with open(index_path, "w", encoding="utf-8", newline='') as f:
            writer = csv.writer(f)
            writer.writerows(index_rows)
        
        # 5. Write docking configuration file for result loading
        import json
        config_data = {
            "engine": "DSDP",
            "docking_type": "rigid",  # DSDP only supports rigid docking
            "flexible_residues": "",
            "scoring": self._scoring,
            "exhaustiveness": self._exhaustiveness,
            "num_modes": self._num_modes
        }
        config_path = job_dir / "docking_config.json"
        with open(config_path, "w", encoding="utf-8") as f:
            json.dump(config_data, f, indent=2)
            
        return str(job_dir)

    def generate_gnina_job_package(self) -> str:
        """
        Generate a job package for Gnina to be run on Linux.
        
        Creates a directory containing:
        - The receptor PDBQT.
        - All ligand PDBQT files.
        - A shell script (run_gnina.sh) to run the docking.
        - An index file (index.csv) mapping filenames to SMILES.
        
        Returns:
            str: Path to the generated job directory.
        """
        import shutil
        import csv
        
        # Create job directory
        job_dir = self._output_dir / "gnina_job"
        job_dir.mkdir(parents=True, exist_ok=True)
        
        # Create ligands subdirectory
        ligands_dir = job_dir / "ligands"
        ligands_dir.mkdir(exist_ok=True)
        
        # Create results subdirectory
        results_dir = job_dir / "results"
        results_dir.mkdir(exist_ok=True)
        
        # 1. Copy protein file
        new_protein_name = "receptor.pdbqt"
        new_protein_path = job_dir / new_protein_name
        shutil.copy2(self._protein_pdbqt, new_protein_path)
        
        # 2. Prepare script and index
        script_lines = [
            "#!/bin/bash",
            "# Generated by BFEE-Docking for Gnina",
            "# Make sure 'gnina' is in your PATH",
            "",
            "mkdir -p results",
            ""
        ]
        
        index_rows = [["Filename", "OriginalPath", "SMILES"]]
        
        center, size = self._docking_range
        
        # Calculate padding width
        total_ligands = len(self._ligand_pdbqts)
        pad_width = len(str(total_ligands))
        pad_width = max(pad_width, 4)
        
        # Check for flexible residues
        is_flexible = self._flexible_residues != ""
        
        for idx, (ligand_path, smiles) in enumerate(zip(self._ligand_pdbqts, self._ligand_smiles)):
            safe_name = f"ligand_{idx + 1:0{pad_width}d}_{ligand_path.name}"
            dest_path = ligands_dir / safe_name
            shutil.copy2(ligand_path, dest_path)
            
            # Record in index
            index_rows.append([safe_name, str(ligand_path), smiles])
            
            # Construct filenames
            stem_name = pathlib.Path(safe_name).stem
            out_pdbqt = f"results/{stem_name}_out.pdbqt"
            out_log = f"results/{stem_name}_log.txt"
            
            # Gnina command construction (same as Smina)
            # Using relative paths
            cmd_parts = [
                "gnina",
                "--receptor", f"\"{new_protein_name}\"",
                "--ligand", f"\"ligands/{safe_name}\"",
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(size[0]),
                "--size_y", str(size[1]),
                "--size_z", str(size[2]),
                "--out", f"\"{out_pdbqt}\"",
                "--exhaustiveness", str(self._exhaustiveness),
                "--num_modes", str(self._num_modes),
                # "--energy_range", str(self._energy_range), # Removed as per requirement
                "--scoring", self._scoring,
                "--cnn_scoring", self._cnn_scoring
            ]
            
            if is_flexible:
                out_flex = f"results/{stem_name}_out_flex_res.pdbqt"
                cmd_parts.extend(["--flexres", self._flexible_residues])
                cmd_parts.extend(["--out_flex", f"\"{out_flex}\""])
            
            # Redirect stdout to log file (since gnina/smina prints to stdout)
            cmd = " ".join(cmd_parts) + f" > \"{out_log}\""
            
            script_lines.append(f"echo 'Processing {safe_name}...'")
            script_lines.append(cmd)
        
        # 3. Write Shell Script
        script_path = job_dir / "run_gnina.sh"
        with open(script_path, "w", encoding="utf-8", newline='\n') as f:
            f.write("\n".join(script_lines))
            
        # 4. Write Index CSV
        index_path = job_dir / "index.csv"
        with open(index_path, "w", encoding="utf-8", newline='') as f:
            writer = csv.writer(f)
            writer.writerows(index_rows)
        
        # 5. Write docking configuration file for result loading
        import json
        config_data = {
            "engine": "gnina",
            "docking_type": "flexible" if is_flexible else "rigid",
            "flexible_residues": self._flexible_residues if is_flexible else "",
            "scoring": self._scoring,
            "cnn_scoring": self._cnn_scoring,
            "exhaustiveness": self._exhaustiveness,
            "num_modes": self._num_modes
        }
        config_path = job_dir / "docking_config.json"
        with open(config_path, "w", encoding="utf-8") as f:
            json.dump(config_data, f, indent=2)
            
        return str(job_dir)

    def load_dsdp_results(self, result_dir: str | pathlib.Path) -> list[dict]:
        """
        Load docking results from a DSDP or Gnina job directory.
        
        This method is used for loading results from server-based engines
        (DSDP, Gnina) after the job has been completed remotely.
        
        Args:
            result_dir: Path to the directory containing results (the folder containing index.csv).
        
        Returns:
            list[dict]: List of docking results.
        """
        import csv
        result_dir = pathlib.Path(result_dir)
        index_path = result_dir / "index.csv"
        
        print(f"Loading DSDP results from: {result_dir}")
        if not index_path.exists():
            print(f"Error: index.csv not found at {index_path}")
            raise FileNotFoundError(f"index.csv not found in {result_dir}. Cannot Map results safely.")
        
        # Load docking configuration if available
        import json
        config_path = result_dir / "docking_config.json"
        config_data = None
        if config_path.exists():
            try:
                with open(config_path, "r", encoding="utf-8") as f:
                    config_data = json.load(f)
                print(f"Docking config loaded: engine={config_data.get('engine')}, "
                      f"type={config_data.get('docking_type')}, "
                      f"flexible_residues={config_data.get('flexible_residues', '')}")
                
                # Update internal settings from config
                if config_data.get('flexible_residues'):
                    self._flexible_residues = config_data['flexible_residues']
                    print(f"Flexible residues restored: {self._flexible_residues}")
                if config_data.get('scoring'):
                    self._scoring = config_data['scoring']
            except Exception as e:
                print(f"Warning: Failed to load docking config: {e}")
        else:
            print("No docking_config.json found, using defaults")
            
        # Read index to map filenames to original context
        # We store multiple keys to help matching: 
        # 1. Full original path (normalized)
        # 2. Filename only
        index_map_path = {} 
        index_map_name = {}
        
        with open(index_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Normalize path for better matching
                orig_path_str = str(pathlib.Path(row["OriginalPath"])).lower().replace('\\', '/')
                entry_data = {
                    "safe_filename": row["Filename"],
                    "smiles": row["SMILES"]
                }
                index_map_path[orig_path_str] = entry_data
                
                # Also map by just the filename for fallback
                orig_name = pathlib.Path(row["OriginalPath"]).name
                index_map_name[orig_name] = entry_data

        print(f"Index loaded with {len(index_map_path)} entries.")
        
        loaded_results = []
        affinity_results = []
        
        # We iterate over the ligands currently configured in this Docking instance
        for ligand_path in self._ligand_pdbqts:
            # Prepare lookup keys
            ligand_path_str = str(ligand_path).lower().replace('\\', '/')
            ligand_name = ligand_path.name
            
            # Try exact path match first, then filename match
            entry = index_map_path.get(ligand_path_str)
            if not entry:
                entry = index_map_name.get(ligand_name)
                if entry:
                    print(f"Match found by filename for {ligand_name}")
                else:
                    print(f"Warning: No index entry found for {ligand_name} (path: {ligand_path})")
            
            result_data = {
                'ligand': ligand_path,
                'output': None,
                'log': None,
                'return_code': 0,
                'success': False,
                'type': 'rigid'
            }
            
            affinities = []
            
            if entry:
                safe_filename = entry["safe_filename"]
                # safe_filename usually looks like "ligand_0001_name.pdbqt"
                # We need to construct the output filenames carefully
                stem = pathlib.Path(safe_filename).stem
                out_pdbqt_name = f"{stem}_out.pdbqt"
                out_log_name = f"{stem}_log.txt"
                
                # Check possibilities: result_dir/results/FILE or result_dir/FILE
                file_locs = [
                    result_dir / "results" / out_pdbqt_name,
                    result_dir / out_pdbqt_name
                ]
                
                found_out = None
                for loc in file_locs:
                    if loc.exists():
                        found_out = loc
                        break
                        
                log_locs = [
                    result_dir / "results" / out_log_name,
                    result_dir / out_log_name
                ]
                
                found_log = None
                for loc in log_locs:
                    if loc.exists():
                        found_log = loc
                        break
                
                if found_out and found_log:
                    result_data['output'] = found_out
                    result_data['log'] = found_log
                    result_data['success'] = True
                    affinities = self._parse_affinity_from_log(found_log)
                    if not affinities:
                         print(f"Warning: Log found for {ligand_name} but no affinities parsed. path: {found_log}")
                    
                    # Check for flexible residues output
                    out_flex_name = f"{stem}_out_flex_res.pdbqt"
                    flex_locs = [
                        result_dir / "results" / out_flex_name,
                        result_dir / out_flex_name
                    ]
                    
                    found_flex = None
                    for loc in flex_locs:
                        if loc.exists():
                            found_flex = loc
                            break
                    
                    if found_flex:
                        result_data['type'] = 'flexible'
                        result_data['flexible_residue_output'] = found_flex
                        print(f"Found flexible residues for {ligand_name}")
                else:
                    print(f"Warning: Output/Log file missing for {ligand_name}. Expected at {result_dir / 'results' / out_log_name}")
            
            loaded_results.append(result_data)
            affinity_results.append(affinities)
            
        self._docking_results = loaded_results
        self._affinity_results = affinity_results
        
        # Update last docking type based on config or detected flex files
        if config_data and config_data.get('docking_type'):
            self._last_docking_type = config_data['docking_type']
        else:
            # Fallback: check if any result has flexible type
            has_flex = any(r.get('type') == 'flexible' for r in loaded_results)
            self._last_docking_type = 'flexible' if has_flex else 'rigid'
        print(f"Docking type set to: {self._last_docking_type}")
        
        self._update_combined_results()
        
        return loaded_results

    # Alias for clarity when loading gnina results
    load_gnina_results = load_dsdp_results

    def _extract_flex_residues_pose(
        self, 
        flex_res_file: pathlib.Path, 
        ligand_index: int,
        pose_index: int, 
        output_path: pathlib.Path
    ) -> pathlib.Path:
        """
        Extract a specific pose from the flexible residues output file.
        
        Smina's flexible residue output has a unique format where each pose
        has multiple MODEL blocks (one per flexible residue), all sharing the
        same MODEL number. This method extracts all MODEL blocks for the
        specified pose and combines them into a single file.
        
        Args:
            flex_res_file: Path to the flexible residues PDBQT file.
            ligand_index: Index of the ligand (0-based), used for output filename.
            pose_index: Index of the pose to extract (0-based).
            output_path: Directory to save the extracted pose file.
        
        Returns:
            pathlib.Path: Path to the extracted pose PDB file.
        
        Raises:
            IndexError: If pose_index is out of range.
            FileNotFoundError: If flex_res_file doesn't exist.
        """
        # Use shared parsing function
        filtered_lines = self._parse_pdbqt_model_atoms(flex_res_file, pose_index)
        
        # Combine all filtered lines for this pose
        extracted_content = ''.join(filtered_lines)
        
        # Write to a PDB file (0-based indices for both ligand and pose)
        extracted_file = output_path / f"flex_residues_ligand{ligand_index}_pose{pose_index}.pdb"
        with open(extracted_file, 'w') as f:
            f.write(extracted_content)
        
        return extracted_file


    def output_MD_prepper_files(
        self, 
        ligand_index: int, 
        pose_index: int, 
        output_dir: str | pathlib.Path
    ) -> dict[str, pathlib.Path]:
        """
        Output files required for MD preparation (molecular dynamics simulation setup).
        
        This method exports the protein, ligand, and optionally flexible residues
        to the specified directory in formats suitable for MD preprocessing tools.
        
        Args:
            ligand_index: Index of the ligand (0-based, corresponding to ligand_pdbqts).
            pose_index: Index of the pose/binding mode for the ligand (0-based).
            output_dir: Directory to save the output files.
        
        Returns:
            dict[str, pathlib.Path]: Dictionary containing paths to the generated files:
                - 'rigid_protein': Path to rigid_protein.pdb
                - 'ligand': Path to ligand.pdb
                - 'ligand_smi': Path to ligand.smi
                - 'flexible_residue': Path to flexible_residue.pdb (only if flexible docking)
        
        Raises:
            IndexError: If ligand_index or pose_index is out of range.
            RuntimeError: If run_docking() hasn't been called yet.
            FileNotFoundError: If required files don't exist.
        
        Note:
            - rigid_protein.pdb: Copy of the protein PDBQT file
            - ligand.pdb: Extracted ATOM/HETATM lines for the specified pose
            - ligand.smi: SMILES string of the ligand
            - flexible_residue.pdb: Flexible residue atoms (only for flexible docking)
        """
        # Check if docking has been run
        if not self._docking_results:
            raise RuntimeError("run_docking() must be called before outputting MD prepper files")
        
        # Validate ligand index
        if ligand_index < 0 or ligand_index >= len(self._docking_results):
            raise IndexError(
                f"ligand_index {ligand_index} out of range [0, {len(self._docking_results) - 1}]"
            )
        
        # Validate pose index
        if pose_index < 0:
            raise IndexError(f"pose_index must be non-negative, got {pose_index}")
        
        # Check if pose_index is within the valid range for this ligand
        if self._affinity_results and ligand_index < len(self._affinity_results):
            num_poses = len(self._affinity_results[ligand_index])
            if num_poses > 0 and pose_index >= num_poses:
                raise IndexError(
                    f"pose_index {pose_index} out of range for ligand {ligand_index}. "
                    f"Valid range: [0, {num_poses - 1}] ({num_poses} poses available)"
                )
        
        # Prepare output directory
        output_path = pathlib.Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize result dictionary
        result_files = {}
        
        # 1. Copy protein PDBQT as rigid_protein.pdb
        rigid_protein_path = output_path / "rigid_protein.pdb"
        with open(self._protein_pdbqt, 'r') as src:
            with open(rigid_protein_path, 'w') as dst:
                dst.write(src.read())
        result_files['rigid_protein'] = rigid_protein_path
        
        # 2. Extract ligand PDB using _parse_pdbqt_model_atoms
        docking_result = self._docking_results[ligand_index]
        output_pdbqt = docking_result['output']
        
        ligand_atom_lines = self._parse_pdbqt_model_atoms(output_pdbqt, pose_index)
        ligand_pdb_path = output_path / "ligand.pdb"
        with open(ligand_pdb_path, 'w') as f:
            f.write(''.join(ligand_atom_lines))
        result_files['ligand'] = ligand_pdb_path
        
        # 3. Output ligand SMILES
        ligand_smi_path = output_path / "ligand.smi"
        smiles = self._ligand_smiles[ligand_index]
        with open(ligand_smi_path, 'w') as f:
            f.write(smiles)
        result_files['ligand_smi'] = ligand_smi_path
        
        # 4. For flexible docking, extract flexible residue PDB
        docking_type = docking_result.get('type', 'rigid')
        flex_res_output = docking_result.get('flexible_residue_output')
        
        if docking_type == "flexible" and flex_res_output and pathlib.Path(flex_res_output).exists():
            # Check if the flex residue file has valid ATOM/HETATM lines
            try:
                with open(flex_res_output, 'r') as f:
                    content = f.read()
                has_atoms = any(
                    line.startswith(('ATOM', 'HETATM')) 
                    for line in content.splitlines()
                )
                
                if has_atoms:
                    # Use _parse_pdbqt_model_atoms to extract flexible residue atoms
                    flex_atom_lines = self._parse_pdbqt_model_atoms(
                        pathlib.Path(flex_res_output), 
                        pose_index
                    )
                    
                    flexible_residue_path = output_path / "flexible_residues.pdb"
                    with open(flexible_residue_path, 'w') as f:
                        f.write(''.join(flex_atom_lines))
                    result_files['flexible_residue'] = flexible_residue_path
            except Exception as e:
                print(f"Warning: Could not extract flexible residues: {e}")
        
        return result_files

    def open_results(self, ligand_index: int, pose_index: int) -> None:
        """
        Open docking results in PyMOL with orthoscopic view.
        
        This method opens the protein PDBQT and a specific ligand pose in PyMOL.
        The ligand is displayed as sticks for the specified binding mode.
        
        For flexible docking results, the flexible residues are also loaded and
        displayed, while the corresponding residues are removed from the protein
        to avoid duplication.
        
        Args:
            ligand_index: Index of the ligand (0-based, corresponding to ligand_pdbqts).
            pose_index: Index of the pose/binding mode for the ligand (0-based).
        
        Raises:
            IndexError: If ligand_index or pose_index is out of range.
            FileNotFoundError: If output PDBQT files don't exist.
            RuntimeError: If run_docking() hasn't been called yet.
        
        Note:
            - Both parameters are 0-based in Python convention
            - PyMOL states are 1-based, so pose_index is converted internally
            - PyMOL must be installed and available in system PATH
            - For flexible docking, a temporary protein file is created with
              flexible residues removed
        
        Example:
            docking.open_results(0, 2)  # Opens first ligand, third pose (state 3 in PyMOL)
        """
        # Check if docking has been run
        if not self._docking_results:
            raise RuntimeError("run_docking() must be called before opening results")
        
        # Validate ligand index
        if ligand_index < 0 or ligand_index >= len(self._docking_results):
            raise IndexError(
                f"ligand_index {ligand_index} out of range [0, {len(self._docking_results) - 1}]"
            )
        
        # Get docking result for the specified ligand
        docking_result = self._docking_results[ligand_index]
        output_pdbqt = docking_result['output']
        docking_type = docking_result.get('type', 'rigid')
        flex_res_output = docking_result.get('flexible_residue_output')
        
        # Check if output file exists
        if not output_pdbqt.exists():
            raise FileNotFoundError(f"Output PDBQT file not found: {output_pdbqt}")
        
        # Convert pose_index from 0-based to 1-based for PyMOL state
        pymol_state = pose_index + 1
        
        # Validate pose index
        if pose_index < 0:
            raise IndexError(f"pose_index must be non-negative, got {pose_index}")
        
        # Check if pose_index is within the valid range for this ligand
        if self._affinity_results and ligand_index < len(self._affinity_results):
            num_poses = len(self._affinity_results[ligand_index])
            if num_poses > 0 and pose_index >= num_poses:
                raise IndexError(
                    f"pose_index {pose_index} out of range for ligand {ligand_index}. "
                    f"Valid range: [0, {num_poses - 1}] ({num_poses} poses available)"
                )
        
        # Get ligand name for PyMOL object naming
        ligand_name = output_pdbqt.stem
        
        # Determine which protein to load and build PyMOL command
        # Check if flexible residue file has valid content (non-empty with ATOM/HETATM lines)
        use_flexible_mode = False
        if docking_type == "flexible" and flex_res_output and flex_res_output.exists():
            # Check if the flex residue file has valid ATOM/HETATM lines
            try:
                with open(flex_res_output, 'r') as f:
                    content = f.read()
                has_atoms = any(
                    line.startswith(('ATOM', 'HETATM')) 
                    for line in content.splitlines()
                )
                if has_atoms:
                    use_flexible_mode = True
                else:
                    print(f"Note: Flexible residue file is empty or has no valid atoms. "
                          f"Falling back to rigid docking visualization.")
            except Exception as e:
                print(f"Warning: Could not read flexible residue file: {e}. "
                      f"Falling back to rigid docking visualization.")
        
        if use_flexible_mode:
            # Flexible docking: create protein without flex residues and extract flex residues for specific pose
            protein_to_load = self._create_protein_without_flex_residues(self._output_dir)
            
            # Extract the specific pose from flexible residue file to ensure all residues are loaded
            # This solves the issue where PyMOL only displays the first residue from smina's output
            flex_res_extracted = self._extract_flex_residues_pose(
                pathlib.Path(flex_res_output), 
                ligand_index,
                pose_index, 
                self._output_dir
            )
            flex_res_name = "flex_residues"
            
            # Construct PyMOL command for flexible docking
            # Note: The extracted file contains only the specified pose, so we load it directly
            # without needing to select a specific state
            # Use 'frame' to switch to specific state and 'set all_states, off' to hide other states
            pymol_cmd = [
                "pymol",
                "-d",
                (
                    f"set orthoscopic, on; "
                    f"load {protein_to_load}, protein; "
                    f"load {output_pdbqt}, {ligand_name}; "
                    f"load {flex_res_extracted}, {flex_res_name}; "
                    f"set all_states, off; "
                    f"frame {pymol_state}; "
                    f"hide everything, {ligand_name}; "
                    f"show sticks, {ligand_name}; "
                    f"show sticks, {flex_res_name}; "
                    f"zoom"
                )
            ]
            
            print(f"Opening PyMOL GUI for FLEXIBLE docking result")
            print(f"Ligand {ligand_index} (SMILES: {self._ligand_smiles[ligand_index]})")
            print(f"Binding mode: {pose_index} (PyMOL state {pymol_state})")
            print(f"Protein (flex residues removed): {protein_to_load}")
            print(f"Ligand output: {output_pdbqt}")
            print(f"Flexible residues (extracted): {flex_res_extracted}")
        else:
            # Rigid docking: use original protein
            # Use 'frame' to switch to specific state and 'set all_states, off' to hide other states
            pymol_cmd = [
                "pymol",
                "-d",
                (
                    f"set orthoscopic, on; "
                    f"load {self._protein_pdbqt}, protein; "
                    f"load {output_pdbqt}, {ligand_name}; "
                    f"set all_states, off; "
                    f"frame {pymol_state}; "
                    f"hide everything, {ligand_name}; "
                    f"show sticks, {ligand_name}; "
                    f"zoom"
                )
            ]
            
            print(f"Opening PyMOL GUI for RIGID docking result")
            print(f"Ligand {ligand_index} (SMILES: {self._ligand_smiles[ligand_index]})")
            print(f"Binding mode: {pose_index} (PyMOL state {pymol_state})")
            print(f"Protein: {self._protein_pdbqt}")
            print(f"Ligand output: {output_pdbqt}")
        
        # Run PyMOL in GUI mode
        try:
            subprocess.run(pymol_cmd, check=True)
            
        except FileNotFoundError:
            raise FileNotFoundError(
                "PyMOL executable not found. Please ensure PyMOL is installed and in your PATH."
            )

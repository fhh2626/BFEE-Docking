# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Haohao Fu <fhh2626@nankai.edu.cn, fhh2626@gmail.com>

# Standard library imports
import collections.abc
import concurrent.futures
import os
import pathlib
import subprocess
import tempfile

# Third-party library imports
import dimorphite_dl
import rdkit.Chem
import rdkit.Chem.AllChem

import third_party_tools


class DrugFileParser:
    def __init__(self, input_data: list[str], output_dir: str = "."):
        """
        Initialize the DrugFileParser.
        
        Args:
            input_data (list[str]): List of file paths. Each file can be pdb, mol2, sdf, or smi format.
                                   Each file can contain one or more molecules.
            output_dir (str): Directory where all output files will be saved. Defaults to current directory.
        """
        self.input_files = input_data if isinstance(input_data, list) else [input_data]
        self.output_dir = pathlib.Path(output_dir).resolve()
        
        # Store molecule data with source information
        # Each entry: {'mol': rdkit.Chem.Mol, 'file': str, 'file_index': int, 'mol_index': int}
        self._molecules = []
        
        # Store processed data for each molecule
        # Each entry: {'smiles': str, 'protonated_smiles': list[str], 'pdbqt_files': list[pathlib.Path]}
        self._processed_data = []
        
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Get absolute path to obabel executable
        self._obabel_executable = third_party_tools.get_obabel_executable()
        
        # Load all molecules from all files
        self._load_all_inputs()
        
        if not self._molecules:
            raise ValueError("No valid molecules found in any of the input files.")
        
        print(f"Successfully loaded {len(self._molecules)} molecule(s) from {len(self.input_files)} file(s).")

    @staticmethod
    def _convert_task_to_pdbqt(task: dict) -> dict:
        """
        Worker function to convert a task (SMILES or Mol) to PDBQT file.
        Defined as static method to support multiprocessing serialization.
        
        Args:
            task (dict): Dictionary containing:
                - 'smiles': The SMILES string to convert (optional if 'mol' is provided)
                - 'mol': RDKit Mol object (optional if 'smiles' is provided)
                - 'mol_idx': Molecule index
                - 'prot_idx': Protonation state index
                - 'base_name': Base name for output file
                - 'output_dir': Output directory path
                - 'obabel_path': Path to obabel executable
                - 'task_id': Unique task identifier
        
        Returns:
            dict: Result dictionary containing:
                - 'success': Whether conversion was successful
                - 'pdbqt_path': Path to generated PDBQT file (if successful)
                - 'smiles': The SMILES string (generated if not provided)
                - 'task_id': The task identifier
                - 'error': Error message (if failed)
        """
        mol_idx = task['mol_idx']
        prot_idx = task['prot_idx']
        base_name = task['base_name']
        output_dir = pathlib.Path(task['output_dir'])
        obabel_path = task['obabel_path']
        task_id = task['task_id']
        
        smiles = task.get('smiles')
        
        result = {
            'success': False,
            'pdbqt_path': None,
            'smiles': smiles,
            'task_id': task_id,
            'error': None
        }
        
        try:
            def generate_3d_and_optimize(mol_obj):
                try:
                    mol_obj = rdkit.Chem.AddHs(mol_obj)
                    
                    # Try 1: Default EmbedMolecule with randomSeed
                    embed_res = rdkit.Chem.AllChem.EmbedMolecule(mol_obj, randomSeed=42)
                    
                    # Try 2: If failed, try with useRandomCoords=True
                    if embed_res != 0:
                        embed_res = rdkit.Chem.AllChem.EmbedMolecule(mol_obj, randomSeed=42, useRandomCoords=True)
                    
                    # Try 3: If still failed, try with ETKDGv3
                    if embed_res != 0:
                        params = rdkit.Chem.AllChem.ETKDGv3()
                        params.randomSeed = 42
                        embed_res = rdkit.Chem.AllChem.EmbedMolecule(mol_obj, params)
                    
                    if embed_res != 0:
                        return False, "Could not generate 3D coordinates"
                    
                    rdkit.Chem.AllChem.MMFFOptimizeMolecule(mol_obj)
                    return True, mol_obj
                except Exception as ex:
                    return False, str(ex)

            mol = None
            if 'mol' in task:
                mol = task['mol']
                # If naive pipeline, check if we need to generate 3D coords
                if mol.GetNumConformers() == 0:
                    success, res = generate_3d_and_optimize(mol)
                    if not success:
                        result['error'] = f"Molecule {mol_idx}: {res}"
                        return result
                    mol = res
                
                if not smiles:
                    try:
                        result['smiles'] = rdkit.Chem.MolToSmiles(mol)
                    except:
                        result['smiles'] = ""
            
            elif smiles:
                # Convert SMILES to RDKit molecule
                mol = rdkit.Chem.MolFromSmiles(smiles)
                if mol is None:
                    result['error'] = f"Could not parse SMILES: {smiles}"
                    return result
                
                success, res = generate_3d_and_optimize(mol)
                if not success:
                    result['error'] = f"{smiles}: {res}"
                    return result
                mol = res
            
            else:
                result['error'] = "No molecule or SMILES provided in task"
                return result

            # Generate output filename
            output_name = f"{base_name}_mol{mol_idx}_prot{prot_idx}"
            output_pdbqt = output_dir / f"{output_name}.pdbqt"
            
            # Write PDB to a temporary file, then convert to PDBQT
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_pdb:
                tmp_pdb_path = tmp_pdb.name
                rdkit.Chem.MolToPDBFile(mol, tmp_pdb_path)
            
            try:
                # Convert PDB to PDBQT using obabel
                subprocess.run(
                    [obabel_path, tmp_pdb_path, "-O", str(output_pdbqt), "-xh"],
                    check=True,
                    capture_output=True,
                    text=True
                )
                result['success'] = True
                result['pdbqt_path'] = output_pdbqt
            except subprocess.CalledProcessError as e:
                result['error'] = f"obabel error: {e.stderr}"
            except FileNotFoundError:
                result['error'] = f"obabel not found at {obabel_path}"
            finally:
                # Clean up temporary PDB file
                if os.path.exists(tmp_pdb_path):
                    os.unlink(tmp_pdb_path)
        
        except Exception as e:
            result['error'] = str(e)
        
        return result

    def _load_all_inputs(self):
        """
        Internal method to load molecules from all input files.
        """
        for file_index, file_path in enumerate(self.input_files):
            if not os.path.isfile(file_path):
                print(f"Warning: File not found: {file_path}")
                continue
                
            mols = self._load_molecules_from_file(file_path)
            
            for mol_index, mol in enumerate(mols):
                self._molecules.append({
                    'mol': mol,
                    'file': file_path,
                    'file_index': file_index,
                    'mol_index': mol_index
                })
            
            print(f"Loaded {len(mols)} molecule(s) from {file_path}")

    def _load_molecules_from_file(self, file_path: str) -> list:
        """
        Load molecules from a single file.
        
        Args:
            file_path (str): Path to the file to load.
            
        Returns:
            list: List of RDKit molecule objects.
        """
        mols = []
        ext = os.path.splitext(file_path)[1].lower()
        
        try:
            if ext == '.pdb':
                mol = rdkit.Chem.MolFromPDBFile(file_path, removeHs=False)
                if mol:
                    mols.append(mol)
            
            elif ext == '.sdf':
                suppl = rdkit.Chem.SDMolSupplier(file_path, removeHs=False)
                for mol in suppl:
                    if mol:
                        mols.append(mol)
            
            elif ext == '.mol2':
                mol = rdkit.Chem.MolFromMol2File(file_path, removeHs=False)
                if mol:
                    mols.append(mol)
            
            elif ext == '.smi':
                # Read SMI file line by line
                with open(file_path, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            # SMI files can have format: "SMILES name" or just "SMILES"
                            smiles = line.split()[0]
                            mol = rdkit.Chem.MolFromSmiles(smiles)
                            if mol:
                                mols.append(mol)
            
            else:
                print(f"Warning: Unsupported file extension '{ext}' for file {file_path}")
                
        except Exception as e:
            print(f"Error loading file {file_path}: {e}")
        
        return mols
    
    def _add_hydrogens(
        self, 
        ph: float, 
        progress_callback: collections.abc.Callable[[int, int, str], None] | None = None
    ) -> None:
        """
        Internal method to add hydrogens to all molecules based on the specified pH value.
        Uses dimorphite_dl to generate protonation states.
        
        Args:
            ph (float): The pH value for protonation.
            progress_callback (callable, optional): A callback function to report progress.
                The callback will be called with (current, total, description) where:
                - current (int): Current molecule index being processed
                - total (int): Total number of molecules
                - description (str): Description of current step
        """
        print(f"\nProcessing {len(self._molecules)} molecules at pH {ph}...")
        total_molecules = len(self._molecules)
        
        for idx, mol_data in enumerate(self._molecules):
            mol = mol_data['mol']
            file_path = mol_data['file']
            mol_index = mol_data['mol_index']
            
            # Report progress
            if progress_callback:
                progress_callback(idx, total_molecules, f"Adding hydrogens ({idx+1}/{total_molecules})...")
            
            # Convert molecule to SMILES
            try:
                smiles = rdkit.Chem.MolToSmiles(mol)
            except Exception as e:
                print(f"Warning: Could not convert molecule {mol_index} from {file_path} to SMILES: {e}")
                self._processed_data.append({
                    'smiles': None,
                    'protonated_smiles': [],
                    'pdbqt_files': [],
                    'generated_smiles': []
                })
                continue
            
            # Use dimorphite_dl to protonate the molecule
            ph_range = 0.2
            try:
                protonated_smiles = dimorphite_dl.protonate_smiles(
                    smiles,
                    ph_min=ph - ph_range,
                    ph_max=ph + ph_range
                )
            except Exception as e:
                print(f"Warning: Could not protonate molecule {mol_index} from {file_path}: {e}")
                protonated_smiles = []
            
            self._processed_data.append({
                'smiles': smiles,
                'protonated_smiles': protonated_smiles,
                'pdbqt_files': [],
                'generated_smiles': []
            })
            
            print(f"  [{idx+1}/{len(self._molecules)}] File: {os.path.basename(file_path)}, "
                  f"Molecule: {mol_index}, Protonation states: {len(protonated_smiles)}")
        
        # Report completion
        if progress_callback:
            progress_callback(total_molecules, total_molecules, f"All hydrogens added ({total_molecules}/{total_molecules})")
    
    def _execute_conversion_tasks(
        self,
        tasks: list[dict],
        max_workers: int | None,
        progress_callback: collections.abc.Callable[[int, int, str], None] | None
    ) -> list[dict]:
        """
        Helper method to execute conversion tasks in parallel.
        
        Args:
            tasks (list[dict]): List of task dictionaries.
            max_workers (int, optional): Max workers.
            progress_callback (callable, optional): Progress callback.
            
        Returns:
            list[dict]: List of result dictionaries.
        """
        total_tasks = len(tasks)
        if total_tasks == 0:
            print("No tasks to convert.")
            if progress_callback:
                progress_callback(0, 0, "No molecules to convert")
            return []
            
        print(f"\nConverting {total_tasks} molecules/states to PDBQT format using multiprocessing...")
        
        results = []
        completed_count = 0
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {executor.submit(DrugFileParser._convert_task_to_pdbqt, task): task 
                              for task in tasks}
            
            for future in concurrent.futures.as_completed(future_to_task):
                completed_count += 1
                try:
                    res = future.result()
                    results.append(res)
                    if not res['success']:
                        print(f"    Warning: {res['error']}")
                except Exception as e:
                    print(f"    Error processing task: {e}")
                
                if progress_callback:
                    progress_callback(
                        completed_count, 
                        total_tasks, 
                        f"Converting to PDBQT ({completed_count}/{total_tasks})..."
                    )
        return results

    def _process_conversion_results(self, results: list[dict], task_to_mol_map: dict) -> list[pathlib.Path]:
        """
        Helper to process conversion results and update internal state.
        Sorts results by protonation index and updates _processed_data.
        
        Args:
            results: List of result dicts from _execute_conversion_tasks
            task_to_mol_map: Dict mapping task_id to (mol_idx, prot_idx) or just mol_idx
            
        Returns:
            List of successfully generated pdbqt paths
        """
        temp_results = {idx: [] for idx in range(len(self._molecules))}
        all_pdbqt_files = []
        
        for result in results:
            if result['success']:
                idx_info = task_to_mol_map[result['task_id']]
                mol_idx = idx_info[0] if isinstance(idx_info, tuple) else idx_info
                prot_idx = idx_info[1] if isinstance(idx_info, tuple) else 0
                
                temp_results[mol_idx].append({
                    'prot_idx': prot_idx,
                    'pdbqt_path': result['pdbqt_path'],
                    'smiles': result['smiles']
                })
                all_pdbqt_files.append(result['pdbqt_path'])
        
        for idx in range(len(self._molecules)):
            items = temp_results[idx]
            # Ensure deterministic order
            items.sort(key=lambda x: x['prot_idx'])
            
            pdbqt_files = [x['pdbqt_path'] for x in items]
            generated_smiles = [x['smiles'] for x in items]
            
            self._processed_data[idx]['pdbqt_files'] = pdbqt_files
            self._processed_data[idx]['generated_smiles'] = generated_smiles
            
            if pdbqt_files:
                mol_data = self._molecules[idx]
                print(f"  Molecule {mol_data['mol_index']} from {os.path.basename(mol_data['file'])}: "
                      f"Generated {len(pdbqt_files)} PDBQT file(s)")
        
        print(f"\nTotal: Generated {len(all_pdbqt_files)} PDBQT file(s)")
        return all_pdbqt_files

    def _convert_to_pdbqt(
        self, 
        progress_callback: collections.abc.Callable[[int, int, str], None] | None = None,
        max_workers: int | None = None
    ) -> list[pathlib.Path]:
        """
        Internal method to convert all protonated SMILES to PDBQT files using RDKit and Open Babel.
        """
        if not self._processed_data:
            raise ValueError("No processed data available. Please run _add_hydrogens() first.")
        
        tasks = []
        task_id = 0
        task_to_mol_map = {}
        
        for idx, (mol_data, proc_data) in enumerate(zip(self._molecules, self._processed_data)):
            file_path = mol_data['file']
            mol_index = mol_data['mol_index']
            protonated_smiles = proc_data['protonated_smiles']
            
            if not protonated_smiles:
                print(f"  Skipping molecule {mol_index} from {os.path.basename(file_path)} (no protonated SMILES)")
                continue
            
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            
            for prot_idx, smiles in enumerate(protonated_smiles):
                task = {
                    'smiles': smiles,
                    'mol_idx': mol_index,
                    'prot_idx': prot_idx,
                    'base_name': base_name,
                    'output_dir': str(self.output_dir),
                    'obabel_path': str(self._obabel_executable),
                    'task_id': task_id
                }
                tasks.append(task)
                task_to_mol_map[task_id] = (idx, prot_idx)
                task_id += 1
        
        results = self._execute_conversion_tasks(tasks, max_workers, progress_callback)
        return self._process_conversion_results(results, task_to_mol_map)
    
    def process_ligand_naive_pipeline(
        self,
        progress_callback: collections.abc.Callable[[int, int, str], None] | None = None,
        max_workers: int | None = None
    ) -> list[pathlib.Path]:
        """
        Run the naive ligand processing pipeline:
        1. Directly convert existing internal molecules to PDBQT files
        2. Skip pH adjustment and dimorphite_dl steps
        
        Args:
            progress_callback (callable, optional): A callback function to report progress.
            max_workers (int, optional): Maximum number of worker processes.
        
        Returns:
            list[pathlib.Path]: List of all generated PDBQT file paths.
        """
        # Clear/Init processed data
        self._processed_data = [{
            'smiles': None, 
            'protonated_smiles': [], 
            'pdbqt_files': [], 
            'generated_smiles': []
        } for _ in self._molecules]
        
        tasks = []
        task_id = 0
        task_to_mol_map = {}
        
        for idx, mol_data in enumerate(self._molecules):
            file_path = mol_data['file']
            mol_index = mol_data['mol_index']
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            
            # For naive pipeline, we just use the mol object directly
            task = {
                'mol': mol_data['mol'],
                'mol_idx': mol_index,
                'prot_idx': 0, # Only one state per molecule in naive mode
                'base_name': base_name,
                'output_dir': str(self.output_dir),
                'obabel_path': str(self._obabel_executable),
                'task_id': task_id
            }
            tasks.append(task)
            task_to_mol_map[task_id] = (idx, 0)
            task_id += 1
            
        results = self._execute_conversion_tasks(tasks, max_workers, progress_callback)
        return self._process_conversion_results(results, task_to_mol_map)

    def process_ligand_pipeline(
        self,
        ph: float,
        progress_callback: collections.abc.Callable[[int, int, str], None] | None = None,
        max_workers: int | None = None
    ) -> list[pathlib.Path]:
        """
        Run the complete ligand processing pipeline:
        1. Add hydrogens / protonate molecules (using pH)
        2. Convert to PDBQT files
        
        Args:
            ph (float): The pH value for protonation.
            progress_callback (callable, optional): A callback function to report progress.
                The callback will be called with (current_step, total_steps, description).
                Total steps = number of molecules + number of protonation states.
            max_workers (int, optional): Maximum number of worker processes for PDBQT conversion.
                Defaults to None (uses number of CPUs).
        
        Returns:
            list[pathlib.Path]: List of all generated PDBQT file paths.
        """
        total_molecules = len(self._molecules)
        
        # Step 1: Add hydrogens
        # During this phase, we don't know the total protonation states yet,
        # so we report progress as molecules completed out of molecules * 2 (estimate)
        def add_h_callback(current, total, description):
            if progress_callback:
                # Use estimated total (2x molecules) during add_hydrogens phase
                progress_callback(current, total_molecules * 2, description)
        
        self._add_hydrogens(ph, progress_callback=add_h_callback)
        
        # After add_hydrogens, calculate actual total protonation states
        total_protonation_states = sum(
            len(proc_data['protonated_smiles']) 
            for proc_data in self._processed_data
        )
        
        # Recalculate total_steps: molecules + protonation_states
        total_steps = total_molecules + total_protonation_states
        
        # Step 2: Convert to PDBQT
        def convert_callback(current, total, description):
            if progress_callback:
                # current is the number of protonation states completed
                progress_callback(total_molecules + current, total_steps, description)
        
        return self._convert_to_pdbqt(progress_callback=convert_callback, max_workers=max_workers)

    def get_generated_pdbqt_paths(self) -> list[pathlib.Path]:
        """
        Get the list of all generated PDBQT file paths.
        
        Returns:
            list[pathlib.Path]: List of all PDBQT file paths (as Path objects) that were generated by _convert_to_pdbqt().
                                Returns an empty list if _convert_to_pdbqt() has not been called yet.
        """
        all_pdbqt_files = []
        for proc_data in self._processed_data:
            all_pdbqt_files.extend(proc_data['pdbqt_files'])
        return all_pdbqt_files
    
    def get_generated_smiles(self) -> list[str]:
        """
        Get the list of all generated protonated SMILES strings that have corresponding PDBQT files.
        
        Returns:
            list[str]: List of all protonated SMILES strings that were generated by _add_hydrogens()
                      AND successfully converted to PDBQT.
                      Returns an empty list if conversion hasn't been run.
        """
        all_smiles = []
        for proc_data in self._processed_data:
            all_smiles.extend(proc_data['generated_smiles'])
        return all_smiles
    
    def get_molecules(self) -> list[dict]:
        """
        Get the list of loaded molecules with their metadata.
        
        Returns:
            list[dict]: List of molecule dictionaries. Each dictionary contains:
                - 'mol': rdkit.Chem.Mol object
                - 'file': source file path
                - 'file_index': index of the source file
                - 'mol_index': index of the molecule within the file
        """
        return self._molecules
    
    def get_processed_data(self) -> list[dict]:
        """
        Get the list of processed data for each molecule.
        
        Returns:
            list[dict]: List of processed data dictionaries. Each dictionary contains:
                - 'smiles': original SMILES string
                - 'protonated_smiles': list of protonated SMILES strings
                - 'pdbqt_files': list of generated PDBQT file paths
                Returns an empty list if _add_hydrogens() has not been called yet.
        """
        return self._processed_data

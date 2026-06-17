# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Haohao Fu <fhh2626@nankai.edu.cn, fhh2626@gmail.com>

"""
Third-party tools path wrapper module.

This module provides functions to get paths to bundled third-party executables
and force field files. All paths are resolved relative to this file's location,
ensuring consistent behavior regardless of the current working directory.

It also includes a wrapper for pdb2pqr with monkey patches for:
1. propka: Fix the '__annotations__' AttributeError bug
2. pdb2pqr.ligand.mol2: Fix the 'formal_charge' ValueError/KeyError bugs
3. pdb2pqr.ligand.peoe: Add compatibility aliases for unsupported ligand atom types
"""

# ========== Standard Library Imports ==========
import os
import pathlib
import shutil
import subprocess
import sys

# ========== Third-Party Imports ==========
# NOTE: Import order matters here! Monkey patches must be applied before
# dependent modules are imported.

# 1. Import propka.parameters first for monkey patching
import propka.parameters as _propka_parameters


# ========== Propka Monkey Patch ==========
# Manually define the type mappings that were supposed to be in __annotations__
# These are extracted from reading the source code of propka/parameters.py
_PROPKA_NUMBER_DICTS = {
    'VanDerWaalsVolume', 'charge', 'model_pkas', 'ions', 
    'valence_electrons', 'custom_model_pkas'
}
_PROPKA_LIST_DICTS = {'backbone_NH_hydrogen_bond', 'backbone_CO_hydrogen_bond'}
_PROPKA_STRING_DICTS = {'protein_group_mapping'}
_PROPKA_STRING_LISTS = {
    'ignore_residues', 'angular_dependent_sidechain_interactions', 'acid_list', 
    'base_list', 'exclude_sidechain_interactions', 'backbone_reorganisation_list', 
    'write_out_order'
}
_PROPKA_MATRICES = {'interaction_matrix', 'sidechain_cutoffs'}
_PROPKA_STRINGS = {'version', 'output_file_tag', 'ligand_typing', 'pH', 'reference'}


def _patched_propka_parse_line(self, line):
    """
    Monkey patched parse_line to avoid usage of self.__annotations__.
    This fixes AttributeError: 'Parameters' object has no attribute '__annotations__'
    """
    # first, remove comments
    comment_pos = line.find('#')
    if comment_pos != -1:
        line = line[:comment_pos]
    # split the line into words
    words = line.split()
    if len(words) == 0:
        return
    
    # Logic to determine type based on manually defined sets
    key = words[0]
    
    if key in _PROPKA_NUMBER_DICTS:
        self.parse_to_number_dictionary(words)
    elif key in _PROPKA_STRING_LISTS:
        self.parse_to_string_list(words)
    elif key in _PROPKA_STRINGS:
        self.parse_string(words)
    elif key in _PROPKA_LIST_DICTS:
        self.parse_to_list_dictionary(words)
    elif key in _PROPKA_MATRICES:
        self.parse_to_matrix(words)
    elif key in _PROPKA_STRING_DICTS:
        self.parse_to_string_dictionary(words)
    else:
        # Fallback to parameter (simple float/int)
        self.parse_parameter(words)


# Apply the propka monkey patch
_propka_parameters.Parameters.parse_line = _patched_propka_parse_line


# ========== PDB2PQR MOL2 Monkey Patch ==========
# We need to patch the formal_charge property of the Atom class in pdb2pqr.ligand.mol2
# This is done at runtime in _apply_pdb2pqr_mol2_patch() because the Atom class
# is defined inside the module and not directly exported.

_pdb2pqr_mol2_patched = False


def _apply_pdb2pqr_mol2_patch():
    """
    Apply monkey patches for pdb2pqr.ligand.mol2 at runtime.
    
    This fixes the ValueError bug where the code assumes certain oxygen atoms
    must be bonded to phosphorus, causing ValueError when P is not found.
    It also adds missing Sybyl atom type aliases to pdb2pqr's nonbonded
    electron table. Open Babel may write carbocation/guanidinium carbons as
    C.cat, but pdb2pqr 3.6.1 omits that key even though its PEOE table supports
    the corresponding C.CAT type.
    """
    global _pdb2pqr_mol2_patched
    if _pdb2pqr_mol2_patched:
        return
    
    import pdb2pqr.ligand as ligand_module
    import pdb2pqr.ligand.mol2 as mol2_module

    if "C.cat" not in ligand_module.NONBONDED_BY_TYPE:
        ligand_module.NONBONDED_BY_TYPE["C.cat"] = ligand_module.NONBONDED_BY_TYPE["C.3"]
    
    # Get the Atom class from the module's namespace
    Atom = getattr(mol2_module, 'Atom', None)
    if Atom is None:
        # Try to find it in module globals
        for name, obj in vars(mol2_module).items():
            if isinstance(obj, type) and 'Atom' in name:
                Atom = obj
                break
    
    if Atom is None:
        print("Warning: Could not find Atom class in pdb2pqr.ligand.mol2 for patching")
        return
    
    # Store the original formal_charge property
    original_formal_charge = Atom.formal_charge
    
    @property
    def patched_formal_charge(self):
        """Return formal charge, with error handling for the P-atom bug."""
        try:
            # Try the original logic
            return original_formal_charge.fget(self)
        except (ValueError, IndexError):
            # If the P-atom lookup fails, return 0 as default charge.
            return 0
    
    # Apply the patch
    Atom.formal_charge = patched_formal_charge
    _pdb2pqr_mol2_patched = True


_pdb2pqr_peoe_patched = False


def _apply_pdb2pqr_peoe_patch():
    """
    Apply monkey patch to pdb2pqr.ligand.peoe.POLY_TERMS at runtime.

    pdb2pqr accepts sulfur-oxide MOL2 atom types in its ligand parser but does
    not provide matching PEOE polynomial terms for S.O. Add a compatibility
    alias so ligand charge assignment can proceed for sulfoxide-like atoms.
    """
    global _pdb2pqr_peoe_patched
    if _pdb2pqr_peoe_patched:
        return

    import pdb2pqr.ligand.peoe as peoe_module

    if "S.O" not in peoe_module.POLY_TERMS:
        peoe_module.POLY_TERMS["S.O"] = peoe_module.POLY_TERMS["S.O2"]

    _pdb2pqr_peoe_patched = True


# ========== PDB2PQR Main Import (After All Patches Applied) ==========
import pdb2pqr.main as _pdb2pqr_main


# ========== Module Constants ==========
# Get the directory containing this file (used as base for all third-party paths)
_SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
_THIRD_PARTY_DIR = _SCRIPT_DIR / "third_party"


def get_third_party_dir() -> pathlib.Path:
    """
    Get the path to the third_party directory.
    
    Returns:
        pathlib.Path: Absolute path to the third_party directory.
    """
    return _THIRD_PARTY_DIR


def add_third_party_to_path() -> None:
    """
    Add the third_party directory to Python's sys.path.
    
    This allows importing bundled Python packages (e.g., pocketeer)
    from the third_party directory.
    """
    third_party_str = str(_THIRD_PARTY_DIR)
    if third_party_str not in sys.path:
        sys.path.insert(0, third_party_str)


# ========== Executable Paths ==========

def get_obabel_executable() -> pathlib.Path:
    """
    Get the path to the Open Babel (obabel) executable.
    
    Returns:
        pathlib.Path: Absolute path to the obabel executable.
                      On Windows: third_party/obabel/obabel.exe
                      On non-Windows: prefer obabel from the current environment,
                      then fall back to third_party/obabel/obabel
    """
    if os.name != 'nt':
        system_obabel = shutil.which("obabel")
        if system_obabel:
            return pathlib.Path(system_obabel).resolve()

    obabel_name = "obabel.exe" if os.name == 'nt' else "obabel"
    return _THIRD_PARTY_DIR / "obabel" / obabel_name


def get_vmd_executable(vmd_path: str) -> pathlib.Path:
    """
    Get the path to the VMD executable.
    
    Args:
        vmd_path: Path to the VMD executable provided by the user.
    
    Returns:
        pathlib.Path: Absolute path to the VMD executable.
    
    Raises:
        FileNotFoundError: If the VMD executable doesn't exist at the given path.
    """
    vmd_executable = pathlib.Path(vmd_path).resolve()
    if not vmd_executable.exists():
        raise FileNotFoundError(
            f"VMD executable not found at: {vmd_executable}\n"
            f"Please provide a valid path to the VMD executable."
        )
    return vmd_executable


def get_vina_executable() -> pathlib.Path:
    """
    Get the path to the AutoDock Vina executable.
    
    Returns:
        pathlib.Path: Absolute path to the vina executable.
                      On Windows: third_party/vina/vina.exe
                      On Unix: third_party/vina/vina
    """
    vina_name = "vina.exe" if os.name == 'nt' else "vina"
    return _THIRD_PARTY_DIR / "vina" / vina_name


def get_smina_executable() -> pathlib.Path:
    """
    Get the path to the smina executable.
    
    Returns:
        pathlib.Path: Absolute path to the smina executable.
                      On Windows: third_party/smina/smina.exe
                      On Unix: third_party/smina/smina
    """
    smina_name = "smina.exe" if os.name == 'nt' else "smina"
    return _THIRD_PARTY_DIR / "smina" / smina_name


def _get_bundled_or_system_executable(
    bundled_dir: pathlib.Path,
    executable_name: str,
    command_name: str,
) -> pathlib.Path:
    """
    Prefer a bundled executable, then fall back to the current environment.
    """
    executable = bundled_dir / executable_name
    if executable.exists():
        return executable

    system_executable = shutil.which(command_name)
    if system_executable:
        return pathlib.Path(system_executable).resolve()

    return pathlib.Path(command_name)


def get_qvina_executable(engine: str) -> pathlib.Path:
    """
    Get the path to a qvina-family executable.

    Args:
        engine: One of "qvina2", "qvinaw", or "vina-classic".

    Returns:
        pathlib.Path: Bundled executable path if available, otherwise a PATH
                      executable path or command name.
    """
    executable_stems = {
        "qvina2": "qvina2",
        "qvinaw": "qvinaw",
        "vina-classic": "vina",
    }
    command_names = {
        "qvina2": "qvina2",
        "qvinaw": "qvinaw",
        "vina-classic": "vina",
    }

    stem = executable_stems[engine]
    executable_name = f"{stem}.exe" if os.name == 'nt' else stem
    return _get_bundled_or_system_executable(
        _THIRD_PARTY_DIR / "qvina",
        executable_name,
        command_names[engine],
    )


def get_docking_executable(engine: str) -> pathlib.Path:
    """
    Get the path to a docking engine executable.
    
    Args:
        engine: Name of the docking engine. Supported values:
                - "vina-new": AutoDock Vina (bundled)
                - "smina": smina (bundled)
                - "vina-classic", "qvina2", "qvinaw": Prefer bundled qvina,
                  then fall back to the current environment

    Returns:
        pathlib.Path: Absolute path to the docking executable.
    
    Raises:
        ValueError: If the engine is not recognized.
        FileNotFoundError: If the bundled executable doesn't exist.
    """
    engine_map = {
        "vina-new": get_vina_executable,
        "smina": get_smina_executable,
    }
    
    if engine in engine_map:
        executable = engine_map[engine]()
        if not executable.exists():
            raise FileNotFoundError(
                f"Docking executable not found at: {executable}\n"
                f"Please ensure the {engine} executable is installed in third_party/"
            )
        return executable
    elif engine in ["vina-classic", "qvina2", "qvinaw"]:
        return get_qvina_executable(engine)
    else:
        raise ValueError(f"Unknown docking engine: {engine}")


def get_docking_subprocess_env(executable: str | pathlib.Path) -> dict[str, str] | None:
    """
    Get subprocess environment overrides for bundled docking executables.

    Bundled qvina binaries ship with same-directory DLL/SO dependencies. When
    such a binary is used, expose its directory to the platform dynamic loader.
    """
    executable_path = pathlib.Path(executable)
    qvina_dir = _THIRD_PARTY_DIR / "qvina"

    try:
        resolved_executable = executable_path.resolve()
        resolved_qvina_dir = qvina_dir.resolve()
    except OSError:
        return None

    if not executable_path.exists() or resolved_executable.parent != resolved_qvina_dir:
        return None

    env = os.environ.copy()
    qvina_dir_str = str(resolved_qvina_dir)

    if os.name == 'nt':
        current_path = env.get("PATH", "")
        env["PATH"] = qvina_dir_str if not current_path else f"{qvina_dir_str}{os.pathsep}{current_path}"
    else:
        current_ld_library_path = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = (
            qvina_dir_str
            if not current_ld_library_path
            else f"{qvina_dir_str}{os.pathsep}{current_ld_library_path}"
        )

    return env


# ========== Force Field Files ==========

def get_protein_topology_file() -> pathlib.Path:
    """
    Get the path to the CHARMM protein topology file.
    
    Returns:
        pathlib.Path: Absolute path to top_all36_prot.rtf
    """
    return _THIRD_PARTY_DIR / "force_field" / "top_all36_prot.rtf"


def get_cgenff_topology_file() -> pathlib.Path:
    """
    Get the path to the CGenFF (ligand) topology file.
    
    Returns:
        pathlib.Path: Absolute path to top_all36_cgenff.rtf
    """
    return _THIRD_PARTY_DIR / "force_field" / "top_all36_cgenff.rtf"


def get_force_field_dir() -> pathlib.Path:
    """
    Get the path to the force field directory.
    
    Returns:
        pathlib.Path: Absolute path to the force_field directory.
    """
    return _THIRD_PARTY_DIR / "force_field"


# ========== Validation Functions ==========

def validate_obabel_exists() -> bool:
    """
    Check if the obabel executable exists.
    
    Returns:
        bool: True if obabel exists, False otherwise.
    """
    return get_obabel_executable().exists()





def validate_docking_engine_exists(engine: str) -> bool:
    """
    Check if a docking engine executable exists.
    
    Args:
        engine: Name of the docking engine.
    
    Returns:
        bool: True if the executable exists, False otherwise.
    """
    try:
        executable = get_docking_executable(engine)
        return executable.exists()
    except (ValueError, FileNotFoundError):
        return False


def validate_force_field_files_exist() -> bool:
    """
    Check if the required force field files exist.
    
    Returns:
        bool: True if all required force field files exist, False otherwise.
    """
    return (
        get_protein_topology_file().exists() and
        get_cgenff_topology_file().exists()
    )


# ========== Utility Functions ==========

def get_all_tool_paths() -> dict:
    """
    Get all tool paths as a dictionary for debugging.
    
    Returns:
        dict: Dictionary containing all tool paths and their existence status.
    """
    return {
        "third_party_dir": str(_THIRD_PARTY_DIR),
        "obabel": {
            "path": str(get_obabel_executable()),
            "exists": validate_obabel_exists()
        },

        "vina": {
            "path": str(get_vina_executable()),
            "exists": get_vina_executable().exists()
        },
        "smina": {
            "path": str(get_smina_executable()),
            "exists": get_smina_executable().exists()
        },
        "qvina2": {
            "path": str(get_docking_executable("qvina2")),
            "exists": get_docking_executable("qvina2").exists()
        },
        "qvinaw": {
            "path": str(get_docking_executable("qvinaw")),
            "exists": get_docking_executable("qvinaw").exists()
        },
        "vina_classic": {
            "path": str(get_docking_executable("vina-classic")),
            "exists": get_docking_executable("vina-classic").exists()
        },
        "protein_topology": {
            "path": str(get_protein_topology_file()),
            "exists": get_protein_topology_file().exists()
        },
        "cgenff_topology": {
            "path": str(get_cgenff_topology_file()),
            "exists": get_cgenff_topology_file().exists()
        }
    }


# ========== PDB2PQR Wrapper ==========

class PDB2PQRError(RuntimeError):
    """Raised when pdb2pqr fails during processing."""


def run_pdb2pqr(
    input_pdb: str,
    output_pqr: str,
    ph: float = 7.0,
    ff: str = "CHARMM",
    titration_state_method: str = "propka",
    keep_chain: bool = True,
    ligand_mol2: str | None = None
) -> None:
    """
    Run pdb2pqr with the propka monkey patch applied.
    
    This function runs pdb2pqr programmatically with the necessary monkey patch
    for propka already applied (to fix the __annotations__ AttributeError).
    
    Args:
        input_pdb: Path to input PDB file.
        output_pqr: Path to output PQR file.
        ph: pH value for protonation (default: 7.0).
        ff: Force field to use (default: "CHARMM"). 
            Options: "CHARMM", "AMBER", "CHARMM_WITH_ION".
            CHARMM_WITH_ION uses custom force field files with metal ion support.
        titration_state_method: Method for titration state (default: "propka").
        keep_chain: Whether to keep chain IDs (default: True).
        ligand_mol2: Path to ligand MOL2 file for parameterization (default: None).
                     When provided, pdb2pqr will use --ligand to generate 
                     parameters for the small molecule.
    
    Raises:
        PDB2PQRError: If pdb2pqr fails.
        ValueError: If ff is not one of the allowed values.
    
    Note:
        The propka monkey patch is applied when this module is imported,
        so this function can safely use pdb2pqr without the __annotations__ bug.
        Hydrogen optimization is enabled by default (noopt=False).
    """
    # Apply ligand patches only when pdb2pqr will parameterize a small molecule.
    if ligand_mol2:
        _apply_pdb2pqr_mol2_patch()
        _apply_pdb2pqr_peoe_patch()
    
    # Validate ff parameter
    allowed_ff = {"CHARMM", "AMBER", "CHARMM_WITH_ION"}
    ff_upper = ff.upper()
    if ff_upper not in allowed_ff:
        raise ValueError(
            f"Invalid force field '{ff}'. "
            f"Allowed values: {', '.join(sorted(allowed_ff))}"
        )
    
    # Build command line arguments as a list
    argv = [
        'pdb2pqr',
        '--titration-state-method', titration_state_method,
        f'--with-ph={ph}'
    ]
    
    # Handle force field selection
    if ff_upper == "CHARMM_WITH_ION":
        # Use custom force field files with metal ion support
        userff_path = _THIRD_PARTY_DIR / "pdb2pqr_ff" / "CHARMM.DAT"
        usernames_path = _THIRD_PARTY_DIR / "pdb2pqr_ff" / "CHARMM.names"
        argv.extend([
            '--userff', str(userff_path.resolve()),
            '--usernames', str(usernames_path.resolve())
        ])
    else:
        # Use built-in force field (CHARMM or AMBER)
        argv.append(f'--ff={ff_upper}')
    
    if keep_chain:
        argv.append('--keep-chain')
    
    # Add ligand parameter if provided
    if ligand_mol2:
        argv.extend(['--ligand', str(ligand_mol2)])
    
    argv.extend([str(input_pdb), str(output_pqr)])
    

    
    # Save original sys.argv and restore after execution
    original_argv = sys.argv
    try:
        sys.argv = argv
        # pdb2pqr.main IS the main function itself (not a module with a main function)
        _pdb2pqr_main()
    except SystemExit as e:
        # pdb2pqr calls sys.exit(0) on success, sys.exit(1) on failure
        if e.code != 0:
            raise PDB2PQRError(f"pdb2pqr failed with exit code {e.code}") from e
    except Exception as e:
        raise PDB2PQRError(f"pdb2pqr failed: {e}") from e
    finally:
        sys.argv = original_argv


def run_vmd(
    vmd_executable: pathlib.Path,
    tcl_script: pathlib.Path,
    cwd: pathlib.Path,
    dispdev: str = "text"
) -> str:
    """
    Execute a TCL script using VMD.

    Args:
        vmd_executable: Path to the VMD executable.
        tcl_script: Path to the TCL script to execute.
        cwd: Working directory for execution.
        dispdev: VMD display device mode (default: "text").

    Returns:
        The standard output from VMD.

    Raises:
        subprocess.CalledProcessError: If VMD execution fails.
    """
    cmd = [str(vmd_executable), "-dispdev", dispdev, "-e", str(tcl_script)]
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            cwd=str(cwd)
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        # Re-raise with captured output for better debugging
        # We assume the caller handles printing the error/output if needed, 
        # or we could print it here. Given the existing code style, printing here ensures visibility.
        print(f"Error running VMD: {e}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        raise

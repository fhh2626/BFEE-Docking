from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout, exit, stderr
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations as trans

# Import custom TABF classes from WTMtABF.py
from WTMtABF import TABFIntegrator, CustomVariableReporter, PMFReporter, MTDReporter, CKPTReporter

# ==============================================================================
# Custom PDB Restart Reporter (Overwrites restart.pdb every N steps)
# ==============================================================================
class PDBRestartReporter:
    """
    A custom reporter that overwrites a PDB file with the current structure
    at regular intervals. Unlike DCDReporter, this overwrites the file each time
    instead of appending to it.
    
    Uses MDAnalysis to write PDB files, preserving original residue names
    from the PSF topology (e.g., TIP3 instead of HOH).
    
    Applies proper PBC handling:
    1. Unwrap molecules to make them whole
    2. Center on protein
    3. Wrap all molecules back into the box
    """
    def __init__(self, file, reportInterval, psf_file, pdb_file, center_selection='protein'):
        self._file = file
        self._reportInterval = reportInterval
        self._center_selection = center_selection
        # Create MDAnalysis Universe from PSF/PDB to preserve CHARMM naming
        self._universe = mda.Universe(psf_file, pdb_file)
        
        # Pre-select the centering group (e.g., protein)
        self._center_group = self._universe.select_atoms(center_selection)
        if len(self._center_group) == 0:
            # Fallback to all atoms if selection is empty
            self._center_group = self._universe.atoms

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False, False)

    def report(self, simulation, state):
        # Get positions in Angstrom (MDAnalysis uses Angstrom)
        positions = state.getPositions(asNumpy=True)
        positions_angstrom = positions.value_in_unit(angstrom)
        
        # Update positions in MDAnalysis Universe
        self._universe.atoms.positions = positions_angstrom
        
        # Get box vectors and set dimensions
        box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
        box_vectors_angstrom = box_vectors.value_in_unit(angstrom)
        # MDAnalysis uses [lx, ly, lz, alpha, beta, gamma] format
        lx = box_vectors_angstrom[0, 0]
        ly = box_vectors_angstrom[1, 1]
        lz = box_vectors_angstrom[2, 2]
        self._universe.dimensions = [lx, ly, lz, 90.0, 90.0, 90.0]
        
        # Apply transformations for proper PBC handling:
        # 1. Unwrap: Make molecules whole (fix molecules split across PBC)
        trans.unwrap(self._universe.atoms)(self._universe.trajectory.ts)
        
        # 2. Center: Move the center group (protein) to the center of the box
        trans.center_in_box(self._center_group, center='geometry')(self._universe.trajectory.ts)
        
        # 3. Wrap: Wrap all atoms back into the unit cell, keeping molecules together
        trans.wrap(self._universe.atoms, compound='fragments')(self._universe.trajectory.ts)
        
        # Write PDB using MDAnalysis (preserves original residue names)
        self._universe.atoms.write(self._file)


# ==============================================================================
# 1. Hardware Platform Configuration (Flexible style matching equilibration.py)
# ==============================================================================
try:
    platform = Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed'}
    print("Using CUDA platform")
except Exception:
    platform = Platform.getPlatformByName('CPU')
    properties = {}
    print("CUDA not found, falling back to CPU")

# ==============================================================================
# 2. Load CHARMM Topology and Parameters
# ==============================================================================
psf_file = 'ionized.psf'
pdb_file = 'ionized.pdb'  # Structure file generated from equilibration output
params_files = [
    'par_all36m_prot.prm',
    'par_all36_cgenff.prm',
    'toppar_water_ions.str'
]

# Check if optional ligand parameters exist
if os.path.exists('ligand.str'):
    params_files.append('ligand.str')

print("Loading CHARMM parameters and PSF...")
psf = CharmmPsfFile(psf_file)
pdb = PDBFile(pdb_file)
params = CharmmParameterSet(*params_files)

# --- Dynamic Box Measurement (Resolves "non-periodic system" error) ---
# This mirrors the logic in equilibration.py to define the periodic box
coords = pdb.positions
pos_nm = np.array([p.value_in_unit(nanometer) for p in coords])
box_size = (np.max(pos_nm, axis=0) - np.min(pos_nm, axis=0)) * nanometer
psf.setBox(box_size[0], box_size[1], box_size[2])
print(f"Measured Box Size from {pdb_file}: {box_size}")

# ==============================================================================
# 3. Create System
# ==============================================================================
system = psf.createSystem(
    params,
    nonbondedMethod=PME,        # Now valid because psf.setBox was called
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005
)

# ==============================================================================
# 4. Initialize TABF Integrator
# ==============================================================================
# Parameters: dt, malpha, T_alpha, friction_alpha, hill_h, hill_s, cutoff, wtm_T, T, friction, T_range, bins, nfull
integrator = TABFIntegrator(
    2*femtoseconds, 
    1.5e4*atomic_mass_unit, 
    300*kelvin, 
    15/picosecond, 
    0.1*kilojoule_per_mole, 
    1, 2, 
    3000*kelvin, 
    300*kelvin, 
    1/picosecond, 
    [300*kelvin, 400*kelvin], 
    100, 1000
)

simulation = Simulation(psf.topology, system, integrator, platform, properties)

# ==============================================================================
# 5. Restart Logic: Load State from Equilibration
# ==============================================================================
# Checkpoint contains precise box vectors, positions, and velocities
if not os.path.exists('state.xml'):
    print("Error: state.xml not found!", file=stderr)
    print("This script requires a checkpoint file from equilibration.py to continue.", file=stderr)
    exit(1)

print("Loading checkpoint from equilibration.py...")
simulation.loadState('state.xml')

# ==============================================================================
# 6. Reporters and Output Configuration
# ==============================================================================
# Native OpenMM DCDReporter (Standard trajectory)
simulation.reporters.append(DCDReporter('wtm-tabf.dcd', 10000))
simulation.reporters.append(PDBRestartReporter('wtm-tabf_restart.pdb', 500000, psf_file, pdb_file))

# Standard log file
simulation.reporters.append(StateDataReporter(
    stdout,
    5000,
    step=True,
    potentialEnergy=True,
    temperature=True,
    volume=True,
    speed=True,
    remainingTime=True,
    totalSteps=50000000))
    
# WTM-TABF specific reporters
simulation.reporters.append(CustomVariableReporter('wtm-tabf.colvars', 100, ['alpha', 'f_alpha', 'Fbias_temp','w']))
simulation.reporters.append(PMFReporter('wtm-tabf.pmf', 10000, 100, 2))
simulation.reporters.append(MTDReporter('wtm-tabf_mtd_part.hist', 10000, 100, 2))
simulation.reporters.append(CKPTReporter('wtm-tabf_checkpoint.npy', 20000, 100, 2))

# ==============================================================================
# 7. Execution
# ==============================================================================
print("Starting WTM-TABF production MD...")
simulation.step(50000000)
import sys
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations as trans
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

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
# Custom XML State Reporter (Replaces buggy CheckpointReporter)
# ==============================================================================
class XmlStateReporter:
    """
    A custom reporter that saves the simulation state using XmlSerializer.
    
    Advantages over CheckpointReporter:
    1. Uses XML format which is more portable and readable
    2. Separately stores box vectors in a companion JSON file
    3. Avoids CheckpointReporter's known issues with state restoration
    
    Files created:
    - <base_name>.xml: Contains serialized State (positions, velocities)
    - <base_name>_box.json: Contains periodic box vectors (in nm)
    """
    def __init__(self, file, reportInterval):
        """
        Parameters
        ----------
        file : str
            Path to the output XML file (e.g., 'output/state.xml')
        reportInterval : int
            Number of steps between each save
        """
        self._file = file
        self._reportInterval = reportInterval
        # Derive the box info file path from the main file path
        base_name = os.path.splitext(file)[0]
        self._box_file = f"{base_name}_box.json"
    
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        # Request positions, velocities, and periodic box vectors
        return (steps, True, True, False, False, False)
    
    def report(self, simulation, state):
        import json
        
        # Get the full state with all required information
        full_state = simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            enforcePeriodicBox=True
        )
        
        # 1. Save state to XML file using XmlSerializer
        state_xml = XmlSerializer.serialize(full_state)
        with open(self._file, 'w') as f:
            f.write(state_xml)
        
        # 2. Save box vectors to companion JSON file
        box_vectors = full_state.getPeriodicBoxVectors(asNumpy=True)
        box_vectors_nm = box_vectors.value_in_unit(nanometer)
        
        box_info = {
            'step': simulation.currentStep,
            'time_ps': full_state.getTime().value_in_unit(picosecond),
            'box_vectors_nm': box_vectors_nm.tolist(),
            'box_lengths_nm': [
                float(box_vectors_nm[0, 0]),
                float(box_vectors_nm[1, 1]),
                float(box_vectors_nm[2, 2])
            ],
            'box_lengths_angstrom': [
                float(box_vectors_nm[0, 0] * 10),
                float(box_vectors_nm[1, 1] * 10),
                float(box_vectors_nm[2, 2] * 10)
            ]
        }
        
        with open(self._box_file, 'w') as f:
            json.dump(box_info, f, indent=2)


# ==============================================================================
# 1. Input File Configuration
# ==============================================================================
psf_file = 'ionized.psf'
pdb_file = 'ionized.pdb'
params_files = [
    'par_all36m_prot.prm',
    'par_all36_cgenff.prm',
    'toppar_water_ions.str'
]

if os.path.exists('ligand.str'):
    params_files.append('ligand.str')

# Define the ligand residue name (used for MDAnalysis selection)
ligand_resname = 'UNL'

# ==============================================================================
# 2. Read Structure and Parameters
# ==============================================================================
print(f"Loading PSF file: {psf_file} ...")
psf = CharmmPsfFile(psf_file)

print(f"Loading PDB file: {pdb_file} ...")
pdb = PDBFile(pdb_file)

print(f"Loading parameter files ...")
params = CharmmParameterSet(*params_files)

# ==============================================================================
# 3. Measure and Set Box Size (Calculated from Coordinates)
# ==============================================================================
print("Calculating box dimensions from atomic coordinate extremes...")
positions = pdb.positions
pos_nm = np.array([p.value_in_unit(nanometer) for p in positions])

min_crds = np.min(pos_nm, axis=0)
max_crds = np.max(pos_nm, axis=0)

padding = 0.0 * nanometer
box_size = (max_crds - min_crds) * nanometer + padding

psf.setBox(box_size[0], box_size[1], box_size[2])
print(f"Setting box dimensions: {box_size}")

# ==============================================================================
# 4. Create System
# ==============================================================================
print("Creating System ...")
system = psf.createSystem(
    params,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005
)

system.addForce(MonteCarloBarostat(1.01325*bar, 300*kelvin, 25))

# ==============================================================================
# 5. Add Position Restraints (Protein Backbone + Ligand Only)
# ==============================================================================
print("Adding position restraints to protein backbone and ligand atoms ...")

# Use MDAnalysis to select atoms:
# - "backbone" selects protein backbone atoms (N, CA, C, O)
# - "resname <ligand_resname>" selects all ligand atoms
u = mda.Universe(psf_file, pdb_file)
selection_string = f"backbone or resname {ligand_resname}"
selected_atoms = u.select_atoms(selection_string)
print(f"  MDAnalysis selection: '{selection_string}'")
print(f"  Number of atoms selected: {len(selected_atoms)}")

# Get the set of selected atom indices for fast lookup
selected_indices = set(selected_atoms.indices)

# Create the restraint force
restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
restraint.addGlobalParameter('k', 10.0 * kilocalories_per_mole/angstrom**2) 
restraint.addPerParticleParameter('x0')
restraint.addPerParticleParameter('y0')
restraint.addPerParticleParameter('z0')

restrained_atom_count = 0

for atom in psf.topology.atoms():
    if atom.index in selected_indices:
        pos = pdb.positions[atom.index]
        restraint.addParticle(atom.index, [pos[0], pos[1], pos[2]])
        restrained_atom_count += 1

system.addForce(restraint)
print(f"Number of atoms with restraints added: {restrained_atom_count}")

# ==============================================================================
# 6. Set Simulation Environment
# ==============================================================================
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

try:
    platform = Platform.getPlatformByName('CUDA')
    prop = {'Precision': 'mixed'}
except Exception:
    platform = Platform.getPlatformByName('CPU')
    prop = {}

simulation = Simulation(psf.topology, system, integrator, platform, prop)
simulation.context.setPositions(pdb.positions)

print(f"Using platform: {simulation.context.getPlatform().getName()}")

# ==============================================================================
# 7. Energy Minimization
# ==============================================================================
print("Starting energy minimization ...")
simulation.minimizeEnergy(maxIterations=1000)
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('minimized.pdb', 'w'))

# ==============================================================================
# 8. Restrained MD (Equilibration)
# ==============================================================================
simulation.context.setVelocitiesToTemperature(300*kelvin)

simulation.reporters.append(StateDataReporter(
    stdout, 5000, step=True, potentialEnergy=True, temperature=True, 
    volume=True, density=True, remainingTime=True, speed=True, totalSteps=50500000
))
simulation.reporters.append(DCDReporter('trajectory.dcd', 5000))
simulation.reporters.append(XmlStateReporter('state.xml', 5000))
simulation.reporters.append(PDBRestartReporter('restart.pdb', 500000, psf_file, pdb_file))

print("Starting restrained MD (1 ns) - protein/ligand restrained ...")
simulation.step(500000)

# ==============================================================================
# 9. Production MD (Unrestrained)
# ==============================================================================
print("Releasing all restraints (k=0) ...")
simulation.context.setParameter('k', 0.0)

print("Starting production MD (100 ns) - free simulation ...")
simulation.step(50000000)

print("Simulation completed!")

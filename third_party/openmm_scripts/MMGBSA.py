import os
import sys
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align
import parmed as pmd
from parmed.charmm import CharmmParameterSet as PmdCharmmParameterSet
from openmm.app import *
from openmm import *
from openmm.unit import *

from MDAnalysis.transformations import unwrap, center_in_box, wrap

# ==============================================================================
# Configuration
# ==============================================================================
PSF_FILE = 'ionized.psf'
DCD_FILE = 'trajectory.dcd'
LIGAND_RESNAME = 'UNL'
OUTPUT_CSV = 'mmgbsa_results.csv'
START_FRAME = 1000  # Start analysis from this frame (0-indexed)
STRIDE = 1  # Analyze every Nth frame

# Parameter files (Must match those used in equilibration.py)
PARAM_FILES = [
    'par_all36m_prot.prm',
    'par_all36_cgenff.prm',
    'toppar_water_ions.str'
]

if os.path.exists('ligand.str'):
    PARAM_FILES.append('ligand.str')

# ==============================================================================
# Helper Functions
# ==============================================================================

def create_subset_psf_pdb(u, pmd_struct, selection_string, prefix):
    """
    Selects atoms using MDAnalysis and saves them to PSF/PDB files using ParmEd.
    MDAnalysis doesn't support writing PSF files, so we use ParmEd for topology.
    
    Args:
        u: MDAnalysis Universe
        pmd_struct: ParmEd Structure loaded from the same PSF
        selection_string: MDAnalysis selection string
        prefix: Output file prefix
    
    Returns: (psf_filename, pdb_filename)
    """
    selection = u.select_atoms(selection_string)
    if len(selection) == 0:
        raise ValueError(f"Selection '{selection_string}' matched no atoms!")
    
    psf_out = f"{prefix}.psf"
    pdb_out = f"{prefix}.pdb"
    
    print(f"Creating subset for '{prefix}' with {len(selection)} atoms...")
    
    # Get atom indices from MDAnalysis selection (0-indexed)
    atom_indices = selection.indices.tolist()
    
    # Create a subset of the ParmEd structure using the same indices
    # ParmEd's atom indexing matches MDAnalysis when loaded from the same PSF
    subset_struct = pmd_struct[atom_indices]
    
    # Get coordinates from MDAnalysis selection (current frame)
    # PSF files don't contain coordinates, so we must set them from the trajectory
    coords = selection.positions  # Angstroms
    
    # Set coordinates for each atom in the ParmEd subset structure
    for i, atom in enumerate(subset_struct.atoms):
        atom.xx = coords[i, 0]
        atom.xy = coords[i, 1]
        atom.xz = coords[i, 2]
    
    # Write PSF and PDB using ParmEd
    subset_struct.save(psf_out, overwrite=True)
    subset_struct.save(pdb_out, overwrite=True)
    
    return psf_out, pdb_out

def setup_openmm_system(psf_path, params):
    """
    Loads a PSF and creates an OpenMM System with Implicit Solvent (OBC2).
    """
    print(f"Loading topology: {psf_path}")
    psf = CharmmPsfFile(psf_path)
    
    # Create system with Implicit Solvent (OBC2)
    # NoCutoff is used because we want the vacuum/implicit energy without PBC artifacts
    print(f"Creating OpenMM System (OBC2 Implicit Solvent)...")
    system = psf.createSystem(
        params,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
        rigidWater=True,
        implicitSolvent=OBC2
    )
    
    # Set up integrator (required for Simulation object, even if just analyzing)
    integrator = VerletIntegrator(0.001*picoseconds)
    
    # Select Platform
    try:
        platform = Platform.getPlatformByName('CUDA')
        prop = {'Precision': 'mixed'}
    except Exception:
        try:
            platform = Platform.getPlatformByName('OpenCL')
            prop = {'Precision': 'mixed'}
        except Exception:
            platform = Platform.getPlatformByName('CPU')
            prop = {}
            
    print(f"Using platform: {platform.getName()}")
    simulation = Simulation(psf.topology, system, integrator, platform, prop)
    return simulation

def calculate_potential_energy(simulation, positions_nm):
    """
    Updates the simulation positions and returns the potential energy in kcal/mol.
    """
    simulation.context.setPositions(positions_nm)
    state = simulation.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)

# ==============================================================================
# Main Execution
# ==============================================================================

def main():
    # 1. Load the Trajectory with MDAnalysis
    if not os.path.exists(PSF_FILE) or not os.path.exists(DCD_FILE):
        print(f"Error: Input files {PSF_FILE} or {DCD_FILE} not found.")
        sys.exit(1)

    print(f"Loading Universe: {PSF_FILE} | {DCD_FILE}")
    u = mda.Universe(PSF_FILE, DCD_FILE)
    
    # ---------------------------------------------------------
    # PBC Handling: Unwrap and Center
    # ---------------------------------------------------------
    # Custom transformation to force ligand to the same periodic image as protein
    # This solves the issue where wrap() might place ligand in a valid box but far from centered protein
    # and handles versions of MDAnalysis where wrap() does not support center=argument.
    
    class ForceLigandNearProtein:
        def __init__(self, protein_indices, ligand_indices):
            self.p_idx = protein_indices
            self.l_idx = ligand_indices
        
        def __call__(self, ts):
            if ts.dimensions is None: return ts
            box = ts.dimensions[:3]
            
            # Calculate geometric centers
            # Note: direct access to ts.positions is faster than atomgroup methods inside loop
            p_pos = ts.positions[self.p_idx]
            l_pos = ts.positions[self.l_idx]
            
            p_cen = np.mean(p_pos, axis=0)
            l_cen = np.mean(l_pos, axis=0)
            
            # Find shift to minimize distance
            diff = l_cen - p_cen
            shift = -np.round(diff / box) * box
            
            if np.any(shift != 0):
                ts.positions[self.l_idx] += shift
            
            return ts

    # 1. Unwrap all atoms (try to assemble continuous molecules)
    # 2. Center the protein in the box
    # 3. Wrap everything specific fragments (optional, but good for visualization)
    # 4. Force Ligand to be near Protein (Critical for energy)
    protein = u.select_atoms("protein")
    ligand_atoms = u.select_atoms(f"resname {LIGAND_RESNAME}")
    
    if len(protein) > 0 and len(ligand_atoms) > 0:
        transform_workflow = [
            unwrap(u.atoms),
            center_in_box(protein, center='geometry'),
            ForceLigandNearProtein(protein.indices, ligand_atoms.indices)
        ]
        u.trajectory.add_transformations(*transform_workflow)
        print("  Added PBC transformations: Unwrap -> Center(Protein) -> ForceLigandNearProtein")
    elif len(protein) > 0:
        transform_workflow = [
            unwrap(u.atoms),
            center_in_box(protein, center='geometry'),
            wrap(u.atoms, compound='fragments')
        ]
        u.trajectory.add_transformations(*transform_workflow)
        print("  Added PBC transformations: Unwrap -> Center(Protein) -> Wrap")
    else:
        print("  Warning: No protein found. Skipping PBC transformations.")

    n_frames = len(u.trajectory)
    print(f"Total frames in trajectory: {n_frames}")

    # Load ParmEd structure for PSF subset extraction
    # MDAnalysis doesn't support writing PSF files, so we use ParmEd
    print(f"Loading ParmEd Structure: {PSF_FILE}")
    pmd_struct = pmd.load_file(PSF_FILE)

    # Load parameters into ParmEd structure to ensure Bond types (req) are available 
    # for LonePair generation when saving subsets.
    print("Loading ParmEd Parameters to fix LonePair saving issue...")
    try:
        # Note: We must use ParmEd's CharmmParameterSet, not the OpenMM one imported via *
        pmd_params = PmdCharmmParameterSet(*PARAM_FILES)
        
        # Apply parameters to the structure (populates bond.type, etc.)
        pmd_struct.load_parameters(pmd_params)
        print("  ParmEd parameters loaded successfully.")
    except Exception as e:
        print(f"  Warning: Failed to load parameters into ParmEd structure: {e}")

    # 2. Define Selections
    complex_sel = f"protein or resname {LIGAND_RESNAME}"
    receptor_sel = "protein"
    ligand_sel = f"resname {LIGAND_RESNAME}"

    # Verify selections
    if len(u.select_atoms(complex_sel)) == 0:
        print("Error: Complex selection found 0 atoms. Check ligand resname.")
        sys.exit(1)

    # 3. Create Temporary Subset Topologies (PSF/PDB)
    # We need these to construct valid OpenMM Systems with the correct CHARMM types.
    run_prefix = "temp_mmgbsa"
    try:
        comp_psf, _ = create_subset_psf_pdb(u, pmd_struct, complex_sel, f"{run_prefix}_complex")
        rec_psf, _ = create_subset_psf_pdb(u, pmd_struct, receptor_sel, f"{run_prefix}_receptor")
        lig_psf, _ = create_subset_psf_pdb(u, pmd_struct, ligand_sel, f"{run_prefix}_ligand")

        # 4. Load CHARMM Parameters
        print("Loading Parameter Sets...")
        params = CharmmParameterSet(*PARAM_FILES)

        # 5. Setup OpenMM Simulations (Contexts)
        # We need three separate contexts to calculate G_complex, G_receptor, G_ligand
        sim_complex = setup_openmm_system(comp_psf, params)
        sim_receptor = setup_openmm_system(rec_psf, params)
        sim_ligand = setup_openmm_system(lig_psf, params)

        # 6. MM-GBSA Energy Analysis Loop
        # NOTE: OpenMM with implicit solvent (OBC2) returns: E = E_MM + G_GB + G_SA
        #       where G_GB is the polar solvation free energy (Generalized Born)
        #       and G_SA is the nonpolar solvation free energy (Surface Area)
        #       This is often called "effective energy" or "quasi-free energy" in literature.
        
        # Pre-select atom groups to extract coordinates efficiently
        ag_complex = u.select_atoms(complex_sel)
        ag_receptor = u.select_atoms(receptor_sel)
        ag_ligand = u.select_atoms(ligand_sel)

        # Create a reference Universe for alignment (to handle PBC/drift/wrapping)
        # We align every frame to this reference to ensure the complex is intact
        # and not split across periodic boundaries before energy calculation.
        print("Loading reference frame for alignment...")
        ref_u = mda.Universe(PSF_FILE, DCD_FILE)
        
        
        # Apply the same transformations to the reference universe
        if len(protein) > 0 and len(ligand_atoms) > 0:
             # Re-create transform for ref_u with its own indices (though they should be same)
            ref_prot = ref_u.select_atoms("protein")
            ref_lig = ref_u.select_atoms(f"resname {LIGAND_RESNAME}")
            
            ref_workflow = [
                unwrap(ref_u.atoms),
                center_in_box(ref_prot, center='geometry'),
                ForceLigandNearProtein(ref_prot.indices, ref_lig.indices)
            ]
            ref_u.trajectory.add_transformations(*ref_workflow)
            print("  Added PBC transformations to reference universe")
        elif len(protein) > 0:
            ref_u.trajectory.add_transformations(unwrap(ref_u.atoms), center_in_box(ref_u.select_atoms("protein"), center='geometry'), wrap(ref_u.atoms, compound='fragments'))
            print("  Added PBC transformations to reference universe")
            
        ref_u.trajectory[START_FRAME]
        # Align on Protein CA atoms (backbone)


        results = []
        frames_to_process = list(range(START_FRAME, n_frames, STRIDE))
        print(f"\nStarting MM-GBSA energy calculation on {len(frames_to_process)} frames (Start={START_FRAME}, Stride={STRIDE})...")

        for i, frame_idx in enumerate(frames_to_process):
            u.trajectory[frame_idx]
            
            # --- Sanity Check: Ligand-Protein Distance ---
            # Detect PBC wrapping artifacts where ligand wraps to the other side of the box
            prot_center = ag_receptor.center_of_geometry()
            lig_center = ag_ligand.center_of_geometry()
            dist = np.linalg.norm(prot_center - lig_center)
            
            # Heuristic: If distance is > 45% of the box dimension, it's likely a wrapping artifact
            # (Assuming protein is centered, max valid distance is < 50% of box)
            dims = u.dimensions
            if dims is not None:
                min_box_dim = min(dims[:3])
                cutoff_dist = min_box_dim * 0.45 
                if dist > cutoff_dist:
                    print(f"  WARINING: Skipping frame {frame_idx} due to large Ligand-Protein distance ({dist:.1f} A > {cutoff_dist:.1f} A). Likely PBC artifact.")
                    continue
            
            # Perform alignment to reference frame
            # This centers the complex and corrects relative positions if they were 
            # wrapped disjointly but consistent with the reference frame's topology
            align.alignto(u, ref_u, select="protein and name CA")

            time_ps = u.trajectory.time

            # Extract coordinates (MDAnalysis uses Angstroms, OpenMM uses Nanometers)
            # Conversion factor: 0.1
            # Add explicit units (nanometer) for OpenMM safety
            pos_complex = (ag_complex.positions * 0.1) * nanometer
            pos_receptor = (ag_receptor.positions * 0.1) * nanometer
            pos_ligand = (ag_ligand.positions * 0.1) * nanometer

            # Calculate MM-GBSA Energies (E_MM + G_solvation)
            G_complex = calculate_potential_energy(sim_complex, pos_complex)
            G_receptor = calculate_potential_energy(sim_receptor, pos_receptor)
            G_ligand = calculate_potential_energy(sim_ligand, pos_ligand)

            # ΔG_MMGBSA = G_complex - (G_receptor + G_ligand)
            # This is the MM-GBSA binding energy
            delta_G_mmgbsa = G_complex - (G_receptor + G_ligand)

            results.append({
                'Frame': frame_idx,
                'Time_ps': time_ps,
                'G_Complex': G_complex,
                'G_Receptor': G_receptor,
                'G_Ligand': G_ligand,
                'Delta_G_MMGBSA': delta_G_mmgbsa
            })

            # Progress output
            if i % 10 == 0:
                print(f"Frame {frame_idx} ({i+1}/{len(frames_to_process)}): ΔG_MMGBSA = {delta_G_mmgbsa:.4f} kcal/mol")

        # 7. Calculate Mean MM-GBSA Binding Energy
        df = pd.DataFrame(results)
        mean_delta_G_mmgbsa = df['Delta_G_MMGBSA'].mean()
        std_delta_G_mmgbsa = df['Delta_G_MMGBSA'].std()

        # 8. Save Detailed Results
        df.to_csv(OUTPUT_CSV, index=False)
        
        # Save summary
        summary_file = OUTPUT_CSV.replace('.csv', '_summary.txt')
        with open(summary_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write("MMGBSA Binding Free Energy Analysis Results\n")
            f.write("="*70 + "\n\n")
            
            f.write(f"Configuration:\n")
            f.write(f"  Frames analyzed: {len(results)}\n")
            f.write(f"  Start frame: {START_FRAME}\n")
            f.write(f"  Stride: {STRIDE}\n\n")
            
            f.write("-"*70 + "\n")
            f.write("RESULTS (all values in kcal/mol)\n")
            f.write("-"*70 + "\n\n")
            
            f.write(f"MM-GBSA Binding Energy (ΔG_MMGBSA = E_MM + G_GB + G_SA):\n")
            f.write(f"   Mean ΔG_MMGBSA = {mean_delta_G_mmgbsa:.4f} ± {std_delta_G_mmgbsa:.4f}\n")
            f.write(f"   Min  ΔG_MMGBSA = {df['Delta_G_MMGBSA'].min():.4f}\n")
            f.write(f"   Max  ΔG_MMGBSA = {df['Delta_G_MMGBSA'].max():.4f}\n\n")
            
            f.write("="*70 + "\n")
        
        # 9. Print Summary to Console
        print("\n" + "="*70)
        print("MMGBSA Binding Free Energy Analysis Results")
        print("="*70)
        
        print(f"\n{'RESULTS':^70}")
        print("-"*70)
        print(f"\n  MM-GBSA Binding Energy (E_MM + G_solvation):")
        print(f"     Mean ΔG_MMGBSA = {mean_delta_G_mmgbsa:.4f} ± {std_delta_G_mmgbsa:.4f} kcal/mol")
        print(f"     Min  ΔG_MMGBSA = {df['Delta_G_MMGBSA'].min():.4f} kcal/mol")
        print(f"     Max  ΔG_MMGBSA = {df['Delta_G_MMGBSA'].max():.4f} kcal/mol")
        
        print("\n" + "="*70)
        print(f"\nDetailed frame-by-frame results saved to: {os.path.abspath(OUTPUT_CSV)}")
        print(f"Summary saved to: {os.path.abspath(summary_file)}")

    finally:
        # Cleanup temporary files
        print("\nCleaning up temporary files...")
        temp_files = [
            f"{run_prefix}_complex.psf", f"{run_prefix}_complex.pdb",
            f"{run_prefix}_receptor.psf", f"{run_prefix}_receptor.pdb",
            f"{run_prefix}_ligand.psf", f"{run_prefix}_ligand.pdb"
        ]
        for f in temp_files:
            if os.path.exists(f):
                os.remove(f)

if __name__ == "__main__":
    main()

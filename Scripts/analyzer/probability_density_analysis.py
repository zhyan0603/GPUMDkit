import sys
import numpy as np
from pymatgen.analysis.diffusion.aimd.pathway import ProbabilityDensityAnalysis
from pymatgen.core import Structure
from ase.io import read

# Check if required arguments are provided
if len(sys.argv) != 5:
    print("Usage: python script.py <ref_struct> <trajectory_file> <species> <interval>")
    sys.exit(1)

# Parse command-line arguments
structure_file = sys.argv[1]
trajectory_file = sys.argv[2]
species = (sys.argv[3],)  # Convert species to tuple
interval = float(sys.argv[4])

# Load structure from current directory
structure = Structure.from_file(structure_file)

# Read trajectory from dump.xyz (non-ULM format)
trajectory = read(trajectory_file, index=":")
trajectory_npy = [atom.get_scaled_positions() for atom in trajectory]

# Perform probability density analysis
pda = ProbabilityDensityAnalysis(
    structure=structure,
    trajectories=trajectory_npy,
    interval=interval,
    species=species
)

# Save probability density to CHGCAR file
pda.to_chgcar(filename=f"probability_density_{interval}.vasp")
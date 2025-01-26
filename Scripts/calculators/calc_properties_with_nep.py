import numpy as np
from ase import Atoms
from ase.io import read, write
from calorine.calculators import CPUNEP
import sys

# Specify input and output files
input_file = sys.argv[1]  # Input file
output_file = sys.argv[2]  # Output file

# Read input file containing all structures
structures = read(input_file, index=':')  # Read all frames

# Ensure structures is a list
if not isinstance(structures, list):
    structures = [structures]

# Process each structure
for atoms in structures:
    # Create a new calculator instance for each atoms object
    calc = CPUNEP(sys.argv[3])
    atoms.set_calculator(calc)
    
    # Calculate properties
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    stress = atoms.get_stress()
    #volume = atoms.get_volume()
    #virial = -stress * volume  # Virial calculation
    
    # Store virial in info
    #atoms.info['virial'] = virial.tolist()  # Convert to list for extxyz compatibility
    #atoms.info['virial'] = ' '.join(map(str, virial.tolist()))
    #atoms.stress = None

# Write the updated structures to the output file
write(output_file, structures)
print(f"Results saved to {output_file}")

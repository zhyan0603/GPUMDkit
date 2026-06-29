"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     calc_properties_with_nep.py
Category:   Calculator Scripts
Purpose:    Calculate energy, forces, and stress for structures using a NEP
            model via calorine, and write results to an output extxyz file.
Usage:      python calc_properties_with_nep.py <input.xyz> <output.xyz> <nep.txt>
Arguments:
  input.xyz   Input extxyz file
  output.xyz  Output extxyz file with computed properties
  nep.txt     Path to the NEP model file
Output:
  <output.xyz>  (structures with NEP-computed energy, forces, stress)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
import numpy as np
from ase import Atoms
from ase.io import read, write
from calorine.calculators import CPUNEP


def print_dependency_notice():
    print(" This function requires the calorine package.")
    print(" If you use this function, we recommend citing:")
    print(" Lindgren et al., J. Open Source Softw. 9, 6264 (2024).")
    print(" https://doi.org/10.21105/joss.06264")


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end='\r')
    if iteration == total:
        print()

# Specify input and output files
input_file = sys.argv[1]  # Input file
output_file = sys.argv[2]  # Output file

print_dependency_notice()

# Read input file containing all structures
structures = read(input_file, index=':')  # Read all frames

# Ensure structures is a list
if not isinstance(structures, list):
    structures = [structures]

total_structures = len(structures)

calc = CPUNEP(sys.argv[3])
# Process each structure
for i, atoms in enumerate(structures):
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
    print_progress_bar(i + 1, total_structures, prefix=' Processing', suffix='Complete', length=50)

# Write the updated structures to the output file
write(output_file, structures)
print(f" Results saved to {output_file}")

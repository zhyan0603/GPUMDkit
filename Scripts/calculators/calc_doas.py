"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     calc_doas.py
Category:   Calculator Scripts
Purpose:    Calculate the density of atomistic states (DOAS) proposed by
            Wang et al. by reading and optimizing structures with NEP_CPU
            and computing per-atom energies.
Usage:      python calc_doas.py <input.xyz> <nep.txt> <output.txt>
Arguments:
  input.xyz   Input extxyz file
  nep.txt     Path to the NEP model file
  output.txt  Output text file for grouped atomic energies
Output:
  <output.txt>  (per-atom energies grouped by element)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""
import sys
from tqdm import tqdm
import numpy as np
from ase import Atoms
from ase.io import read
from ase.optimize import BFGS
from calorine.calculators import CPUNEP
from collections import defaultdict


def print_dependency_notice():
    print(" This function requires the calorine package.")
    print(" If you use this function, we recommend citing:")
    print(" Lindgren et al., J. Open Source Softw. 9, 6264 (2024).")
    print(" https://doi.org/10.21105/joss.06264")


# Specify input and output files
input_file = sys.argv[1]  # Input extxyz file
model_path = sys.argv[2]  # Path to the model for CPUNEP
output_file = sys.argv[3]  # Output text file for grouped atomic energies

print_dependency_notice()

# Read input file containing all structures
structures = read(input_file, index=':')  # Read all frames

# Ensure structures is a list
if not isinstance(structures, list):
    structures = [structures]

# Dictionary to store atomic energies grouped by element
energies_by_element = defaultdict(list)

# Process each structure with tqdm progress bar
for atoms in tqdm(structures, desc="Optimizing structures", unit="frame"):
    # Create a new calculator instance and assign it
    calc = CPUNEP(model_path)
    atoms.calc = calc  # Updated to avoid FutureWarning

    # Perform structure optimization, suppress console output
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=0.05)  # Convergence criterion: max force 0.05 eV/Å

    # Get per-atom energies using CPUNEP's nepy interface
    energies, forces, virials = atoms.calc.nepy.get_potential_forces_and_virials()
    atomic_energies = np.array(energies)

    # Verify that the sum of per-atom energies matches the total energy
    total_energy = atoms.get_potential_energy()
    if not np.isclose(atomic_energies.sum(), total_energy, rtol=1e-5):
        raise ValueError(f"Sum of per-atom energies ({atomic_energies.sum()}) does not match total energy ({total_energy}).")

    # Group energies by element
    symbols = atoms.get_chemical_symbols()
    for sym, energy in zip(symbols, atomic_energies):
        energies_by_element[sym].append(energy)

# Write the grouped atomic energies to the output text file
with open(output_file, 'w') as f:
    for element, energies in energies_by_element.items():
        f.write(f"# {element}\n")
        for energy in energies:
            f.write(f"{energy}\n")
        f.write("\n")  # Blank line to separate elements

print(f"Results saved to {output_file}")

"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     analyze_chem_species.py
Category:   Analyzer Scripts
Purpose:    Identify all unique chemical elements present in an extxyz
            trajectory.
Usage:      gpumdkit.sh -chem_species <input.xyz>
            python analyze_chem_species.py <input.xyz>
Arguments:
  input.xyz  Input extxyz file
Output:
  List of unique chemical elements sorted by atomic number
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 1 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -chem_species <exyzfile>")
    print("    or: python analyze_chem_species.py <input.xyz>")
    print("")
    print(" Arguments:")
    print("   exyzfile    Input extxyz trajectory file")
    print("")
    print(" Example: gpumdkit.sh -chem_species train.xyz")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

from ase.io import read
from ase.data import chemical_symbols, atomic_numbers

atoms_list = read(sys.argv[1], index=':')

# Identify all unique elements present across all frames
element_set = set()
for atoms in atoms_list:
    element_set.update(atoms.get_chemical_symbols())

# Sort elements according to their atomic numbers (Periodic Table order)
sorted_elements = sorted(element_set, key=lambda x: atomic_numbers[x])

# Display results
print(f'Total unique elements: {len(sorted_elements)}')
print('Element list (Periodic Table order):')
print(' '.join(sorted_elements))

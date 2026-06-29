"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     analyze_chem_species.py
Category:   Analyzer Scripts
Purpose:    Identify all unique chemical elements present in an extxyz
            trajectory.
Usage:      python analyze_chem_species.py <input.xyz>
Arguments:
  input.xyz  Input extxyz file
Output:
  List of unique chemical elements sorted by atomic number
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
from ase.io import read
from ase.data import chemical_symbols, atomic_numbers

# Load all frames from the provided extxyz file
try:
    atoms_list = read(sys.argv[1], index=':')
except IndexError:
    print("Usage: python analyze_chem_species.py <input.xyz>")
    sys.exit(1)

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

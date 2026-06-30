"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     cif2exyz.py
Category:   Format Conversion Scripts
Purpose:    Convert CIF file to extended XYZ format using ASE.
Usage:      gpumdkit.sh -cif2exyz <input.cif> <output.xyz>
            python cif2exyz.py <input.cif> <output.xyz>
Arguments:
  input.cif   Input CIF file
  output.xyz  Output extxyz file
Output:
  <output.xyz>  (converted structure in extxyz format)
Author:     Boyi Situ (situboyi@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
from ase.io import read, write

# Check command line arguments
if len(sys.argv) != 3:
    print(" Usage: gpumdkit.sh -cif2exyz <input.cif> <output.xyz>")
    print("    or: python cif2exyz.py <input.cif> <output.xyz>")
    sys.exit(1)

# Read input and output file paths
cif_file = sys.argv[1]
xyz_file = sys.argv[2]

# Read CIF structure
try:
    atoms = read(cif_file, format='cif')
except Exception as e:
    print(f" Error reading CIF file {cif_file}: {e}")
    sys.exit(1)

# Remove specified fields from atoms.info
fields_to_remove = ['spacegroup', 'unit_cell', 'occupancy']
for field in fields_to_remove:
    atoms.info.pop(field, None)

if hasattr(atoms, 'arrays') and 'spacegroup_kinds' in atoms.arrays:
    del atoms.arrays['spacegroup_kinds']

# Write exerts format
try:
    write(xyz_file, atoms, format='extxyz')
    print(f" Converted {cif_file} to {xyz_file}")
except Exception as e:
    print(f" Error writing extxyz file {xyz_file}: {e}")
    sys.exit(1)
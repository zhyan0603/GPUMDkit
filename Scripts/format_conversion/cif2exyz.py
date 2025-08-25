"""
Convert CIF file to VASP extxyz format using ASE.
Usage:
    python cif2exyz.py input.cif model.xyz
    Author: Boyi Situ (situboyi@westlake.edu.cn)
    Date: 2025-08-21
"""

import sys
from ase.io import read, write

# Check command line arguments
if len(sys.argv) != 3:
    print(" Usage: python cif2exyz.py input.cif output.xyz")
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
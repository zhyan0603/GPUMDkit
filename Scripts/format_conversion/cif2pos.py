"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Convert CIF file to VASP POSCAR format using ASE.
Usage:
    python cif2pos.py input.cif output.vasp
    Author: Boyi Situ (situboyi@westlake.edu.cn)
    Date: 2025-08-21
"""


import sys
from ase.io import read, write

# Check command line arguments
if len(sys.argv) != 3:
    print(" Usage: python cif2poscar.py input.cif output.vasp")
    sys.exit(1)

# Read input and output file paths
cif_file = sys.argv[1]
poscar_file = sys.argv[2]

# Read CIF structure
atoms = read(cif_file)

# Write POSCAR format
write(poscar_file, atoms, format='vasp', vasp5=True, direct=True)

print(f" Converted {cif_file} to {poscar_file}")

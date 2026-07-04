"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     cif2pos.py
Category:   Format Conversion Scripts
Purpose:    Convert CIF file to VASP POSCAR format using ASE.
Usage:      gpumdkit.sh -cif2pos <input.cif> <output.vasp>
            python cif2pos.py <input.cif> <output.vasp>
Arguments:
  input.cif    Input CIF file
  output.vasp  Output VASP POSCAR file
Output:
  <output.vasp>  (converted structure in POSCAR format)
Author:     Boyi Situ (situboyi@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 2 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -cif2pos <input.cif> <output.vasp>")
    print("    or: python cif2pos.py <input.cif> <output.vasp>")
    print("")
    print(" Arguments:")
    print("   input.cif    Input CIF file")
    print("   output.vasp  Output VASP POSCAR file")
    print("")
    print(" Example: gpumdkit.sh -cif2pos structure.cif POSCAR")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

from ase.io import read, write

# Read input and output file paths
cif_file = sys.argv[1]
poscar_file = sys.argv[2]

# Read CIF structure
atoms = read(cif_file)

# Write POSCAR format
write(poscar_file, atoms, format='vasp', vasp5=True, direct=True)

print(f" Converted {cif_file} to {poscar_file}")

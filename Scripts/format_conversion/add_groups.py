"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     add_groups.py
Category:   Format Conversion Scripts
Purpose:    Add group information to atoms in a structure file based on
            element types, outputting the result to model.xyz.
Usage:      gpumdkit.sh -addgroup <input.xyz> <element1> <element2> ...
            gpumdkit.sh -addlabel <input.xyz> <element1> <element2> ...
            python add_groups.py <input.xyz> <element1> <element2> ...
Arguments:
  input.xyz   Input structure file
  elementX    Element symbols to assign group indices
Output:
  model.xyz  Structure with group information
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 2 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -addgroup <input.xyz> <element1> <element2> ...")
    print("        gpumdkit.sh -addlabel <input.xyz> <element1> <element2> ...")
    print("    or: python add_groups.py <input.xyz> <element1> <element2> ...")
    print("")
    print(" Arguments:")
    print("   input.xyz   Input structure file")
    print("   elementX    Element symbols to assign group indices (in order)")
    print("")
    print(" Output:")
    print("   model.xyz   Structure with group information")
    print("")
    print(" Example: gpumdkit.sh -addgroup POSCAR Li Ti O")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

from ase.io import read, write
import numpy as np

# Read the file name from command line arguments
file_name = sys.argv[1]

# Read the elements order from command line arguments
elements = sys.argv[2:]

# Read the atomic data from the file
atoms = read(file_name)

# Create groups corresponding to the elements
groups = []
for element in atoms.get_chemical_symbols():
    if element in elements:
        group = elements.index(element)
    else:
        raise ValueError(f"Element {element} not found in the provided elements list")
    groups.append(group)

# Add the group information to the atomic data
groups_array = np.array(groups)
atoms.new_array("group", groups_array)

# Write the output to a file
write("model.xyz", atoms)

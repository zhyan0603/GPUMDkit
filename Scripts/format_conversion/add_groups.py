"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     add_groups.py
Category:   Format Conversion Scripts
Purpose:    Add group information to atoms in a structure file based on
            element types, outputting the result to stdout.
Usage:      python add_groups.py <input.xyz> <element1> <element2> ...
Arguments:
  input.xyz   Input structure file
  elementX    Element symbols to assign group indices
Output:
  Modified structure with group information (printed to stdout)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
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

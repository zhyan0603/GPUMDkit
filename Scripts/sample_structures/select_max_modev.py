"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     select_max_modev.py
Category:   Sample Structure Scripts
Purpose:    Extract the top N structures with the highest max force
            deviation from active learning output (active.out + active.xyz).
Usage:      gpumdkit.sh
            choose 205) Select max force deviation structs
            python select_max_modev.py <top_n> <min_deviation>
Arguments:
  top_n          Number of top structures to extract
  min_deviation  Minimum force deviation threshold for filtering
Input files:
  active.out     Max force deviation per structure
  active.xyz     Corresponding structures in extxyz format
Output:
  selected.xyz
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
import numpy as np
from ase.io import read, write


def print_usage():
    print(" Usage: gpumdkit.sh")
    print("        choose 205) Select max force deviation structs")
    print("    or: python select_max_modev.py <top_n> <min_deviation>")
    print("")
    print(" Arguments:")
    print("   top_n          Number of top structures to extract")
    print("   min_deviation  Minimum force deviation threshold for filtering")
    print("")
    print(" Input files:")
    print("   active.out     Max force deviation per structure")
    print("   active.xyz     Corresponding structures in extxyz format")
    print("")
    print(" Output:")
    print("   selected.xyz")
    print("")
    print(" Example: in interactive mode, enter: 200 0.15")
    print("          python select_max_modev.py 200 0.15")
    print("")


args = sys.argv[1:]
if len(args) != 2 or args[0] in ("-h", "--help"):
    print_usage()
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

try:
    top_n = int(args[0])
    min_deviation = float(args[1])
except ValueError:
    print(" Error: <top_n> must be an integer and <min_deviation> must be a number.")
    print("")
    print_usage()
    sys.exit(1)

# Define file paths
input_file = 'active.out'  # Input file for max force deviation
#output_file = 'filtered_max_force_deviations.txt'  # Output file for filtered max force deviations
xyz_file = 'active.xyz'  # Input .extxyz file
output_xyz_file = 'selected.xyz'  # Output file for extracted structures

# Lists to store filtered max force deviations and corresponding indices
filtered_deviations = []
filtered_indices = []

# Read and filter data based on max force deviation
with open(input_file, 'r') as f:
    lines = f.readlines()
    
    for idx, line in enumerate(lines):
        time, max_force_deviation = line.split()
        max_force_deviation = float(max_force_deviation)

        if max_force_deviation > min_deviation:
            filtered_deviations.append(max_force_deviation)
            filtered_indices.append(idx)  # Store the index in the original file

# Output the number of filtered structures
print(f' Number of structures in active.xyz: {len(filtered_deviations)}')

# Save filtered max force deviations to file
#with open(output_file, 'w') as f:
#    for deviation in filtered_deviations:
#        f.write(f'{deviation}\n')

# Get the indices of the top N max force deviations (sorted in descending order)
top_indices_in_filtered = sorted(range(len(filtered_deviations)), key=lambda i: filtered_deviations[i], reverse=True)[:top_n]
top_deviations = [filtered_deviations[i] for i in top_indices_in_filtered]

# Output the indices of the top N max force deviations in active.xyz
print(f' Indices of the top {top_n} max force deviations in active.xyz:\n {top_indices_in_filtered}')

# Function to extract structures from .extxyz file based on given indices
def extract_xyz_structures(xyz_file, indices, output_xyz_file):
    # Read all structures from the .extxyz file using ASE
    atoms_list = read(xyz_file, index=':')  # Read all structures as a list of ASE atoms objects

    # Extract the corresponding structures based on the indices from the filtered list
    selected_atoms = [atoms_list[i] for i in indices]
    
    # Write the selected structures to the output .xyz file
    write(output_xyz_file, selected_atoms)

    print(f' Extracted structures saved to {output_xyz_file}')

# Extract the top N structures and save them to the output .xyz file
extract_xyz_structures(xyz_file, top_indices_in_filtered, output_xyz_file)

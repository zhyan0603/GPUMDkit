"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     filter_dist_range.py
Category:   Analyzer Scripts
Purpose:    Filter structures by minimum interatomic distance between two
            elements falling within a specified range.
Usage:      gpumdkit.sh -filter_range <input.xyz> <element1> <element2> <min_dist> <max_dist>
            python filter_dist_range.py <input.xyz> <element1> <element2> <min_dist> <max_dist>
Arguments:
  input.xyz  Input extxyz file
  element1   First element symbol (e.g., Li)
  element2   Second element symbol (e.g., Cl)
  min_dist   Minimum distance (Angstrom)
  max_dist   Maximum distance (Angstrom)
Output:
  filtered_<elem1>_<elem2>_<min>_<max>.xyz
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 5 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -filter_range <exyzfile> <element1> <element2> <min_dist> <max_dist>")
    print("    or: python filter_dist_range.py <input.xyz> <element1> <element2> <min_dist> <max_dist>")
    print("")
    print(" Arguments:")
    print("   exyzfile    Input extxyz trajectory file")
    print("   element1    First element symbol (e.g., Li)")
    print("   element2    Second element symbol (e.g., Cl)")
    print("   min_dist    Minimum distance (Angstrom)")
    print("   max_dist    Maximum distance (Angstrom)")
    print("")
    print(" Output:")
    print("   filtered_<element1>_<element2>_<min_dist>_<max_dist>.xyz")
    print("")
    print(" Example: gpumdkit.sh -filter_range train.xyz Li Cl 1.5 3.0")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

import numpy as np
from ase.io import read, write

def get_min_distance(atoms, symbol1, symbol2):
    """Calculate the minimum distance between atoms of symbol1 and symbol2."""
    # Get indices of atoms for both symbols
    indices1 = [i for i, symbol in enumerate(atoms.get_chemical_symbols()) if symbol == symbol1]
    indices2 = [i for i, symbol in enumerate(atoms.get_chemical_symbols()) if symbol == symbol2]
    
    if len(indices1) < 1 or len(indices2) < 1:  # Need at least one atom of each type
        return float('inf')
    
    # Get positions of relevant atoms
    positions1 = atoms.get_positions()[indices1]
    positions2 = atoms.get_positions()[indices2]
    
    # Compute pairwise distances using vectorized operations
    diff = positions1[:, np.newaxis, :] - positions2[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff**2, axis=-1))
    
    # If same symbol, exclude self-distances
    if symbol1 == symbol2:
        np.fill_diagonal(distances, np.inf)
    
    # Return minimum distance
    return np.min(distances)

def main():
    input_file = sys.argv[1]
    atom1 = sys.argv[2]
    atom2 = sys.argv[3]
    try:
        min_dist = float(sys.argv[4])
        max_dist = float(sys.argv[5])
    except ValueError:
        print(" Error: min_dist and max_dist must be numbers")
        sys.exit(1)
    
    # Set default output file name
    output_file = f"filtered_{atom1}_{atom2}_{min_dist}_{max_dist}.xyz"
    
    # Read all frames
    try:
        structures = read(input_file, index=':', format='extxyz')
    except FileNotFoundError:
        print(f" Error: Input file {input_file} not found")
        sys.exit(1)
    
    # List to store filtered structures
    filtered_structures = []
    
    # Process each frame
    for idx, atoms in enumerate(structures):
        min_distance = get_min_distance(atoms, atom1, atom2)
        if min_dist <= min_distance <= max_dist:
            print(f" Frame {idx}: Minimum {atom1}-{atom2} distance = {min_distance:.3f} Å")
            filtered_structures.append(atoms)
    
    # Write filtered structures to output file
    if filtered_structures:
        write(output_file, filtered_structures, format='extxyz')
        print(f" Filtered structures written to {output_file}")
    else:
        print(f" No structures found with {atom1}-{atom2} distance between {min_dist} and {max_dist} Å")

if __name__ == "__main__":
    main()
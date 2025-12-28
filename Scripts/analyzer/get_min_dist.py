"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Calculate minimum distance between atoms

Usage:
    python get_min_dist.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

import sys
import numpy as np
from ase.io import read
from scipy.spatial.distance import pdist, squareform

# Read the file name from command line arguments
file_name = sys.argv[1]

# Read all frames from the extxyz file
frames = read(file_name, index=':')

# Get the original order of atom types from the first frame
first_frame_symbols = frames[0].get_chemical_symbols()
# Get unique symbols while preserving order
unique_symbols = []
[unique_symbols.append(s) for s in first_frame_symbols if s not in unique_symbols]

# Dictionary to store minimum distances between atom pairs
min_distances = {}
# Variable to track overall minimum distance
overall_min_distance = float('inf')

# Iterate over each frame
for frame in frames:
    # Get atomic positions and chemical symbols
    positions = frame.get_positions()
    symbols = frame.get_chemical_symbols()

    # Compute pairwise distances
    distances = squareform(pdist(positions))

    # Calculate minimum distance for each atom pair type
    for i, sym1 in enumerate(unique_symbols):
        for sym2 in unique_symbols[i:]:  # Start from i to include same-type pairs
            # Get indices of atoms for both types
            idx1 = [i for i, s in enumerate(symbols) if s == sym1]
            idx2 = [i for i, s in enumerate(symbols) if s == sym2]

            # Get all distances between these atom types
            pair_distances = distances[np.ix_(idx1, idx2)]

            # Exclude zero distances (same atom) if sym1 == sym2
            if sym1 == sym2:
                pair_distances = pair_distances[pair_distances > 0]

            # Update minimum distance if we have valid distances
            if len(pair_distances) > 0:
                min_dist = np.min(pair_distances)
                pair_key = f"{sym1}-{sym2}"
                if pair_key not in min_distances or min_dist < min_distances[pair_key]:
                    min_distances[pair_key] = min_dist
                # Update overall minimum distance
                overall_min_distance = min(overall_min_distance, min_dist)

# Print results in table format
print(" +---------------------------+")
print(" |   PBC ignored for speed   |")
print(" | use -min_dist_pbc for PBC |")
print(" +---------------------------+")
print(" Minimum interatomic distances:")
print(" +---------------------------+")
print(" | Atom Pair |  Distance (Å) |")
print(" +---------------------------+")
for i, sym1 in enumerate(unique_symbols):
    for sym2 in unique_symbols[i:]:
        pair = f"{sym1}-{sym2}"
        distance = min_distances[pair]
        print(f" |   {pair:<6}  |     {distance:>5.3f}     |")
print(" +---------------------------+")
print(f" Overall min_distance: {overall_min_distance:.3f} Å")
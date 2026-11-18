"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     get_min_dist_pbc.py
Category:   Analyzer Scripts
Purpose:    Compute the minimum interatomic distance for each element pair
            across all frames in an extxyz file, accounting for periodic
            boundary conditions.
Usage:      gpumdkit.sh -min_dist_pbc <file_name>
            python get_min_dist_pbc.py <file_name>
Example:    gpumdkit.sh -min_dist_pbc train.xyz
Arguments:
  file_name  Input extxyz file
Output:
  Table of minimum distances for each element pair and overall minimum
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-11-10
=============================================================================
"""

import os
import sys

args = sys.argv[1:]
if len(args) < 1 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -min_dist_pbc <file_name>")
    print("    or: python get_min_dist_pbc.py <file_name>")
    print("")
    print(" Arguments:")
    print("   file_name   Input extxyz trajectory file")
    print("")
    print(" Output:")
    print("   Table of minimum distances (with PBC) for each element pair")
    print("")
    print(" Example: gpumdkit.sh -min_dist_pbc train.xyz")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

import numpy as np
from ase.geometry import get_distances
from ase.io import read

# Read the file name from command line arguments
file_name = sys.argv[1]

if not os.path.isfile(file_name):
    print(f" Error: file '{file_name}' does not exist.")
    sys.exit(1)

# Read all frames from the extxyz file
frames = read(file_name, index=':')

# Get unique symbols across all frames while preserving first appearance order.
unique_symbols = []
for frame in frames:
    for symbol in frame.get_chemical_symbols():
        if symbol not in unique_symbols:
            unique_symbols.append(symbol)

# Dictionary to store minimum distances between atom pairs
min_distances = {}
# Variable to track overall minimum distance
overall_min_distance = float('inf')

for frame in frames:
    symbols = np.asarray(frame.get_chemical_symbols())
    positions = frame.positions
    indices = {symbol: np.flatnonzero(symbols == symbol) for symbol in unique_symbols}
    for i, sym1 in enumerate(unique_symbols):
        idx1 = indices[sym1]
        for sym2 in unique_symbols[i:]:
            idx2 = indices[sym2]
            if not len(idx1) or not len(idx2):
                continue
            block_size = max(1, 8192 // len(idx2))
            pair_min = np.inf
            for start in range(0, len(idx1), block_size):
                rows = idx1[start:start + block_size]
                _, distances = get_distances(positions[rows], positions[idx2], frame.cell, frame.pbc)
                distances[rows[:, None] == idx2] = np.inf
                pair_min = min(pair_min, distances.min())
            if np.isfinite(pair_min):
                key = f"{sym1}-{sym2}"
                min_distances[key] = min(min_distances.get(key, np.inf), pair_min)
                overall_min_distance = min(overall_min_distance, pair_min)

# Print results in table format
print(" Minimum interatomic distances (with PBC):")
print(" +---------------------------+")
print(" | Atom Pair |  Distance (Å) |")
print(" +---------------------------+")
for i, sym1 in enumerate(unique_symbols):
    for sym2 in unique_symbols[i:]:
        pair = f"{sym1}-{sym2}"
        if pair in min_distances:
            distance = min_distances[pair]
            print(f" |   {pair:<6}  |     {distance:>5.3f}     |")
        else:
            print(f" |   {pair:<6}  |       N/A     |")
print(" +---------------------------+")
if overall_min_distance < float('inf'):
    print(f" Overall min_distance: {overall_min_distance:.3f} Å")
else:
    print(" Overall min_distance: N/A")

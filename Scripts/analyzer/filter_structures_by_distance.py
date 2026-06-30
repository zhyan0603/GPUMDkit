"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     filter_structures_by_distance.py
Category:   Analyzer Scripts
Purpose:    Filter structures by minimum interatomic distance without
            periodic boundary conditions (no PBC). Faster than the PBC
            version for large systems.
Usage:      gpumdkit.sh -filter_dist <file_name> [distance_threshold]
            python filter_structures_by_distance.py <file_name> [distance_threshold]
Example:    gpumdkit.sh -filter_dist train.xyz 1.5
Arguments:
  file_name          Input extxyz file
  distance_threshold (optional) Minimum allowed interatomic distance
Output:
  filtered_<file_name>       Structures that pass the filter
  filtered_out_<file_name>   Structures that fail the filter
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
import numpy as np
from ase.io import read, write
from scipy.spatial.distance import pdist

def print_progress_bar(iteration, total, length=50):
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = '█' * filled_length + '-' * (length - filled_length)
    print(f'\r Progress: |{bar}| {percent}% Complete', end='\r')
    # Print New Line on Complete
    if iteration == total:
        print()

# Check command-line arguments
if len(sys.argv) < 2:
    print(" Usage: gpumdkit.sh -filter_dist <file_name> [distance_threshold]")
    print("    or: python filter_structures_by_distance.py <file_name> [distance_threshold]")
    sys.exit(1)

# Read the file name from command line arguments
file_name = sys.argv[1]

# Check if the distance threshold is provided
if len(sys.argv) > 2:
    distance_threshold = float(sys.argv[2])
else:
    distance_threshold = None

# Read all frames from the extxyz file
frames = read(file_name, index=':')

total_frames = len(frames)
filtered_frames = []
filtered_out_frames = []

# Iterate over each frame
for i, frame in enumerate(frames):
    # Get atomic positions
    positions = frame.get_positions()
    # Compute distances between all pairs of atoms
    distances = pdist(positions)
    # Find the minimum distance in the current frame
    min_distance = np.min(distances)
    
    # Check if the minimum distance meets the threshold
    if distance_threshold is None or min_distance >= distance_threshold:
        filtered_frames.append(frame)
    else:
        filtered_out_frames.append(frame)
    
    # Print progress bar
    print_progress_bar(i + 1, total_frames)

filtered_count = len(filtered_frames)
filtered_out_count = total_frames - filtered_count

# Output the filtered frames to a new XYZ file
output_file_name = 'filtered_' + file_name
write(output_file_name, filtered_frames)

# Output the filtered-out frames to a new XYZ file
filtered_out_file_name = 'filtered_out_' + file_name
if filtered_out_count > 0:  # Only write if there are filtered-out structures
   write(filtered_out_file_name, filtered_out_frames)

# Print summary of filtering results
print(f' Total structures processed: {total_frames}')
print(f' Structures filtered out: {filtered_out_count}')
print(f' Structures retained: {filtered_count}')
print(f' Filtered structures saved to {output_file_name}')
if filtered_out_count > 0:
    print(f' Filtered-out structures saved to: {filtered_out_file_name}')
else:
    print(' No structures were filtered out')
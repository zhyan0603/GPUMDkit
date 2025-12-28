"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Filter structures by minimum distance with PBC

Usage:
    python filter_structures_by_distance_pbc.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
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
    if iteration == total:
        print()

# The following function is a part of NepTrainKit (https://github.com/aboys-cb/NepTrainKit)
# Thanks to Chengbing Chen for providing this function.
def calculate_pairwise_distances(lattice_params, atom_coords, fractional=True):
    """
    Calculate the distances between all pairs of atoms in a crystal, considering periodic boundary conditions.

    Parameters:
    lattice_params: Lattice parameters, a 3x3 numpy array representing the lattice vectors (a, b, c).
    atom_coords: Atomic coordinates, an Nx3 numpy array.
    fractional: Whether the coordinates are fractional (True) or Cartesian (False).

    Returns:
    distances: An NxN numpy array containing the distances between all pairs of atoms.
    """
    if fractional:
        atom_coords = np.dot(atom_coords, lattice_params)

    diff = atom_coords[np.newaxis, :, :] - atom_coords[:, np.newaxis, :]
    shifts = np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).T.reshape(-1, 3)
    lattice_shifts = np.dot(shifts, lattice_params)
    all_diffs = diff[:, :, np.newaxis, :] + lattice_shifts[np.newaxis, np.newaxis, :, :]
    all_distances = np.sqrt(np.sum(all_diffs ** 2, axis=-1))
    distances = np.min(all_distances, axis=-1)
    np.fill_diagonal(distances, 0)
    return distances

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
    # Get cell parameters and atomic positions
    lattice_params = frame.get_cell()
    atom_coords = frame.get_positions()
    # Calculate distances using custom function
    distances = calculate_pairwise_distances(lattice_params, atom_coords, fractional=False)
    np.fill_diagonal(distances, float('inf'))
    min_distance = np.min(distances)
    
    # Check if the minimum distance meets the threshold
    if distance_threshold is None or min_distance >= distance_threshold:
        filtered_frames.append(frame)
    else:
        filtered_out_frames.append(frame)

    # Print progress bar
    print_progress_bar(i + 1, total_frames)

filtered_count = len(filtered_frames)
filtered_out_count = len(filtered_out_frames)

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
print(f' Filtered structures saved to: {output_file_name}')
if filtered_out_count > 0:
    print(f' Filtered-out structures saved to: {filtered_out_file_name}')
else:
    print(' No structures were filtered out')
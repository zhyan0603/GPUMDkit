#!/usr/bin/env python3
"""
Get minimum interatomic distances from extended XYZ file with PBC.
Developer: Zihan YAN

Usage: get_min_dist_pbc.py <input.xyz>

Example:
    python get_min_dist_pbc.py train.xyz
    python get_min_dist_pbc.py dump.xyz

Note: This script considers periodic boundary conditions (PBC).
      Use get_min_dist.py for faster calculations without PBC.
"""

import sys
import os
import numpy as np
from ase.io import read


def print_error(message, usage=None, example=None):
    """Print formatted error message to stderr."""
    print(f"Error: {message}", file=sys.stderr)
    
    if usage:
        print(f"Usage: {usage}", file=sys.stderr)
    
    if example:
        print(f"Example: {example}", file=sys.stderr)


def check_file_exists(filepath, description="File"):
    """Check if a file exists."""
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"{description} '{filepath}' does not exist")
    return filepath


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


def validate_arguments():
    """Validate command line arguments."""
    if len(sys.argv) < 2:
        print_error(
            "Missing required argument: input file",
            usage="get_min_dist_pbc.py <input.xyz>",
            example="python get_min_dist_pbc.py train.xyz"
        )
        sys.exit(1)
    
    if sys.argv[1] in ['-h', '--help']:
        print(__doc__)
        sys.exit(0)
    
    input_file = sys.argv[1]
    
    try:
        check_file_exists(input_file, "Input file")
    except FileNotFoundError as e:
        print_error(
            str(e),
            usage="get_min_dist_pbc.py <input.xyz>",
            example="python get_min_dist_pbc.py train.xyz"
        )
        sys.exit(1)
    
    return input_file


def main():
    """Main function to calculate minimum interatomic distances with PBC."""
    # Validate arguments
    file_name = validate_arguments()
    
    # Read all frames from the extxyz file
    try:
        frames = read(file_name, index=':')
    except Exception as e:
        print_error(
            f"Failed to read file '{file_name}': {str(e)}",
            usage="get_min_dist_pbc.py <input.xyz>"
        )
        sys.exit(1)
    
    if len(frames) == 0:
        print_error(
            f"No frames found in '{file_name}'",
            usage="get_min_dist_pbc.py <input.xyz>"
        )
        sys.exit(1)
    
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
        symbols = frame.get_chemical_symbols()
        # Get cell parameters and atomic positions
        lattice_params = frame.get_cell()
        atom_coords = frame.get_positions()
        # Calculate distances using custom function
        distances = calculate_pairwise_distances(lattice_params, atom_coords, fractional=False)
        
        # Calculate minimum distance for each atom pair type
        for i, sym1 in enumerate(unique_symbols):
            for sym2 in unique_symbols[i:]:
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
    print(" Minimum interatomic distances (with PBC):")
    print(" +---------------------------+")
    print(" | Atom Pair |  Distance (A) |")
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
        print(f" Overall min_distance: {overall_min_distance:.3f} A")
    else:
        print(" Overall min_distance: N/A")


if __name__ == "__main__":
    main()

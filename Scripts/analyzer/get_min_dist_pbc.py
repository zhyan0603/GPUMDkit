import sys
import numpy as np
from ase.io import read

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
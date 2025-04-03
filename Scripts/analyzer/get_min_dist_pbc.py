import sys
import numpy as np
from ase.io import read

# Read the file name from command line arguments
file_name = sys.argv[1]

# Read all frames from the extxyz file
frames = read(file_name, index=':')

min_distance = float('inf')

# Iterate over each frame
for frame in frames:
    # Get all distances between atoms considering PBC
    distances = frame.get_all_distances(mic=True)
    # Exclude self-distances (diagonal elements) by setting them to a large value
    np.fill_diagonal(distances, float('inf'))
    # Update the minimum distance
    min_distance = min(min_distance, np.min(distances))

print(f' Minimum interatomic distance: {min_distance:.3f} Å')
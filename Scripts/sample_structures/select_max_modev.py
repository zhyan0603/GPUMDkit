import sys
import numpy as np
from ase.io import read, write

# Define file paths
input_file = 'active.out'  # Input file for max force deviation
#output_file = 'filtered_max_force_deviations.txt'  # Output file for filtered max force deviations
xyz_file = 'active.xyz'  # Input .extxyz file
output_xyz_file = 'selected.xyz'  # Output file for extracted structures

# Number of top structures to extract
top_n = int(sys.argv[1])  # Specify the number of top structures to extract (e.g., top 100)

# Lists to store filtered max force deviations and corresponding indices
filtered_deviations = []
filtered_indices = []

# Read and filter data based on max force deviation
with open(input_file, 'r') as f:
    lines = f.readlines()
    
    for idx, line in enumerate(lines):
        time, max_force_deviation = line.split()
        max_force_deviation = float(max_force_deviation)

        if max_force_deviation > float(sys.argv[2]):
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


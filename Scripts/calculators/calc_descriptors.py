"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Calculate structural descriptors

Usage:
    python calc_descriptors.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

import numpy as np
from ase.io import read
from calorine.nep import get_descriptors
import sys
import os

def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end='\r')
    if iteration == total:
        print()

# Check command-line arguments
if len(sys.argv) != 5:
    print(" Usage: python calc_descriptors.py input.xyz output.npy nep.txt element")
    sys.exit(1)

input_file, output_file, model_file, target_element = sys.argv[1:5]

# Check if input and model files exist
if not os.path.isfile(input_file) or not os.path.isfile(model_file):
    print(" Error: Input file or NEP model file not found.")
    sys.exit(1)

# Read structures from XYZ file
atoms_list = read(input_file, index=':')
if not atoms_list:
    print(" Error: No structures found in input file.")
    sys.exit(1)

# Store descriptors for target element
all_descriptors = []
num_structures = len(atoms_list)

# Process each structure with progress bar
for i, atoms in enumerate(atoms_list, 1):   
    # Get indices of target element
    symbols = np.array(atoms.get_chemical_symbols())
    element_indices = np.where(symbols == target_element)[0]
    
    if len(element_indices) == 0:
        continue
    
    # Calculate descriptors
    descriptors = get_descriptors(atoms, model_filename=model_file)
    # Update progress bar
    print_progress_bar(i, num_structures, prefix=' Processing:', suffix='Complete')
    
    # Extract descriptors for target element
    selected_descriptors = descriptors[element_indices, :]

    # if you want to output the mean results
    # Uncomment the next two lines and comment the original 'all_descriptors'

    # mean_descriptors = np.mean(selected_descriptors, axis=0)
    # all_descriptors.append(mean_descriptors)

    all_descriptors.append(selected_descriptors)

# Check if any descriptors were collected
if not all_descriptors:
    print(f" Error: No {target_element} atoms found in any structure.")
    sys.exit(1)

# Combine and save descriptors
all_descriptors = np.concatenate(all_descriptors, axis=0)
np.save(output_file, all_descriptors)

# Print output summary
print(f" Saved descriptors to '{output_file}'")
print(f" Output shape: {all_descriptors.shape} ({all_descriptors.shape[0]} {target_element} atoms)")
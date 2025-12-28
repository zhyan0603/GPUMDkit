"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Compare structural descriptors

Usage:
    python plt_descriptors_compare.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

import numpy as np
from ase.io import read
from calorine.nep import get_descriptors
import sys
import os
import matplotlib.pyplot as plt


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end='\r')
    if iteration == total:
        print()

# Check command-line arguments
if len(sys.argv) < 5:
    print(" Usage: python calc_descriptors.py method nep.txt element train1.xyz [train2.xyz ...]")
    print(" method: 'pca' or 'umap'")
    print(" element: chemical symbol (e.g., Li) or 'total'")
    sys.exit(1)

method, model_file, target_element = sys.argv[1:4]
xyz_files = sys.argv[4:]

# Validate method
if method.lower() not in ['pca', 'umap']:
    print(" Error: Method must be 'pca' or 'umap'.")
    sys.exit(1)

# Check if model file exists
if not os.path.isfile(model_file):
    print(f" Error: NEP model file '{model_file}' not found.")
    sys.exit(1)

# Process each XYZ file
all_descriptors = []
file_labels = []
descriptors_by_file = []

for xyz_file in xyz_files:
    # Check if XYZ file exists
    if not os.path.isfile(xyz_file):
        print(f" Error: XYZ file '{xyz_file}' not found.")
        sys.exit(1)
    
    # Read structures
    print(f" Processing '{xyz_file}'...")
    atoms_list = read(xyz_file, index=':')
    if not atoms_list:
        print(f" Error: No structures found in '{xyz_file}'.")
        sys.exit(1)
    
    # Compute mean descriptors for each structure
    file_descriptors = []
    num_structures = len(atoms_list)
    
    for i, atoms in enumerate(atoms_list, 1):
        # Calculate descriptors
        descriptors = get_descriptors(atoms, model_filename=model_file)
        
        # Update progress bar
        print_progress_bar(i, num_structures, prefix=f'  Structures in {os.path.basename(xyz_file)}:', suffix='Complete')
        
        # Select descriptors
        if target_element.lower() == 'total':
            selected_descriptors = descriptors
        else:
            symbols = np.array(atoms.get_chemical_symbols())
            element_indices = np.where(symbols == target_element)[0]
            if len(element_indices) == 0:
                continue
            selected_descriptors = descriptors[element_indices, :]
        
        # Compute mean descriptors
        mean_descriptors = np.mean(selected_descriptors, axis=0)
        file_descriptors.append(mean_descriptors)
    
    if not file_descriptors:
        print(f" Warning: No descriptors collected from '{xyz_file}' (check if '{target_element}' exists).")
        continue
    
    # Store descriptors and label
    file_descriptors = np.array(file_descriptors)
    all_descriptors.append(file_descriptors)
    descriptors_by_file.append(file_descriptors)
    file_labels.append(os.path.splitext(os.path.basename(xyz_file))[0])
    
    # Save descriptors
    # output_file = f"{os.path.splitext(xyz_file)[0]}_descriptors.npy"
    # np.save(output_file, file_descriptors)
    # print(f"Saved descriptors to '{output_file}' (shape: {file_descriptors.shape})")

# Check if any descriptors were collected
if not all_descriptors:
    print(f" Error: No descriptors collected for '{target_element}' across all files.")
    sys.exit(1)

# Combine all descriptors for reduction
combined_descriptors = np.concatenate(all_descriptors, axis=0)

# Apply dimensionality reduction
print(f"Applying {method.upper()} reduction...")
if method.lower() == 'pca':
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    scaled_descriptors = scaler.fit_transform(combined_descriptors)
    reducer = PCA(n_components=2)
    reduced = reducer.fit_transform(scaled_descriptors)
elif method.lower() == 'umap':
    import umap
    reducer = umap.UMAP(n_components=2, random_state=42)
    reduced = reducer.fit_transform(combined_descriptors)
else:  # umap
    print("Error: Method must be 'pca' or 'umap'.")
    sys.exit(1)

# Split reduced descriptors back by file
reduced_by_file = []
start = 0
for descriptors in all_descriptors:
    end = start + descriptors.shape[0]
    reduced_by_file.append(reduced[start:end])
    start = end

# Plot reduced descriptors
plt.figure(figsize=(6, 4.5), dpi=100)
colors = ['#9BBBE1', '#EAB883', '#A9CA70', '#DD7C4F', '#F09BA0', '#B58C9A'] 
for reduced, label, color in zip(reduced_by_file, file_labels, colors):
    plt.scatter(reduced[:, 0], reduced[:, 1], label=label, color=color, alpha=0.6, s=15)

#plt.title(f'{method.upper()} Analysis')
plt.xlabel(f'{method.upper()} Component 1')
plt.ylabel(f'{method.upper()} Component 2')
plt.legend()
plt.tight_layout()

from matplotlib import get_backend
if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
    print("Unable to display the plot due to the non-interactive backend.")
    print("The plot has been automatically saved as 'descriptors.png'.")
    plt.savefig('descriptors.png', dpi=300)
else:
    plt.show()
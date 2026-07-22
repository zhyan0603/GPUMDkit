"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_charge.py
Category:   Plot Scripts
Purpose:    Plot charge distribution from charge_train.out for qNEP model
            analysis.
Usage:      gpumdkit.sh -plt charge [save]
            python plt_charge.py [save]
Arguments:
  save      Save the plot as 'charge.png' instead of displaying it
Input files (required in working directory):
  train.xyz         Training structures in extxyz format
  charge_train.out  Charge data from qNEP model
Output:
  charge.png   Charge distribution plot (saved or displayed)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ase.io import read
from collections import defaultdict

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
})

def load_xyz_and_charges(xyz_file, charge_file):
    """
    Load atomic charges from charge.out and match with atom symbols from an extxyz file.
    
    Args:
        xyz_file (str): Path to the extxyz file (e.g., train.xyz).
        charge_file (str): Path to the charge file (e.g., charge.out).
    
    Returns:
        dict: Dictionary with element symbols as keys and lists of corresponding charges as values.
    """
    # Read extxyz file
    structures = read(xyz_file, index=':', format='extxyz')
    
    # Read charge file
    charges = np.loadtxt(charge_file, ndmin=1)
    
    # Verify that the number of charges matches the total number of atoms
    total_atoms = sum(len(struct) for struct in structures)
    if len(charges) != total_atoms:
        raise ValueError(f"Charge data length ({len(charges)}) does not match total atoms in XYZ ({total_atoms})")
    
    # Collect charges by element
    charges_by_element = defaultdict(list)
    charge_index = 0
    for struct in structures:
        symbols = struct.get_chemical_symbols()  # Get element symbols for the structure
        for symbol in symbols:
            charges_by_element[symbol].append(charges[charge_index])
            charge_index += 1
    
    return charges_by_element

# File paths
xyz_file = "train.xyz"
charge_file = "charge_train.out"

# Load data and plot
charges_by_element = load_xyz_and_charges(xyz_file, charge_file)

plt.figure(figsize=(5, 3.5), dpi=150)
# Plot charge distribution for each element
for element, charge_values in charges_by_element.items():
    sns.kdeplot(charge_values, fill=True, label=element, alpha=0.6)

plt.xlabel("Charge (e)")
plt.ylabel("Density")
plt.legend()
#plt.grid()
plt.tight_layout()

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('charge.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'charge.png'.")
        plt.savefig('charge.png', dpi=300)
    else:
        plt.show()
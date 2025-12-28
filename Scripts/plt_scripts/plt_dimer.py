"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Plot dimer potential curve

Usage:
    python plt_dimer.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from calorine.calculators import CPUNEP
import sys

# Check the number of command-line arguments
if len(sys.argv) < 4:
    print("Usage: python plt_dimer.py <atom_symbol1> <atom_symbol2> <nep_dir>")
    print("Example: python plt_dimer.py Li Li ./nep.txt")
    sys.exit(1)

# Get command-line inputs
symbol1 = sys.argv[1]  # Type of the first atom
symbol2 = sys.argv[2]  # Type of the second atom
nep_path = sys.argv[3]  # Directory of the NEP potential file

# Set the range and step size for distances
start_distance = 0.1  # Starting distance in Ångström
end_distance = 6.0    # Ending distance in Ångström
step_size = 0.01      # Step size in Ångström

# Initialize arrays for results
distances = np.arange(start_distance, end_distance + step_size, step_size)
energies = []
forces = []

# Create CPUNEP calculator
calc = CPUNEP(nep_path)

# Define the unit cell
cell = [[30, 0, 0],
        [0, 30, 0],
        [0, 0, 30]]

# Calculate energy and force for each distance
for distance in distances:
    # Create a structure with two atoms and set the unit cell
    structure = Atoms([symbol1, symbol2], 
                      positions=[[0, 0, 0], [distance, 0, 0]], 
                      cell=cell, 
                      pbc=[True, True, True])
    structure.calc = calc

    # Calculate potential energy
    energy = structure.get_potential_energy()
    energies.append(energy)

    # Calculate force (x-component on the second atom)
    force = structure.get_forces()[1, 0]
    forces.append(force)

# Convert lists to NumPy arrays
energies = np.array(energies)
forces = np.array(forces)

# Set reference energy as the last energy value and shift energies
reference_energy = energies[-1]  # Use the energy at the largest distance as reference
energies_shifted = energies - reference_energy

# Save data to file after calculations
with open(f'dimer_{symbol1}_{symbol2}.txt', 'w') as f:
    # Write the header
    f.write('Distance (Å)  Energy (eV)  Force (eV/Å)\n')
    
    # Write shifted energy and force data
    for i, distance in enumerate(distances):
        f.write(f"{distance:.2f} {energies_shifted[i]:.6f} {forces[i]:.6f}\n")

# Create two subplots (top and bottom)
fig, axs = plt.subplots(2, 1, figsize=(6, 6), dpi=150)

# Plot energy difference vs. distance
axs[0].plot(distances, energies_shifted, marker='o')
axs[0].set_xlabel('Dimer Distance (Å)')
axs[0].set_ylabel(r'$\Delta$E (eV)')
title = f"{symbol1}-{symbol2} Dimer Interaction"
axs[0].set_title(title)

# Plot force vs. distance
axs[1].plot(distances, forces, marker='o', color='C1')
axs[1].set_xlabel('Dimer Distance (Å)')
axs[1].set_ylabel('Fx (eV/Å)')

# Adjust layout and display the plot
plt.tight_layout()

if len(sys.argv) > 4 and sys.argv[4] == 'save':
    plt.savefig('dimer.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'dimer.png'.")
        plt.savefig(f'dimer-{symbol1}-{symbol2}.png', dpi=300)
    else:
        plt.show()
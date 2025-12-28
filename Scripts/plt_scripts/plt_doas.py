"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Plot density of atomic states

Usage:
    python plt_doas.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Specify the input file and element
input_file = sys.argv[1]  # Input file (e.g., doas.out)
element = sys.argv[2]     # Element to plot (e.g., Li)

# Function to read energies for a specific element
def read_energies(filename, element):
    energies = []
    with open(filename, 'r') as f:
        in_section = False
        for line in f:
            line = line.strip()
            if line == f"# {element}":
                in_section = True
                continue
            if line == "" and in_section:
                break
            if in_section and line:
                energies.append(float(line))
    if not energies:
        raise ValueError(f"No energies found for element {element} in {filename}")
    return np.array(energies)

# Read energies for the specified element
energies = read_energies(input_file, element)

# Set up the plot style
plt.figure(figsize=(6, 4))

# Plot KDE
sns.kdeplot(energies, label=f'{element}', color='C0', linewidth=2, alpha=0.2, fill=True)

# Customize the plot
plt.xlabel('Atomistic Energy (eV)')
plt.ylabel('Density')
plt.title(f'Density of Atomistic States of {element}')
plt.legend()
plt.tight_layout()

if len(sys.argv) > 3 and sys.argv[3] == 'save':
    plt.savefig('doas.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'doas.png'.")
        plt.savefig('doas.png', dpi=300)
    else:
        plt.show()
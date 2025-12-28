"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    plt_msd_all.py
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# Function to read data from the given file
def read_data(file_name):
    data = np.loadtxt(file_name)
    num_columns = data.shape[1]
    num_groups = (num_columns - 1) // 6  # Each group has 6 columns (3 MSD + 3 SDC)
    time = data[:, 0]
    msd_data = [data[:, 1 + i*6:4 + i*6] for i in range(num_groups)]  # Extract MSD (x, y, z) for each group
    return time, msd_data, num_groups

# Determine the input file and elements to plot
input_file = sys.argv[1] 
elements_to_plot = sys.argv[2:] if len(sys.argv) > 2 else None

# Read the data from the file
time, msd_data, num_groups = read_data(input_file)

# If no elements specified, default to all groups with generic labels
if elements_to_plot is None:
    elements_to_plot = [f'Group {i}' for i in range(num_groups)]
    groups_to_plot = list(range(num_groups))
else:
    # Map elements to group indices (assume order matches file)
    if len(elements_to_plot) > num_groups:
        print(f"Error: Specified {len(elements_to_plot)} elements, but file contains only {num_groups} groups.")
        sys.exit(1)
    groups_to_plot = list(range(len(elements_to_plot)))

# Create a plot for the MSD data
plt.figure(figsize=(4.5, 3.5), dpi=150)
colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']  # Color cycle for groups

# Plot total MSD for each specified group
for group_idx, element in zip(groups_to_plot, elements_to_plot):
    msd_x, msd_y, msd_z = msd_data[group_idx].T
    msd_total = msd_x + msd_y + msd_z
    color = colors[group_idx % len(colors)]
    plt.plot(time, msd_total, label=f'{element} total-MSD', color=color)

# Add titles and labels
plt.xlabel('dt (ps)')
plt.ylabel(r'MSD ($\AA^2$)')
plt.legend()
plt.tight_layout()

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('msd_all.png', dpi=300)
else:
    # Handle saving or displaying the plot
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'msd_all.png'.")
        plt.savefig('msd_all.png', dpi=300)
    else:
        plt.show()
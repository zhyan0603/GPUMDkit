"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Check MSD convergence

Usage:
    python plt_msd_convergence_check.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Function to read data from the given file
def read_data(file_name):
    data = np.loadtxt(file_name)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# Function to calculate the slope of the 40% to 80% portion of the data
def calculate_slope(x, y):
    len_data = len(x)
    start_idx = int(0.4 * len_data)  # 40% of the data
    end_idx = int(0.8 * len_data)    # 80% of the data
    x_range, y_range = x[start_idx:end_idx], y[start_idx:end_idx]
    coeffs = np.polyfit(x_range, y_range, 1)  # Linear fit
    return coeffs[0]  # Slope of the linear fit

# Function to read time_step from run.in file
def get_time_step():
    default_time_step = 1.0  # Default: 1 fs
    try:
        with open('run.in', 'r') as file:
            for line in file:
                if 'time_step' in line.lower():
                    # Extract the number after 'time_step'
                    parts = line.split()
                    for part in parts:
                        try:
                            return float(part)
                        except ValueError:
                            continue
    except FileNotFoundError:
        pass
    return default_time_step

# Find all msd_step*.out files and extract step numbers
file_pattern = './msd_step*.out'
files = glob.glob(file_pattern)
steps = []
for file in files:
    # Extract step number from filename (between 'msd_step' and '.out')
    step_str = os.path.basename(file).replace('msd_step', '').replace('.out', '')
    try:
        steps.append(int(step_str))
    except ValueError:
        continue
steps = np.sort(steps)  # Sort steps in ascending order

# Convert steps to ns using time_step (1 step = time_step fs)
time_step = get_time_step()  # Get time_step in fs
time_points = steps * time_step / 1e6  # Convert steps to ns (1 fs = 1e-6 ns)

# Initialize lists to store slopes for each direction and total
slopes_x = []
slopes_y = []
slopes_z = []
slopes_total = []

# Process each file
for step in steps:
    input_file = f'./msd_step{step}.out'
    try:
        # Read the data from the file
        time, msd_x, msd_y, msd_z = read_data(input_file)
        msd_total = msd_x + msd_y + msd_z  # Calculate total MSD
        
        # Calculate slopes for each direction and total
        slopes_x.append(calculate_slope(time, msd_x))
        slopes_y.append(calculate_slope(time, msd_y))
        slopes_z.append(calculate_slope(time, msd_z))
        slopes_total.append(calculate_slope(time, msd_total))
        
    except FileNotFoundError:
        print(f"File {input_file} not found, skipping.")
        slopes_x.append(np.nan)
        slopes_y.append(np.nan)
        slopes_z.append(np.nan)
        slopes_total.append(np.nan)

# Define custom colors
color_total = '#1f77b4'  # Blue for total MSD
color_x = '#ff7f0e'      # Orange for MSD X
color_y = '#2ca02c'      # Green for MSD Y
color_z = '#d62728'      # Red for MSD Z

# Create the plot with two y-axes
fig, ax1 = plt.subplots(figsize=(5.5, 4), dpi=150)

# Plot total MSD on the left y-axis
ax1.plot(time_points, slopes_total, marker='o', linestyle='-', color=color_total, 
         fillstyle='none', label='Total')
ax1.set_xlabel('Simulation Time (ns)')
ax1.set_ylabel(r'Total Diffusion Rate ($\AA^2$/ps)', color=color_total)
ax1.tick_params(axis='y', labelcolor=color_total)

# Create second y-axis for individual components
ax2 = ax1.twinx()
ax2.plot(time_points, slopes_x, marker='o', linestyle='-', color=color_x, 
         fillstyle='none', label='x')
ax2.plot(time_points, slopes_y, marker='o', linestyle='-', color=color_y, 
         fillstyle='none', label='y')
ax2.plot(time_points, slopes_z, marker='o', linestyle='-', color=color_z, 
         fillstyle='none', label='z')
ax2.set_ylabel(r'Component Diffusion Rate ($\AA^2$/ps)', color='gray')
ax2.tick_params(axis='y', labelcolor='gray')

# Combine legends from both axes
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')

plt.tight_layout()

# Save or show the plot
if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('msd_convergence.png', dpi=300)
else:
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been saved as 'msd_convergence.png'.")
        plt.savefig('msd_convergence.png', dpi=300)
    else:
        plt.show()
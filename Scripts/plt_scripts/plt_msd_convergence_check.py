import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Function to read data from the given file
def read_data(file_name):
    data = np.loadtxt(file_name)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# Function to calculate diffusion coefficient from the 40% to 80% portion of the data
def calculate_diffusion_coefficient(x, y, dimension='total'):
    len_data = len(x)
    start_idx = int(0.4 * len_data)
    end_idx = int(0.8 * len_data)
    x_range = x[start_idx:end_idx]
    y_range = y[start_idx:end_idx]
    
    if len(x_range) < 2:
        return np.nan
    
    coeffs = np.polyfit(x_range, y_range, 1)
    slope = coeffs[0]
    
    if dimension in ['x', 'y', 'z']:
        return slope / 2.0
    else:  # total
        return slope / 6.0

# Function to read time_step from run.in file
def get_time_step():
    default_time_step = 1.0  # Default: 1 fs
    try:
        with open('run.in', 'r') as file:
            for line in file:
                if 'time_step' in line.lower():
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
    step_str = os.path.basename(file).replace('msd_step', '').replace('.out', '')
    try:
        steps.append(int(step_str))
    except ValueError:
        continue
steps = np.sort(steps)

# Convert steps to ns using time_step (1 step = time_step fs)
time_step = get_time_step()
time_points = steps * time_step / 1e6  # steps → ns

# Lists to store diffusion coefficients
D_x_list = []
D_y_list = []
D_z_list = []
D_total_list = []

# Process each file
for step in steps:
    input_file = f'./msd_step{step}.out'
    try:
        time, msd_x, msd_y, msd_z = read_data(input_file)
        msd_total = msd_x + msd_y + msd_z
        
        D_x = calculate_diffusion_coefficient(time, msd_x, 'x')
        D_y = calculate_diffusion_coefficient(time, msd_y, 'y')
        D_z = calculate_diffusion_coefficient(time, msd_z, 'z')
        D_total = calculate_diffusion_coefficient(time, msd_total, 'total')
        
        D_x_list.append(D_x)
        D_y_list.append(D_y)
        D_z_list.append(D_z)
        D_total_list.append(D_total)
        
    except Exception:
        print(f" File {input_file} not found or invalid, skipping.")
        D_x_list.append(np.nan)
        D_y_list.append(np.nan)
        D_z_list.append(np.nan)
        D_total_list.append(np.nan)

# Define custom colors
color_total = '#1f77b4'  # Blue for total
color_x     = '#ff7f0e'  # Orange for x
color_y     = '#2ca02c'  # Green for y
color_z     = '#d62728'  # Red for z

# Create the plot (single y-axis recommended)
fig, ax = plt.subplots(figsize=(4.5, 3.6), dpi=150)

ax.plot(time_points, D_total_list, marker='o', linestyle='-', linewidth=2.0,
        color=color_total, label='Total', zorder=10)

ax.plot(time_points, D_x_list, marker='s', linestyle='-', linewidth=1.4,
        color=color_x, label='x', alpha=0.9)
ax.plot(time_points, D_y_list, marker='^', linestyle='-', linewidth=1.4,
        color=color_y, label='y', alpha=0.9)
ax.plot(time_points, D_z_list, marker='D', linestyle='-', linewidth=1.4,
        color=color_z, label='z', alpha=0.9)

ax.set_xlabel('Simulation Time (ns)')
ax.set_ylabel(r'Diffusion Coefficient ($\AA^2$/ps)')

ax.grid(True, alpha=0.25, linestyle='--')
ax.legend(loc='best', frameon=True)

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
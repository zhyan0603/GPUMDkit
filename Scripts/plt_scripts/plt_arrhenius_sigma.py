"""
This script processes GPUMD outputs to plot the Arrhenius relationship for conductivity (sigma) as a function of temperature. 
It reads MSD and thermo data from specified temperature folders, calculates diffusivity and conductivity, and fits the data 
to extract activation energy. The plot follows PRL journal style with professional academic formatting.

Author: Modified from original by Zihan YAN (yanzihan@westlake.edu.cn)
"""

import re
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import get_backend
from matplotlib.lines import Line2D
import scipy.constants as consts
from ase.io import read
from scipy.stats import linregress
import matplotlib.ticker as ticker
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],  
    "font.size": 11,
    "axes.labelsize": 11.5,
    "axes.titlesize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.figsize": (4.3, 3.8),  
    "axes.spines.top": True,
    "axes.spines.right": True,
    "axes.linewidth": 0.8,
})

k = 8.617333262145e-5

def fit_arrhenius(T, sigma_T):
    """Fit Arrhenius equation to calculate activation energy"""
    mask = ~np.isnan(sigma_T) & ~np.isnan(T) & (sigma_T > 0)
    T_valid, val_valid = T[mask], sigma_T[mask]
    if len(T_valid) < 2: 
        return None, None, None, None
    
    x = 1000 / T_valid
    y = np.log10(val_valid) 
    res = linregress(x, y)
    Ea = -res.slope * k * 1000 * np.log(10)  
    return Ea, res.intercept, res.slope, res.rvalue

def tick_function(X):
    """Convert inverse temperature back to temperature for top axis"""
    V = 1000 / X  
    return ["%.0f" % z for z in V]

# Get current directory and temperature folders (e.g., 400K, 500K)
base_dir = os.getcwd()
temp_folders = [f for f in os.listdir(base_dir) if re.match(r'\d+K$', f) and os.path.isdir(f)]
temp_folders = sorted(temp_folders, key=lambda x: int(x[:-1]))  # Sort by temperature

# Optional: Uncomment to manually specify temperatures
# temp_list = [400, 500, 600, 700]
# temp_folders = [f"{t}K" for t in temp_list if os.path.isdir(os.path.join(base_dir, f"{t}K"))]

if not temp_folders:
    raise ValueError("No valid temperature folders (e.g., 400K, 500K) found in current directory")

# Define colors following PRL style
default_colors = ['#457B9D', '#D62828', '#2A9D8F', '#E9C46A']
color = default_colors[0]  # Using first PRL-style color

# Get structure label from current directory name
cell = os.path.basename(base_dir)

# Read model.xyz from the first temperature folder
xyz_file = os.path.join(base_dir, temp_folders[0], 'model.xyz')
if not os.path.exists(xyz_file):
    raise FileNotFoundError(f"Missing {temp_folders[0]}/model.xyz file")

# Count Li or Na ions
atoms = read(xyz_file, format='extxyz')
num_ions = sum(1 for atom in atoms if atom.symbol in ['Li', 'Na'])

# Check for replicate parameter in run.in
run_in_path = os.path.join(base_dir, 'run.in')
rep_factor = 1
if os.path.exists(run_in_path):
    with open(run_in_path, 'r') as f:
        for line in f:
            if line.strip().startswith('replicate'):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        rep_x = int(parts[1])
                        rep_y = int(parts[2])
                        rep_z = int(parts[3])
                        rep_factor = rep_x * rep_y * rep_z
                        print(f"[Info] Detected replicate {rep_x} {rep_y} {rep_z}, multiplying ion count by {rep_factor}")
                        num_ions *= rep_factor
                        break
                    except ValueError:
                        print("[Warning] Invalid replicate format, assuming no replication")
else:
    print("[Note] No run.in file found, assuming no replication")

# Initialize lists for plotting and data storage
raw_group_Ts = []
raw_group_sigmaTs = []
raw_group_sigmas = []

# Process each temperature folder
for temp_folder in temp_folders:
    temp = int(temp_folder[:-1])  # Extract temperature value
    folder_path = os.path.join(base_dir, temp_folder)
    thermo_file = os.path.join(folder_path, 'thermo.out')
    msd_file = os.path.join(folder_path, 'msd.out')

    if not os.path.exists(thermo_file) or not os.path.exists(msd_file):
        print(f"[Warning] Missing thermo.out or msd.out in {temp_folder}, skipping")
        continue

    try:
        # Extract volume from thermo.out
        data = np.loadtxt(thermo_file)
        num_columns = data.shape[1]
        if num_columns == 12:
            vol = data[:, 9] * data[:, 10] * data[:, 11]
        elif num_columns == 18:
            a = np.column_stack((data[:, 9], data[:, 10], data[:, 11]))
            b = np.column_stack((data[:, 12], data[:, 13], data[:, 14]))
            c = np.column_stack((data[:, 15], data[:, 16], data[:, 17]))
            vol = np.abs(np.einsum('ij,ij->i', a, np.cross(b, c)))
        else:
            raise ValueError("Unsupported number of columns in thermo.out")
        volume = np.mean(vol)

        # Read MSD data
        data = np.loadtxt(msd_file)
        dts = data[:, 0]
        msds = np.sum(data[:, 1:4], axis=1)
        start = int(len(dts) * 0.4)
        end = int(len(dts) * 0.8)

        # Calculate diffusivity
        k_diff, _ = np.polyfit(dts[start:end], msds[start:end], 1)
        D = k_diff * 1e-4 / 6  # Convert A/ps to cm/s

        # Calculate conductivity
        z = 1  # Species charge
        vol_cm3 = volume * 1e-24  # Convert volume to cm^3
        conversion_factor = (num_ions / (vol_cm3 * consts.N_A) * z**2 * (consts.N_A * consts.e)**2 / (consts.R * temp))
        conductivity = conversion_factor * D
        sigma_T = conductivity * temp

        if sigma_T <= 0 or np.isnan(sigma_T):
            print(f"[Skip] {temp_folder}: Invalid sigma_T = {sigma_T}")
            continue

        # Store data
        raw_group_Ts.append(temp)
        raw_group_sigmaTs.append(sigma_T)
        raw_group_sigmas.append(conductivity)

    except Exception as e:
        print(f"[Error] Processing {temp_folder}: {e}")
        continue

# Sort by temperature
if len(raw_group_Ts) > 0:
    sorted_indices = np.argsort(raw_group_Ts)
    group_Ts = np.array(raw_group_Ts)[sorted_indices]
    group_sigmaTs = np.array(raw_group_sigmaTs)[sorted_indices]
    group_sigmas = np.array(raw_group_sigmas)[sorted_indices]

    # Print header for conductivity data
    w_t, w_sigma, w_sigmaT = 10, 16, 18
    line = f"+{'-'*(w_t+2)}-{'-'*(w_sigma+2)}-{'-'*(w_sigmaT+2)}+"
    print(line)
    print(f"| {'Conductivity (unit: S/cm, Sigma*T: K S/cm)':^{len(line)-4}} |")
    print(line)
    print(f"| {'T (K)':^{w_t}} | {'Sigma':^{w_sigma}} | {'Sigma*T':^{w_sigmaT}} |")
    print(line)

    # Print conductivity data in sorted order
    for temp, conductivity, sigma_T in zip(group_Ts, group_sigmas, group_sigmaTs):
        print(f"| {temp:^{w_t}d} | {conductivity:^{w_sigma}.3e} | {sigma_T:^{w_sigmaT}.3e} |")

    print(line)

    if len(group_Ts) < 2:
        raise ValueError("Insufficient data points for fitting, need at least two temperatures")

    # Create figure with PRL style
    fig, ax = plt.subplots(figsize=(4.3, 3.8))

    # Calculate activation energy using the new function
    Ea, intercept, slope, r_value = fit_arrhenius(group_Ts, group_sigmaTs)

    # Prepare data for plotting
    x_fit_range = np.linspace(1000/np.max(group_Ts), 1000/np.min(group_Ts), 100)
    y_fit_ln = (slope * x_fit_range + intercept) * np.log(10)

    # Plot fitted line
    ax.plot(x_fit_range, y_fit_ln, '--', c=color, linewidth=1.5)

    # Plot data points
    ax.plot(1000/group_Ts, np.log(group_sigmaTs), 'o', label=f'{cell} ({Ea:.3f} eV)', 
            c=color, markersize=8, markerfacecolor='none', markeredgewidth=1.5)

    # Set axis labels with proper formatting
    ax.set_xlabel('1000/T (1/K)', labelpad=7)
    ax.set_ylabel(r'ln($\sigma$T) (S$\cdot$K/cm)')

    # Add legend
    ax.legend(loc='lower left', frameon=False, fontsize=11)

    # Add secondary x-axis on top showing actual temperature
    ax_top = ax.secondary_xaxis('top')
    ax_top.set_xlabel('Temperature (K)', labelpad=7)
    xticks = ax.get_xticks()
    ax_top.set_xticks(xticks)
    ax_top.set_xticklabels([f'{int(1000/xt)}' if xt>0 else '' for xt in xticks])

    # Format y-axis
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

    # Print results
    print(f"\n{cell}, Ea: {Ea:.3f} eV")

    # Calculate conductivity at 300K
    target_T = 300
    sig_300 = (10**(slope * (1000/target_T) + intercept)) / target_T
    print(f"at {target_T}K, {cell}: Sigma = {sig_300:.3e} S/cm")

    # Generate temperature and conductivity arrays for export
    temp_array = group_Ts.tolist()
    conductivity_array = [float(f'{val:.3e}') for val in group_sigmas.tolist()]
    
    print("\nExportable arrays for Python:")
    print(f"temperatures = {temp_array}")
    print(f"conductivity_values = {conductivity_array}")

    # Adjust layout and save
    plt.tight_layout()

    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        plt.savefig('Arrhenius_sigma.png', dpi=300, bbox_inches='tight')
    else:
        if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
            print("Unable to display the plot due to the non-interactive backend.")
            print("The plot has been automatically saved as 'Arrhenius_sigma.png'.")
            plt.savefig('Arrhenius_sigma.png', dpi=300, bbox_inches='tight')
        else:
            plt.show()

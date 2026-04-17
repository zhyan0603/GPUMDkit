"""
This script processes MSD and thermo data from multiple temperature folders (e.g., 400K, 500K) 
to calculate ionic conductivity and its temperature dependence for different directions. 
It extracts volume from thermo.out, calculates diffusivity from msd.out for x, y, z directions, 
computes conductivity using the Nernst-Einstein relation, 
and plots ln(sigma*T) vs 1000/T to determine activation energy for each direction. 
The script also extrapolates conductivity to 300K based on the fitted Arrhenius behavior.

Author: Modified from original by Zihan YAN (yanzihan@westlake.edu.cn)
"""

import re
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import get_backend
import scipy.constants as consts
from ase.io import read
from scipy.stats import linregress
import matplotlib.ticker as ticker
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
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
default_colors = ['#2A9D8F', '#D62828', '#457B9D', '#E9C46A']

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
group_Ts = []
group_sigmaTs_all = []
group_sigmaTs_x = []
group_sigmaTs_y = []
group_sigmaTs_z = []
group_sigmas_all = []
group_sigmas_x = []
group_sigmas_y = []
group_sigmas_z = []

# Create figure with PRL style
fig, ax = plt.subplots(figsize=(4.3, 3.8))

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
        msd_x = data[:, 1]
        msd_y = data[:, 2]
        msd_z = data[:, 3]
        msd_all = msd_x + msd_y + msd_z
        start = int(len(dts) * 0.4)
        end = int(len(dts) * 0.8)

        # Calculate diffusivity
        k_all, _ = np.polyfit(dts[start:end], msd_all[start:end], 1)
        k_x, _ = np.polyfit(dts[start:end], msd_x[start:end], 1)
        k_y, _ = np.polyfit(dts[start:end], msd_y[start:end], 1)
        k_z, _ = np.polyfit(dts[start:end], msd_z[start:end], 1)

        D_all = k_all * 1e-4 / 6  # Convert A/ps to cm/s
        D_x = k_x * 1e-4 / 2
        D_y = k_y * 1e-4 / 2
        D_z = k_z * 1e-4 / 2

        # Calculate conductivity
        z = 1  # Species charge
        vol_cm3 = volume * 1e-24  # Convert volume to cm^3
        conversion_factor = (num_ions / (vol_cm3 * consts.N_A) * z**2 * (consts.N_A * consts.e)**2 / (consts.R * temp))
        sigma_all = conversion_factor * D_all
        sigma_x = conversion_factor * D_x
        sigma_y = conversion_factor * D_y
        sigma_z = conversion_factor * D_z
        sigmaT_all = sigma_all * temp
        sigmaT_x = sigma_x * temp
        sigmaT_y = sigma_y * temp
        sigmaT_z = sigma_z * temp

        if any(v <= 0 or np.isnan(v) for v in [sigmaT_all, sigmaT_x, sigmaT_y, sigmaT_z]):
            print(f"[Skip] {temp_folder}: Invalid sigma_T values")
            continue

        group_Ts.append(temp)
        group_sigmaTs_all.append(sigmaT_all)
        group_sigmaTs_x.append(sigmaT_x)
        group_sigmaTs_y.append(sigmaT_y)
        group_sigmaTs_z.append(sigmaT_z)
        group_sigmas_all.append(sigma_all)
        group_sigmas_x.append(sigma_x)
        group_sigmas_y.append(sigma_y)
        group_sigmas_z.append(sigma_z)

    except Exception as e:
        print(f"[Error] Processing {temp_folder}: {e}")
        continue

if len(group_Ts) < 2:
    raise ValueError("Insufficient data points for fitting, need at least two temperatures")

# Sort by temperature
sorted_indices = np.argsort(group_Ts)
group_Ts = np.array(group_Ts)[sorted_indices]
group_sigmaTs_all = np.array(group_sigmaTs_all)[sorted_indices]
group_sigmaTs_x = np.array(group_sigmaTs_x)[sorted_indices]
group_sigmaTs_y = np.array(group_sigmaTs_y)[sorted_indices]
group_sigmaTs_z = np.array(group_sigmaTs_z)[sorted_indices]
group_sigmas_all = np.array(group_sigmas_all)[sorted_indices]
group_sigmas_x = np.array(group_sigmas_x)[sorted_indices]
group_sigmas_y = np.array(group_sigmas_y)[sorted_indices]
group_sigmas_z = np.array(group_sigmas_z)[sorted_indices]

# Prepare data for plotting
x_min, x_max = 1000/np.max(group_Ts), 1000/np.min(group_Ts)
x_fit_range = np.linspace(x_min, x_max, 100)

# Calculate activation energies for each direction
Ea_all, intercept_all, slope_all, r_all = fit_arrhenius(group_Ts, group_sigmaTs_all)
Ea_x, intercept_x, slope_x, r_x = fit_arrhenius(group_Ts, group_sigmaTs_x)
Ea_y, intercept_y, slope_y, r_y = fit_arrhenius(group_Ts, group_sigmaTs_y)
Ea_z, intercept_z, slope_z, r_z = fit_arrhenius(group_Ts, group_sigmaTs_z)

# Prepare labels with activation energies
labels_with_ea = [
    (f'Total ({Ea_all:.3f} eV)', group_sigmaTs_all, intercept_all, slope_all, '#457B9D', 'o'),
    (f'X ({Ea_x:.3f} eV)', group_sigmaTs_x, intercept_x, slope_x, '#D62828', 's'),
    (f'Y ({Ea_y:.3f} eV)', group_sigmaTs_y, intercept_y, slope_y, '#2A9D8F', '^'),
    (f'Z ({Ea_z:.3f} eV)', group_sigmaTs_z, intercept_z, slope_z, '#E9C46A', 'd')
]

# Print conductivity tables
w_t, w_all, w_xyz = 10, 16, 12
line = f"+{'-'*(w_t+2)}-{'-'*(w_all+2)}-{'-'*(w_xyz+2)}-{'-'*(w_xyz+2)}-{'-'*(w_xyz+2)}+"

print(line)
print(f"| {'Conductivity (unit: S/cm)':^{len(line)-4}} |")
print(line)
print(f"| {'T (K)':^{w_t}} | {'Sigma_all':^{w_all}} | {'Sigma_x':^{w_xyz}} | {'Sigma_y':^{w_xyz}} | {'Sigma_z':^{w_xyz}} |")
print(line)
for temp, s_all, s_x, s_y, s_z in zip(group_Ts, group_sigmas_all, group_sigmas_x, group_sigmas_y, group_sigmas_z):
    print(
        f"| {temp:^{w_t}d} | {s_all:^{w_all}.3e} | {s_x:^{w_xyz}.3e} | {s_y:^{w_xyz}.3e} | {s_z:^{w_xyz}.3e} |"
    )
print(line)

print(line)
print(f"| {'Sigma*T (unit: K S/cm)':^{len(line)-4}} |")
print(line)
print(f"| {'T (K)':^{w_t}} | {'SigmaT_all':^{w_all}} | {'SigmaT_x':^{w_xyz}} | {'SigmaT_y':^{w_xyz}} | {'SigmaT_z':^{w_xyz}} |")
print(line)
for temp, s_all, s_x, s_y, s_z in zip(group_Ts, group_sigmaTs_all, group_sigmaTs_x, group_sigmaTs_y, group_sigmaTs_z):
    print(
        f"| {temp:^{w_t}d} | {s_all:^{w_all}.3e} | {s_x:^{w_xyz}.3e} | {s_y:^{w_xyz}.3e} | {s_z:^{w_xyz}.3e} |"
    )
print(line)

# Print activation energies
print(f"\n{cell}, Ea_total: {Ea_all:.3f} eV")
print(f"{cell}, Ea_x: {Ea_x:.3f} eV")
print(f"{cell}, Ea_y: {Ea_y:.3f} eV")
print(f"{cell}, Ea_z: {Ea_z:.3f} eV")

# Plot each direction
for label, sigmaT_data, intercept, slope, color, marker in labels_with_ea:
    # Fitted line
    y_fit_ln = (slope * x_fit_range + intercept) * np.log(10)
    ax.plot(x_fit_range, y_fit_ln, '--', c=color, linewidth=1.5)
    # Data points
    ax.plot(1000/group_Ts, np.log(sigmaT_data), marker, label=label, 
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

# Calculate conductivity at 300K for each direction
target_T = 300
sig_300_all = (10**(slope_all * (1000/target_T) + intercept_all)) / target_T
sig_300_x = (10**(slope_x * (1000/target_T) + intercept_x)) / target_T
sig_300_y = (10**(slope_y * (1000/target_T) + intercept_y)) / target_T
sig_300_z = (10**(slope_z * (1000/target_T) + intercept_z)) / target_T

print(f"at {target_T}K, {cell}: Sigma_total = {sig_300_all:.3e} S/cm")
print(f"at {target_T}K, {cell}: Sigma_x = {sig_300_x:.3e} S/cm")
print(f"at {target_T}K, {cell}: Sigma_y = {sig_300_y:.3e} S/cm")
print(f"at {target_T}K, {cell}: Sigma_z = {sig_300_z:.3e} S/cm")

# Adjust layout and save
plt.tight_layout()

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('Arrhenius_sigma_xyz.png', dpi=300, bbox_inches='tight')
else:
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'Arrhenius_sigma_xyz.png'.")
        plt.savefig('Arrhenius_sigma_xyz.png', dpi=300, bbox_inches='tight')
    else:
        plt.show()

"""
This script reads MSD data from multiple temperature folders (e.g., 600K, 700K, etc.), 
calculates diffusion coefficients (D) for all directions and individual x, y, z components, 
and plots the Arrhenius relationship (log10(D) vs. 1000/T). 
It also fits the data to extract activation energies (Ea) for each component. 

Author: Modified from original by Zihan YAN (yanzihan@westlake.edu.cn)

Usage:
    python plt_arrhenius_d_xyz.py [save]
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as consts
from scipy.stats import linregress
import matplotlib.ticker as ticker
import warnings
warnings.filterwarnings('ignore')

# PRL 风格基础设置
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

def fit_arrhenius(T, D):
    """Fit Arrhenius equation to calculate activation energy for diffusion"""
    mask = ~np.isnan(D) & ~np.isnan(T) & (D > 0)
    T_valid, D_valid = T[mask], D[mask]
    if len(T_valid) < 2: 
        return None, None, None, None
    
    x = 1000 / T_valid
    y = np.log10(D_valid) 
    res = linregress(x, y)
    Ea = -res.slope * k * 1000 * np.log(10)  
    return Ea, res.intercept, res.slope, res.rvalue

def read_msd_file(msd_file):
    data = np.loadtxt(msd_file)
    dts = data[:, 0]
    msd_x = data[:, 1]
    msd_y = data[:, 2]
    msd_z = data[:, 3]
    msd_all = msd_x + msd_y + msd_z
    return dts, msd_all, msd_x, msd_y, msd_z

def tick_function(X):
    """Convert inverse temperature back to temperature for top axis"""
    V = 1000 / X  
    return ["%.0f" % z for z in V]

# Plot setup with PRL style
fig, ax = plt.subplots(figsize=(4.3, 3.8))

# Automatically find *K folders
raw_temps = []
raw_D_alls = []
raw_D_xs = []
raw_D_ys = []
raw_D_zs = []

current_dir = os.getcwd()
for folder in os.listdir(current_dir):
    if folder.endswith('K') and os.path.isdir(folder):
        try:
            temp = int(folder[:-1])  # Remove 'K' and convert to int
            msd_file = os.path.join(folder, 'msd.out')
            if os.path.exists(msd_file):
                dts, msd_all, msd_x, msd_y, msd_z = read_msd_file(msd_file)

                start = int(len(dts) * 0.4)
                end = int(len(dts) * 0.8)
                k_all, _ = np.polyfit(dts[start:end], msd_all[start:end], 1)
                k_x, _ = np.polyfit(dts[start:end], msd_x[start:end], 1)
                k_y, _ = np.polyfit(dts[start:end], msd_y[start:end], 1)
                k_z, _ = np.polyfit(dts[start:end], msd_z[start:end], 1)

                D_all = k_all * 1e-4 / 6
                D_x = k_x * 1e-4 / 2
                D_y = k_y * 1e-4 / 2
                D_z = k_z * 1e-4 / 2

                raw_temps.append(temp)
                raw_D_alls.append(D_all)
                raw_D_xs.append(D_x)
                raw_D_ys.append(D_y)
                raw_D_zs.append(D_z)
                
        except ValueError:
            continue  # Skip if folder name doesn't convert to int properly

# Sort by temperature
if len(raw_temps) > 0:
    sorted_indices = np.argsort(raw_temps)
    temps = np.array(raw_temps)[sorted_indices]
    D_alls = np.array(raw_D_alls)[sorted_indices]
    D_xs = np.array(raw_D_xs)[sorted_indices]
    D_ys = np.array(raw_D_ys)[sorted_indices]
    D_zs = np.array(raw_D_zs)[sorted_indices]

    # Calculate activation energies for each direction
    Ea_all, intercept_all, slope_all, r_all = fit_arrhenius(temps, D_alls)
    Ea_x, intercept_x, slope_x, r_x = fit_arrhenius(temps, D_xs)
    Ea_y, intercept_y, slope_y, r_y = fit_arrhenius(temps, D_ys)
    Ea_z, intercept_z, slope_z, r_z = fit_arrhenius(temps, D_zs)

    # Prepare labels with activation energies
    labels_with_ea = [
        (f'Total ({Ea_all:.3f} eV)', D_alls, intercept_all, slope_all, '#2A9D8F', 'o'),
        (f'X ({Ea_x:.3f} eV)', D_xs, intercept_x, slope_x, '#D62828', 's'),
        (f'Y ({Ea_y:.3f} eV)', D_ys, intercept_y, slope_y, '#457B9D', '^'),
        (f'Z ({Ea_z:.3f} eV)', D_zs, intercept_z, slope_z, '#E9C46A', 'd')
    ]

    # Print activation energies
    print(f"Ea_all: {Ea_all:.3f} eV")
    print(f"Ea_x: {Ea_x:.3f} eV")
    print(f"Ea_y: {Ea_y:.3f} eV")
    print(f"Ea_z: {Ea_z:.3f} eV")

    # Prepare data for plotting
    x_min, x_max = 1000/np.max(temps), 1000/np.min(temps)
    x_fit_range = np.linspace(x_min, x_max, 100)

    # Plot each direction
    for label, D_data, intercept, slope, color, marker in labels_with_ea:
        # Fitted line
        y_fit_log10 = slope * x_fit_range + intercept
        ax.plot(x_fit_range, y_fit_log10, '--', c=color, linewidth=1.5)
        # Data points
        ax.plot(1000/temps, np.log10(D_data), marker, label=label, 
                c=color, markersize=8, markerfacecolor='none', markeredgewidth=1.5)

    # Generate temperature and diffusion coefficient arrays for export
    # Using total diffusion coefficient (D_all) for the export
    temp_array = temps.tolist()
    diff_coeff_array = D_alls.tolist()
    
    print("\nExportable arrays for Python:")
    print(f"temperatures = {temp_array}")
    print(f"diffusion_coeffs = [{', '.join([f'{val:.3e}' for val in diff_coeff_array])}]\n")

# Print diffusion coefficients table
w_t, w_all, w_xyz = 10, 15, 12
line = f"+{'-'*(w_t+2)}-{'-'*(w_all+2)}-{'-'*(w_xyz+2)}-{'-'*(w_xyz+2)}-{'-'*(w_xyz+2)}+"

print(line)
print(f"| {'Diffusion coefficients (unit: cm2/s)':^{len(line)-4}} |")
print(line)
print(f"| {'T (K)':^{w_t}} | {'D_all':^{w_all}} | {'D_x':^{w_xyz}} | {'D_y':^{w_xyz}} | {'D_z':^{w_xyz}} |")
print(line)

for temp, D_all, D_x, D_y, D_z in zip(temps, D_alls, D_xs, D_ys, D_zs):
    print(
        f"| {temp:^{w_t}d} | {D_all:^{w_all}.3e} | {D_x:^{w_xyz}.3e} | {D_y:^{w_xyz}.3e} | {D_z:^{w_xyz}.3e} |"
    )
print(line)

# Axes setup with PRL styling
ax.set_xlabel('1000/T (1/K)', labelpad=7)
ax.set_ylabel(r'log10(D) (cm$^2$/s)')

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

plt.tight_layout()

# Check if last argument is 'save' and handle accordingly
should_save = len(sys.argv) > 1 and sys.argv[-1] == 'save'
if should_save:
    plt.savefig('Arrhenius_D_xyz.png', dpi=300, bbox_inches='tight')
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'Arrhenius_D_xyz.png'.")
        plt.savefig('Arrhenius_D_xyz.png', dpi=300, bbox_inches='tight')
    else:
        plt.show()

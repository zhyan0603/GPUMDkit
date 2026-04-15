"""
This script reads msd.out files from temperature-specific folders (e.g., 600K, 700K, etc.), 
calculates the diffusion coefficients, and plots them in an Arrhenius plot. 
It also fits the data to extract the activation energy. 

Author: Modified from original by Zihan YAN (yanzihan@westlake.edu.cn)

Usage:
    python plt_arrhenius_d.py [save]
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as consts
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
    msds_all = np.sum(data[:, 1:4], axis=1)
    return dts, msds_all

def tick_function(X):
    """Convert inverse temperature back to temperature for top axis"""
    V = 1000 / X  
    return ["%.0f" % z for z in V]

# Plot setup with PRL style
fig, ax = plt.subplots(figsize=(4.3, 3.8))

# Automatically find *K folders
raw_temps = []
raw_D_alls = []

current_dir = os.getcwd()
for folder in os.listdir(current_dir):
    if folder.endswith('K') and os.path.isdir(folder):
        try:
            temp = int(folder[:-1])  # Remove 'K' and convert to int
            msd_file = os.path.join(folder, 'msd.out')
            if os.path.exists(msd_file):
                dts, msds_all = read_msd_file(msd_file)
                
                start = int(len(dts) * 0.4)
                end = int(len(dts) * 0.8)
                k_all, _ = np.polyfit(dts[start:end], msds_all[start:end], 1)
                D_all = k_all * 1e-4 / 6
                raw_temps.append(temp)
                raw_D_alls.append(D_all)
        except ValueError:
            continue  # Skip if folder name doesn't convert to int properly

# Sort by temperature
if len(raw_temps) > 0:
    sorted_indices = np.argsort(raw_temps)
    temps = np.array(raw_temps)[sorted_indices]
    D_alls = np.array(raw_D_alls)[sorted_indices]

    # Fit and calculate activation energy
    Ea_all, intercept_all, slope_all, r_all = fit_arrhenius(temps, D_alls)
    
    if Ea_all is not None:
        print(f"\nActivation Energy: {Ea_all:.3f} eV")
        
        # Prepare data for plotting
        x_min, x_max = 1000/np.max(temps), 1000/np.min(temps)
        x_fit_range = np.linspace(x_min, x_max, 100)
        y_fit_log10 = slope_all * x_fit_range + intercept_all
        
        # Plot fitted line
        ax.plot(x_fit_range, y_fit_log10, '--', c='#D62828', linewidth=1.5)
        
        # Plot data points
        ax.plot(1000/temps, np.log10(D_alls), 'o', label=f'D_total ({Ea_all:.3f} eV)', 
                c='#D62828', markersize=8, markerfacecolor='none', markeredgewidth=1.5)
    else:
        # If fitting failed, just plot the data points
        ax.plot(1000/temps, np.log10(D_alls), 'o', label='D_total', 
                c='#D62828', markersize=8, markerfacecolor='none', markeredgewidth=1.5)

    # Generate temperature and diffusion coefficient arrays for export
    temp_array = temps.tolist()
    diff_coeff_array = [float(f'{val:.3e}') for val in D_alls.tolist()]
    
    print("\nExportable arrays for Python:")
    print(f"temperatures = {temp_array}")
    print(f"diffusion_coeffs = {diff_coeff_array}\n")

# Print diffusion coefficients table
w_t, w_all = 10, 15
line = f"+{'-'*(w_t+2)}-{'-'*(w_all+2)}+"

print(line)
print(f"| {'  Diffusivity (unit: cm2/s)':^{len(line)-4}} |")
print(line)
print(f"| {'T (K)':^{w_t}} | {'D_total':^{w_all}} |")
print(line)

for temp, D_all in zip(temps, D_alls):
    print(f"| {temp:^{w_t}d} | {D_all:^{w_all}.3e} |")
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
    plt.savefig('Arrhenius_D.png', dpi=300, bbox_inches='tight')
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'Arrhenius_D.png'.")
        plt.savefig('Arrhenius_D.png', dpi=300, bbox_inches='tight')
    else:
        plt.show()

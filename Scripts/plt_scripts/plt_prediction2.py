"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_prediction2.py
Category:   Plot Scripts
Purpose:    Alternative NEP prediction parity plot (energy, forces, stresses)
            with different layout from plt_prediction.py.
Usage:      gpumdkit.sh -plt prediction2
            python plt_prediction2.py
Output:
    prediction.png  (if save is used, or if backend is non-interactive)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-21
=============================================================================
"""

import sys
import matplotlib.pyplot as plt
import numpy as np

# Plot style (self-contained; keep in sync across plt_*.py scripts)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "mathtext.fontset": "dejavusans",
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.linewidth": 0.9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.minor.size": 2.0,
    "ytick.minor.size": 2.0,
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "legend.frameon": False,
    "legend.fontsize": 9,
    "figure.dpi": 100,
})
COLORS = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
REF_STYLE = dict(color='0.55', linestyle='--', linewidth=1.2, zorder=0)
PT_STYLE = dict(marker='o', linestyle='none', markersize=4.5, markeredgewidth=0, rasterized=True)

# Load data
energy_data = np.loadtxt('energy_train.out')
force_data = np.loadtxt('force_train.out')
#virial_data = np.loadtxt('virial_train.out')
stress_data = np.loadtxt('stress_train.out')

# Filter out rows with invalid stress data
valid_rows = ~np.any(np.abs(stress_data[:, :12]) >= 1e6, axis=1)
stress_data = stress_data[valid_rows]

# Function to calculate RMSE
def calculate_rmse(pred, actual):
    return np.sqrt(np.mean((pred - actual) ** 2))

# Function to calculate MAE
def calculate_mae(pred, actual):
    return np.mean(np.abs(pred - actual))

# Function to calculate R²
def calculate_r2(pred, actual):
    ss_tot = np.sum((actual - np.mean(actual)) ** 2)
    ss_res = np.sum((pred - actual) ** 2)
    return 1 - ss_res / ss_tot if ss_tot != 0 else 1.0

# Function to calculate dynamic axis limits
def calculate_limits(train_data, padding=0.08):
    data_min = np.min(train_data)
    data_max = np.max(train_data)
    data_range = data_max - data_min
    return data_min - padding * data_range, data_max + padding * data_range

def parity_panel(ax, actual, pred, labels, colors, lim):
    """Draw a square parity panel with a reference diagonal and legend."""
    ax.plot(lim, lim, **REF_STYLE)
    for i, (label, color) in enumerate(zip(labels, colors)):
        ax.plot(actual[:, i], pred[:, i], color=color, label=label, **PT_STYLE)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.legend(loc='upper left', ncol=2 if len(labels) > 3 else 1, markerscale=1.0, handletextpad=0.2, columnspacing=1.0)

def metrics_text(ax, pred, actual, scale, unit, ndigits):
    """Write R2/MAE/RMSE in the lower-right corner of a parity panel."""
    pred, actual = pred.flatten(), actual.flatten()
    r2 = calculate_r2(pred, actual)
    mae = calculate_mae(pred, actual) * scale
    rmse = calculate_rmse(pred, actual) * scale
    ax.text(0.96, 0.05,
            r'$R^2$ = ' + f'{r2:.4f}\nMAE = {mae:.{ndigits}f} {unit}\nRMSE = {rmse:.{ndigits}f} {unit}',
            transform=ax.transAxes, fontsize=9, va='bottom', ha='right', linespacing=1.5)

# Create a subplot with 1 row and 3 columns
fig, axs = plt.subplots(1, 3, figsize=(10.8, 3.2))

# Energy panel
energy_lim = calculate_limits(energy_data[:, 1])
parity_panel(axs[0], energy_data[:, [1]], energy_data[:, [0]], ['Energy'], COLORS[:1], energy_lim)
axs[0].set_xlabel('DFT Energy (eV/atom)')
axs[0].set_ylabel('NEP Energy (eV/atom)')
metrics_text(axs[0], energy_data[:, 0], energy_data[:, 1], 1000.0, 'meV/atom', 2)

# Force panel
force_lim = calculate_limits(force_data[:, 3:6].reshape(-1))
parity_panel(axs[1], force_data[:, 3:6], force_data[:, 0:3], ['fx', 'fy', 'fz'], COLORS[:3], force_lim)
axs[1].set_xlabel(r'DFT Force (eV/$\mathrm{\AA}$)')
axs[1].set_ylabel(r'NEP Force (eV/$\mathrm{\AA}$)')
metrics_text(axs[1], force_data[:, 0:3], force_data[:, 3:6], 1000.0, r'meV/$\mathrm{\AA}$', 2)

# Stress panel
stress_lim = calculate_limits(stress_data[:, 6:12].reshape(-1))
parity_panel(axs[2], stress_data[:, 6:12], stress_data[:, 0:6], ['xx', 'yy', 'zz', 'xy', 'yz', 'zx'], COLORS, stress_lim)
axs[2].set_xlabel('DFT Stress (GPa)')
axs[2].set_ylabel('NEP Stress (GPa)')
metrics_text(axs[2], stress_data[:, 0:6], stress_data[:, 6:12], 1.0, 'GPa', 4)

# Subplot geometry (figure fractions); tune these to adjust margins and spacing
fig.subplots_adjust(left=0.065, right=0.985, top=0.96, bottom=0.18, wspace=0.34)

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('prediction.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'prediction.png'.")
        plt.savefig('prediction.png', dpi=300)
    else:
        plt.show()

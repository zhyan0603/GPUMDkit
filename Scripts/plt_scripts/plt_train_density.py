"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_train_density.py
Category:   Plot Scripts
Purpose:    Generate density-based parity plots for NEP training results
            (energy, forces, stress or virial), useful for large datasets.
            Stress is preferred when both tensor files exist.
Usage:      gpumdkit.sh -plt train_density [save]
            python plt_train_density.py [save]
Arguments:
  save      Save the plot as 'train_density.png' instead of displaying it
Output:
  train_density.png  (if save is used, or if backend is non-interactive)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-21
=============================================================================
"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, ListedColormap

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
def density_cmap(name, start=0.35):
    """Colormap without its near-white low end, so sparse bins stay visible."""
    return ListedColormap(plt.get_cmap(name)(np.linspace(start, 1.0, 256)))

def density_bins(n_points, low=50, high=90):
    """Fewer bins for small data sets so that each bin holds enough points."""
    return int(np.clip(np.sqrt(n_points) / 2, low, high))

REF_STYLE = dict(color='0.55', linestyle='--', linewidth=1.2, zorder=2)

def panel_label(ax, text, x=-0.14, y=1.0):
    """Place a panel label outside the upper-left corner of an axes."""
    ax.text(x, y, text, transform=ax.transAxes, fontsize=12, va='bottom', ha='right')

def metrics_text(ax, r2, mae, rmse, unit, ndigits):
    """Write R2/MAE/RMSE in the lower-right corner of a parity panel."""
    ax.text(0.96, 0.05,
            r'$R^2$ = ' + f'{r2:.4f}\nMAE = {mae:.{ndigits}f} {unit}\nRMSE = {rmse:.{ndigits}f} {unit}',
            transform=ax.transAxes, fontsize=9, va='bottom', ha='right', linespacing=1.5)

# Load data
loss = np.loadtxt('loss.out')
energy_data = np.loadtxt('energy_train.out')
force_data = np.loadtxt('force_train.out')
tensor_name = "Stress"
tensor_axis_unit = "GPa"
tensor_metric_unit = "GPa"
tensor_metric_scale = 1.0
tensor_data = None
if os.path.isfile('stress_train.out'):
    tensor_data = np.atleast_2d(np.loadtxt('stress_train.out'))
elif os.path.isfile('virial_train.out'):
    tensor_name = "Virial"
    tensor_axis_unit = "eV/atom"
    tensor_metric_unit = "meV/atom"
    tensor_metric_scale = 1000.0
    tensor_data = np.atleast_2d(np.loadtxt('virial_train.out'))

# Filter out rows with invalid stress/virial target data
if tensor_data is not None:
    valid_rows = ~np.any(np.abs(tensor_data[:, :12]) >= 1e6, axis=1)
    tensor_data = tensor_data[valid_rows]

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

def print_parity_metrics(energy, force, tensor, tensor_label, tensor_unit_label):
    """Print energy, force, and tensor parity metrics as a three-line table."""
    tensor_r2, tensor_mae, tensor_rmse = tensor if tensor is not None else (None, None, None)
    tensor_r2_text = f"{tensor_r2:.4f}" if tensor_r2 is not None else "N/A"
    tensor_mae_text = f"{tensor_mae:.4f}" if tensor_mae is not None else "N/A"
    tensor_rmse_text = f"{tensor_rmse:.4f}" if tensor_rmse is not None else "N/A"
    line = " " + "-" * 48
    print(" Parity metrics (training set)")
    print(f" Energy: meV/atom, Force: meV/Ang, {tensor_label}: {tensor_unit_label}")
    print(line)
    print(f" {'Metric':<8}{'Energy':>12}{'Force':>12}{tensor_label:>12}")
    print(line)
    print(f" {'R^2':<8}{energy[0]:>12.4f}{force[0]:>12.4f}{tensor_r2_text:>12}")
    print(f" {'MAE':<8}{energy[1]:>12.2f}{force[1]:>12.2f}{tensor_mae_text:>12}")
    print(f" {'RMSE':<8}{energy[2]:>12.2f}{force[2]:>12.2f}{tensor_rmse_text:>12}")
    print(line)

# Function to calculate dynamic axis limits
def calculate_limits(train_data, padding=0.08):
    data_min = np.min(train_data)
    data_max = np.max(train_data)
    data_range = data_max - data_min
    return data_min - padding * data_range, data_max + padding * data_range

# Create a subplot with 2 rows and 2 columns 
fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.1))

if loss[0, 0] == 100:
    xlabel = 'Generation/100'
    plot_cols = slice(1, 7)
    legend_labels = ['Total', 'L1', 'L2', 'Energy', 'Force', 'Virial']
elif loss[0, 0] == 1:
    xlabel = 'Epoch'
    plot_cols = slice(1, 5)
    legend_labels = ['Total', 'Energy', 'Force', 'Virial']
else:
    raise ValueError("Unexpected loss data format.")


for i, label in enumerate(legend_labels):
    axs[0, 0].loglog(loss[:, plot_cols][:, i], '-', linewidth=1.4, color=COLORS[i % len(COLORS)], label=label)
axs[0, 0].set_xlabel(xlabel)
axs[0, 0].set_ylabel('Loss functions')
axs[0, 0].margins(x=0.02)
axs[0, 0].legend(loc='lower left', ncol=2, fontsize=8.5, handlelength=1.4, handletextpad=0.5, columnspacing=1.0, labelspacing=0.25)
panel_label(axs[0, 0], '(a)', x=-0.11)

xmin_energy, xmax_energy = calculate_limits(energy_data[:, 0])
axs[0, 1].set_xlim(xmin_energy, xmax_energy)
axs[0, 1].set_ylim(xmin_energy, xmax_energy)

axs[0, 1].hist2d(energy_data[:, 1], energy_data[:, 0], 
                 bins=density_bins(len(energy_data)), cmap=density_cmap('Blues'), cmin=1, norm=LogNorm(),
                 range=[[xmin_energy, xmax_energy], [xmin_energy, xmax_energy]])

axs[0, 1].plot([xmin_energy, xmax_energy], [xmin_energy, xmax_energy], **REF_STYLE)
axs[0, 1].set_xlabel('DFT Energy (eV/atom)')
axs[0, 1].set_ylabel('NEP Energy (eV/atom)')

energy_rmse = calculate_rmse(energy_data[:, 1], energy_data[:, 0]) * 1000
energy_mae = calculate_mae(energy_data[:, 1], energy_data[:, 0]) * 1000
energy_r2 = calculate_r2(energy_data[:, 1], energy_data[:, 0])
metrics_text(axs[0, 1], energy_r2, energy_mae, energy_rmse, 'meV/atom', 2)
panel_label(axs[0, 1], '(b)')


f_pred = force_data[:, 0:3].flatten()
f_target = force_data[:, 3:6].flatten()

xmin_force, xmax_force = calculate_limits(f_target)
axs[1, 0].set_xlim(xmin_force, xmax_force)
axs[1, 0].set_ylim(xmin_force, xmax_force)

axs[1, 0].hist2d(f_target, f_pred, 
                 bins=density_bins(len(f_target)), cmap=density_cmap('Oranges'), cmin=1, norm=LogNorm(),
                 range=[[xmin_force, xmax_force], [xmin_force, xmax_force]])

axs[1, 0].plot([xmin_force, xmax_force], [xmin_force, xmax_force], **REF_STYLE)
axs[1, 0].set_xlabel(r'DFT Force (eV/$\mathrm{\AA}$)')
axs[1, 0].set_ylabel(r'NEP Force (eV/$\mathrm{\AA}$)')

force_rmse = calculate_rmse(f_pred, f_target) * 1000
force_mae = calculate_mae(f_pred, f_target) * 1000
force_r2 = calculate_r2(f_pred, f_target)
metrics_text(axs[1, 0], force_r2, force_mae, force_rmse, r'meV/$\mathrm{\AA}$', 2)
panel_label(axs[1, 0], '(c)', x=-0.11)


if tensor_data is None or tensor_data.shape[0] == 0:
    axs[1, 1].axis('off')
    tensor_metrics = None
else:
    tensor_pred = tensor_data[:, 0:6].flatten()
    tensor_target = tensor_data[:, 6:12].flatten()

    xmin_tensor, xmax_tensor = calculate_limits(tensor_target)
    axs[1, 1].set_xlim(xmin_tensor, xmax_tensor)
    axs[1, 1].set_ylim(xmin_tensor, xmax_tensor)
    
    axs[1, 1].hist2d(tensor_target, tensor_pred,
                     bins=density_bins(len(tensor_target)), cmap=density_cmap('Greens'), cmin=1, norm=LogNorm(),
                     range=[[xmin_tensor, xmax_tensor], [xmin_tensor, xmax_tensor]])
        
    axs[1, 1].plot([xmin_tensor, xmax_tensor], [xmin_tensor, xmax_tensor], **REF_STYLE)
    axs[1, 1].set_xlabel(f'DFT {tensor_name} ({tensor_axis_unit})')
    axs[1, 1].set_ylabel(f'NEP {tensor_name} ({tensor_axis_unit})')

    tensor_rmse = calculate_rmse(tensor_pred, tensor_target) * tensor_metric_scale
    tensor_mae = calculate_mae(tensor_pred, tensor_target) * tensor_metric_scale
    tensor_r2 = calculate_r2(tensor_pred, tensor_target)
    tensor_metrics = (tensor_r2, tensor_mae, tensor_rmse)
    
    metrics_text(axs[1, 1], tensor_r2, tensor_mae, tensor_rmse, tensor_metric_unit, 4)
    panel_label(axs[1, 1], '(d)')

print_parity_metrics(
    (energy_r2, energy_mae, energy_rmse),
    (force_r2, force_mae, force_rmse),
    tensor_metrics,
    tensor_name,
    tensor_metric_unit,
)

# Subplot geometry (figure fractions); tune these to adjust margins and spacing
fig.subplots_adjust(left=0.11, right=0.98, top=0.94, bottom=0.09, hspace=0.30, wspace=0.34)

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('train_density.png', dpi=300)
else:
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'train_density.png'.")
        plt.savefig('train_density.png', dpi=300)
    else:
        plt.show()

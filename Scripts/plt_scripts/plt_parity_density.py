"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_parity_density.py
Category:   Plot Scripts
Purpose:    Generate density-based parity plots for energies, forces, and
            stresses, especially useful for large NEP training datasets.
Usage:      gpumdkit.sh -plt parity_density
            python plt_parity_density.py
Output:
  parity_density_plot.png  (if save is used, or if backend is non-interactive)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-07-11
=============================================================================
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext, LogLocator
from matplotlib.gridspec import GridSpec
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
})

# Configuration
FONT_SIZE = 10
BINS = 200
COLORMAP = 'viridis' 
# preferred colormaps: 'rainbow', 'plasma', 'inferno', 'magma', 'cividis', 'viridis'

def print_parity_metrics(energy, force, stress):
    """Print energy, force, and stress parity metrics as a three-line table."""
    stress_r2, stress_mae, stress_rmse = stress if stress is not None else (None, None, None)
    stress_r2_text = f"{stress_r2:.4f}" if stress_r2 is not None else "N/A"
    stress_mae_text = f"{stress_mae:.4f}" if stress_mae is not None else "N/A"
    stress_rmse_text = f"{stress_rmse:.4f}" if stress_rmse is not None else "N/A"
    line = " " + "-" * 48
    print(" Parity metrics (training set)")
    print(" Energy: meV/atom, Force: meV/Ang, Stress: GPa")
    print(line)
    print(f" {'Metric':<8}{'Energy':>12}{'Force':>12}{'Stress':>12}")
    print(line)
    print(f" {'R^2':<8}{energy[0]:>12.4f}{force[0]:>12.4f}{stress_r2_text:>12}")
    print(f" {'MAE':<8}{energy[1]:>12.2f}{force[1]:>12.2f}{stress_mae_text:>12}")
    print(f" {'RMSE':<8}{energy[2]:>12.2f}{force[2]:>12.2f}{stress_rmse_text:>12}")
    print(line)

# Load data
energy_data = np.loadtxt('energy_train.out')
force_data = np.loadtxt('force_train.out')
stress_data = np.loadtxt('stress_train.out')

# Filter out rows with invalid stress data
valid_rows = ~np.any(np.abs(stress_data[:, :12]) > 1e6, axis=1)
stress_data = stress_data[valid_rows]

# === Layout Setup ===
fig = plt.figure(figsize=(12, 4.2), dpi=100)
gs = GridSpec(2, 3, height_ratios=[4, 0.2], hspace=0.35)
cmap = plt.get_cmap(COLORMAP)

# === Energy ===
ax1 = fig.add_subplot(gs[0, 0])
cb1_ax = fig.add_subplot(gs[1, 0])
xmin, xmax = np.min([energy_data[:, 1], energy_data[:, 0]]), np.max([energy_data[:, 1], energy_data[:, 0]])
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(xmin, xmax)
hb1 = ax1.hexbin(energy_data[:, 1], energy_data[:, 0], gridsize=BINS, cmap=cmap,
                 norm=LogNorm(), mincnt=1, linewidths=0)
ax1.plot([xmin, xmax], [xmin, xmax], color='gray', linestyle='-', linewidth=1)
ax1.set_xlabel("DFT energy (eV/atom)", fontsize=FONT_SIZE)
ax1.set_ylabel("NEP energy (eV/atom)", fontsize=FONT_SIZE)
ax1.tick_params(labelsize=FONT_SIZE)
energy_rmse = np.sqrt(mean_squared_error(energy_data[:, 0], energy_data[:, 1])) * 1000
energy_mae = mean_absolute_error(energy_data[:, 0], energy_data[:, 1]) * 1000
energy_r2 = r2_score(energy_data[:, 0], energy_data[:, 1])
ax1.text(0.05, 0.95,
         f"RMSE = {energy_rmse:.2f} meV/atom\n"
         f"MAE = {energy_mae:.2f} meV/atom\n"
         r"$\mathrm{{R^2}}$" + f" = {energy_r2:.5f}",
         transform=ax1.transAxes,
         fontsize=FONT_SIZE,
         va='top', ha='left')
cb1 = fig.colorbar(hb1, cax=cb1_ax, orientation='horizontal')
cb1.set_label("Data density", fontsize=FONT_SIZE)
cb1.ax.tick_params(labelsize=FONT_SIZE)
cb1.locator = LogLocator(base=10.0)
cb1.formatter = LogFormatterMathtext(base=10, labelOnlyBase=True)
cb1.update_ticks()

# === Force ===
ax2 = fig.add_subplot(gs[0, 1])
cb2_ax = fig.add_subplot(gs[1, 1])
x_force = force_data[:, 3:6].reshape(-1)
y_force = force_data[:, 0:3].reshape(-1)
xmin, xmax = np.min([x_force, y_force]), np.max([x_force, y_force])
ax2.set_xlim(xmin, xmax)
ax2.set_ylim(xmin, xmax)
hb2 = ax2.hexbin(x_force, y_force, gridsize=BINS, cmap=cmap,
                 norm=LogNorm(), mincnt=1, linewidths=0)
ax2.plot([xmin, xmax], [xmin, xmax], color='gray', linestyle='-', linewidth=1)
ax2.set_xlabel(r"DFT force (eV/$\mathrm{{\AA}}$)", fontsize=FONT_SIZE)
ax2.set_ylabel(r"NEP force (eV/$\mathrm{{\AA}}$)", fontsize=FONT_SIZE)
ax2.tick_params(labelsize=FONT_SIZE)
force_rmse = np.sqrt(mean_squared_error(y_force, x_force)) * 1000
force_mae = mean_absolute_error(y_force, x_force) * 1000
force_r2 = r2_score(y_force, x_force)
ax2.text(0.05, 0.95,
         f'RMSE = {force_rmse:.2f} meV/'+r'$\mathrm{{\AA}}$'+f'\n'
         f'MAE = {force_mae:.2f} meV/'+r'$\mathrm{{\AA}}$'+f'\n'
         r"$\mathrm{{R^2}}$" + f' = {force_r2:.5f}',
         transform=ax2.transAxes,
         fontsize=FONT_SIZE,
         va='top', ha='left')
cb2 = fig.colorbar(hb2, cax=cb2_ax, orientation='horizontal')
cb2.set_label("Data density", fontsize=FONT_SIZE)
cb2.ax.tick_params(labelsize=FONT_SIZE)
cb2.locator = LogLocator(base=10.0)
cb2.formatter = LogFormatterMathtext(base=10, labelOnlyBase=True)
cb2.update_ticks()

# === Stress ===
ax3 = fig.add_subplot(gs[0, 2])
cb3_ax = fig.add_subplot(gs[1, 2])
if stress_data.shape[0] == 0:
    ax3.axis('off')
    cb3_ax.axis('off')
    stress_metrics = None
else:
    x_stress = stress_data[:, 6:12].reshape(-1)
    y_stress = stress_data[:, 0:6].reshape(-1)
    xmin, xmax = np.min([x_stress, y_stress]), np.max([x_stress, y_stress])
    ax3.set_xlim(xmin, xmax)
    ax3.set_ylim(xmin, xmax)
    hb3 = ax3.hexbin(x_stress, y_stress, gridsize=BINS, cmap=cmap,
                     norm=LogNorm(), mincnt=1, linewidths=0)
    ax3.plot([xmin, xmax], [xmin, xmax], color='gray', linestyle='-', linewidth=1)
    ax3.set_xlabel("DFT stress (GPa)", fontsize=FONT_SIZE)
    ax3.set_ylabel("NEP stress (GPa)", fontsize=FONT_SIZE)
    ax3.tick_params(labelsize=FONT_SIZE)
    stress_rmse = np.sqrt(mean_squared_error(y_stress, x_stress))
    stress_mae = mean_absolute_error(y_stress, x_stress)
    stress_r2 = r2_score(y_stress, x_stress)
    stress_metrics = (stress_r2, stress_mae, stress_rmse)
    ax3.text(0.05, 0.95,
             f"RMSE = {stress_rmse:.4f} GPa\n"
             f"MAE = {stress_mae:.4f} GPa\n"
             r"$\mathrm{{R^2}}$" + f" = {stress_r2:.5f}",
             transform=ax3.transAxes,
             fontsize=FONT_SIZE,
             va='top', ha='left')
    cb3 = fig.colorbar(hb3, cax=cb3_ax, orientation='horizontal')
    cb3.set_label("Data density", fontsize=FONT_SIZE)
    cb3.ax.tick_params(labelsize=FONT_SIZE)
    cb3.locator = LogLocator(base=10.0)
    cb3.formatter = LogFormatterMathtext(base=10, labelOnlyBase=True)
    cb3.update_ticks()

print_parity_metrics(
    (energy_r2, energy_mae, energy_rmse),
    (force_r2, force_mae, force_rmse),
    stress_metrics,
)

plt.subplots_adjust(top=0.968, bottom=0.122, left=0.073, right=0.983, hspace=0.2, wspace=0.286)

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('parity_density_plot.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'parity_density_plot.png'.")
        plt.savefig('parity_density_plot.png', dpi=300)
    else:
        plt.show()

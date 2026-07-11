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
Last-modified: 2026-07-11
=============================================================================
"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm  

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
})  

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
fig, axs = plt.subplots(2, 2, figsize=(9, 7), dpi=100)

if loss[0, 0] == 100:
    xlabel = 'Generation/100'
    plot_cols = slice(1, 7)
    legend_labels = ['Total', 'L1-Reg', 'L2-Reg', 'Energy-train', 'Force-train', 'Virial-train']
elif loss[0, 0] == 1:
    xlabel = 'Epoch'
    plot_cols = slice(1, 5)
    legend_labels = ['Total', 'Energy-train', 'Force-train', 'Virial-train']
else:
    raise ValueError("Unexpected loss data format.")


axs[0, 0].loglog(loss[:, plot_cols], '-', linewidth=2)
axs[0, 0].set_xlabel(xlabel, fontsize=10)
axs[0, 0].set_ylabel('Loss functions', fontsize=10)
axs[0, 0].tick_params(axis='both', labelsize=10)
axs[0, 0].legend(legend_labels, prop={'size': 8}, loc='lower left', frameon=False)
axs[0, 0].axis('tight')


xmin_energy, xmax_energy = calculate_limits(energy_data[:, 0])
axs[0, 1].set_xlim(xmin_energy, xmax_energy)
axs[0, 1].set_ylim(xmin_energy, xmax_energy)

axs[0, 1].hist2d(energy_data[:, 1], energy_data[:, 0], 
                 bins=100, cmap='Blues', cmin=1, norm=LogNorm(), 
                 range=[[xmin_energy, xmax_energy], [xmin_energy, xmax_energy]])

axs[0, 1].plot([xmin_energy, xmax_energy], [xmin_energy, xmax_energy], linewidth=1.5, color='grey', linestyle='--')
axs[0, 1].set_xlabel('DFT energy (eV/atom)', fontsize=10)
axs[0, 1].set_ylabel('NEP energy (eV/atom)', fontsize=10)
axs[0, 1].tick_params(axis='both', labelsize=10)

energy_rmse = calculate_rmse(energy_data[:, 1], energy_data[:, 0]) * 1000
energy_mae = calculate_mae(energy_data[:, 1], energy_data[:, 0]) * 1000
energy_r2 = calculate_r2(energy_data[:, 1], energy_data[:, 0])
axs[0, 1].text(0.7, 0.12, r'R$^2$'+f': {energy_r2:.4f}\nMAE: {energy_mae:.2f} meV/atom\nRMSE: {energy_rmse:.2f} meV/atom', 
               transform=axs[0, 1].transAxes, fontsize=10, verticalalignment='center', horizontalalignment='center')


f_pred = force_data[:, 0:3].flatten()
f_target = force_data[:, 3:6].flatten()

xmin_force, xmax_force = calculate_limits(f_target)
axs[1, 0].set_xlim(xmin_force, xmax_force)
axs[1, 0].set_ylim(xmin_force, xmax_force)

axs[1, 0].hist2d(f_target, f_pred, 
                 bins=150, cmap='Oranges', cmin=1, norm=LogNorm(),
                 range=[[xmin_force, xmax_force], [xmin_force, xmax_force]])

axs[1, 0].plot([xmin_force, xmax_force], [xmin_force, xmax_force], linewidth=1.5, color='grey', linestyle='--')
axs[1, 0].set_xlabel(r'DFT force (eV/$\mathrm{\AA}$)', fontsize=10)
axs[1, 0].set_ylabel(r'NEP force (eV/$\mathrm{\AA}$)', fontsize=10)
axs[1, 0].tick_params(axis='both', labelsize=10)

force_rmse = calculate_rmse(f_pred, f_target) * 1000
force_mae = calculate_mae(f_pred, f_target) * 1000
force_r2 = calculate_r2(f_pred, f_target)
axs[1, 0].text(0.7, 0.12, 
               r'R$^2$'+f': {force_r2:.4f}\nMAE: {force_mae:.2f} meV/'+r'$\mathrm{{\AA}}$'
               +f'\nRMSE: {force_rmse:.2f} meV/'+r'$\mathrm{{\AA}}$', 
               transform=axs[1, 0].transAxes, fontsize=10, verticalalignment='center', horizontalalignment='center')


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
                     bins=150, cmap='Greens', cmin=1, norm=LogNorm(),
                     range=[[xmin_tensor, xmax_tensor], [xmin_tensor, xmax_tensor]])
        
    axs[1, 1].plot([xmin_tensor, xmax_tensor], [xmin_tensor, xmax_tensor], linewidth=1.5, color='grey', linestyle='--')
    axs[1, 1].set_xlabel(f'DFT {tensor_name.lower()} ({tensor_axis_unit})', fontsize=10)
    axs[1, 1].set_ylabel(f'NEP {tensor_name.lower()} ({tensor_axis_unit})', fontsize=10)
    axs[1, 1].tick_params(axis='both', labelsize=10)

    tensor_rmse = calculate_rmse(tensor_pred, tensor_target) * tensor_metric_scale
    tensor_mae = calculate_mae(tensor_pred, tensor_target) * tensor_metric_scale
    tensor_r2 = calculate_r2(tensor_pred, tensor_target)
    tensor_metrics = (tensor_r2, tensor_mae, tensor_rmse)
    
    axs[1, 1].text(0.7, 0.12, r'R$^2$'+f': {tensor_r2:.4f}\nMAE: {tensor_mae:.4f} {tensor_metric_unit}\nRMSE: {tensor_rmse:.4f} {tensor_metric_unit}',
                transform=axs[1, 1].transAxes, fontsize=10, verticalalignment='center', horizontalalignment='center')

print_parity_metrics(
    (energy_r2, energy_mae, energy_rmse),
    (force_r2, force_mae, force_rmse),
    tensor_metrics,
    tensor_name,
    tensor_metric_unit,
)

# plt.tight_layout()
fig.subplots_adjust(top=0.968,bottom=0.088,left=0.086,right=0.983,hspace=0.22,wspace=0.24)

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

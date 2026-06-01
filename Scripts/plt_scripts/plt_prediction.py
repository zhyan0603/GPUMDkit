"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_prediction.py
Category:   Plot Scripts
Purpose:    Visualize NEP prediction results with parity plots for energy,
            forces, and stresses, including marginal distributions and
            residual histograms.
Usage:      gpumdkit.sh -plt prediction [save]
            python plt_prediction.py [save]
Arguments:
  save      Save the plot as 'prediction.png' instead of displaying it
Output:
  prediction.png  (if save is used, or if backend is non-interactive)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# If running on a server without GUI, keep this commented unless needed
# matplotlib.use('Agg')

# =========================
# Global style
# =========================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False

plt.rcParams['axes.labelsize'] = 13
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['text.color'] = 'black'

# =========================
# Load data - No longer loading loss.out
# =========================
energy_data = np.loadtxt('energy_train.out')
force_data = np.loadtxt('force_train.out')
stress_data = np.loadtxt('stress_train.out')

# Filter invalid stress rows
valid_rows = ~np.any(np.abs(stress_data[:, :12]) >= 1e6, axis=1)
stress_data = stress_data[valid_rows]

# =========================
# Utility functions
# =========================
def rmse(pred, true):
    return np.sqrt(np.mean((pred - true) ** 2))

def mae(pred, true):
    return np.mean(np.abs(pred - true))

def r2_score_np(true, pred):
    true = np.asarray(true).reshape(-1)
    pred = np.asarray(pred).reshape(-1)
    ss_res = np.sum((true - pred) ** 2)
    ss_tot = np.sum((true - np.mean(true)) ** 2)
    if ss_tot == 0:
        return np.nan
    return 1.0 - ss_res / ss_tot

def get_limits(true, pred, padding=0.05):
    true = np.asarray(true).reshape(-1)
    pred = np.asarray(pred).reshape(-1)
    data_min = min(np.min(true), np.min(pred))
    data_max = max(np.max(true), np.max(pred))
    data_range = data_max - data_min
    if data_range == 0:
        data_range = 1.0
    pad = padding * data_range
    return data_min - pad, data_max + pad

def beautify_axes(ax):
    # ax.grid(True, linestyle='-', linewidth=0.8, alpha=0.5, color='lightgray')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', colors='black')

def add_marginal_distributions(ax, true, pred, color, bins=35):
    divider = make_axes_locatable(ax)

    ax_top = divider.append_axes("top", size="18%", pad=0.06, sharex=ax)
    ax_right = divider.append_axes("right", size="18%", pad=0.06, sharey=ax)

    ax_top.hist(true, bins=bins, color=color, alpha=0.55, edgecolor='gray', linewidth=0.6)
    ax_right.hist(pred, bins=bins, orientation='horizontal',
                  color=color, alpha=0.55, edgecolor='gray', linewidth=0.6)

    # Top marginal
    ax_top.tick_params(axis='x', labelbottom=False, bottom=False)
    ax_top.tick_params(axis='y', left=False, labelleft=False)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['left'].set_visible(False)
    ax_top.grid(False)

    # Right marginal
    ax_right.tick_params(axis='y', labelleft=False, left=False)
    ax_right.tick_params(axis='x', bottom=False, labelbottom=False)
    ax_right.spines['top'].set_visible(False)
    ax_right.spines['right'].set_visible(False)
    ax_right.spines['bottom'].set_visible(False)
    ax_right.grid(False)

def add_residual_inset(ax, residuals, color):
    # Moved upward to avoid overlapping with main x-axis labels
    ax_inset = ax.inset_axes([0.60, 0.20, 0.28, 0.18])

    ax_inset.hist(residuals, bins=28, color=color, alpha=0.70, edgecolor='gray', linewidth=0.6)

    mean_res = np.mean(residuals)
    ax_inset.axvline(mean_res, color='k', linestyle='--', linewidth=1.0)

    x_left, x_right = ax_inset.get_xlim()
    y_bottom, y_top = ax_inset.get_ylim()
    x_shift = 0.04 * (x_right - x_left)

    ax_inset.text(mean_res + x_shift, y_top * 0.82, f'{mean_res:.2f}',
                  color=color, fontsize=9, fontweight='bold')

    ax_inset.set_xlabel('Residual', fontsize=9, labelpad=1)
    ax_inset.set_ylabel('')
    ax_inset.set_yticks([])
    ax_inset.tick_params(axis='x', labelsize=8, pad=1)

    ax_inset.spines['top'].set_visible(False)
    ax_inset.spines['right'].set_visible(False)
    ax_inset.spines['left'].set_visible(False)
    ax_inset.patch.set_alpha(0.0)

def plot_parity_with_marginals(ax, true, pred, xlabel, ylabel, color, mae_text, rmse_text):
    true = np.asarray(true).reshape(-1)
    pred = np.asarray(pred).reshape(-1)

    xmin, xmax = get_limits(true, pred, padding=0.05)

    ax.scatter(true, pred, s=35, c=color, alpha=0.35, edgecolors='none', rasterized=True)
    ax.plot([xmin, xmax], [xmin, xmax], color='grey', linestyle='--', linewidth=2.0)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(xmin, xmax)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    beautify_axes(ax)

    r2_val = r2_score_np(true, pred)
    ax.text(0.05, 0.96, f'$R^2 = {r2_val:.4f}$',
            transform=ax.transAxes, color='black', fontsize=11, va='top')
    ax.text(0.05, 0.86, mae_text,
            transform=ax.transAxes, color='black', fontsize=11, va='top')
    ax.text(0.05, 0.76, rmse_text,
            transform=ax.transAxes, color='black', fontsize=11, va='top')

    add_marginal_distributions(ax, true, pred, color, bins=35)

    residuals = pred - true
    add_residual_inset(ax, residuals, color)

# =========================
# Figure and colors - Changed to 1x3 layout
# =========================
fig, axs = plt.subplots(1, 3, figsize=(12, 4.0), dpi=100)

energy_color = '#1f77b4'   # blue
force_color = '#2ca02c'    # green
stress_color = '#ff7f0e'   # orange

# =========================
# Energy panel
# energy_data[:, 0] -> NEP
# energy_data[:, 1] -> DFT
# =========================
energy_true = energy_data[:, 1]
energy_pred = energy_data[:, 0]
energy_rmse = rmse(energy_pred, energy_true) * 1000.0
energy_mae = mae(energy_pred, energy_true) * 1000.0

plot_parity_with_marginals(
    axs[0],
    energy_true,
    energy_pred,
    'DFT energy (eV/atom)',
    'NEP energy (eV/atom)',
    energy_color,
    f'MAE = {energy_mae:.2f} meV/atom', 
    f'RMSE = {energy_rmse:.2f} meV/atom'
)

# =========================
# Force panel
# force_data[:, 0:3] -> NEP
# force_data[:, 3:6] -> DFT
# =========================
force_true = force_data[:, 3:6]
force_pred = force_data[:, 0:3]

force_true_flat = force_true.reshape(-1)
force_pred_flat = force_pred.reshape(-1)

force_rmse = rmse(force_pred_flat, force_true_flat) * 1000.0
force_mae = mae(force_pred_flat, force_true_flat) * 1000.0

plot_parity_with_marginals(
    axs[1],
    force_true_flat,
    force_pred_flat,
    r'DFT force (eV/$\mathrm{\AA}$)',
    r'NEP force (eV/$\mathrm{\AA}$)',
    force_color,
    rf'MAE = {force_mae:.2f} meV/$\mathrm{{\AA}}$',
    rf'RMSE = {force_rmse:.2f} meV/$\mathrm{{\AA}}$'
)

# =========================
# Stress panel
# stress_data[:, 0:6] -> NEP
# stress_data[:, 6:12] -> DFT
# =========================
stress_true = stress_data[:, 6:12]
stress_pred = stress_data[:, 0:6]

stress_true_flat = stress_true.reshape(-1)
stress_pred_flat = stress_pred.reshape(-1)

stress_rmse = rmse(stress_pred_flat, stress_true_flat)
stress_mae = mae(stress_pred_flat, stress_true_flat)

plot_parity_with_marginals(
    axs[2],
    stress_true_flat,
    stress_pred_flat,
    'DFT stress (GPa)',
    'NEP stress (GPa)',
    stress_color,
    f'MAE = {stress_mae:.4f} GPa',
    f'RMSE = {stress_rmse:.4f} GPa'
)

# =========================
# Layout and output
# =========================
abc = [f"({chr(97 + i)})" for i in range(3)]

for i, ax in enumerate(axs.flat):
    ax.text(-0.15, 1.05, abc[i],
            transform=ax.transAxes,
            fontsize=15,
            ha='left',
            va='bottom')

plt.tight_layout()
plt.subplots_adjust(wspace=0.2)

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('prediction.png', dpi=300, bbox_inches='tight')
else:
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Non-interactive backend detected. Plot saved as 'prediction.png'.")
        plt.savefig('prediction.png', dpi=300, bbox_inches='tight')
    else:
        plt.show()

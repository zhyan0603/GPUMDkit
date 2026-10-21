"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_train_test.py
Category:   Plot Scripts
Purpose:    Combined train/test parity plots for energy, force, and stress.
Usage:      gpumdkit.sh -plt train_test [save]
            python plt_train_test.py [save]
Arguments:
  save      Save the plot as 'train_test.png' instead of displaying it
Output:
  train_test.png  (if save is used, or if backend is non-interactive)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-20
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
TRAIN_COLOR = '#237B9F'
TEST_COLOR = '#EC817E'
REF_STYLE = dict(color='0.55', linestyle='--', linewidth=1.2, zorder=0)
KDE_LINE_CUTOFF = 0.005  # outline is drawn only where the KDE exceeds this fraction of its peak
MARGINAL_RATIO = 0.12  # marginal histogram size relative to the parity panel
PT_STYLE = dict(marker='o', linestyle='none', markersize=4.5, markeredgewidth=0, rasterized=True)

def panel_label(ax, text, x=-0.14, y=1.0):
    """Place a panel label outside the upper-left corner of an axes."""
    ax.text(x, y, text, transform=ax.transAxes, fontsize=12, va='bottom', ha='right')


def rmse(pred, true):
    return np.sqrt(np.mean((pred - true) ** 2))


def mae(pred, true):
    return np.mean(np.abs(pred - true))


def r2_score(true, pred):
    true = np.asarray(true).reshape(-1)
    pred = np.asarray(pred).reshape(-1)
    ss_res = np.sum((true - pred) ** 2)
    ss_tot = np.sum((true - np.mean(true)) ** 2)
    if ss_tot == 0:
        return np.nan
    return 1.0 - ss_res / ss_tot


def get_limits(*arrays, padding=0.05):
    values = np.concatenate([np.asarray(array).reshape(-1) for array in arrays])
    data_min = np.min(values)
    data_max = np.max(values)
    data_range = data_max - data_min
    if data_range == 0:
        data_range = 1.0
    pad = padding * data_range
    return data_min - pad, data_max + pad


def gaussian_kde_1d(values, grid):
    """Gaussian KDE (Scott bandwidth) on a uniform grid, in sample counts per grid bin.

    The data are binned onto the grid and smoothed with a Gaussian kernel, which
    keeps the cost independent of the number of samples. Returns None when a KDE
    is not meaningful (fewer than 3 finite samples or zero spread).
    """
    values = values[np.isfinite(values)]
    if values.size < 3 or np.std(values) == 0:
        return None
    dx = grid[1] - grid[0]
    bandwidth = 1.06 * np.std(values) * values.size ** (-0.2)
    sigma = max(bandwidth / dx, 1.0)
    half = int(4 * sigma) + 1
    kernel = np.exp(-0.5 * (np.arange(-half, half + 1) / sigma) ** 2)
    edges = np.append(grid - dx / 2, grid[-1] + dx / 2)
    counts, _ = np.histogram(values, bins=edges)
    return np.convolve(counts, kernel / kernel.sum(), mode="same")


def add_marginal_distributions(ax_top, ax_right, train_true, train_pred, test_true, test_pred, grid):
    """Draw KDE curves (histogram if a KDE is not possible) flush against the parity panel."""
    for ax_m, horizontal, datasets in (
        (ax_top, False, ((train_true, TRAIN_COLOR), (test_true, TEST_COLOR))),
        (ax_right, True, ((train_pred, TRAIN_COLOR), (test_pred, TEST_COLOR))),
    ):
        for values, color in datasets:
            density = gaussian_kde_1d(values, grid)
            if density is None:  # fall back to a histogram
                ax_m.hist(values, bins=np.linspace(grid[0], grid[-1], 36), color=color,
                          alpha=0.7, edgecolor="none",
                          orientation="horizontal" if horizontal else "vertical")
            else:
                # Hide the outline where the density is negligible so it does not trace the axes edge
                outline = np.where(density > KDE_LINE_CUTOFF * density.max(), density, np.nan)
                if horizontal:
                    ax_m.fill_betweenx(grid, density, color=color, alpha=0.35, linewidth=0)
                    ax_m.plot(outline, grid, color=color, linewidth=1.2)
                else:
                    ax_m.fill_between(grid, density, color=color, alpha=0.35, linewidth=0)
                    ax_m.plot(grid, outline, color=color, linewidth=1.2)

    ax_top.set_ylim(bottom=0)
    ax_right.set_xlim(left=0)
    for marginal in (ax_top, ax_right):
        marginal.axis("off")


def plot_parity(cell, train_true, train_pred, test_true, test_pred,
                xlabel, ylabel, metric_scale, metric_unit):
    """Draw one parity panel with marginal histograms inside the given subplot cell."""
    inner = cell.subgridspec(2, 2, width_ratios=[1, MARGINAL_RATIO],
                             height_ratios=[MARGINAL_RATIO, 1], hspace=0, wspace=0)
    ax = fig.add_subplot(inner[1, 0])
    ax_top = fig.add_subplot(inner[0, 0], sharex=ax)
    ax_right = fig.add_subplot(inner[1, 1], sharey=ax)
    train_true = np.asarray(train_true).reshape(-1)
    train_pred = np.asarray(train_pred).reshape(-1)
    test_true = np.asarray(test_true).reshape(-1)
    test_pred = np.asarray(test_pred).reshape(-1)

    xmin, xmax = get_limits(train_true, train_pred, test_true, test_pred)
    grid = np.linspace(xmin, xmax, 512)

    ax.plot([xmin, xmax], [xmin, xmax], **REF_STYLE)
    ax.plot(train_true, train_pred, color=TRAIN_COLOR, label="Train", **PT_STYLE)
    ax.plot(test_true, test_pred, color=TEST_COLOR, label="Test", **PT_STYLE)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(xmin, xmax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc="upper left", markerscale=1.0, handletextpad=0.2)

    train_r2 = r2_score(train_true, train_pred)
    test_r2 = r2_score(test_true, test_pred)
    train_mae = mae(train_pred, train_true) * metric_scale
    test_mae = mae(test_pred, test_true) * metric_scale
    train_rmse = rmse(train_pred, train_true) * metric_scale
    test_rmse = rmse(test_pred, test_true) * metric_scale

    metrics = (
        "Train / Test\n"
        f"$R^2$ = {train_r2:.3f} / {test_r2:.3f}\n"
        f"MAE = {train_mae:.2f} / {test_mae:.2f} {metric_unit}\n"
        f"RMSE = {train_rmse:.2f} / {test_rmse:.2f} {metric_unit}"
    )
    ax.text(0.96, 0.05, metrics, transform=ax.transAxes, fontsize=9,
            ha="right", va="bottom", linespacing=1.5)

    add_marginal_distributions(
        ax_top, ax_right, train_true, train_pred, test_true, test_pred, grid
    )
    return ax


energy_train = np.atleast_2d(np.loadtxt("energy_train.out"))
energy_test = np.atleast_2d(np.loadtxt("energy_test.out"))
force_train = np.atleast_2d(np.loadtxt("force_train.out"))
force_test = np.atleast_2d(np.loadtxt("force_test.out"))
stress_train = np.atleast_2d(np.loadtxt("stress_train.out"))
stress_test = np.atleast_2d(np.loadtxt("stress_test.out"))

valid_rows_train = (
    np.isfinite(stress_train[:, :12]).all(axis=1)
    & np.all(np.abs(stress_train[:, :12]) < 1e6, axis=1)
)
valid_rows_test = (
    np.isfinite(stress_test[:, :12]).all(axis=1)
    & np.all(np.abs(stress_test[:, :12]) < 1e6, axis=1)
)
stress_train = stress_train[valid_rows_train]
stress_test = stress_test[valid_rows_test]

fig = plt.figure(figsize=(12, 3.5))
grid = fig.add_gridspec(1, 3)
axs = []

axs.append(plot_parity(
    grid[0],
    energy_train[:, 1], energy_train[:, 0],
    energy_test[:, 1], energy_test[:, 0],
    "DFT Energy (eV/atom)", "NEP Energy (eV/atom)",
    1000.0, "meV/atom",
))
axs.append(plot_parity(
    grid[1],
    force_train[:, 3:6].reshape(-1), force_train[:, 0:3].reshape(-1),
    force_test[:, 3:6].reshape(-1), force_test[:, 0:3].reshape(-1),
    r"DFT Force (eV/$\mathrm{\AA}$)",
    r"NEP Force (eV/$\mathrm{\AA}$)",
    1000.0, "meV/Å",
))
axs.append(plot_parity(
    grid[2],
    stress_train[:, 6:12].reshape(-1), stress_train[:, 0:6].reshape(-1),
    stress_test[:, 6:12].reshape(-1), stress_test[:, 0:6].reshape(-1),
    "DFT Stress (GPa)", "NEP Stress (GPa)",
    1.0, "GPa",
))

panel_labels = ["(a)", "(b)", "(c)"]
for label, ax in zip(panel_labels, axs):
    panel_label(ax, label, y=1.0)

# Subplot geometry (figure fractions); tune these to adjust margins and spacing
fig.subplots_adjust(left=0.06, right=0.985, top=0.96, bottom=0.16, wspace=0.30)

if len(sys.argv) > 1 and sys.argv[1] == "save":
    plt.savefig("train_test.png", dpi=300)
else:
    from matplotlib import get_backend

    if get_backend().lower() in ["agg", "cairo", "pdf", "ps", "svg"]:
        print("Non-interactive backend detected. Plot saved as 'train_test.png'.")
        plt.savefig("train_test.png", dpi=300)
    else:
        plt.show()

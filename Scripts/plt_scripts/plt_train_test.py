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
Last-modified: 2026-09-30
=============================================================================
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "axes.unicode_minus": False,
    "axes.labelsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.linewidth": 1.5,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.labelcolor": "black",
    "text.color": "black",
})

TRAIN_COLOR = "#237B9F"
TEST_COLOR = "#EC817E"


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


def beautify_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(top=False, right=False)
    ax.tick_params(axis="both", colors="black")


def add_marginal_distributions(ax, train_true, train_pred, test_true, test_pred, bins):
    divider = make_axes_locatable(ax)
    ax_top = divider.append_axes("top", size="18%", pad=0, sharex=ax)
    ax_right = divider.append_axes("right", size="18%", pad=0, sharey=ax)

    for values, color in ((train_true, TRAIN_COLOR), (test_true, TEST_COLOR)):
        ax_top.hist(values, bins=bins, color=color, alpha=0.45,
                    edgecolor="gray", linewidth=0.5)
    for values, color in ((train_pred, TRAIN_COLOR), (test_pred, TEST_COLOR)):
        ax_right.hist(values, bins=bins, orientation="horizontal", color=color,
                      alpha=0.45, edgecolor="gray", linewidth=0.5)

    ax_top.tick_params(axis="both", which="both", length=0,
                       labelbottom=False, labelleft=False)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.spines["left"].set_visible(False)
    ax_top.grid(False)

    ax_right.tick_params(axis="both", which="both", length=0,
                         labelbottom=False, labelleft=False)
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.spines["bottom"].set_visible(False)
    ax_right.grid(False)


def plot_parity(ax, train_true, train_pred, test_true, test_pred,
                xlabel, ylabel, metric_scale, metric_unit):
    train_true = np.asarray(train_true).reshape(-1)
    train_pred = np.asarray(train_pred).reshape(-1)
    test_true = np.asarray(test_true).reshape(-1)
    test_pred = np.asarray(test_pred).reshape(-1)

    xmin, xmax = get_limits(train_true, train_pred, test_true, test_pred)
    bins = np.linspace(xmin, xmax, 36)

    ax.scatter(train_true, train_pred, s=35, c=TRAIN_COLOR, alpha=0.35,
               edgecolors="none", rasterized=True, label="Train")
    ax.scatter(test_true, test_pred, s=35, c=TEST_COLOR, alpha=0.35,
               edgecolors="none", rasterized=True, label="Test")
    ax.plot([xmin, xmax], [xmin, xmax], color="grey",
            linestyle="--", linewidth=2.0)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(xmin, xmax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    beautify_axes(ax)
    ax.legend(frameon=False, loc="upper left", fontsize=10)

    train_r2 = r2_score(train_true, train_pred)
    test_r2 = r2_score(test_true, test_pred)
    train_mae = mae(train_pred, train_true) * metric_scale
    test_mae = mae(test_pred, test_true) * metric_scale
    train_rmse = rmse(train_pred, train_true) * metric_scale
    test_rmse = rmse(test_pred, test_true) * metric_scale

    metrics = (
        "Train / Test\n"
        f"$R^2$: {train_r2:.3f} / {test_r2:.3f}\n"
        f"MAE: {train_mae:.2f} / {test_mae:.2f} {metric_unit}\n"
        f"RMSE: {train_rmse:.2f} / {test_rmse:.2f} {metric_unit}"
    )
    ax.text(0.95, 0.04, metrics, transform=ax.transAxes, fontsize=10,
            ha="right", va="bottom", linespacing=1.25)

    add_marginal_distributions(
        ax, train_true, train_pred, test_true, test_pred, bins
    )


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

fig, axs = plt.subplots(1, 3, figsize=(12, 4.0), dpi=100)

plot_parity(
    axs[0],
    energy_train[:, 1], energy_train[:, 0],
    energy_test[:, 1], energy_test[:, 0],
    "DFT energy (eV/atom)", "NEP energy (eV/atom)",
    1000.0, "meV/atom",
)
plot_parity(
    axs[1],
    force_train[:, 3:6].reshape(-1), force_train[:, 0:3].reshape(-1),
    force_test[:, 3:6].reshape(-1), force_test[:, 0:3].reshape(-1),
    r"DFT force (eV/$\mathrm{\AA}$)",
    r"NEP force (eV/$\mathrm{\AA}$)",
    1000.0, "meV/Å",
)
plot_parity(
    axs[2],
    stress_train[:, 6:12].reshape(-1), stress_train[:, 0:6].reshape(-1),
    stress_test[:, 6:12].reshape(-1), stress_test[:, 0:6].reshape(-1),
    "DFT stress (GPa)", "NEP stress (GPa)",
    1.0, "GPa",
)

panel_labels = ["(a)", "(b)", "(c)"]
for label, ax in zip(panel_labels, axs.flat):
    ax.text(-0.15, 1.05, label, transform=ax.transAxes,
            fontsize=15, ha="left", va="bottom")

plt.tight_layout()
plt.subplots_adjust(wspace=0.2)

if len(sys.argv) > 1 and sys.argv[1] == "save":
    plt.savefig("train_test.png", dpi=300, bbox_inches="tight")
else:
    from matplotlib import get_backend

    if get_backend().lower() in ["agg", "cairo", "pdf", "ps", "svg"]:
        print("Non-interactive backend detected. Plot saved as 'train_test.png'.")
        plt.savefig("train_test.png", dpi=300, bbox_inches="tight")
    else:
        plt.show()

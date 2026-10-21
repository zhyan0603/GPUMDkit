"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_born_charge.py
Category:   Plot Scripts
Purpose:    Train/test parity plot for Born effective charges (BEC).
            Structures with all-zero reference BEC are filtered out.
Usage:      gpumdkit.sh -plt born_charge
            python plt_born_charge.py [save]
Arguments:
  save      Save the plot as 'bec.png' instead of displaying it
Output:
  bec.png  (if save is used, or if backend is non-interactive)
Author:     Denan LI (lidenan@westlake.edu.cn)
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


def calculate_rmse(pred, true):
    return np.sqrt(np.mean((pred - true) ** 2))


def calculate_mae(pred, true):
    return np.mean(np.abs(pred - true))


def calculate_r2(true, pred):
    true = np.asarray(true).reshape(-1)
    pred = np.asarray(pred).reshape(-1)
    ss_res = np.sum((true - pred) ** 2)
    ss_tot = np.sum((true - np.mean(true)) ** 2)
    if ss_tot == 0:
        return np.nan
    return 1.0 - ss_res / ss_tot


def calculate_limits(*arrays, padding=0.05):
    values = np.concatenate([np.asarray(array).reshape(-1) for array in arrays])
    data_min = np.min(values)
    data_max = np.max(values)
    data_range = data_max - data_min
    if data_range == 0:
        data_range = 1.0
    pad = padding * data_range
    return data_min - pad, data_max + pad


def load_bec(filename):
    data = np.loadtxt(filename)
    pred = data[:, :9]
    target = data[:, 9:]
    mask = ~np.all(target == 0, axis=1)
    return pred[mask], target[mask]


def beautify_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(top=False, right=False)
    ax.tick_params(axis="both", colors="black")


def add_marginal_distributions(ax, datasets, bins):
    divider = make_axes_locatable(ax)
    ax_top = divider.append_axes("top", size="18%", pad=0, sharex=ax)
    ax_right = divider.append_axes("right", size="18%", pad=0, sharey=ax)

    for _, true, _, color in datasets:
        ax_top.hist(true, bins=bins, color=color, alpha=0.45,
                    edgecolor="gray", linewidth=0.5)
    for _, _, pred, color in datasets:
        ax_right.hist(pred, bins=bins, orientation="horizontal", color=color,
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


def calculate_metrics(true, pred):
    return (
        calculate_r2(true, pred),
        calculate_mae(pred, true),
        calculate_rmse(pred, true),
    )


def format_metrics(datasets):
    if len(datasets) == 2:
        train_metrics = calculate_metrics(datasets[0][1], datasets[0][2])
        test_metrics = calculate_metrics(datasets[1][1], datasets[1][2])
        return (
            "Train / Test\n"
            f"$R^2$: {train_metrics[0]:.3f} / {test_metrics[0]:.3f}\n"
            f"MAE: {train_metrics[1]:.2f} / {test_metrics[1]:.2f} e\n"
            f"RMSE: {train_metrics[2]:.2f} / {test_metrics[2]:.2f} e"
        )

    train_metrics = calculate_metrics(datasets[0][1], datasets[0][2])
    return (
        "Train\n"
        f"$R^2$: {train_metrics[0]:.3f}\n"
        f"MAE: {train_metrics[1]:.2f} e\n"
        f"RMSE: {train_metrics[2]:.2f} e"
    )


pred_train, tgt_train = load_bec("bec_train.out")
datasets = [
    ("Train", tgt_train.ravel(), pred_train.ravel(), TRAIN_COLOR),
]
try:
    pred_test, tgt_test = load_bec("bec_test.out")
except (OSError, IOError):
    pass
else:
    datasets.append(("Test", tgt_test.ravel(), pred_test.ravel(), TEST_COLOR))

xmin, xmax = calculate_limits(*[values for _, true, pred, _ in datasets
                                for values in (true, pred)])
bins = np.linspace(xmin, xmax, 36)

fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
for label, true, pred, color in datasets:
    ax.scatter(true, pred, s=35, c=color, alpha=0.35, edgecolors="none",
               rasterized=True, label=label)

ax.plot([xmin, xmax], [xmin, xmax], color="grey", linestyle="--", linewidth=2.0)
ax.set_xlim(xmin, xmax)
ax.set_ylim(xmin, xmax)
ax.set_xlabel("DFT BEC (e)")
ax.set_ylabel("NEP BEC (e)")
beautify_axes(ax)
ax.legend(frameon=False, loc="upper left", fontsize=10)
ax.text(0.95, 0.04, format_metrics(datasets), transform=ax.transAxes,
        fontsize=10, ha="right", va="bottom", linespacing=1.25)

add_marginal_distributions(ax, datasets, bins)
plt.tight_layout()

if len(sys.argv) > 1 and sys.argv[1] == "save":
    plt.savefig("bec.png", dpi=300, bbox_inches="tight")
else:
    from matplotlib import get_backend

    if get_backend().lower() in ["agg", "cairo", "pdf", "ps", "svg"]:
        print("Non-interactive backend detected. Plot saved as 'bec.png'.")
        plt.savefig("bec.png", dpi=300, bbox_inches="tight")
    else:
        plt.show()

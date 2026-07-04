"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_born_charge.py
Category:   Plot Scripts
Purpose:    Parity plots for Born effective charges (BEC) on training and
            testing datasets. Structures with all-zero reference BEC are
            filtered out.
Usage:      gpumdkit.sh -plt born_charge
            python plt_born_charge.py
Output:
  Display of BEC parity plot
Author:     Denan LI (lidenan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
})

train_color = '#237B9F'
test_color = '#EC817E'


def calculate_rmse(pred, actual):
    return np.sqrt(np.mean((pred - actual) ** 2))


def calculate_limits(train_data, test_data, padding=0.08):
    data_min = np.min(np.concatenate((train_data, test_data)))
    data_max = np.max(np.concatenate((train_data, test_data)))
    data_range = data_max - data_min
    return data_min - padding * data_range, data_max + padding * data_range


def load_bec(filename):
    data = np.loadtxt(filename)
    pred = data[:, :9]
    target = data[:, 9:]
    mask = ~np.all(target == 0, axis=1)
    return pred[mask], target[mask]


pred_train, tgt_train = load_bec('bec_train.out')
try:
    pred_test, tgt_test = load_bec('bec_test.out')
    has_test = True
except (OSError, IOError):
    has_test = False

fig, ax = plt.subplots(figsize=(5, 5), dpi=100)

x_train = tgt_train.ravel()
y_train = pred_train.ravel()

if has_test:
    x_test = tgt_test.ravel()
    y_test = pred_test.ravel()
    xmin, xmax = calculate_limits(x_train, x_test)
else:
    xmin, xmax = calculate_limits(x_train, x_train)

ax.set_xlim(xmin, xmax)
ax.set_ylim(xmin, xmax)

ax.plot(x_train, y_train, '.', markersize=10, label='Train', color=train_color)
if has_test:
    ax.plot(x_test, y_test, '.', markersize=10, label='Test', color=test_color)
ax.plot([xmin, xmax], [xmin, xmax], linewidth=1, color='red')

ax.set_xlabel('DFT BEC (e)', fontsize=10)
ax.set_ylabel('NEP BEC (e)', fontsize=10)
ax.tick_params(axis='both', labelsize=10)
ax.legend(frameon=False, loc='upper left')

train_rmse = calculate_rmse(y_train, x_train)
ax.text(0.35, 0.15, f'RMSE (Train): {train_rmse:.4f} e', transform=ax.transAxes, fontsize=10, verticalalignment='center')
if has_test:
    test_rmse = calculate_rmse(y_test, x_test)
    ax.text(0.35, 0.07, f'RMSE (Test): {test_rmse:.4f} e', transform=ax.transAxes, fontsize=10, verticalalignment='center')

plt.tight_layout()

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('bec.png', dpi=300)
else:
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'bec_parity.png'.")
        plt.savefig('bec.png', dpi=300)
    else:
        plt.show()

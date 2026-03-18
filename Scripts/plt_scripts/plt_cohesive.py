"""
Plot GPUMD cohesive.out: isotropic scaling factor vs cohesive energy.

File format:
- column 1: scaling factor (dimensionless)
- column 2: potential energy (eV) after minimization

Author: Qilin Guo (guoqilin@buaa.edu.cn)
Modified by Zihan Yan (yanzihan@westlake.edu.cn)
Last modified: 2026-03-18

Usage:
  python plot_cohesive.py [save]
"""

import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# Backend Configuration
# ---------------------------------------------------------
backend = matplotlib.get_backend()
if not os.environ.get("DISPLAY") and "agg" not in str(backend).lower():
    matplotlib.use("Agg", force=True)

# ---------------------------------------------------------
# Data Loading
# ---------------------------------------------------------
in_path = "cohesive.out"

if not os.path.exists(in_path):
    print(f"Missing input file: {in_path}", file=sys.stderr)
    sys.exit(1)

data = np.loadtxt(in_path)
if data.ndim == 1:
    data = data.reshape(1, -1)
if data.shape[1] < 2:
    raise ValueError(f"{in_path} must have at least 2 columns")

scale = data[:, 0]
energy = data[:, 1]

# ---------------------------------------------------------
# PRL-style Plotting Parameters
# ---------------------------------------------------------
plt.rcParams.update({
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 10,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "axes.linewidth": 1.0,
})

# Fit a parabola E(a) ~ a*s^2 + b*s + c near the minimum
x = scale
y = energy
x0_theory = 1.0
mask = np.isfinite(x) & np.isfinite(y)
x_fit = x[mask]
y_fit = y[mask]

if x_fit.size >= 5:
    order = np.argsort(np.abs(x_fit - x0_theory))
    k = max(5, x_fit.size // 2)
    idx = order[:k]
    x_fit = x_fit[idx]
    y_fit = y_fit[idx]

x_min_fit = None
if x_fit.size >= 3:
    coeffs = np.polyfit(x_fit, y_fit, deg=2)
    a_coeff, b_coeff, c_coeff = coeffs
    if a_coeff != 0.0:
        x_min_fit = -b_coeff / (2.0 * a_coeff)

# 4:3 aspect ratio, pure white background
fig, ax = plt.subplots(figsize=(4, 3))

# Raw data: C0 hollow circles
ax.plot(scale, energy, marker="o", ms=4.5, mfc='none', mec='C0', mew=1.0, 
        linestyle='none', label="Data")

# Fitted curve: C1 solid line
if x_fit.size >= 3 and x_min_fit is not None:
    x_dense = np.linspace(x.min(), x.max(), 200)
    y_dense = np.polyval(coeffs, x_dense)
    ax.plot(x_dense, y_dense, color="C1", lw=1.5, ls="-", label="Fit")

    # Display minimum position in the bottom right corner
    ax.text(0.95, 0.05, rf"$a_{{\rm min}} = {x_min_fit:.4f}$", 
            transform=ax.transAxes, ha='right', va='bottom', fontsize=10)

# Axes labels
ax.set_xlabel(r"Isotropic scaling factor $a$")
ax.set_ylabel(r"Potential energy $E$ (eV)")

# Set Y-axis to scientific notation
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
ax.yaxis.get_offset_text().set_fontsize(9)

# Simple legend
ax.legend(loc="upper center", frameon=False, handlelength=1.5)

plt.tight_layout(pad=0.5)

# ---------------------------------------------------------
# Output Logic
# ---------------------------------------------------------
if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('cohesive.png', dpi=300, facecolor='white')
    print("The plot has been saved as 'cohesive.png'.")
else:
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'cohesive.png'.")
        plt.savefig('cohesive.png', dpi=300, facecolor='white')
    else:
        plt.show()
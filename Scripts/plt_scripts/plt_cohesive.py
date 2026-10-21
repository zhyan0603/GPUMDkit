"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_cohesive.py
Category:   Plot Scripts
Purpose:    Plot cohesive energy curve from GPUMD cohesive.out: isotropic
            scaling factor vs cohesive energy, with the minimum marked.
Usage:      gpumdkit.sh -plt cohesive [save]
            python plt_cohesive.py [save]
Arguments:
  save      Save the plot as 'cohesive.png' instead of displaying it
Output:
  cohesive.png  (if save is used)
Author:     Qilin Guo (guoqilin@buaa.edu.cn), Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-21
=============================================================================
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

# Position of the energy minimum (the minimum of the sampled curve)
mask = np.isfinite(scale) & np.isfinite(energy)
i_min = np.nanargmin(np.where(mask, energy, np.nan))
x_min = scale[i_min]

fig, ax = plt.subplots(figsize=(4.0, 3.2))

# Raw data: C0 hollow circles
ax.plot(scale, energy, marker="o", ms=4.5, mfc='none', mec='C0', mew=1.0,
        linestyle='none')

# Vertical dashed line and label at the minimum
ax.axvline(x_min, color='0.55', linestyle='--', linewidth=1.0, zorder=0)
ax.annotate(rf"$a_{{\rm min}} = {x_min:.4f}$", xy=(x_min, 0.95),
            xycoords=ax.get_xaxis_transform(), xytext=(6, 0), textcoords='offset points',
            ha='left', va='top', fontsize=9)

# Axes labels
ax.set_xlabel(r"Isotropic scaling factor $a$")
ax.set_ylabel(r"Potential energy $E$ (eV)")

# Set Y-axis to scientific notation
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
ax.yaxis.get_offset_text().set_fontsize(9)

# Subplot geometry (figure fractions); tune these to adjust margins
fig.subplots_adjust(left=0.19, right=0.97, top=0.93, bottom=0.17)

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

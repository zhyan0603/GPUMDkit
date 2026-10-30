"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_nep_restart.py
Category:   Plot Scripts
Purpose:    Visualize parameters stored in nep.restart file from NEP training.
Usage:      gpumdkit.sh -plt restart
            python plt_nep_restart.py
Output:
  Display of NEP restart parameters plot
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-28
=============================================================================
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

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

# nep.restart stores, for every parameter, the mean (mu) and the standard deviation (sigma)
# of the SNES search distribution (column 1 and column 2, respectively)
# Function to read two-column data from a file
def read_data(file_path):
    with open(file_path) as f:
        skiprows = 3 if f.readline().strip() == "nep_restart" else 0
    data = np.loadtxt(file_path, skiprows=skiprows, ndmin=2)
    return data[:, 0], data[:, 1]

# Read data from files
col1, col2 = read_data('./nep.restart')

# Create subplots
fig, axes = plt.subplots(2, 1, figsize=(5.6, 4.0), sharex=True)

# First subplot
axes[0].plot(range(len(col1)), col1, label='Mean', color='C0', linewidth=0.8)
axes[0].hlines(0, 0, len(col1), color='0.55', linestyle='--', linewidth=1.0, zorder=0)
axes[0].set_ylabel(r'$\mu$')
axes[0].legend(loc='upper left', handlelength=1.4)
axes[0].set_title('Parameters in nep.restart')
axes[0].margins(x=0.01)
#axes[0].grid(True)

# Second subplot
axes[1].plot(range(len(col2)), col2, label='Std. dev.', color='C1', linewidth=0.8)
axes[1].set_ylabel(r'$\sigma$')
axes[1].set_xlabel('Index')
axes[1].legend(loc='upper left', handlelength=1.4)
#axes[1].grid(True)

# Adjust layout
# Subplot geometry (figure fractions); tune these to adjust margins and spacing
fig.subplots_adjust(left=0.125, right=0.975, top=0.93, bottom=0.12, hspace=0.0)
fig.align_ylabels(axes)

# Show the plot
if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('nep_restart.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'nep_restart.png'.")
        plt.savefig('nep_restart.png', dpi=300)
    else:
        plt.show()

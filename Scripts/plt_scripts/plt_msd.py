"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_msd.py
Category:   Plot Scripts
Purpose:    Plot mean square displacement (MSD) from msd.out for all
            directions (x, y, z) with slope annotations.
Usage:      gpumdkit.sh -plt msd [save]
            python plt_msd.py [save]
Arguments:
  save      Save the plot as 'msd.png' instead of displaying it
Output:
  msd.png   (if save is used, or if backend is non-interactive)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-20
=============================================================================
"""

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
COLORS = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']

# Function to read data from the given file
def read_data(file_name):
    data = np.loadtxt(file_name)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# Function to calculate the slope over the middle 40%-80% of the data
def calculate_slope(x, y):
    len_data = len(x)
    start_idx = int(0.4 * len_data)  # 40% of the data
    end_idx = int(0.8 * len_data)    # 80% of the data
    x_range, y_range = x[start_idx:end_idx], y[start_idx:end_idx]
    coeffs = np.polyfit(x_range, y_range, 1)  # Linear fit
    return coeffs[0]  # Slope of the linear fit

# Determine the input file
#input_file = sys.argv[1] if len(sys.argv) > 1 else './msd.out'
input_file = './msd.out'

# Read the data from the file
time, msd_x, msd_y, msd_z = read_data(input_file)

slope_x = calculate_slope(time, msd_x)
slope_y = calculate_slope(time, msd_y)
slope_z = calculate_slope(time, msd_z)

# Calculate the total MSD and mean MSD
#msd_total = msd_x + msd_y + msd_z
#msd_mean = (msd_x + msd_y + msd_z)/3

# Create a plot for the MSD data
fig, ax = plt.subplots(figsize=(4.2, 3.5))
n_fit = len(time)
fit_start, fit_end = time[int(0.4 * n_fit)], time[int(0.8 * n_fit) - 1]  # slope-fit window (40%-80%)
ax.axvspan(fit_start, fit_end, color='0.5', alpha=0.08, linewidth=0)
for t_edge, text in ((fit_start, '40%'), (fit_end, '80%')):
    ax.axvline(t_edge, color='0.6', linestyle='--', linewidth=0.8, zorder=0)
    ax.text(t_edge, 0.98, text + ' ', transform=ax.get_xaxis_transform(), fontsize=8, color='0.4', ha='right', va='top')
for label, msd, slope, color in (('x', msd_x, slope_x, COLORS[0]),
                                 ('y', msd_y, slope_y, COLORS[1]),
                                 ('z', msd_z, slope_z, COLORS[2])):
    ax.plot(time, msd, color=color, linewidth=1.6, label=f'{label}: {slope:.3f}')
#ax.plot(time, msd_mean, label='mean', color='C4')
#ax.plot(time, msd_total, label='total', color='C5')

# Add labels
ax.set_xlabel('dt (ps)')
ax.set_ylabel(r'MSD ($\mathrm{\AA}^2$)')
ax.margins(x=0.02)
ax.legend(loc='upper left', title=r'Slope ($\mathrm{\AA}^2$/ps)', title_fontsize=9)

# Subplot geometry (figure fractions); tune these to adjust margins and spacing
fig.subplots_adjust(left=0.135, right=0.975, top=0.97, bottom=0.135)

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('msd.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'msd.png'.")
        plt.savefig('msd.png', dpi=300)
    else:
        plt.show()

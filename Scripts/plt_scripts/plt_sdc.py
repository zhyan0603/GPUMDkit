"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_sdc.py
Category:   Plot Scripts
Purpose:    Plot self-diffusion coefficient (SDC) vs time from msd.out with
            an inset showing the last 80% of data with moving average overlay.
Usage:      gpumdkit.sh -plt sdc
            python plt_sdc.py [save]
Arguments:
  save      Save the plot as 'SDC.png' instead of displaying it
Output:
  SDC.png   (if save is used, or if backend is non-interactive)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-21
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
COLORS = ['C0', 'C1', 'C2']

# Function to read data from file
def read_data(file_name):
    data = np.loadtxt(file_name)
    return data[:, 0], data[:, 4], data[:, 5], data[:, 6]

input_file = './msd.out'

# Read the data
time, sdc_x, sdc_y, sdc_z = read_data(input_file)

# Create main figure
fig, ax = plt.subplots(figsize=(4.2, 3.5))

# Main plot: full range, original data
ax.plot(time, sdc_x, label='x', color='C0', linewidth=1.0)
ax.plot(time, sdc_y, label='y', color='C1', linewidth=1.0)
ax.plot(time, sdc_z, label='z', color='C2', linewidth=1.0)

# Legend: one row at the top, where the curves leave room
ax.legend(ncol=3, loc='upper right', handlelength=1.4, columnspacing=1.0, handletextpad=0.5)

ax.set_xlabel('dt (ps)')
ax.set_ylabel(r'SDC ($\mathrm{\AA}^2$/ps)')

# ────────────────────────────────────────────────
# Inset: last 80% of the data + smoothing overlay
# ────────────────────────────────────────────────

n = len(time)
zoom_start = int(0.2 * n)           # start index for last 80%
time_zoom = time[zoom_start:]

sdc_x_zoom = sdc_x[zoom_start:]
sdc_y_zoom = sdc_y[zoom_start:]
sdc_z_zoom = sdc_z[zoom_start:]

# Create inset axes - as requested
ax_inset = ax.inset_axes([0.25, 0.25, 0.65, 0.65])

# Moving average parameters
window_frac = 0.01                  # ~1% of total data length
window_size = max(5, int(window_frac * n))

# Helper function to compute moving average with proper edge padding
def moving_average(data, window):
    if window <= 1:
        return data.copy()
    pad_width = window // 2
    padded = np.pad(data, (pad_width, pad_width), mode='edge')
    kernel = np.ones(window) / window
    smoothed = np.convolve(padded, kernel, mode='valid')
    # Ensure output length matches input (in case of rounding issues)
    return smoothed[:len(data)]

# Compute smoothed versions for inset only
sdc_x_smooth = moving_average(sdc_x_zoom, window_size)
sdc_y_smooth = moving_average(sdc_y_zoom, window_size)
sdc_z_smooth = moving_average(sdc_z_zoom, window_size)

# Plot raw data (transparent) in inset
ax_inset.plot(time_zoom, sdc_x_zoom, color='C0', alpha=0.35, lw=0.8)
ax_inset.plot(time_zoom, sdc_y_zoom, color='C1', alpha=0.35, lw=0.8)
ax_inset.plot(time_zoom, sdc_z_zoom, color='C2', alpha=0.35, lw=0.8)

# Overlay smoothed lines (solid, opaque)
ax_inset.plot(time_zoom, sdc_x_smooth, color='C0', lw=1.4, label='x smooth')
ax_inset.plot(time_zoom, sdc_y_smooth, color='C1', lw=1.4, label='y smooth')
ax_inset.plot(time_zoom, sdc_z_smooth, color='C2', lw=1.4, label='z smooth')

# Subplot geometry (figure fractions); tune these to adjust margins
fig.subplots_adjust(left=0.15, right=0.975, top=0.97, bottom=0.135)

# Save or show
if len(sys.argv) > 1 and sys.argv[1] == 'save':
    fig.savefig('sdc.png', dpi=300)
else:
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to non-interactive backend.")
        print("Plot saved automatically as 'sdc.png'.")
        fig.savefig('sdc.png', dpi=300)
    else:
        plt.show()
"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_msd_sdc.py
Category:   Plot Scripts
Purpose:    Plot MSD and SDC data side-by-side from msd.out with slope
            annotations and moving average overlay.
Usage:      gpumdkit.sh -plt msd_sdc [save]
            python plt_msd_sdc.py [save]
Arguments:
  save      Save the plot as 'MSD_SDC.png' instead of displaying it
Output:
  MSD_SDC.png  (if save is used, or if backend is non-interactive)
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

# Function to read all required data from file
def read_data(file_name):
    data = np.loadtxt(file_name)
    # time, msd(x,y,z), sdc(x,y,z)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6]

# Function to calculate the slope over the middle 40%-80% of the data
def calculate_slope(x, y):
    len_data = len(x)
    start_idx = int(0.4 * len_data)  # 40% of the data
    end_idx = int(0.8 * len_data)    # 80% of the data
    x_range, y_range = x[start_idx:end_idx], y[start_idx:end_idx]
    coeffs = np.polyfit(x_range, y_range, 1)  # Linear fit
    return coeffs[0]  # Slope of the linear fit

# Helper function to compute moving average with proper edge padding
def moving_average(data, window):
    if window <= 1:
        return data.copy()
    pad_width = window // 2
    padded = np.pad(data, (pad_width, pad_width), mode='edge')
    kernel = np.ones(window) / window
    smoothed = np.convolve(padded, kernel, mode='valid')
    return smoothed[:len(data)]

input_file = './msd.out'

# Read the data
time, msd_x, msd_y, msd_z, sdc_x, sdc_y, sdc_z = read_data(input_file)

# Create a figure with 1 row and 2 columns
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.3))

# ==========================================
# Left Subplot (ax1): MSD vs dt
# ==========================================
slope_x = calculate_slope(time, msd_x)
slope_y = calculate_slope(time, msd_y)
slope_z = calculate_slope(time, msd_z)

n_fit = len(time)
fit_start, fit_end = time[int(0.4 * n_fit)], time[int(0.8 * n_fit) - 1]  # slope-fit window (40%-80%)
ax1.axvspan(fit_start, fit_end, color='0.5', alpha=0.08, linewidth=0)
for t_edge, text in ((fit_start, '40%'), (fit_end, '80%')):
    ax1.axvline(t_edge, color='0.6', linestyle='--', linewidth=0.8, zorder=0)
    ax1.text(t_edge, 0.98, text + ' ', transform=ax1.get_xaxis_transform(), fontsize=8, color='0.4', ha='right', va='top')
for label, msd, slope, color in (('x', msd_x, slope_x, COLORS[0]),
                                 ('y', msd_y, slope_y, COLORS[1]),
                                 ('z', msd_z, slope_z, COLORS[2])):
    ax1.plot(time, msd, color=color, linewidth=1.6, label=f'{label}: {slope:.3f}')

ax1.set_xlabel('dt (ps)')
ax1.set_ylabel(r'MSD ($\mathrm{\AA}^2$)')
ax1.margins(x=0.02)
ax1.legend(loc='upper left', title=r'Slope ($\mathrm{\AA}^2$/ps)', title_fontsize=9)

# ==========================================
# Right Subplot (ax2): SDC vs dt
# ==========================================
ax2.plot(time, sdc_x, label='x', color='C0', linewidth=1.0)
ax2.plot(time, sdc_y, label='y', color='C1', linewidth=1.0)
ax2.plot(time, sdc_z, label='z', color='C2', linewidth=1.0)

ax2.set_xlabel('dt (ps)')
ax2.set_ylabel(r'SDC ($\mathrm{\AA}^2$/ps)')
ax2.legend(ncol=3, loc='upper right', handlelength=1.4, columnspacing=1.0, handletextpad=0.5)

# ────────────────────────────────────────────────
# Inset (ax2): last 80% of the data + smoothing overlay
# ────────────────────────────────────────────────
n = len(time)
zoom_start = int(0.2 * n)           # start index for last 80%
time_zoom = time[zoom_start:]

sdc_x_zoom = sdc_x[zoom_start:]
sdc_y_zoom = sdc_y[zoom_start:]
sdc_z_zoom = sdc_z[zoom_start:]

# Create inset axes inside ax2
ax_inset = ax2.inset_axes([0.25, 0.25, 0.65, 0.65])

# Moving average parameters
window_frac = 0.01                  # ~1% of total data length
window_size = max(5, int(window_frac * n))

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

# ==========================================
# Final Layout and Save/Show logic
# ==========================================
# Subplot geometry (figure fractions); tune these to adjust margins and spacing
fig.subplots_adjust(left=0.085, right=0.985, top=0.97, bottom=0.15, wspace=0.28)

save_name = 'msd_sdc.png'

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    fig.savefig(save_name, dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print(f"The plot has been automatically saved as '{save_name}'.")
        fig.savefig(save_name, dpi=300)
    else:
        plt.show()

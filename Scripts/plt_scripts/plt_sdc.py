"""
This script reads SDC data from 'msd.out', plots the full range of SDC vs time, 
and includes an inset showing the last 80% of the data with a moving average overlay. 

Author: Zihan YAN (yanzihan@westlake.edu.cn)
Last modified: 2026-03-19

Usage:
    python plt_sdc.py [save]
"""    

import sys
import numpy as np
import matplotlib.pyplot as plt

# Function to read data from file
def read_data(file_name):
    data = np.loadtxt(file_name)
    return data[:, 0], data[:, 4], data[:, 5], data[:, 6]

input_file = './msd.out'

# Read the data
time, sdc_x, sdc_y, sdc_z = read_data(input_file)

# Create main figure
fig, ax = plt.subplots(figsize=(4.5, 3.5), dpi=150)

# Main plot: full range, original data
ax.plot(time, sdc_x, label='x', color='C0')
ax.plot(time, sdc_y, label='y', color='C1')
ax.plot(time, sdc_z, label='z', color='C2')

# Legend: one row, bottom right, no frame border
ax.legend(ncol=3, loc='lower right', frameon=False)

ax.set_xlabel('dt (ps)')
ax.set_ylabel(r'SDC ($\AA^2$/ps)')

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
ax_inset.plot(time_zoom, sdc_x_zoom, color='C0', alpha=0.35, lw=1.0)
ax_inset.plot(time_zoom, sdc_y_zoom, color='C1', alpha=0.35, lw=1.0)
ax_inset.plot(time_zoom, sdc_z_zoom, color='C2', alpha=0.35, lw=1.0)

# Overlay smoothed lines (solid, opaque)
ax_inset.plot(time_zoom, sdc_x_smooth, color='C0', lw=1.6, label='x smooth')
ax_inset.plot(time_zoom, sdc_y_smooth, color='C1', lw=1.6, label='y smooth')
ax_inset.plot(time_zoom, sdc_z_smooth, color='C2', lw=1.6, label='z smooth')

plt.tight_layout()

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
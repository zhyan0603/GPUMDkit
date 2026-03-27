"""
This script reads data from 'msd.out' and plots both MSD and SDC data 
side-by-side. The left plot shows MSD vs time with slope annotations. 
The right plot shows SDC vs time with an inset showing the last 80% 
of the data with a moving average overlay.

Usage:
    python plt_msd_sdc.py [save]
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# Function to read all required data from file
def read_data(file_name):
    data = np.loadtxt(file_name)
    # time, msd(x,y,z), sdc(x,y,z)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6]

# Function to calculate the slope of the first 50% of the data
def calculate_slope(x, y):
    len_data = len(x)
    start_idx = int(0.1 * len_data)  # 10% of the data
    end_idx = int(0.3 * len_data)    # 30% of the data
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
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5), dpi=150)

# ==========================================
# Left Subplot (ax1): MSD vs dt
# ==========================================
slope_x = calculate_slope(time, msd_x)
slope_y = calculate_slope(time, msd_y)
slope_z = calculate_slope(time, msd_z)

ax1.plot(time, msd_x, label='x', color='C0')
ax1.plot(time, msd_y, label='y', color='C1')
ax1.plot(time, msd_z, label='z', color='C2')

ax1.text(0.5, 0.3, f'Slope (x): {slope_x:.3f} Å$^2$/ps', fontsize=10, color='C0', verticalalignment='top', transform=ax1.transAxes)
ax1.text(0.5, 0.21, f'Slope (y): {slope_y:.3f} Å$^2$/ps', fontsize=10, color='C1', verticalalignment='top', transform=ax1.transAxes)
ax1.text(0.5, 0.12, f'Slope (z): {slope_z:.3f} Å$^2$/ps', fontsize=10, color='C2', verticalalignment='top', transform=ax1.transAxes)

ax1.set_xlabel('dt (ps)')
ax1.set_ylabel(r'MSD ($\AA^2$)')
ax1.legend()

# ==========================================
# Right Subplot (ax2): SDC vs dt
# ==========================================
ax2.plot(time, sdc_x, label='x', color='C0')
ax2.plot(time, sdc_y, label='y', color='C1')
ax2.plot(time, sdc_z, label='z', color='C2')

# Legend: one row, bottom right, no frame border
# ax2.legend(ncol=3, loc='lower right', frameon=False)
ax2.set_xlabel('dt (ps)')
ax2.set_ylabel(r'SDC ($\AA^2$/ps)')

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
ax_inset.plot(time_zoom, sdc_x_zoom, color='C0', alpha=0.35, lw=1.0)
ax_inset.plot(time_zoom, sdc_y_zoom, color='C1', alpha=0.35, lw=1.0)
ax_inset.plot(time_zoom, sdc_z_zoom, color='C2', alpha=0.35, lw=1.0)

# Overlay smoothed lines (solid, opaque)
ax_inset.plot(time_zoom, sdc_x_smooth, color='C0', lw=1.6, label='x smooth')
ax_inset.plot(time_zoom, sdc_y_smooth, color='C1', lw=1.6, label='y smooth')
ax_inset.plot(time_zoom, sdc_z_smooth, color='C2', lw=1.6, label='z smooth')

# ==========================================
# Final Layout and Save/Show logic
# ==========================================
plt.tight_layout()

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
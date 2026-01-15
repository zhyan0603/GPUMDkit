import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import ScalarFormatter
from collections import Counter
from scipy import interpolate

# Load data from 'force_train.out' file
data = np.loadtxt('force_train.out', skiprows=1)  # Assuming the first row is the header

# Extract NEP and DFT force components
nep_forces = data[:, 0:3]  # First three columns: NEP forces (Fx_NEP, Fy_NEP, Fz_NEP)
dft_forces = data[:, 3:6]  # Last three columns: DFT forces (Fx_DFT, Fy_DFT, Fz_DFT)

# Calculate force magnitudes
nep_f_mag = np.linalg.norm(nep_forces, axis=1)  # Magnitude of NEP forces
dft_f_mag = np.linalg.norm(dft_forces, axis=1)  # Magnitude of DFT forces

# Calculate the required values
dft_f = dft_f_mag           # DFT force magnitude
delta_f = nep_f_mag - dft_f_mag  # Difference in force magnitudes

# Calculate the angle between NEP and DFT forces (in degrees)
dot_product = np.sum(nep_forces * dft_forces, axis=1)  # Dot product
with np.errstate(divide='ignore', invalid='ignore'):  # Handle division by zero or invalid values
    cos_theta = dot_product / (nep_f_mag * dft_f_mag)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # Ensure cos_theta is within [-1, 1]
    delta_angle = np.arccos(cos_theta) * (180.0 / np.pi)  # Convert to degrees
delta_angle[np.isnan(delta_angle)] = 0.0  # Set NaN (zero vector cases) to 0

# Set up the figure and grid for plotting
fig = plt.figure(figsize=(12, 6), dpi=100)
gs = GridSpec(2, 3)

# Subplot 1: Scatter plot of delta_f vs dft_f
ax1 = plt.subplot(gs[0, 0])
ax1.scatter(dft_f, delta_f, s=10, c='g', marker='o')
ax1.set_xlabel(r'|$F_{DFT}$| (eV/$\mathrm{{\AA}}$)')
ax1.set_ylabel(r'$\delta_F$ (eV/$\mathrm{{\AA}}$)')
ax1.set_title('The errors of $\delta_F$')

# Subplot 2: Scatter plot of delta_angle vs dft_f
ax2 = plt.subplot(gs[1, 0])
ax2.scatter(dft_f, delta_angle, s=10, c='orange', marker='s')
ax2.set_xlabel(r'|$F_{DFT}$| (eV/$\mathrm{{\AA}}$)')
ax2.set_ylabel(r'$\delta_{\theta}$ ($\degree$)')
ax2.set_title(r'The errors of $\delta_{\theta}$')

# Subplot 3: Distribution of delta_f
ax3 = plt.subplot(gs[0, 1])
min_f = np.min(delta_f)
max_f = np.max(delta_f)
interval_width = (max_f - min_f) / 20
interval_counts = Counter()
for value in delta_f:
    interval = int((value - min_f) / interval_width)
    interval_counts[interval] += 1

delta_f_vals = []
counts_f = []
for interval, count in sorted(interval_counts.items()):
    interval_start = min_f + interval * interval_width
    interval_end = interval_start + interval_width
    delta_f_vals.append((interval_start + interval_end) / 2)
    counts_f.append(count)

ax3.plot(delta_f_vals, counts_f, marker='o', linestyle='-', c='g')
ax3.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax3.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax3.set_xlabel(r'$\delta_F$ (eV/$\mathrm{{\AA}}$)')
ax3.set_ylabel('Counts')
ax3.set_title('The distribution of $\delta_F$')
ax3.grid(False)

# Subplot 4: Distribution of delta_angle
ax4 = plt.subplot(gs[1, 1])
min_angle = np.min(delta_angle)
max_angle = np.max(delta_angle)
interval_width = (max_angle - min_angle) / 20
interval_counts = Counter()
for value in delta_angle:
    interval = int((value - min_angle) / interval_width)
    interval_counts[interval] += 1

delta_angle_vals = []
counts_angle = []
for interval, count in sorted(interval_counts.items()):
    interval_start = min_angle + interval * interval_width
    interval_end = interval_start + interval_width
    delta_angle_vals.append((interval_start + interval_end) / 2)
    counts_angle.append(count)

ax4.plot(delta_angle_vals, counts_angle, marker='s', linestyle='-', c='orange')
ax4.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax4.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax4.set_xlabel(r'$\delta_{\theta}$ ($\degree$)')
ax4.set_ylabel('Counts')
ax4.set_title(r'The distribution of $\delta_{\theta}$')
ax4.grid(False)

# Subplot 5: Cumulative Distribution Function (CDF) of delta_f
ax5 = plt.subplot(gs[0, 2])
sorted_df = np.sort(np.abs(delta_f))
min_f = sorted_df[0]
max_f = sorted_df[-1]
xs = np.arange(min_f, max_f, (max_f - min_f) / 20)
cdf = np.arange(len(sorted_df)) / float(len(sorted_df))
interp = interpolate.interp1d(sorted_df, cdf)
ys = interp(xs) * 100
ax5.plot(xs, ys, marker='o', linestyle='-', c='g')
ax5.set_xlabel(r'|$\delta_F$| (eV/$\mathrm{{\AA}}$)')
ax5.set_ylabel('Probability (%)')
ax5.set_title(r'The CDF of $\delta_F$')

# Subplot 6: Cumulative Distribution Function (CDF) of delta_angle
ax6 = plt.subplot(gs[1, 2])
sorted_da = np.sort(np.abs(delta_angle))
min_a = sorted_da[0]
max_a = sorted_da[-1]
xs_a = np.arange(min_a, max_a, (max_a - min_a) / 20)
cdf_a = np.arange(len(sorted_da)) / float(len(sorted_da))
interp_a = interpolate.interp1d(sorted_da, cdf_a)
ys_a = interp_a(xs_a) * 100
ax6.plot(xs_a, ys_a, marker='s', linestyle='-', c='orange')
ax6.set_xlabel(r'$\delta_{\theta}$ ($\degree$)')
ax6.set_ylabel('Probability (%)')
ax6.set_title(r'The CDF of $\delta_{\theta}$')

# Adjust layout and save the figure
plt.tight_layout(pad=1.0)

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('force_errors.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'force_errors.png'.")
        plt.savefig('force_errors.png', dpi=300)
    else:
        plt.show()
"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_vac.py
Category:   Plot Scripts
Purpose:    Plot velocity autocorrelation function (VAC) from SDC output data.
Usage:      gpumdkit.sh -plt vac [save]
            python plt_vac.py [save]
Arguments:
  save      Save the plot as PNG instead of displaying it
Output:
  vac.png  (if save is used, or if backend is non-interactive)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
})

# Function to read data from the given file
def read_data(file_name):
    data = np.loadtxt(file_name)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# Determine the input file
#input_file = sys.argv[1] if len(sys.argv) > 1 else './sdc.out'
input_file = './sdc.out'

# Read the data from the file
time, vac_x, vac_y, vac_z = read_data(input_file)

# Calculate the total SDC
vac_mean = (vac_x + vac_y + vac_z)/3

# Create a plot for the MSD data
plt.figure(figsize=(6, 4.5), dpi=100)
plt.plot(time, vac_x, label='x')
plt.plot(time, vac_y, label='y')
plt.plot(time, vac_z, label='z')
plt.plot(time, vac_mean, label='mean', color='C4')

# Add titles and labels
#plt.title('VAC vs dt')
plt.xlabel('dt (ps)')
plt.ylabel(r'VAC ($\AA^2/ps^2$)')
plt.legend()
plt.tight_layout()
#plt.grid(True)

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('vac.png')
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'vac.png'.")
        plt.savefig('vac.png', dpi=300)
    else:
        plt.show()
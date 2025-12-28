"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Plot NEP restart information

Usage:
    python plt_nep_restart.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Function to read two-column data from a file
def read_data(file_path):
    data = np.loadtxt(file_path)
    return data[:, 0], data[:, 1]

# Read data from files
col1, col2 = read_data('./nep.restart')

# Create subplots
fig, axes = plt.subplots(2, 1, figsize=(6, 4.2), sharex=True)

# First subplot
axes[0].plot(range(len(col1)), col1, label='Column 1', color='C0')
axes[0].hlines(0, 0, len(col1), color='grey', linestyle='--')
axes[0].set_ylabel('Column 1')
axes[0].legend()
axes[0].set_title('Parameters in nep.restart')
#axes[0].grid(True)

# Second subplot
axes[1].plot(range(len(col2)), col2, label='Column 2', color='C1')
axes[1].set_ylabel('Column 2')
axes[1].set_xlabel('Index')
axes[1].legend()
#axes[1].grid(True)

# Adjust layout
#plt.tight_layout()
fig.subplots_adjust(top=0.918,bottom=0.142,left=0.139,right=0.975,hspace=0.0,wspace=0.24)

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

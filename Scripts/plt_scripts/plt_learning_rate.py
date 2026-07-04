"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_learning_rate.py
Category:   Plot Scripts
Purpose:    Visualize learning rate evolution during GNEP training process.
Usage:      gpumdkit.sh -plt lr
            python plt_learning_rate.py [save]
Arguments:
  save      Save the plot as 'learning_rate.png' instead of displaying it
Output:
  learning_rate.png  (if save is used)
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

data = np.loadtxt('loss.out')

epochs = np.arange(1, len(data) + 1) 
learning_rates = data[:, -2]  

print(f"Ploting the learning rate during GNEP training process...")

plt.figure(figsize=(5, 3.5))
plt.plot(epochs, learning_rates, c='C0', label='Learning Rate')
plt.xlabel('Epoch')
plt.ylabel('Learning Rate')
plt.legend()
plt.tight_layout()

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('learning_rate.png', dpi=300)
else:
    # Handle saving or displaying the plot
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'learning_rate.png'.")
        plt.savefig('learning_rate.png', dpi=300)
    else:
        plt.show()
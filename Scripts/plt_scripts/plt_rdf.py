"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_rdf.py
Category:   Plot Scripts
Purpose:    Plot radial distribution function (RDF) from rdf.out showing
            pair correlations for all atom pairs.
Usage:      gpumdkit.sh -plt rdf [save]
            python plt_rdf.py [save]
Arguments:
  save     Save the plot as 'rdf.png' instead of displaying it
Output:
  rdf.png  RDF plot (saved or displayed)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-21
=============================================================================
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import math

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
COLORS = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Plot RDF data from rdf.out file.')
parser.add_argument(
    'action',
    nargs='?',
    choices=['save'],
    default=None,
    help='Use "save" to write rdf.png; otherwise, show the plot.',
)
args = parser.parse_args()

# Read the rdf.out file
with open('rdf.out', 'r') as f:
    lines = f.readlines()

# Extract header from the first line (strip # and split)
header = lines[0].strip().lstrip('#').split()

# Read data lines, skipping empty or comment lines
data_lines = [line.split() for line in lines[1:] if line.strip() and not line.startswith('#')]

# Create DataFrame
df = pd.DataFrame(data_lines, columns=header).astype(float)

# Extract radius (first column)
radius = df[header[0]]

# RDF columns are from the second column onward
rdf_columns = header[1:]

# Number of subplots equals number of RDF columns
num_subplots = len(rdf_columns)

# Calculate nrows and ncols: max 3 columns
ncols = min(3, num_subplots)
nrows = math.ceil(num_subplots / ncols)

# Create subplots
fig_width, fig_height = 3.4 * ncols, 2.2 * nrows
fig, axs = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))

# Flatten axs for easy iteration
axs = axs.flat if num_subplots > 1 else [axs]

# Plot each RDF in its own subplot
for i, col in enumerate(rdf_columns):
    axs[i].plot(radius, df[col], color=COLORS[i % len(COLORS)], linewidth=1.4)
    axs[i].set_title(col, fontsize=9, pad=3)
    axs[i].set_xlabel(r'$r$ ($\mathrm{\AA}$)')
    axs[i].set_ylabel(r'$g(r)$')
    axs[i].set_xlim(radius.iloc[0], radius.iloc[-1])

# Hide unused subplots if any
for i in range(num_subplots, nrows * ncols):
    axs[i].axis('off')

# Subplot geometry (margins are given in inches and converted to figure fractions);
# tune these numbers to adjust margins and spacing
fig.subplots_adjust(left=0.52 / fig_width, right=1 - 0.15 / fig_width,
                    top=1 - 0.25 / fig_height, bottom=0.48 / fig_height,
                    wspace=0.2, hspace=0.5)

# Save or show based on argument
if args.action == 'save':
    plt.savefig('rdf.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'rdf.png'.")
        plt.savefig('rdf.png', dpi=300)
    else:
        plt.show()

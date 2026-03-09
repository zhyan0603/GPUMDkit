import pandas as pd
import matplotlib.pyplot as plt
import argparse
import math

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Plot RDF data from rdf.out file.')
parser.add_argument('action', nargs='?', default=None, help='If "save", save the plot as PNG; otherwise, show the plot.')
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
fig, axs = plt.subplots(nrows, ncols, figsize=(3.3 * ncols, 2.3 * nrows)) 

# Flatten axs for easy iteration
axs = axs.flat if num_subplots > 1 else [axs]

# Plot each RDF in its own subplot
for i, col in enumerate(rdf_columns):
    axs[i].plot(radius, df[col], label=col)
    axs[i].legend(loc='upper left')
    axs[i].set_xlabel('Radius')
    axs[i].set_ylabel('g(r)')

# Hide unused subplots if any
for i in range(num_subplots, nrows * ncols):
    axs[i].axis('off')

# Adjust layout
plt.tight_layout()

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
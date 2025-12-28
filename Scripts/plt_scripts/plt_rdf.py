"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Plot radial distribution function

Usage:
    python plt_rdf.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as plt

# Function to plot RDF data
def plot_rdf(file_path, column_index=None, save=False):
    # Load data from the file
    data = np.loadtxt(file_path)
    
    # Extract distances and RDF values
    distances = data[:, 0]
    rdf_values = data[:, 1:]

    # Generate column names (assumes columns after the first one are types)
    column_names = ['distance'] + [f'column{int(i)}' for i in range(1, rdf_values.shape[1] + 1)]
    
    # If a column index is provided, plot only that column
    if column_index is not None:
        column_index -= 1  # Adjust to zero-indexing
        if column_index < 0 or column_index >= len(column_names) - 1:
            print(f"Invalid column index. Please choose a value between 1 and {len(column_names) - 1}")
            return
        
        plt.figure(figsize=(6, 3.5))
        # Plot the selected column
        plt.plot(distances, rdf_values[:, column_index], label=column_names[column_index + 1])
        plt.xlabel(r'Distance ($\AA$)')
        plt.ylabel('g(r)')
        plt.tight_layout()
        plt.legend()
        
        if save:
            plt.savefig(f'rdf.png', dpi=300)  # Save the plot as PNG
        else:
            plt.show()
    else:
        # If no column index is provided, plot all columns
        num_plots = len(column_names) - 1  # Exclude 'total'
        
        # Calculate optimal grid size (rows, columns)
        ncols = int(math.ceil(math.sqrt(num_plots)))  # Determine number of columns
        nrows = int(math.ceil(num_plots / ncols))     # Determine number of rows
        
        # Create subplots with optimal grid layout
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 3.5, nrows * 2.3))
        
        # Flatten axes if there's more than 1 row/column
        axes = axes.flatten() if nrows > 1 else axes
        
        # Plot each column in a subplot
        for i, ax in enumerate(axes):
            if i < num_plots:  # Only plot for the available number of columns
                ax.plot(distances, rdf_values[:, i], label=column_names[i + 1])
                ax.set_xlabel(r'Distance ($\AA$)')
                ax.set_ylabel('g(r)')
                ax.legend()
            else:
                ax.axis('off')  # Turn off unused axes if the grid is larger than the number of columns

        # Adjust layout
        plt.tight_layout()
        
        if save:
            plt.savefig(f'rdf.png', dpi=300)  # Save the plot as PNG
        else:
            plt.show()

# Check if a column index is provided via command line arguments
if len(sys.argv) > 1:
    # Check if the first argument is 'save'
    if sys.argv[1] == 'save':
        column_index = None  # If 'save', we want to plot all columns
        save = True
    else:
        # Otherwise, assume the first argument is the column index
        column_index = int(sys.argv[1])
else:
    column_index = None

# Check if 'save' is in sys.argv, this handles the case where 'save' could be in any position
save = 'save' in sys.argv

file_path = 'rdf.out'

# Call the plot function
plot_rdf(file_path, column_index, save)

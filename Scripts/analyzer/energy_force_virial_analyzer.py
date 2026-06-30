"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     energy_force_virial_analyzer.py
Category:   Analyzer Scripts
Purpose:    Compute the range (min, max) of energy, force, or virial in an
            extxyz file. Optionally plot a histogram.
Usage:      gpumdkit.sh -range <filename> <property> [hist]
            python energy_force_virial_analyzer.py <filename> <property> [hist]
Arguments:
  filename   Input extxyz file
  property   Property to analyze: energy, force, or virial
  hist       (optional) If "hist", plot a histogram
Output:
  Printed min/max values. PNG histogram if requested.
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 2 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -range <exyzfile> <property> [hist]")
    print("    or: python energy_force_virial_analyzer.py <filename> <property> [hist]")
    print("")
    print(" Arguments:")
    print("   exyzfile    Input extxyz trajectory file")
    print("   property    Property to analyze: energy, force, or virial")
    print("   hist        (optional) If 'hist', plot a histogram")
    print("")
    print(" Example: gpumdkit.sh -range train.xyz force hist")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

import numpy as np
from ase.io import read
import matplotlib.pyplot as plt

def calculate_range(frames, property_name):
    property_name = property_name.lower()
    values = []
    
    for frame in frames:
        info_lower = {k.lower(): v for k, v in frame.info.items()}
        
        if property_name == "energy":
            if 'energy' in info_lower:
                values.append(info_lower['energy'])
            else:
                raise ValueError("Energy information not found in frame info.")
        elif property_name in ["force", "forces"]:
            forces = frame.get_forces()
            values.extend(np.linalg.norm(forces, axis=1))
        elif property_name == "virial":
            if 'virial' in info_lower:
                virial = info_lower['virial']
                values.extend(virial)  
            else:
                raise ValueError("Virial information not found in frame info.")
        else:
            raise ValueError("Invalid property. Choose from 'energy', 'force', or 'virial'.")
    
    return np.min(values), np.max(values), values

def plot_histogram(values, property_name):
    plt.figure(figsize=(6,4), dpi=100)
    plt.hist(values, bins=30, edgecolor='black')
    plt.title(f'{property_name.capitalize()} Histogram')
    plt.xlabel(f'{property_name.capitalize()}')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plt.show()
    #plt.savefig(f'range_{property_name.capitalize()}.png')

if __name__ == "__main__":
    filename = sys.argv[1]
    property_name = sys.argv[2]
    plot_hist = len(sys.argv) > 3 and sys.argv[3] == 'hist'
    
    # Read the extxyz file
    frames = read(filename, index=":")
    
    # Calculate the range of the specified property
    min_val, max_val, values = calculate_range(frames, property_name)
    
    # Print the range
    print(f"{property_name.capitalize()} range: {min_val:.3f} to {max_val:.3f}")
    
    # Plot histogram if requested
    if plot_hist:
        plot_histogram(values, property_name)

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
Last-modified: 2026-10-14
=============================================================================
"""

import os
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

def get_property(frame, name):
    """Return a frame property from frame.info or its calculator, ignoring key case.

    ASE 3.23+ moves extxyz keys such as energy= into frame.calc, while
    capitalized keys such as Energy= stay in frame.info.
    """
    sources = [frame.info]
    if frame.calc is not None:
        sources.append(frame.calc.results)
    for source in sources:
        for key, value in source.items():
            if key.lower() == name:
                return value
    return None

def calculate_range(frames, property_name):
    property_name = property_name.lower()
    if property_name == "forces":
        property_name = "force"
    if property_name not in ("energy", "force", "virial"):
        print(" Error: invalid property. Choose from 'energy', 'force', or 'virial'.")
        sys.exit(1)
    values = []
    missing = 0

    for frame in frames:
        value = get_property(frame, "forces" if property_name == "force" else property_name)
        if value is None:
            missing += 1
            continue
        value = np.asarray(value, dtype=float)
        if property_name == "force":
            values.extend(np.linalg.norm(value, axis=1))
        else:
            values.extend(value.ravel())

    if missing == len(frames):
        print(f" Error: {property_name} information not found in any frame.")
        sys.exit(1)
    if missing:
        print(f" Note: {missing} frames have no {property_name} and were skipped.")
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
    if not os.path.isfile(filename):
        print(f" Error: file '{filename}' does not exist.")
        sys.exit(1)
    plot_hist = len(sys.argv) > 3 and sys.argv[3] == 'hist'
    
    # Read the extxyz file
    frames = read(filename, index=":")
    
    # Calculate the range of the specified property
    min_val, max_val, values = calculate_range(frames, property_name)
    
    # Print the range
    print(f" {property_name.capitalize()} range: {min_val:.3f} to {max_val:.3f}")
    
    # Plot histogram if requested
    if plot_hist:
        plot_histogram(values, property_name)

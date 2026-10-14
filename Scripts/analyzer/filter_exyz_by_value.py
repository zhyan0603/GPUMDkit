"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     filter_exyz_by_value.py
Category:   Analyzer Scripts
Purpose:    Filter structures based on energy, force, or virial thresholds.
Usage:      gpumdkit.sh -filter_value <input_file> <property> <threshold>
            python filter_exyz_by_value.py <input_file> <property> <threshold>
Arguments:
  input_file  Input extxyz file
  property    Filtering property: energy, force, or virial
  threshold   Maximum allowed value for the specified property
Output:
  filtered.xyz
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-14
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 3 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -filter_value <input_file> <property> <threshold>")
    print("    or: python filter_exyz_by_value.py <input_file> <property> <threshold>")
    print("")
    print(" Arguments:")
    print("   input_file   Input extxyz trajectory file")
    print("   property     Filtering property: energy, force, or virial")
    print("   threshold    Maximum allowed value for the specified property")
    print("")
    print(" Output:")
    print("   filtered.xyz  Structures that pass the filter")
    print("")
    print(" Example: gpumdkit.sh -filter_value train.xyz energy -3.5")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

import numpy as np
from ase.io import read, write


def get_property(atoms, name):
    """Return a frame property from atoms.info or its calculator, ignoring key case.

    ASE 3.23+ moves extxyz keys such as energy= into atoms.calc, while
    capitalized keys such as Energy= stay in atoms.info.
    """
    sources = [atoms.info]
    if atoms.calc is not None:
        sources.append(atoms.calc.results)
    for source in sources:
        for key, value in source.items():
            if key.lower() == name:
                return value
    return None


def main():
    # Command line arguments
    filename = sys.argv[1]  # Input .extxyz file
    property = sys.argv[2].lower()  # Filtering property: 'energy', 'force', or 'virial'
    if property == 'forces':
        property = 'force'
    if property not in ('energy', 'force', 'virial'):
        print(" Error: unsupported property. Please use 'energy', 'force', or 'virial'.")
        sys.exit(1)
    try:
        threshold = float(sys.argv[3])  # Threshold value for filtering
    except ValueError:
        print(f" Error: threshold must be a number, got '{sys.argv[3]}'.")
        sys.exit(1)

    # Read the structures from the .extxyz file
    images = read(filename, index=':')

    # Keep a structure only when every value is at most the threshold
    # Frames without the property cannot be checked, so they are kept.
    filtered_images = []
    missing = 0
    for atoms in images:
        values = get_property(atoms, 'forces' if property == 'force' else property)
        if values is None:
            missing += 1
            filtered_images.append(atoms)
            continue
        values = np.asarray(values, dtype=float)
        if property == 'force':
            # Each force component is compared by magnitude
            values = np.abs(values)
        if np.all(values <= threshold):
            filtered_images.append(atoms)

    if missing == len(images):
        print(f" Error: no frame contains {property}; filtered.xyz was not written.")
        sys.exit(1)
    if missing:
        print(f" Note: {missing} frames have no {property} and were kept unchecked.")

    # Write the filtered structures back to a new .extxyz file
    write("filtered.xyz", filtered_images)
    print(f" {len(filtered_images)} of {len(images)} structures kept; written to filtered.xyz")


if __name__ == "__main__":
    main()

"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     xdatcar2exyz.py
Category:   Format Conversion Scripts
Purpose:    Convert VASP XDATCAR file to extended XYZ format using ASE.
Usage:      gpumdkit.sh -xdat2exyz <XDATCAR> <output.xyz>
            python xdatcar2exyz.py <XDATCAR> <output.xyz>
Arguments:
  input       Path to the input XDATCAR file
  output      Path to the output .xyz file
Output:
  <output.xyz>  (converted trajectory in extxyz format)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 2 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -xdat2exyz <XDATCAR> <output.xyz>")
    print("    or: python xdatcar2exyz.py <XDATCAR> <output.xyz>")
    print("")
    print(" Arguments:")
    print("   XDATCAR     Path to the input XDATCAR file")
    print("   output.xyz  Path to the output extxyz file")
    print("")
    print(" Example: gpumdkit.sh -xdat2exyz XDATCAR output.xyz")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

from ase.io import read, write

input_file = sys.argv[1]
output_file = sys.argv[2]

# Execute conversion logic
try:
    print(f" Reading XDATCAR: {input_file} ...")
    
    # Read all frames from the trajectory (index=":")
    # ASE automatically parses scaling, lattice, and species
    trajectory = read(input_file, index=":", format="vasp-xdatcar")
    
    num_frames = len(trajectory)
    print(f" Detected {num_frames} frames.")
    
    print(f" Writing to extxyz: {output_file} ...")
    
    # Writing in extxyz format preserves the Lattice matrix in the header
    write(output_file, trajectory, format="extxyz")
    
    print(" Conversion complete.")

except FileNotFoundError:
    print(f" Error: Could not find '{input_file}'. Please check the path.")
    sys.exit(1)
except Exception as e:
    print(f" An unexpected error occurred: {e}")
    sys.exit(1)
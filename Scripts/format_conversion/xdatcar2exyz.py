"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     xdatcar2exyz.py
Category:   Format Conversion Scripts
Purpose:    Convert VASP XDATCAR file to extended XYZ format using ASE.
Usage:      python xdatcar2exyz.py <XDATCAR> <output.xyz>
Arguments:
  input       Path to the input XDATCAR file
  output      Path to the output .xyz file
Output:
  <output.xyz>  (converted trajectory in extxyz format)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import argparse
import sys
from ase.io import read, write

# 1. Setup Argument Parser for positional arguments
parser = argparse.ArgumentParser(
    description=" Convert VASP XDATCAR to extxyz format."
)
parser.add_argument("input", help=" Path to the input XDATCAR file")
parser.add_argument("output", help=" Path to the output .xyz file")

args = parser.parse_args()

# 2. Execute conversion logic
try:
    print(f" Reading XDATCAR: {args.input} ...")
    
    # Read all frames from the trajectory (index=":")
    # ASE automatically parses scaling, lattice, and species
    trajectory = read(args.input, index=":", format="vasp-xdatcar")
    
    num_frames = len(trajectory)
    print(f" Detected {num_frames} frames.")
    
    print(f" Writing to extxyz: {args.output} ...")
    
    # Writing in extxyz format preserves the Lattice matrix in the header
    write(args.output, trajectory, format="extxyz")
    
    print(" Conversion complete.")

except FileNotFoundError:
    print(f" Error: Could not find '{args.input}'. Please check the path.")
    sys.exit(1)
except Exception as e:
    print(f" An unexpected error occurred: {e}")
    sys.exit(1)
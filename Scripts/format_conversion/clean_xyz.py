"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     clean_xyz.py
Category:   Format Conversion Scripts
Purpose:    Remove stress, virial, and force information from an extxyz
            training file, keeping only structural information.
Usage:      gpumdkit.sh -clean_xyz <input.xyz> <output.xyz>
            python clean_xyz.py <input.xyz> <output.xyz>
Arguments:
  input.xyz    Input extxyz training file
  output.xyz   Output cleaned extxyz file
Output:
  Cleaned extxyz file with only structural information
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 2 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -clean_xyz <input.xyz> <output.xyz>")
    print("    or: python clean_xyz.py <input.xyz> <output.xyz>")
    print("")
    print(" Arguments:")
    print("   input.xyz    Input extxyz training file")
    print("   output.xyz   Output cleaned extxyz file")
    print("")
    print(" Output:")
    print("   Cleaned extxyz file with only structural information")
    print("   (stress, virial, and force data are removed)")
    print("")
    print(" Example: gpumdkit.sh -clean_xyz train.xyz cleaned.xyz")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

from ase.io import read, write

# Read all frames from the input train.xyz file
frames = read(sys.argv[1], index=':')

# Iterate over each frame to remove unnecessary information
for frame in frames:
    # Remove stress and virial information from frame.info
    for key in list(frame.info.keys()):  # Use list() to avoid errors when modifying the dictionary
        if 'stress' in key.lower() or 'virial' in key.lower():
            del frame.info[key]
    
    # Remove forces from frame.arrays
    if 'forces' in frame.arrays:
        del frame.arrays['forces']
    
    # Remove calculation results (including stress and forces) from frame.calc
    if hasattr(frame, 'calc') and frame.calc is not None:
        if hasattr(frame.calc, 'results'):
            # Remove stress and forces from calc.results
            for key in list(frame.calc.results.keys()):
                if 'stress' in key.lower() or 'forces' in key.lower():
                    del frame.calc.results[key]
    
    # Ensure that the calculator for the frame is cleared
    frame.calc = None

# Export the cleaned frames to a new extxyz file, containing only lattice and atomic information
write(sys.argv[2], frames, format='extxyz')
print(f" Cleaned {len(frames)} frames, saved to {sys.argv[2]}")

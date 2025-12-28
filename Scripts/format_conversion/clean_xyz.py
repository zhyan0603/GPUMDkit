"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Clean and normalize XYZ files

Usage:
    python clean_xyz.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

from ase.io import read, write
import sys

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

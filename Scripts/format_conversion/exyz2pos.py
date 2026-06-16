"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     exyz2pos.py
Category:   Format Conversion Scripts
Purpose:    Convert extended XYZ file(s) to VASP POSCAR format(s).
Usage:      python exyz2pos.py [extxyz_file]
Arguments:
  extxyz_file Input extxyz file (default: train.xyz)
Output:
  POSCAR file(s) in VASP format
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
from ase.io import read, write

def print_progress_bar(iteration, total, length=50):
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = '█' * filled_length + '-' * (length - filled_length)
    print(f'\r Progress: |{bar}| {percent}% Complete', end='\r')
    # Print New Line on Complete
    if iteration == total:
        print()

# Check arguments
if len(sys.argv) < 2:
    print("Usage: python exyz2pos.py <input.xyz>")
    print("   or: gpumdkit.sh -exyz2pos <input.xyz>")
    sys.exit(1)

input_file = sys.argv[1]

# Read all frames
frames = read(input_file, index=':')
total_frames= len(frames)

# Save to POSCAR
for i, frame in enumerate(frames):
    poscar_filename = f'POSCAR_{i + 1}.vasp'
    write(poscar_filename, frame)
    print_progress_bar(i + 1, total_frames)

print(f' All frames have been converted.')
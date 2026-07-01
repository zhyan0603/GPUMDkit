"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     frame_range.py
Category:   Sample Structure Scripts
Purpose:    Extract a range of frames from an extxyz trajectory by start
            and end fractions (e.g., first 80% or last 50%).
Usage:      gpumdkit.sh -frame_range <input.xyz> <start_fraction> <end_fraction>
            python frame_range.py <input.xyz> <start_fraction> <end_fraction>
Arguments:
  input.xyz       Input extxyz trajectory file
  start_fraction  Start fraction (0.0 to 1.0, e.g., 0 for beginning)
  end_fraction    End fraction (0.0 to 1.0, e.g., 0.8 for first 80%)
Output:
  <input>_<start>_<end>.xyz
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 3 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -frame_range <exyzfile> <start_fraction> <end_fraction>")
    print("    or: python frame_range.py <input.xyz> <start_fraction> <end_fraction>")
    print("")
    print(" Arguments:")
    print("   exyzfile         Input extxyz trajectory file")
    print("   start_fraction   Start fraction (0.0 to 1.0, e.g., 0 for beginning)")
    print("   end_fraction     End fraction (0.0 to 1.0, e.g., 0.8 for first 80%)")
    print("")
    print(" Output:")
    print("   <input>_<start>_<end>.xyz  Extracted frame range")
    print("")
    print(" Example: gpumdkit.sh -frame_range dump.xyz 0.2 0.5")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

from ase.io import read, write

input_file = sys.argv[1]
start_fraction = float(sys.argv[2])
end_fraction = float(sys.argv[3])

# Read the entire trajectory
images = read(input_file, index=':')

# Calculate indices
n_frames = len(images)
start_idx = int(start_fraction * n_frames)
end_idx = int(end_fraction * n_frames)

# Extract subset
subset_images = images[start_idx:end_idx]

# Generate output filename
output_file = input_file.replace('.xyz', f'_{start_fraction:.2f}_{end_fraction:.2f}.xyz')

# Write to new file
write(output_file, subset_images, format='extxyz')

print(f" Extracted {len(subset_images)} frames from fraction {start_fraction} to {end_fraction}")
print(f" Corresponding to frames {start_idx} to {end_idx-1}")
print(f" Saved to {output_file}")
"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    frame_range.py
"""

import sys
from ase.io import read, write

# Check command line arguments
if len(sys.argv) != 4:
    print(" Usage: python frame_range.py input.xyz start_fraction end_fraction")
    print(" Example: python frame_range.py dump.xyz 0 0.8")
    print(" Example: python frame_range.py dump.xyz 0.5 1")
    sys.exit(1)

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
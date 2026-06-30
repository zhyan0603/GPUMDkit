"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     get_frame.py
Category:   Format Conversion Scripts
Purpose:    Extract a single frame from an extxyz trajectory by frame number.
Usage:      gpumdkit.sh -get_frame <input.xyz> <frame_number>
            python get_frame.py <input.xyz> <frame_number>
Arguments:
  input.xyz      Input extxyz trajectory file
  frame_number   Frame number to extract (1-indexed)
Output:
  frame_<N>.xyz  (extracted single frame)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 2 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -get_frame <input.xyz> <frame_number>")
    print("    or: python get_frame.py <input.xyz> <frame_number>")
    print("")
    print(" Arguments:")
    print("   input.xyz      Input extxyz trajectory file")
    print("   frame_number   Frame number to extract (1-indexed)")
    print("")
    print(" Output:")
    print("   frame_<N>.xyz  Extracted single frame")
    print("")
    print(" Example: gpumdkit.sh -get_frame train.xyz 10")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

from ase.io import read, write

def get_frame():
    input_file = sys.argv[1]
    frame_number = int(sys.argv[2])
    frames = read(input_file, index=':')
    
    # Check if the frame number is valid
    if frame_number < 1 or frame_number > len(frames):
        print(f"Error: Frame {frame_number} is out of range. The file contains {len(frames)} frames.")
        sys.exit(1)
    selected_frame = frames[frame_number - 1]
    output_file = f"frame_{frame_number}.xyz"
    write(output_file, selected_frame)
    
    print(f"Frame {frame_number} has been successfully written to {output_file}.")

if __name__ == "__main__":
    get_frame()

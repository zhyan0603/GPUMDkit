"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     filter_exyz_by_box.py
Category:   Analyzer Scripts
Purpose:    Filter structures by a maximum box-edge length. Structures with
            any box edge exceeding the limit are discarded.
Usage:      gpumdkit.sh -filter_box <input_file> <edge_limit>
            python filter_exyz_by_box.py <input_file> <edge_limit>
Arguments:
  input_file   Input extxyz file
  edge_limit   Maximum allowed box-edge length (Angstrom)
Output:
  filtered_by_box.xyz
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 2 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -filter_box <exyzfile> <edge_limit>")
    print("    or: python filter_exyz_by_box.py <input_file> <edge_limit>")
    print("")
    print(" Arguments:")
    print("   exyzfile    Input extxyz trajectory file")
    print("   edge_limit  Maximum allowed box-edge length (Angstrom)")
    print("")
    print(" Output:")
    print("   filtered_by_box.xyz  Structures that pass the filter")
    print("")
    print(" Example: gpumdkit.sh -filter_box train.xyz 20.0")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

from ase.io import read, write
import numpy as np

def filter_frames(input_file, output_file, edge_limit):
    # Read all frames from the input file
    frames = read(input_file, index=":")

    filtered_frames = []
    discarded_count = 0

    for frame in frames:
        box = frame.get_cell()  # Get the box matrix (3x3)
        edges = np.linalg.norm(box, axis=1)  # Compute lengths of box edges
        
        if np.all(edges <= edge_limit):
            filtered_frames.append(frame)
        else:
            discarded_count += 1

    # Write the filtered frames to the output file
    if filtered_frames:
        write(output_file, filtered_frames)
    
    # Print the number of filtered-out structures
    print(f"Number of structures discarded: {discarded_count}")

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = "filtered_by_box.xyz"
    edge_limit = float(sys.argv[2])

    filter_frames(input_file, output_file, edge_limit)

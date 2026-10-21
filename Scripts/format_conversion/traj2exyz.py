"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     traj2exyz.py
Category:   Format Conversion Scripts
Purpose:    Convert an ASE .traj trajectory file to extended XYZ format.
Usage:      gpumdkit.sh -traj2exyz <input.traj> <output.xyz>
            python traj2exyz.py <input.traj> <output.xyz>
Arguments:
  input.traj   Input ASE trajectory file
  output.xyz   Output extxyz file
Output:
  <output.xyz>  (converted trajectory in extxyz format)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 2 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -traj2exyz <input.traj> <output.xyz>")
    print("    or: python traj2exyz.py <input.traj> <output.xyz>")
    print("")
    print(" Arguments:")
    print("   input.traj   Input ASE trajectory file")
    print("   output.xyz   Output extxyz file")
    print("")
    print(" Example: gpumdkit.sh -traj2exyz md.traj output.xyz")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

import os
from ase.io import read, write

def convert_traj_to_extxyz(input_file, output_file):
    """
    Converts an ASE .traj file to .extxyz format.
    """
    if not os.path.exists(input_file):
        print(f" Error: The file '{input_file}' does not exist.")
        return

    try:
        print(f" Reading frames from: {input_file} ...")
        # index=':' ensures all images in the trajectory are loaded
        images = read(input_file, index=':')
        
        print(f" Writing {len(images)} frames to: {output_file} ...")
        write(output_file, images, format='extxyz')
        
        print(" Conversion completed successfully.")
        
    except Exception as e:
        print(f" An unexpected error occurred: {e}")

if __name__ == "__main__":
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    convert_traj_to_extxyz(input_path, output_path)
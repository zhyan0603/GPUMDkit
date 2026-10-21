"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     pos2exyz.py
Category:   Format Conversion Scripts
Purpose:    Convert one or more VASP POSCAR files to extended XYZ format.
Usage:      gpumdkit.sh -pos2exyz <POSCAR> <output.xyz>
            python pos2exyz.py <POSCAR> <output.xyz>
Arguments:
  POSCAR       One or more VASP POSCAR files
  output.xyz   Output extxyz file
Output:
  <output.xyz>  (converted structures in extxyz format)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 2 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -pos2exyz <POSCAR> <output.xyz>")
    print("    or: python pos2exyz.py <POSCAR> <output.xyz>")
    print("")
    print(" Arguments:")
    print("   POSCAR       One or more VASP POSCAR files (supports wildcards)")
    print("   output.xyz   Output extxyz file")
    print("")
    print(" Example: gpumdkit.sh -pos2exyz POSCAR output.xyz")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

import glob
from ase.io import read, write

def convert_poscar_to_extxyz(poscar_filenames, extxyz_filename):

    # Create an empty list to store the frames
    all_frames = []

    # Read the structure from each POSCAR file and append to the list
    for poscar_filename in poscar_filenames:
        frames = read(poscar_filename, format='vasp')
        
        # If 'frames' is a list (multiple frames), we should iterate over it
        if isinstance(frames, list):
            all_frames.extend(frames)  # Add all frames from this file
        else:
            all_frames.append(frames)  # Add the single frame

    # Write all frames to extxyz file (using append=False ensures overwrite)
    write(extxyz_filename, all_frames, format='extxyz', append=False)

if __name__ == '__main__':
    poscar_filename_pattern = sys.argv[1]
    extxyz_filename = sys.argv[2]

    # Use glob to handle file patterns (with wildcards)
    poscar_filenames = glob.glob(poscar_filename_pattern)

    # If no files matched the pattern, print a message and exit
    if not poscar_filenames:
        print(f"No files matched the pattern: {poscar_filename_pattern}")
        sys.exit(1)

    # Call function to convert POSCAR to extxyz
    convert_poscar_to_extxyz(poscar_filenames, extxyz_filename)

    print(f"Conversion completed, structures saved to {extxyz_filename}")

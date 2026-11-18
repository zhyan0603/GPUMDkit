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
Usage:      gpumdkit.sh -pos2exyz <POSCAR> [<POSCAR2> ...] <output.xyz>
            python pos2exyz.py <POSCAR> [<POSCAR2> ...] <output.xyz>
Arguments:
  POSCAR       One or more VASP POSCAR files (wildcards allowed; files are
               read in natural order, e.g. POSCAR_2 before POSCAR_10)
  output.xyz   Output extxyz file (always the last argument)
Output:
  <output.xyz>  (converted structures in extxyz format)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-11-10
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 2 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -pos2exyz <POSCAR> [<POSCAR2> ...] <output.xyz>")
    print("    or: python pos2exyz.py <POSCAR> [<POSCAR2> ...] <output.xyz>")
    print("")
    print(" Arguments:")
    print("   POSCAR       One or more VASP POSCAR files (supports wildcards)")
    print("   output.xyz   Output extxyz file (always the last argument)")
    print("")
    print(" Example: gpumdkit.sh -pos2exyz POSCAR output.xyz")
    print("          gpumdkit.sh -pos2exyz POSCAR_*.vasp output.xyz")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

import glob
import os
import re
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

def natural_key(name):
    """Sort key that orders POSCAR_2 before POSCAR_10."""
    return [int(part) if part.isdigit() else part for part in re.split(r'(\d+)', name)]

if __name__ == '__main__':
    # The last argument is the output; everything before it is an input file
    # or a quoted wildcard pattern. This also covers an unquoted wildcard that
    # the shell has already expanded into many file names.
    poscar_patterns = args[:-1]
    extxyz_filename = args[-1]

    if not extxyz_filename.lower().endswith(('.xyz', '.extxyz')):
        print(f" Error: the last argument '{extxyz_filename}' must be the output file (.xyz or .extxyz).")
        sys.exit(1)

    poscar_filenames = []
    for pattern in poscar_patterns:
        matched = [pattern] if os.path.isfile(pattern) else glob.glob(pattern)
        if not matched:
            print(f" Error: no files matched the pattern: {pattern}")
            sys.exit(1)
        poscar_filenames.extend(matched)
    # Sort after collecting so both quoted and shell-expanded wildcards use natural order.
    poscar_filenames.sort(key=natural_key)

    if os.path.exists(extxyz_filename) and any(
        os.path.samefile(extxyz_filename, name) for name in poscar_filenames
    ):
        print(f" Error: the output file '{extxyz_filename}' is also an input file.")
        print(" Give the output extxyz file as the last argument.")
        sys.exit(1)

    # Call function to convert POSCAR to extxyz
    try:
        convert_poscar_to_extxyz(poscar_filenames, extxyz_filename)
    except Exception as exc:
        print(f" Error: conversion failed: {exc}")
        sys.exit(1)

    print(f" Conversion completed, {len(poscar_filenames)} file(s) saved to {extxyz_filename}")

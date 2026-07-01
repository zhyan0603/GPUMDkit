"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     lmp2exyz.py
Category:   Format Conversion Scripts
Purpose:    Convert LAMMPS dump file to extended XYZ format with proper
            element type mapping.
Usage:      gpumdkit.sh -lmp2exyz <dump_file> <element1> <element2> ...
            python lmp2exyz.py <dump_file> <element1> <element2> ...
Arguments:
  dump_file  Input LAMMPS dump file
  elementX   Chemical element symbols in order of atomic types
Output:
  dump.xyz   Converted structures in extxyz format
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

args = sys.argv[1:]
if len(args) < 2 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -lmp2exyz <dump_file> <element1> <element2> ...")
    print("    or: python lmp2exyz.py <dump_file> <element1> <element2> ...")
    print("")
    print(" Arguments:")
    print("   dump_file   Input LAMMPS dump file")
    print("   elementX    Chemical element symbols in order of atomic types")
    print("")
    print(" Output:")
    print("   dump.xyz    Converted structures in extxyz format")
    print("")
    print(" Example: gpumdkit.sh -lmp2exyz dump.lammps Si O")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

from ase.io import read, write

def lmp2exyz(dump_file, elements):
    # Read LAMMPS dump file
    frames = read(dump_file, format='lammps-dump-text', index=':')
    
    # Ensure all atomic types are within the provided elements range
    for frame in frames:
        unique_types = set(frame.get_atomic_numbers())
        if max(unique_types) > len(elements):
            raise ValueError(f" Found atomic type {max(unique_types)} which is out of the provided elements range.")
    
    # Map atomic types to specified elements
    type_to_element = {i+1: elements[i] for i in range(len(elements))}
    
    # Convert atomic types to specified elements
    for frame in frames:
        new_symbols = [type_to_element[number] for number in frame.get_atomic_numbers()]
        frame.set_chemical_symbols(new_symbols)
    
    # Write to extxyz file
    extxyz_file = 'dump.xyz'
    write(extxyz_file, frames, format='extxyz')

if __name__ == '__main__':
    dump_file = sys.argv[1]
    elements = sys.argv[2:]

    lmp2exyz(dump_file, elements)

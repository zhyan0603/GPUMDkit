"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     exyz2pos.py
Category:   Format Conversion Scripts
Purpose:    Convert extended XYZ file(s) to VASP POSCAR format(s).
Usage:      gpumdkit.sh -exyz2pos <input.xyz>
            python3 exyz2pos.py <input.xyz>
Arguments:
  input.xyz  Input extxyz file
Output:
  POSCAR_N.vasp   POSCAR file(s) in VASP format in the current directory
Notes:
  Elements are grouped using their first-seen order in the input trajectory.
  Velocities are not exported.
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
            Huan Wang (huan.wang@whut.edu.cn)
Last-modified: 2026-08-31
=============================================================================
"""

import os
import sys


def print_help():
    print(" Usage: gpumdkit.sh -exyz2pos <input.xyz>")
    print("    or: python3 exyz2pos.py <input.xyz>")
    print(" Input:")
    print("   input.xyz       Input extxyz trajectory file")
    print(" Output:")
    print("   POSCAR_1.vasp   One POSCAR file per frame in the current directory")
    print(" Notes:")
    print("   Elements are grouped using their first-seen order in the input trajectory.")
    print("   Velocities are not exported.")


def print_progress_bar(iteration, total, length=50):
    percent = f"{100 * iteration / total:.1f}"
    filled_length = int(length * iteration // total)
    bar = "#" * filled_length + "-" * (length - filled_length)
    print(f"\r Progress: |{bar}| {percent}% Complete", end="\r")
    if iteration == total:
        print()


def first_seen_element_order(frames):
    """Return unique element symbols in their first-seen order."""
    element_order = []
    seen = set()
    for atoms in frames:
        for symbol in atoms.get_chemical_symbols():
            if symbol not in seen:
                seen.add(symbol)
                element_order.append(symbol)
    return element_order


def reorder_atoms(atoms, element_order):
    """Group atoms by the supplied element order without changing coordinates."""
    symbols = atoms.get_chemical_symbols()
    indices = [
        index
        for symbol in element_order
        for index, atom_symbol in enumerate(symbols)
        if atom_symbol == symbol
    ]
    reordered = atoms[indices]
    reordered.arrays.pop("velocities", None)
    reordered.arrays.pop("momenta", None)
    return reordered


def main():
    args = sys.argv[1:]
    if args and args[0] in ("-h", "--help"):
        print_help()
        return 0
    if len(args) != 1:
        print_help()
        print(" Error: provide exactly one input extxyz file.")
        return 1

    input_file = args[0]
    if not os.path.isfile(input_file):
        print(f" Error: input file not found: {input_file}")
        return 1

    # Import ASE after handling help so command discovery does not require it.
    try:
        from ase.io import read, write
    except ImportError as exc:
        print(f" Error: ASE is required for exyz2pos: {exc}")
        return 1

    try:
        frames = read(input_file, index=":")
    except Exception as exc:
        print(f" Error: failed to read {input_file}: {exc}")
        return 1

    if not frames:
        print(f" Error: no frames found in {input_file}")
        return 1

    element_order = first_seen_element_order(frames)
    if not element_order:
        print(f" Error: no atoms found in {input_file}")
        return 1
    total_frames = len(frames)

    for index, atoms in enumerate(frames, start=1):
        output_file = f"POSCAR_{index}.vasp"
        reordered = reorder_atoms(atoms, element_order)
        try:
            write(output_file, reordered)
        except Exception as exc:
            print(f" Error: failed to write {output_file}: {exc}")
            return 1
        print_progress_bar(index, total_frames)

    print(" All frames have been converted.")
    print(f" Element order: {' '.join(element_order)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

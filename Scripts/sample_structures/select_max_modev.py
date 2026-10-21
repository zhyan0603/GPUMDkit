"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     select_max_modev.py
Category:   Sample Structure Scripts
Purpose:    Extract the top N structures with the highest max force
            deviation from active learning output (active.out + active.xyz).
Usage:      gpumdkit.sh
            choose 205) Select max force deviation structs
            python select_max_modev.py <top_n> <min_deviation>
Arguments:
  top_n          Number of top structures to extract
  min_deviation  Minimum force deviation threshold for filtering
Input files:
  active.out     Time and max force deviation of every check
  active.xyz     Structures GPUMD saved above its active threshold
Output:
  selected.xyz
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-15
=============================================================================
"""

import os
import sys


def print_usage():
    print(" Usage: gpumdkit.sh")
    print("        choose 205) Select max force deviation structs")
    print("    or: python select_max_modev.py <top_n> <min_deviation>")
    print("")
    print(" Arguments:")
    print("   top_n          Number of top structures to extract")
    print("   min_deviation  Minimum force deviation threshold for filtering")
    print("")
    print(" Input files:")
    print("   active.out     Time and max force deviation of every check")
    print("   active.xyz     Structures GPUMD saved above its active threshold")
    print("")
    print(" Output:")
    print("   selected.xyz")
    print("")
    print(" Example: in interactive mode, enter: 200 0.15")
    print("          python select_max_modev.py 200 0.15")
    print("")


args = sys.argv[1:]
if len(args) != 2 or args[0] in ("-h", "--help"):
    print_usage()
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

try:
    top_n = int(args[0])
    min_deviation = float(args[1])
except ValueError:
    print(" Error: <top_n> must be an integer and <min_deviation> must be a number.")
    print("")
    print_usage()
    sys.exit(1)

from ase.io import read, write

input_file = 'active.out'       # time and max force deviation of every check
xyz_file = 'active.xyz'         # only the checks above the GPUMD active threshold
output_xyz_file = 'selected.xyz'

for required in (input_file, xyz_file):
    if not os.path.isfile(required):
        print(f" Error: file '{required}' does not exist in the current directory.")
        sys.exit(1)


def read_active_out(filename):
    """Return (time, deviation) rows of active.out in file order."""
    rows = []
    with open(filename, 'r') as f:
        for line_number, line in enumerate(f, start=1):
            fields = line.split()
            if not fields:
                continue
            try:
                rows.append((float(fields[0]), float(fields[1])))
            except (IndexError, ValueError):
                print(f" Error: cannot parse line {line_number} of {filename}.")
                sys.exit(1)
    return rows


def frame_deviations(frames, rows):
    """Deviation of every active.xyz frame.

    GPUMD writes a line to active.out at every check but appends a frame to
    active.xyz only above its threshold, so line numbers of the two files do
    not correspond. Each frame carries its own uncertainty=; frames without it
    are matched to active.out by Time= in file order.
    """
    deviations = []
    row = 0
    for index, atoms in enumerate(frames):
        if 'uncertainty' in atoms.info:
            deviations.append(float(atoms.info['uncertainty']))
            continue
        time = atoms.info.get('Time')
        if time is None:
            print(f" Error: frame {index} of {xyz_file} has neither uncertainty= nor Time=.")
            sys.exit(1)
        while row < len(rows) and abs(rows[row][0] - float(time)) > 1e-6 * max(1.0, abs(float(time))):
            row += 1
        if row == len(rows):
            print(f" Error: frame {index} of {xyz_file} (Time={time}) has no matching line in {input_file}.")
            print(" Regenerate both files from the same GPUMD active run.")
            sys.exit(1)
        deviations.append(rows[row][1])
        row += 1
    return deviations


rows = read_active_out(input_file)
frames = read(xyz_file, index=':')
deviations = frame_deviations(frames, rows)

# Keep frames above min_deviation, then take the top N in descending order
candidates = [i for i, d in enumerate(deviations) if d > min_deviation]
selected = sorted(candidates, key=lambda i: deviations[i], reverse=True)[:top_n]

print(f' Checks in {input_file}: {len(rows)}; frames in {xyz_file}: {len(frames)}')
print(f' Frames above {min_deviation}: {len(candidates)}; selected: {len(selected)}')
print(f' Selected frame indices in {xyz_file} (0-based):\n {selected}')

if not selected:
    print(f" Error: no frame exceeds min_deviation; {output_xyz_file} was not written.")
    sys.exit(1)

write(output_xyz_file, [frames[i] for i in selected])
print(f' Extracted structures saved to {output_xyz_file}')

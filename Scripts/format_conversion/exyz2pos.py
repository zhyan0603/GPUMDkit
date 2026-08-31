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
            python exyz2pos.py <input.xyz>
Arguments:
  input.xyz  Input extxyz file
Output:
  POSCAR_*.vasp   POSCAR file(s) in VASP format (one per frame)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
            Huan Wang (huan.wang@whut.edu.cn)
Last-modified: 2026-08-31
=============================================================================
"""

import argparse as ap
from ase.io import read, write
from collections import Counter
from pathlib import Path


def parse_args():
    parser = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,
    usage="gpumdkit.sh -exyz2pos <input.xyz>\n   or: python exyz2pos.py <input.xyz>",
    epilog="output:\n"
        "  POSCAR_*.vasp   POSCAR file(s) in VASP format (one per frame)\n"
        "  (Atoms are sorted alphabetically by element symbol before writing)",
    add_help=True)

    parser.add_argument('input_file',
                        type=Path,
                        default=Path.cwd() / 'train.xyz',
                        help="Input extxyz trajectory file")
    parser.add_argument('-o', '--output_dir',
                        type=Path,
                        default=Path.cwd(),
                        help="Output directory for POSCAR files (default: current directory)",)
    parser.add_argument('-v', '--keep-velocities',
                        action='store_true',
                        help='Keep velocity information in POSCAR and label it as "Velocities"')
    return parser.parse_args()

def print_progress_bar(iteration, total, length=50):
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = '\u2588' * filled_length + '-' * (length - filled_length)
    print(f'\r Progress: |{bar}| {percent}% Complete', end='\r')
    # Print New Line on Complete
    if iteration == total:
        print()

def fix_poscar_velocity_label(filename: Path) -> None:
    """
    Replaceing the second 'Cartesian' to 'Velocities' in the POSCAR file if it contains a velocity block.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    indices = [i for i, line in enumerate(lines) if 'Cartesian' in line]
    if len(indices) >= 2:
        lines[indices[1]] = lines[indices[1]].replace('Cartesian', 'Velocities')
        with open(filename, 'w') as f:
            f.writelines(lines)
        print(f"   (Velocity label fixed in {filename})")
    else:
        print(f"   (No velocity block found in {filename})")

def main():
  args = parse_args()
  input_file = args.input_file
  output_dir = args.output_dir
  keep_vel = args.keep_velocities

  # Read all frames in the extended XYZ file
  frames = read(input_file, index=':')
  total_frames= len(frames)
  width = len(str(total_frames))

  # Save each seperated frame to POSCAR
  for i, atoms in enumerate(frames):
      indices = sorted(range(len(atoms)), key=lambda idx: atoms[idx].symbol)
      atoms = atoms[indices]

      if not keep_vel:
          for key in ['velocities', 'momenta']:
              atoms.arrays.pop(key, None)
          if hasattr(atoms, 'velocities'):
              atoms.set_velocities(None)

      poscar_filename = f'POSCAR_{i+1:0{width}d}.vasp'
      write(output_dir /  poscar_filename, atoms)

      if keep_vel:
          fix_poscar_velocity_label(output_dir / poscar_filename)

      print_progress_bar(i + 1, total_frames)

  print(f' All frames have been converted.')


if __name__ == "__main__":
    main()
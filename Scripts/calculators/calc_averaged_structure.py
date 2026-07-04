"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     calc_averaged_structure.py
Category:   Calculator Scripts
Purpose:    Calculate and save the averaged structure from trajectory frames
            using ASE and NumPy.
Usage:      gpumdkit.sh -calc avg-struct [options]
            python calc_averaged_structure.py -i <movie.xyz> -o <output.xyz> [options]
Example:    gpumdkit.sh -calc avg-struct -i movie.xyz -l 0.2 -o averaged_structure.xyz
Arguments:
  -i, --input   Input trajectory file (default: movie.xyz)
  -s            Start index
  -t            End index
  -p            Step (process every p-th frame)
  -l            Fraction or number of last frames to average
  -o, --output  Output averaged structure file (default: averaged_structure.xyz)
Output:
  Averaged structure in extxyz format
Author:     Denan LI (lidenan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

if __name__ == "__main__" and (
    len(sys.argv) == 1 or (len(sys.argv) > 1 and sys.argv[1] in ("-h", "--help"))
):
    print(" Usage: gpumdkit.sh -calc avg-struct -i <movie.xyz> -o <output.xyz> [options]")
    print("    or: python calc_averaged_structure.py -i <movie.xyz> -o <output.xyz> [options]")
    print("")
    print(" Arguments:")
    print("   -i, --input   Input trajectory file")
    print("   -o, --output  Output averaged structure file")
    print("   -s, --start   Slice start index")
    print("   -t, --stop    Slice stop index")
    print("   -p, --step    Slice step")
    print("   -l, --last    Use the last frames; integer N or 0 < ratio < 1")
    print("")
    print(" Example: gpumdkit.sh -calc avg-struct -i movie.xyz -l 0.2 -o averaged_structure.xyz")
    print("")
    sys.exit(0 if len(sys.argv) > 1 and sys.argv[1] in ("-h", "--help") else 1)

from ase.io import read, write
import argparse
import math
from ase import Atoms
import numpy as np

def __select_traj(traj: list[Atoms] | Atoms, select: list[int] | slice | None = None) -> list[Atoms]:
    nframe = len(traj)

    # 1. default: select last 50% frames
    if select is None:
        select = slice(nframe//2, nframe, 1)

    # 2. select frames
    if isinstance(traj, Atoms):
        selected_traj = [traj]
    else:
        selected_traj: list[Atoms] = traj[select]
    print(f"Number of Selected Frames: {len(selected_traj)}")
    return selected_traj

def calculate_averaged_structure(traj: list[Atoms], select: list[int] | slice | None = None) -> Atoms:
    """Compute the time-averaged atomic structure from an MD trajectory.

    Atomic coordinates are unwrapped with respect to the first selected frame
    before averaging to avoid artefacts from periodic boundary crossings.

    Parameters
    ----------
    traj : list[Atoms]
        Full MD trajectory as a list of ASE Atoms objects.
    select : list[int] | slice | None, optional
        Frame selection. ``None`` selects the last 50 % of frames.
        Defaults to ``None``.

    Returns
    -------
    Atoms
        ASE Atoms object with averaged positions and cell. Element symbols and
        PBC flags are taken from the first frame of the trajectory.
    """

    selected_traj: list[Atoms] = __select_traj(traj, select)
    coords = np.array([atoms.get_positions() for atoms in selected_traj])
    cells = np.array([atoms.get_cell().array for atoms in selected_traj])

    cells_inv = np.linalg.inv(cells)
    coords_frac = np.matmul(coords, cells_inv)
    coords_frac_diff = coords_frac - coords_frac[0]
    coords_frac[coords_frac_diff > 0.5] -= 1
    coords_frac[coords_frac_diff < -0.5] += 1
    coords_unwrapped = np.matmul(coords_frac, cells)

    # 3. update coordinates to account for PBC
    #coords_frac = np.array([np.dot(coords[i], np.linalg.inv(cells[i])) for i in range(len(coords))])
    #coords_frac_diff = coords_frac - coords_frac[0]
    #coords_frac[coords_frac_diff > 0.5] -= 1
    #coords_frac[coords_frac_diff < -0.5] += 1
    #coords = np.array([np.dot(coords_frac[i], cells[i]) for i in range(len(coords))])

    # 4. compute averaged structure
    avg_cell = np.mean(cells, axis=0)
    avg_coords = np.mean(coords_unwrapped, axis=0)
    symbols = [atom.symbol for atom in traj[0]]
    atoms = Atoms(symbols=symbols, positions=avg_coords, cell=avg_cell, pbc=True)
    return atoms


class HelpFormatter(argparse.RawTextHelpFormatter):
    def _get_help_string(self, action):
        help_msg = action.help
        if action.default is not argparse.SUPPRESS and action.default is not None:
            if "%(default)" not in help_msg:
                help_msg += " (default: %(default)s)"
        return help_msg


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate and save averaged structure from trajectory frames.",
        formatter_class=HelpFormatter,
        epilog=(
            "Examples :\n"
            "1. Use all frames:\n"
            ">>  gpumdkit.sh -calc avg-struct -i movie.xyz -o averaged_structure.xyz\n"
            ">>  python calc_averaged_structure.py -i movie.xyz -o averaged_structure.xyz\n"
            "2. Average every 2 frames from index 100 to 500:\n"
            ">>  python calc_averaged_structure.py -i movie.xyz -s 100 -t 500 -p 2 -o avg.xyz\n"
            "3. Use the last 20% frames:\n"
            ">>  python calc_averaged_structure.py -i movie.xyz -l 0.2 -o avg_last20.xyz\n"
            "4. Use the last 100 frames:\n"
            ">>  python calc_averaged_structure.py -i movie.xyz -l 100 -o avg_last100.xyz"
        ),
    )
    parser.add_argument("-i", "--input", default="movie.xyz", help="Input trajectory file.")
    parser.add_argument("-o", "--output", default="averaged_structure.xyz", help="Output structure file.")
    parser.add_argument("-s", "--start", type=int, default=0, help="Slice start index.")
    parser.add_argument("-t", "--stop", type=int, default=None, help="Slice stop index. Default is end.")
    parser.add_argument("-p", "--step", type=int, default=1, help="Slice step.")
    parser.add_argument(
        "-l",
        "--last",
        type=float,
        default=None,
        help="Use the last frames. Integer: last N frames; 0<value<1: last ratio.",
    )
    args = parser.parse_args()
    has_slice_args = any(
        flag in sys.argv[1:]
        for flag in ("-s", "--start", "-t", "--stop", "-p", "--step")
    )
    if args.last is not None and has_slice_args:
        parser.error("--last is mutually exclusive with -s/--start, -t/--stop, -p/--step.")
    return args


def get_nframe(file_name):
    with open(file_name, "r") as f:
        nframe = sum(1 for line in f if "Lattice=" in line)
    print(f"Detected {nframe} frames in {file_name}.")
    return nframe


def apply_last_to_slice(args, nframe):
    if args.last is None:
        return

    if 0 < args.last < 1:
        n_last = max(1, math.ceil(nframe * args.last))
    elif args.last > 0 and args.last.is_integer():
        n_last = int(args.last)
    else:
        raise ValueError("--last must be a positive integer, or 0 < --last < 1.")

    n_last = min(n_last, nframe)
    args.start = max(0, nframe - n_last)
    args.stop = nframe
    args.step = 1


def main():
    args = parse_args()
    nframe = get_nframe(args.input)
    apply_last_to_slice(args, nframe)
    select = slice(args.start, args.stop, args.step)
    atoms = read(args.input, index=":")
    averaged_atoms = calculate_averaged_structure(atoms, select=select)
    write(args.output, averaged_atoms)


if __name__ == "__main__":
    main()

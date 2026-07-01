"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     calc_oct_tilt.py
Category:   Calculator Scripts
Purpose:    Calculate octahedral tilt angles from trajectory/model and B-O
            neighbor list for perovskite systems using ferrodispcalc.
Usage:      python calc_oct_tilt.py -i <input.xyz> -n <nl-B-O.dat> -o <output.dat> [options]
Arguments:
  -i, --input   Input xyz file (default: model.xyz)
  -n            Neighbor list file for B-O pairs (required)
  -o, --output  Output file for tilt data
  -s            Start index for slicing frames
  -t            Stop index for slicing frames
  -p            Step for slicing frames
  -l            Use the last frames (integer for count, float for ratio)
Output:
  Octahedral tilt angles saved to output file
Author:     Denan LI (lidenan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

from ase.io import read
import argparse
import math
import numpy as np
import sys

try:
    from ferrodispcalc.compute import calculate_octahedral_tilt
except ImportError:
    raise ImportError(
        "The 'ferrodispcalc' package is required to run this script.\n"
        "Install using: `pip3 install git+https://github.com/MoseyQAQ/ferrodispcalc.git`"
    )


class HelpFormatter(argparse.RawTextHelpFormatter):
    def _get_help_string(self, action):
        help_msg = action.help
        if action.default is not argparse.SUPPRESS and action.default is not None:
            if "%(default)" not in help_msg:
                help_msg += " (default: %(default)s)"
        return help_msg


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate octahedral tilt from trajectory/model and B-O neighbor list.",
        formatter_class=HelpFormatter,
        epilog=(
            "Examples :\n"
            "1. Single-frame model.xyz:\n"
            ">>  gpumdkit.sh -calc oct-tilt -i model.xyz -n nl-Ti-O.dat -o oct_tilt.dat\n"
            ">>  python calc_oct_tilt.py -i model.xyz -n nl-Ti-O.dat -o oct_tilt.dat\n"
            "2. Multi-frame movie.xyz with slice:\n"
            ">>  python calc_oct_tilt.py -i movie.xyz -n nl-Ti-O.dat -s 100 -t 500 -p 2 -o oct_tilt_slice.dat\n"
            "3. Use the last 20% frames:\n"
            ">>  python calc_oct_tilt.py -i movie.xyz -n nl-Ti-O.dat -l 0.2 -o oct_tilt_last20.dat"
        ),
    )
    parser.add_argument("-i", "--input", default="model.xyz", help="Input xyz file.")
    parser.add_argument(
        "-n",
        "--neighbor-list",
        required=True,
        help="B-O neighbor list file.",
    )
    parser.add_argument("-o", "--output", default="octahedral_tilt.dat", help="Output text file.")
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


def save_tilt(file_name, oct_tilt):
    save_data = oct_tilt
    if oct_tilt.ndim > 2:
        save_data = oct_tilt.reshape(-1, 3)
    np.savetxt(file_name, save_data)
    print(f"Saved octahedral tilt data to {file_name}, shape={save_data.shape}.")


def main():
    args = parse_args()
    nframe = get_nframe(args.input)
    apply_last_to_slice(args, nframe)
    select = slice(args.start, args.stop, args.step)

    nl = np.loadtxt(args.neighbor_list, dtype=int)
    traj = read(args.input, index=":")
    oct_tilt = calculate_octahedral_tilt(traj, nl, select=select)
    save_tilt(args.output, oct_tilt)


if __name__ == "__main__":
    main()

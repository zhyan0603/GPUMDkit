"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     calc_neighbor_list.py
Category:   Calculator Scripts
Purpose:    Build and save neighbor lists between center and neighbor atoms
            using the ferrodispcalc package.
Usage:      python calc_neighbor_list.py -c <cutoff> -n <neighbor_num>
            -C <center_elements> -E <neighbor_elements> [options]
Arguments:
  -i, --input       Input structure file (default: model.xyz)
  -x, --index       Which frame to read from the input file (default: 0)
  -c, --cutoff      Cutoff distance (required)
  -n, --neighbor-num  Number of neighbors (required)
  -C                Center element symbols (e.g., Pb Sr)
  -E                Neighbor element symbols (e.g., O)
  -o, --output      Output file path (default: nl-<center>-<neighbor>.dat)
Output:
  nl-<center>-<neighbor>.dat (neighbor list file)
Author:     Denan LI (lidenan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

from ase.io import read
import argparse

try:
    from ferrodispcalc.neighborlist import build_neighbor_list, save_neighbor_list
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
    default_output_name = "nl-<center>-<neighbor>.dat"
    parser = argparse.ArgumentParser(
        description="Build and save neighbor list.",
        formatter_class=HelpFormatter,
        epilog=(
            "Examples :\n"
            "1. Find the Pb and Sr's nearest 12 oxygen atoms: \n"
            ">>  gpumdkit.sh -calc nlist -c 4 -n 12 -C Pb Sr -E O\n"
            ">>  python calc_neighbor_list.py -c 4 -n 12 -C Pb Sr -E O\n"
        ),
    )
    parser.add_argument("-i", "--input", default="model.xyz", help="Input structure file.")
    parser.add_argument("-x", "--index", default="0", help="Which frame to read from the input file.")
    parser.add_argument("-c", "--cutoff", type=float, required=True, help="Cutoff distance.")
    parser.add_argument(
        "-n",
        "--neighbor-num",
        type=int,
        required=True,
        help="Number of neighbors.",
    )
    parser.add_argument("-d", "--defect", action="store_true", help="Enable defect mode.")
    parser.add_argument(
        "-C",
        "--center-elements",
        nargs="+",
        required=True,
        help="Center element symbols (space separated).",
    )
    parser.add_argument(
        "-E",
        "--neighbor-elements",
        nargs="+",
        required=True,
        help="Neighbor element symbols (space separated).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=default_output_name,
        help="Output file path.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    atoms = read(args.input, index=args.index)
    nl = build_neighbor_list(
        atoms=atoms,
        center_elements=args.center_elements,
        neighbor_elements=args.neighbor_elements,
        cutoff=args.cutoff,
        neighbor_num=args.neighbor_num,
        defect=args.defect,
    )
    if args.output == "nl-<center>-<neighbor>.dat":
        save_file_name = (
            f'nl-{"_".join(args.center_elements)}-{"_".join(args.neighbor_elements)}.dat'
        )
    else:
        save_file_name = args.output
    save_neighbor_list(nl, save_file_name)


if __name__ == "__main__":
    main()

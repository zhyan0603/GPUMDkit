from ase.io import read, write
import argparse
import math
import sys

try:
    from ferrodispcalc.compute import calculate_averaged_structure
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
        description="Calculate and save averaged structure from trajectory frames.",
        formatter_class=HelpFormatter,
        epilog=(
            "Examples :\n"
            "1. Use all frames:\n"
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

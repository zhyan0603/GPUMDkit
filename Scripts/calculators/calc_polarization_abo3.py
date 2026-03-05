from ase.io import read
import argparse
import math
import numpy as np
import sys

try:
    from ferrodispcalc.compute import calculate_polarization
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
        description="Calculate local polarization for ABO3 systems from neighbor lists.",
        formatter_class=HelpFormatter,
        epilog=(
            "Examples :\n"
            "1. Single-frame model.xyz:\n"
            ">>  python calc_polarization_abo3.py -i model.xyz --nl-ba nl-Ti-A.dat --nl-bo nl-Ti-O.dat\n"
            "2. Multi-frame movie.xyz with slice:\n"
            ">>  python calc_polarization_abo3.py -i movie.xyz --nl-ba nl-Ti-A.dat --nl-bo nl-Ti-O.dat -s 100 -t 500 -p 2 -o pol_slice.dat\n"
            "3. Set Born effective charges:\n"
            ">>  python calc_polarization_abo3.py --nl-ba nl-Ti-A.dat --nl-bo nl-Ti-O.dat --bec Pb=2.5 Sr=2.0 Ti=4.0 O=-2.0"
        ),
    )
    parser.add_argument("-i", "--input", default="model.xyz", help="Input xyz file.")
    parser.add_argument("--nl-ba", required=True, help="Neighbor list file for B-A pairs.")
    parser.add_argument("--nl-bo", required=True, help="Neighbor list file for B-O pairs.")
    parser.add_argument("-o", "--output", default="polarization.dat", help="Output text file.")
    parser.add_argument(
        "--bec",
        nargs="+",
        required=True,
        help="Born effective charge terms in the format Element=value.",
    )
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


def parse_bec_list(bec_terms):
    bec = {}
    for term in bec_terms:
        if "=" not in term:
            raise ValueError(
                f"Invalid --bec term '{term}'. Expected format Element=value, e.g. Ti=4.0."
            )
        key, value = term.split("=", 1)
        key = key.strip()
        if not key:
            raise ValueError(f"Invalid --bec term '{term}': empty element symbol.")
        try:
            bec[key] = float(value)
        except ValueError as exc:
            raise ValueError(
                f"Invalid --bec term '{term}': value must be a float."
            ) from exc
    return bec


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


def save_polarization(file_name, polarization):
    save_data = polarization
    if polarization.ndim > 2:
        save_data = polarization.reshape(-1, 3)
    np.savetxt(file_name, save_data)
    print(f"Saved polarization data to {file_name}, shape={save_data.shape}.")


def main():
    args = parse_args()
    nframe = get_nframe(args.input)
    apply_last_to_slice(args, nframe)
    select = slice(args.start, args.stop, args.step)
    bec = parse_bec_list(args.bec)

    nl_ba = np.loadtxt(args.nl_ba, dtype=int)
    nl_bo = np.loadtxt(args.nl_bo, dtype=int)
    traj = read(args.input, index=":")

    polarization = calculate_polarization(
        traj=traj,
        nl_ba=nl_ba,
        nl_bx=nl_bo,
        born_effective_charge=bec,
        select=select,
    )
    save_polarization(args.output, polarization)


if __name__ == "__main__":
    main()

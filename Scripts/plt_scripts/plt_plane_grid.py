from ase.io import read
import argparse
import numpy as np

try:
    from ferrodispcalc.vis import plane_profile, grid_data
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
        description="Plot plane profiles from displacement grid data.",
        formatter_class=HelpFormatter,
        epilog=(
            "Examples :\n"
            "1. Plot single-frame displacement:\n"
            ">>  python plt_plane_grid.py -i model.xyz -d displacements.dat -e Pb Sr\n"
            "2. Plot frame 3 from movie.xyz:\n"
            ">>  python plt_plane_grid.py -i movie.xyz -d displacements.dat -e Pb Sr -f 3"
        ),
    )
    parser.add_argument("-i", "--input", default="model.xyz", help="Input xyz file.")
    parser.add_argument("-d", "--disp", default="displacements.dat", help="Displacement data file.")
    parser.add_argument("-e", "--elements", nargs="+", required=True, help="Center element symbols.")
    parser.add_argument("-m", "--tol", type=float, default=1.0, help="Tolerance for grid mapping.")
    parser.add_argument(
        "-g",
        "--target-size",
        nargs=3,
        type=int,
        default=[10, 10, 10],
        help="Grid size as nx ny nz.",
    )
    parser.add_argument("-o", "--save-dir", default="plot", help="Directory to save figures.")
    parser.add_argument(
        "-f",
        "--frame",
        type=int,
        default=0,
        help="Absolute frame index to plot (only one frame is plotted).",
    )

    parser.add_argument(
        "--select-xy",
        nargs="+",
        type=int,
        default=None,
        help="XY plane indices. Default: plot all XY layers.",
    )
    parser.add_argument(
        "--select-xz",
        nargs="+",
        type=int,
        default=None,
        help="XZ plane indices. Default: plot all XZ layers.",
    )
    parser.add_argument(
        "--select-yz",
        nargs="+",
        type=int,
        default=None,
        help="YZ plane indices. Default: plot all YZ layers.",
    )

    return parser.parse_args()


def get_nframe(file_name):
    with open(file_name, "r") as f:
        nframe = sum(1 for line in f if "Lattice=" in line)
    if nframe == 0:
        nframe = 1
    print(f"Detected {nframe} frames in {file_name}.")
    return nframe


def count_elements(atoms, elements):
    symbols = np.array(atoms.get_chemical_symbols())
    return int(np.sum(np.isin(symbols, elements)))


def load_disp_frame(file_name, n_ele, frame_idx):
    disp = np.loadtxt(file_name)
    assert disp.ndim == 2, f"Expected 2D array, got {disp.ndim}D."

    nline = disp.shape[0]
    if nline % n_ele != 0:
        raise ValueError(f"Invalid disp shape {disp.shape}: nline must be divisible by n_ele={n_ele}.")

    nframe_disp = nline // n_ele
    print(f"Detected {nframe_disp} displacement frames in {file_name}.")
    disp = disp.reshape(nframe_disp, n_ele, -1)

    if nframe_disp == 1:
        return disp[0]

    if frame_idx < 0 or frame_idx >= nframe_disp:
        raise ValueError(
            f"Frame index {frame_idx} is out of range for displacement frames ({nframe_disp})."
        )
    return disp[frame_idx]


def build_select_dict(args):
    select = {}
    if args.select_yz is not None:
        select["x"] = args.select_yz
    if args.select_xz is not None:
        select["y"] = args.select_xz
    if args.select_xy is not None:
        select["z"] = args.select_xy
    return select if select else None


def main():
    args = parse_args()
    nframe_struct = get_nframe(args.input)

    atoms = read(args.input)
    n_ele = count_elements(atoms, args.elements)
    if n_ele == 0:
        raise ValueError(f"No atoms found for elements: {args.elements}")

    disp_frame = load_disp_frame(args.disp, n_ele, args.frame)
    disp_grid = grid_data(
        atoms,
        disp_frame,
        args.elements,
        tol=args.tol,
        target_size=tuple(args.target_size),
    )
    plane_profile(
        data=disp_grid,
        save_dir=args.save_dir,
        select=build_select_dict(args),
    )


if __name__ == "__main__":
    main()

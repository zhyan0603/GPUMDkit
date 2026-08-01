#!/usr/bin/env python3
"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_xrd.py
Category:   Plot Scripts
Purpose:    Plot X-ray diffraction data from an XRD calculator output file.
            The second numeric column is used as the angle and the fourth
            numeric column is used as the intensity.
Usage:      gpumdkit.sh -plt xrd [xrd_file] [save]
            python plt_xrd.py [xrd_file] [save]
Arguments:
  xrd_file  XRD output file; defaults to 'xrd.out'
  save      Save the plot as 'xrd.png' instead of displaying it
Output:
  xrd.png   XRD plot (saved or displayed)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-08-01
=============================================================================
"""

from __future__ import annotations

import sys
from pathlib import Path


def _print_help() -> None:
    """Print the command-line usage without importing plotting packages."""
    print(" Usage: gpumdkit.sh -plt xrd [xrd_file] [save]")
    print("        python plt_xrd.py [xrd_file] [save]")
    print("")
    print(" xrd_file  XRD output file; defaults to 'xrd.out'")
    print(" save      Save the plot as 'xrd.png' instead of displaying it")
    print("")
    print(" Data columns: column 2 is the angle and column 4 is the intensity.")


def _parse_arguments() -> tuple[Path, bool] | None:
    """Parse the small, positional interface used by GPUMDkit plot scripts."""
    arguments = sys.argv[1:]
    if arguments and arguments[0] in {"-h", "--help"}:
        _print_help()
        return None

    if len(arguments) > 2:
        raise ValueError("expected at most an XRD file and the optional 'save'")

    save = False
    if arguments and arguments[-1] == "save":
        save = True
        arguments = arguments[:-1]

    if arguments and arguments[0] == "save":
        raise ValueError("'save' must be the last argument")

    input_file = Path(arguments[0]) if arguments else Path("xrd.out")
    return input_file, save


def _read_xrd_data(input_file: Path):
    """Read angle and intensity from columns 2 and 4 of an XRD output file."""
    import numpy as np

    if not input_file.is_file():
        raise OSError(f"XRD input file not found: {input_file}")

    try:
        data = np.loadtxt(input_file, comments="#", ndmin=2)
    except ValueError as error:
        raise ValueError(f"unable to read numeric XRD data from '{input_file}'") from error

    if data.ndim != 2 or data.shape[0] == 0 or data.shape[1] < 4:
        raise ValueError(
            "XRD data must contain at least four numeric columns "
            "(angle in column 2 and intensity in column 4)"
        )

    angle = data[:, 1]
    intensity = data[:, 3]
    if not np.all(np.isfinite(angle)) or not np.all(np.isfinite(intensity)):
        raise ValueError("XRD angle and intensity columns must contain finite values")

    return angle, intensity


def _plot_xrd(angle, intensity, save: bool) -> None:
    """Display the XRD curve or save it for a non-interactive workflow."""
    import matplotlib.pyplot as plt
    from matplotlib import get_backend

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    })

    plt.figure(figsize=(5.0, 2.5), dpi=150)
    plt.plot(angle, intensity, color="C0", linewidth=1)
    plt.xlabel(r"$2\theta$ (degree)")
    plt.ylabel("Intensity")
    plt.tight_layout()

    backend = get_backend().lower()
    if save or backend in {"agg", "cairo", "pdf", "ps", "svg"}:
        plt.savefig("xrd.png", dpi=300)
        print(" Saved XRD plot to 'xrd.png'.")
    else:
        plt.show()


def main() -> int:
    """Run the XRD plotting command."""
    try:
        parsed = _parse_arguments()
        if parsed is None:
            return 0
        input_file, save = parsed
        angle, intensity = _read_xrd_data(input_file)
        _plot_xrd(angle, intensity, save)
    except (OSError, ValueError) as error:
        print(f" Error: {error}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

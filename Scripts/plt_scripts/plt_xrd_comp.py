#!/usr/bin/env python3
"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_xrd_comp.py
Category:   Plot Scripts
Purpose:    Compare X-ray diffraction curves from several temperature folders.
            The current working directory must contain folders named like
            <temperature>K, each with an xrd.out file produced by calculator
            413.
Usage:      gpumdkit.sh -plt xrd_comp [save]
            python3 plt_xrd_comp.py [save]
Arguments:
  save      Save the comparison figure as 'xrd_comp.png' instead of displaying
            it when an interactive Matplotlib backend is available.
Input:
  <temperature>K/xrd.out  One calculator-413 XRD output for each temperature.
Output:
  xrd_comp.png  Comparison figure when saving is requested or the selected
                Matplotlib backend cannot display figures.
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-21
=============================================================================
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


def _print_help() -> None:
    """Print command-line usage without importing plotting dependencies."""
    print(" Usage: gpumdkit.sh -plt xrd_comp [save]")
    print("    or: python3 plt_xrd_comp.py [save]")
    print("")
    print(" Input:  <temperature>K/xrd.out in the current working directory")
    print(" save:   Save the figure as 'xrd_comp.png'")
    print("")
    print(" Example: cd /path/to/xrd_series && gpumdkit.sh -plt xrd_comp save")


args = sys.argv[1:]
if args and args[0] in ("-h", "--help"):
    _print_help()
    sys.exit(0)
if len(args) > 1 or (args and args[0] != "save"):
    print(" Error: expected only the optional 'save' argument.")
    print("")
    _print_help()
    sys.exit(1)

save = bool(args and args[0] == "save")


def find_temperature_dirs(base_dir: Path) -> list[tuple[int, Path]]:
    """Find integer-temperature folders that contain an XRD output file."""
    if not base_dir.is_dir():
        raise OSError(f"current working directory is not accessible: {base_dir}")

    temperature_dirs = []
    for subdirectory in base_dir.iterdir():
        match = re.fullmatch(r"(\d+)K", subdirectory.name)
        if not subdirectory.is_dir() or match is None:
            continue
        xrd_file = subdirectory / "xrd.out"
        if not xrd_file.is_file():
            print(f" [Warning] Missing xrd.out in {subdirectory.name}; skipping")
            continue
        temperature_dirs.append((int(match.group(1)), subdirectory))

    temperature_dirs.sort(key=lambda item: item[0], reverse=True)
    return temperature_dirs


def read_xrd(xrd_file: Path):
    """Read and validate one calculator-413 XRD output file."""
    from plt_xrd import _read_xrd_data

    try:
        return _read_xrd_data(xrd_file)
    except (OSError, ValueError) as error:
        raise ValueError(f"failed to read '{xrd_file}': {error}") from error


def plot_comparison(temperature_dirs: list[tuple[int, Path]], save_plot: bool) -> None:
    """Display or save XRD curves in one panel, shifted along y with temperature."""
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import get_backend
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    # Plot style (self-contained; keep in sync across plt_*.py scripts)
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "mathtext.fontset": "dejavusans",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.linewidth": 0.9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.minor.size": 2.0,
        "ytick.minor.size": 2.0,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "figure.dpi": 100,
    })
    Y_STEP = 1.0  # vertical shift between neighbouring temperatures, in units of the highest peak (1.0: peaks never reach the curve above)

    datasets = []
    for temperature, directory in temperature_dirs:
        angle, intensity = read_xrd(directory / "xrd.out")
        datasets.append((temperature, angle, intensity))
    datasets.sort(key=lambda item: item[0])  # lowest temperature at the bottom

    temperatures = [temperature for temperature, _, _ in datasets]
    # Colors run from cold (blue) through teal and gold to hot (red)
    cmap = LinearSegmentedColormap.from_list(
        "xrd_temperature", ["#457B9D", "#2A9D8F", "#8AB17D", "#E9C46A", "#F4A261", "#D62828"])
    norm = Normalize(vmin=min(temperatures), vmax=max(temperatures))

    step = Y_STEP * max(float(intensity.max()) for _, _, intensity in datasets)
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    for index, (temperature, angle, intensity) in enumerate(datasets):
        color = cmap(norm(temperature))
        offset = index * step
        ax.plot(angle, intensity + offset, color=color, linewidth=1.2)
        ax.annotate(
            f"{temperature} K",
            xy=(0.02, offset),
            xycoords=("axes fraction", "data"),
            xytext=(0, 2),
            textcoords="offset points",
            ha="left",
            va="bottom",
            fontsize=8,
            color=color,
        )

    angle_min = min(float(angle.min()) for _, angle, _ in datasets)
    angle_max = max(float(angle.max()) for _, angle, _ in datasets)
    if angle_min >= angle_max:
        raise ValueError("XRD angle data must span a non-zero range")
    ax.set_xlim(angle_min, angle_max)
    ax.set_ylim(-0.05 * step, (len(datasets) - 1) * step + 1.05 * step)
    ax.set_yticks([])
    ax.set_xlabel(r"$2\theta$ (degree)")

    # Colorbar as an inset so that it spans exactly the height of the axes
    colorbar_axes = ax.inset_axes([1.03, 0.0, 0.035, 1.0])
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=colorbar_axes)
    colorbar.set_label("Temperature (K)", fontsize=9)
    colorbar.ax.tick_params(labelsize=8, length=3, width=0.8)

    # Subplot geometry (figure fractions); tune these to adjust margins
    fig.subplots_adjust(left=0.035, right=0.84, top=0.97, bottom=0.13)

    if save_plot or get_backend().lower() in {"agg", "cairo", "pdf", "ps", "svg"}:
        fig.savefig("xrd_comp.png", dpi=300)
        print(" Saved XRD comparison plot to 'xrd_comp.png'.")
    else:
        plt.show()
    plt.close(fig)


def main() -> int:
    """Find temperature folders, validate their XRD files, and plot them."""
    try:
        temperature_dirs = find_temperature_dirs(Path.cwd())
        if not temperature_dirs:
            raise ValueError(
                "no '<temperature>K/xrd.out' files found in the current directory"
            )
        plot_comparison(temperature_dirs, save)
    except (ImportError, OSError, ValueError) as error:
        print(f" Error: {error}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

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
Last-modified: 2026-08-29
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
    """Display or save stacked XRD curves ordered from high to low temperature."""
    import matplotlib.pyplot as plt
    from matplotlib import get_backend

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    })

    datasets = []
    for temperature, directory in temperature_dirs:
        angle, intensity = read_xrd(directory / "xrd.out")
        datasets.append((temperature, angle, intensity))

    number_of_temperatures = len(datasets)
    figure_height = max(1.1, number_of_temperatures * 0.65)
    figure, axes = plt.subplots(
        number_of_temperatures,
        1,
        figsize=(6.0, figure_height),
        dpi=150,
        sharex=True,
        squeeze=False,
    )
    axes = axes[:, 0]
    figure.subplots_adjust(
        hspace=0,
        left=0.09,
        right=0.97,
        top=0.97,
        bottom=0.07,
    )

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for index, (temperature, angle, intensity) in enumerate(datasets):
        axis = axes[index]
        color = colors[(number_of_temperatures - 1 - index) % len(colors)]
        axis.plot(angle, intensity, color=color, linewidth=1.2)
        axis.tick_params(labelsize=8, bottom=True, labelbottom=True)
        axis.set_ylim(bottom=0)
        axis.set_yticks([])
        axis.annotate(
            f"{temperature} K",
            xy=(0.98, 0.90),
            xycoords="axes fraction",
            ha="right",
            va="top",
            fontsize=8,
            color=color,
        )

    angle_min = min(float(angle.min()) for _, angle, _ in datasets)
    angle_max = max(float(angle.max()) for _, angle, _ in datasets)
    if angle_min >= angle_max:
        raise ValueError("XRD angle data must span a non-zero range")
    axes[0].set_xlim(angle_min, angle_max)
    axes[-1].set_xlabel(r"$2\theta$ (degree)", fontsize=9)
    axes[number_of_temperatures // 2].set_ylabel("Intensity", fontsize=9)

    if save_plot or get_backend().lower() in {"agg", "cairo", "pdf", "ps", "svg"}:
        figure.savefig("xrd_comp.png", dpi=300)
        print(" Saved stacked XRD plot to 'xrd_comp.png'.")
    else:
        plt.show()
    plt.close(figure)


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

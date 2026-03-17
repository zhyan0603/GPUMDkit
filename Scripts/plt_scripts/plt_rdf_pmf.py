#!/usr/bin/env python3
"""
RDF + PMF analysis for GPUMD rdf.out (single total RDF column).

This script reads GPUMD's `rdf.out` (radius, total RDF, [optional pair columns]),
computes:
  - Smoothed RDF g(r)
  - Potential of mean force (PMF) from g(r)
  - First and second RDF peaks
  - First minimum between the first two peaks
and plots:
  - Top panel: raw and smoothed RDF with peaks/minimum marked
  - Bottom panel: PMF with the bound-state minimum, barrier maximum, and ΔG

Usage:
  conda activate pymatgen_env
  python plot_rdf_pmf.py rdf.out -T 2000

Output:
  - PNG figure next to input file: <rdf_root>_pmf.png
  - JSON with extracted metrics:   <rdf_root>_pmf.json
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import Axes
from scipy.signal import find_peaks, savgol_filter


def configure_backend() -> None:
    """Choose a sensible matplotlib backend (Agg on headless, default otherwise)."""
    backend = matplotlib.get_backend()
    if not os.environ.get("DISPLAY") and "agg" not in str(backend).lower():
        matplotlib.use("Agg", force=True)


configure_backend()


def round_to_nice_bound(x: np.ndarray) -> tuple[float, float, int]:
    """Calculate nice rounded bounds for an array of numbers."""
    x = np.asarray(x)
    x_min, x_max = float(x.min()), float(x.max())
    x_range = x_max - x_min
    if x_range <= 0:
        return x_min, x_max, 0
    round_digit = int(math.floor(math.log10(x_range)) - 1)
    x_lower = math.floor(x_min * 10 ** (-round_digit)) / 10 ** (-round_digit)
    x_upper = math.ceil(x_max * 10 ** (-round_digit)) / 10 ** (-round_digit)
    return x_lower, x_upper, round_digit


def set_nice_ticks(ax: Axes, axis: str = "x") -> None:
    """Set nice rounded ticks on a given axis using round_to_nice_bound."""
    if axis == "x":
        vals = ax.get_xticks()
        lower, upper, rd = round_to_nice_bound(vals)
        ax.set_xlim(lower, upper)
        loc = plt.MaxNLocator(nbins=5)
        ax.xaxis.set_major_locator(loc)
        ax.figure.canvas.draw()
        new_ticks = ax.get_xticks()
        fmt = f"{{:.{-rd}f}}" if rd < 0 else "{:.0f}"
        ax.set_xticks(new_ticks)
        ax.set_xticklabels([fmt.format(t) for t in new_ticks])
    elif axis == "y":
        vals = ax.get_yticks()
        lower, upper, rd = round_to_nice_bound(vals)
        ax.set_ylim(lower, upper)
        loc = plt.MaxNLocator(nbins=5)
        ax.yaxis.set_major_locator(loc)
        ax.figure.canvas.draw()
        new_ticks = ax.get_yticks()
        fmt = f"{{:.{-rd}f}}" if rd < 0 else "{:.0f}"
        ax.set_yticks(new_ticks)
        ax.set_yticklabels([fmt.format(t) for t in new_ticks])


def load_gpumd_rdf(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """
    Load GPUMD rdf.out and return (r, g_all, labels).

    Columns:
      col 0: radius r
      col 1: total RDF
      col 2+: pair-resolved RDFs (if present)
    """
    # Try to read a header for labels; fall back to generic names.
    header_labels: list[str] = []
    with path.open("r", encoding="utf-8") as f:
        first = f.readline().strip()
        if first.startswith("#"):
            header_labels = first.lstrip("#").split()
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError(f"{path} must have at least 2 columns (radius and total RDF)")

    r = data[:, 0]
    g_all = data[:, 1:]

    # Build labels: use header if available and consistent, else generic.
    if header_labels and len(header_labels) - 1 == g_all.shape[1]:
        labels = header_labels[1:]
    else:
        labels = ["total"] + [f"pair_{i}" for i in range(2, g_all.shape[1] + 1)]

    return r, g_all, labels


def compute_pmf_from_rdf(g: np.ndarray, temperature_K: float) -> np.ndarray:
    """
    Compute PMF from RDF:
      PMF(r) = -k_B T ln g(r), converted to kJ/mol as in rdf_pmf_analysis.py.
    """
    boltzmann_constant = 8.617333262145e-5  # eV/K
    ev_to_kj_per_mol = 96.48533212
    energy_factor = boltzmann_constant * temperature_K * ev_to_kj_per_mol
    g_safe = np.clip(g, 1e-8, None)
    pmf = -energy_factor * np.log(g_safe)
    return pmf


def analyze_rdf_pmf(r: np.ndarray, g: np.ndarray, pmf: np.ndarray) -> dict:
    """
    Smooth RDF, find first/second RDF peaks, first minimum between them,
    and associated PMF barrier.
    """
    # Smooth RDF
    win = min(5, g.size if g.size % 2 == 1 else g.size - 1)
    if win < 5:
        g_smooth = g.copy()
    else:
        g_smooth = savgol_filter(x=g, window_length=win, polyorder=3)

    # Find peaks in smoothed RDF
    peak_indices, _ = find_peaks(g_smooth, distance=1.0, prominence=0.05)
    if peak_indices.size < 2:
        raise RuntimeError("Could not find two RDF peaks in smoothed g(r).")

    first_peak_index = int(peak_indices[0])
    second_peak_index = int(peak_indices[1])

    first_max_r = float(r[first_peak_index])
    first_max_g = float(g_smooth[first_peak_index])
    second_max_r = float(r[second_peak_index])
    second_max_g = float(g_smooth[second_peak_index])

    # First minimum between first and second peaks
    search_region = slice(first_peak_index, second_peak_index)
    local_min_relative = int(np.argmin(g_smooth[search_region]))
    first_min_idx = first_peak_index + local_min_relative

    first_min_r = float(r[first_min_idx])
    first_min_g = float(g_smooth[first_min_idx])

    # PMF at first peak and first minimum
    pmf_at_first_peak = float(pmf[first_peak_index])
    pmf_at_first_min = float(pmf[first_min_idx])
    delta_G = pmf_at_first_min - pmf_at_first_peak

    return {
        "first_peak_index": first_peak_index,
        "second_peak_index": second_peak_index,
        "first_min_index": first_min_idx,
        "first_max_r": first_max_r,
        "first_max_g": first_max_g,
        "second_max_r": second_max_r,
        "second_max_g": second_max_g,
        "first_min_r": first_min_r,
        "first_min_g": first_min_g,
        "pmf_at_first_peak": pmf_at_first_peak,
        "pmf_at_first_min": pmf_at_first_min,
        "delta_G": delta_G,
        "g_smooth": g_smooth,
    }


def plot_rdf_pmf(
    r: np.ndarray,
    g: np.ndarray,
    pmf: np.ndarray,
    analysis: dict,
    temperature_K: float,
    output_png: Path,
) -> None:
    """Create the RDF + PMF figure, following the style of rdf_pmf_analysis.py."""
    first_max_r = analysis["first_max_r"]
    first_max_g = analysis["first_max_g"]
    second_max_r = analysis["second_max_r"]
    second_max_g = analysis["second_max_g"]
    first_min_r = analysis["first_min_r"]
    first_min_g = analysis["first_min_g"]
    pmf_at_first_peak = analysis["pmf_at_first_peak"]
    pmf_at_first_min = analysis["pmf_at_first_min"]
    delta_G = analysis["delta_G"]
    g_smooth = analysis["g_smooth"]

    # Colors
    color_rdf_raw = (0.00, 0.45, 0.74)
    color_rdf_smooth = (0.30, 0.75, 0.93)
    color_pmf = (0.85, 0.33, 0.10)
    color_marker = (0.93, 0.69, 0.13)
    color_ref = (0.5, 0.5, 0.5)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    # RDF plot
    ax1.plot(r, g, lw=2, color=color_rdf_raw, label="Raw g(r)")
    ax1.plot(r, g_smooth, lw=2, color=color_rdf_smooth, label="Smoothed g(r)")

    ax1.scatter(first_max_r, first_max_g, s=50, color=color_marker, edgecolor="k", zorder=5)
    ax1.scatter(first_min_r, first_min_g, s=50, color=color_marker, edgecolor="k", zorder=5)
    ax1.scatter(second_max_r, second_max_g, s=50, color=color_marker, edgecolor="k", zorder=5)
    ax1.axhline(1.0, ls="--", color=color_ref, lw=1)

    ax1.set_ylabel("RDF, g(r)")
    ax1.legend(frameon=False)
    ax1.set_ylim(0)
    set_nice_ticks(ax1, axis="y")

    # Annotate RDF
    ax1.text(
        x=first_max_r + 0.2,
        y=first_max_g,
        s=f"(r = {first_max_r:.3f} Å, g = {first_max_g:.3f})",
        verticalalignment="center",
    )
    ax1.text(
        x=first_min_r,
        y=first_min_g - 0.2,
        s=f"(r = {first_min_r:.3f} Å, g = {first_min_g:.3f})",
        verticalalignment="center",
        horizontalalignment="center",
    )
    ax1.text(
        x=second_max_r,
        y=second_max_g + 0.2,
        s=f"(r = {second_max_r:.3f} Å, g = {second_max_g:.3f})",
        verticalalignment="center",
        horizontalalignment="center",
    )

    # PMF plot
    ax2.plot(r, pmf, lw=2, color=color_pmf)

    ax2.scatter(first_max_r, pmf_at_first_peak, s=50, color=color_marker, edgecolor="k", zorder=5)
    ax2.scatter(first_min_r, pmf_at_first_min, s=50, color=color_marker, edgecolor="k", zorder=5)

    ax2.axhline(0, ls="--", color=color_ref, lw=1)
    ax2.set_ylabel("PMF (kJ/mol)")
    ax2.set_xlabel("Radial distance, r (Å)")

    # Annotate PMF barrier
    ax2.vlines(x=first_min_r, ymin=pmf_at_first_peak, ymax=pmf_at_first_min)
    ax2.hlines(xmin=first_min_r - 0.2, xmax=first_min_r + 0.2, y=pmf_at_first_peak)
    ax2.hlines(xmin=first_min_r - 0.2, xmax=first_min_r + 0.2, y=pmf_at_first_min)
    ax2.text(
        x=first_min_r + 0.2,
        y=0.5 * (pmf_at_first_min + pmf_at_first_peak),
        s=f"ΔG = {delta_G:.3f} kJ/mol @ {temperature_K:.1f} K",
        verticalalignment="center",
    )

    # Symmetric y-limit
    y_lim = float(np.ceil(max(abs(pmf_at_first_peak), abs(pmf_at_first_min))))
    ax2.set_ylim(-y_lim, y_lim)
    set_nice_ticks(ax2, axis="y")

    # Reasonable x-range if RDF extends further
    if r.max() - r.min() > 0:
        ax2.set_xlim(r.min(), r.max())

    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    plt.close()


def rdf_pmf_analysis(
    input_data_file: Path, temperature_kelvin: float, column_index: int
) -> None:
    """
    Perform RDF+PMF analysis on a chosen RDF column.

    column_index is 1-based with respect to the data file:
      1 -> radius (ignored)
      2 -> total RDF (default)
      3+ -> specific pair RDFs (if present)
    """
    r, g_all, labels = load_gpumd_rdf(input_data_file)

    ncols = g_all.shape[1]
    # Convert 1-based file column index (2..ncols+1) to 0-based index into g_all.
    if column_index < 2 or column_index > ncols + 1:
        raise ValueError(
            f"Invalid column index {column_index}. "
            f"Valid RDF columns are 2..{ncols+1} (1=radius)."
        )
    g_idx = column_index - 2
    g = g_all[:, g_idx]

    pmf = compute_pmf_from_rdf(g, temperature_kelvin)
    analysis = analyze_rdf_pmf(r, g, pmf)

    file_root = str(input_data_file.with_suffix(""))
    output_image_file = Path(file_root + "_pmf.png")
    output_json_file = Path(file_root + "_pmf.json")

    plot_rdf_pmf(r, g, pmf, analysis, temperature_kelvin, output_image_file)

    # Collect results into dictionary
    rdf_pmf_data_dict = {
        "input_data_file": str(input_data_file),
        "temperature_K": f"{temperature_kelvin:.3f}",
        "rdf_column_index": column_index,
        "rdf_column_label": labels[g_idx] if 0 <= g_idx < len(labels) else "",
        "rdf": {
            "first_peak": {
                "r": f"{analysis['first_max_r']:.3f}",
                "g": f"{analysis['first_max_g']:.3f}",
            },
            "second_peak": {
                "r": f"{analysis['second_max_r']:.3f}",
                "g": f"{analysis['second_max_g']:.3f}",
            },
            "first_minimum": {
                "r": f"{analysis['first_min_r']:.3f}",
                "g": f"{analysis['first_min_g']:.3f}",
            },
        },
        "pmf": {
            "at_rdf_first_peak_kJ_per_mol": f"{analysis['pmf_at_first_peak']:.3f}",
            "at_rdf_first_minimum_kJ_per_mol": f"{analysis['pmf_at_first_min']:.3f}",
            "dissociation_barrier_kJ_per_mol": f"{analysis['delta_G']:.3f}",
        },
    }

    with output_json_file.open("w", encoding="utf-8") as f:
        json.dump(rdf_pmf_data_dict, f, indent=4)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="RDF + PMF analysis for GPUMD rdf.out (total RDF only)."
    )
    parser.add_argument(
        "input",
        type=str,
        help="Input RDF data file (e.g. rdf.out)",
    )
    parser.add_argument(
        "-T",
        "--temperature",
        type=float,
        required=True,
        help="Temperature in Kelvin (e.g. 2000)",
    )
    parser.add_argument(
        "-c",
        "--column",
        type=int,
        default=2,
        help=(
            "RDF column index to analyze (1=radius, 2=total, 3+=pair RDF). "
            "Default: 2 (total RDF)."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    rdf_pmf_analysis(Path(args.input), temperature_kelvin=args.temperature, column_index=args.column)


if __name__ == "__main__":
    main()


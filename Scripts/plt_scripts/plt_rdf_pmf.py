"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_rdf_pmf.py
Category:   Plot Scripts
Purpose:    RDF + PMF visualization for GPUMD rdf.out with optional
            temperature-dependent PMF calculation.
Usage:      gpumdkit.sh -plt rdf_pmf [temperature] [column_index] [save]
            python plt_rdf_pmf.py [temperature] [column_index] [save]
Arguments:
  temperature   Temperature in K for PMF calculation (optional)
  column_index  Column index for specific RDF pair (optional)
  save          Save the plot as PNG instead of displaying it
Output:
  RDF + PMF plot displayed or saved
Author:     Qilin Guo (guoqilin@buaa.edu.cn), Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-21
=============================================================================
"""

import json
import math
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import Axes
from scipy.signal import find_peaks, savgol_filter

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
RDF_COLOR = 'C0'
PMF_COLOR = 'C3'
REF_STYLE = dict(color='0.6', linestyle='--', linewidth=0.8, zorder=0)

# Use Agg backend in headless environments
if "DISPLAY" not in os.environ:
    plt.switch_backend("Agg")


def load_gpumd_rdf(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Load GPUMD rdf.out: r + total + optional pair RDFs"""
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    r = data[:, 0]
    g_all = data[:, 1:]
    labels = ["total"] + [f"pair_{i+1}" for i in range(g_all.shape[1] - 1)]
    return r, g_all, labels


def compute_pmf(g: np.ndarray, temperature: float) -> np.ndarray:
    """PMF(r) = -k_B T ln(g(r)) in kJ/mol"""
    kB = 8.617333262145e-5  # eV/K
    ev_to_kJ_mol = 96.48533212
    factor = kB * temperature * ev_to_kJ_mol
    g_safe = np.clip(g, 1e-10, None)
    return -factor * np.log(g_safe)


def analyze_peaks_min(r: np.ndarray, g: np.ndarray, pmf: np.ndarray) -> dict:
    """Find first/second peak and first minimum (r > 1.0 Å filter)"""
    win = min(9, len(g) if len(g) % 2 == 1 else len(g) - 1)
    g_smooth = savgol_filter(g, win, 3) if win >= 5 else g.copy()

    mask = r > 1.0
    peaks, _ = find_peaks(g_smooth[mask], distance=4, prominence=0.25)
    if len(peaks) < 2:
        raise RuntimeError("Could not find two clear peaks after r > 1 Å")

    idx = np.where(mask)[0][peaks]
    i1 = int(idx[0])
    i2 = int(idx[1])

    i_min = i1 + np.argmin(g_smooth[i1:i2])

    return {
        "r_peak1": float(r[i1]),
        "r_min":   float(r[i_min]),
        "r_peak2": float(r[i2]),
        "pmf_peak1": float(pmf[i1]),
        "pmf_min":   float(pmf[i_min]),
        "delta_G":   float(pmf[i_min] - pmf[i1]),
    }


def plot_rdf_pmf(
    r: np.ndarray,
    g: np.ndarray,
    pmf: np.ndarray,
    info: dict,
    temperature: float,
    save: bool = False
) -> None:
    """Generate clean RDF + PMF plot with adjustable margins"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.4, 4.6), sharex=True)
    marks = [info["r_peak1"], info["r_min"], info["r_peak2"]]

    # ─── RDF panel ───────────────────────────────────────
    ax1.plot(r, g, linewidth=1.4, color=RDF_COLOR)
    for pos in marks:
        ax1.axvline(pos, **REF_STYLE)
    ax1.axhline(1.0, **REF_STYLE)

    # Key positions in top-right corner (multi-line, no box)
    pos_text = (
        f"1st peak: {info['r_peak1']:.2f} Å\n"
        f"1st min : {info['r_min']:.2f} Å\n"
        f"2nd peak: {info['r_peak2']:.2f} Å"
    )
    ax1.text(0.98, 0.96, pos_text, transform=ax1.transAxes,
             ha="right", va="top", fontsize=9, color="0.2", linespacing=1.4)
    ax1.set_ylabel(r"$g(r)$")

    # ─── PMF panel ───────────────────────────────────────
    ax2.plot(r, pmf, linewidth=1.4, color=PMF_COLOR)
    for pos in marks:
        ax2.axvline(pos, **REF_STYLE)
    ax2.axhline(0, **REF_STYLE)

    # ΔG bracket
    ax2.vlines(info["r_min"], info["pmf_peak1"], info["pmf_min"], color="0.35", linewidth=1.2)
    ax2.hlines([info["pmf_peak1"], info["pmf_min"]],
               info["r_min"] - 0.20, info["r_min"] + 0.20, color="0.35", linewidth=1.2)

    # ΔG text in top-right (no box)
    ax2.text(0.98, 0.96,
             f"ΔG = {info['delta_G']:.2f} kJ/mol\n(T = {temperature:.0f} K)",
             transform=ax2.transAxes, ha="right", va="top",
             fontsize=9, color=PMF_COLOR, linespacing=1.4)
    ax2.set_ylabel("PMF (kJ/mol)")
    ax2.set_xlabel(r"$r$ ($\mathrm{\AA}$)")

    # View range for the PMF: from the first RDF peak outwards, so that the steep
    # repulsive wall (g(r) ~ 0 is clipped to a huge PMF) does not set the y scale
    view_mask = (r >= info["r_peak1"]) & (r < info["r_peak2"] + 5)
    y_vals = pmf[view_mask]
    y_min, y_max = np.min(y_vals), np.max(y_vals)
    margin = (y_max - y_min) * 0.12
    ax2.set_ylim(y_min - margin, y_max + margin * 1.5)
    ax2.set_xlim(r[0], r[-1])

    # Subplot geometry (figure fractions); tune these to adjust margins and spacing
    fig.subplots_adjust(left=0.16, right=0.97, top=0.97, bottom=0.11, hspace=0.08)

    if save:
        plt.savefig("rdf_pmf.png", dpi=300)
        print("Figure saved → rdf_pmf.png")
    else:
        plt.show()

    plt.close()


def main():
    if len(sys.argv) == 1:
        print("1. temperature is your simulation temperature in K\n2. [column_index] is the RDF column to analyze\n3. use 'save' to save the figure instead of showing it")
        print("----------------------------------------")
        print("Enter: temperature [column_index] [save]")
        line = input("> ").strip()
        parts = line.split()
        T = float(parts[0])
        col = int(parts[1]) if len(parts) > 1 else 2
        save_flag = len(parts) > 2 and parts[-1].lower().startswith(("s", "y"))
    else:
        T = float(sys.argv[1])
        col = int(sys.argv[2]) if len(sys.argv) > 2 else 2
        save_flag = len(sys.argv) > 3 and sys.argv[3].lower().startswith(("s", "y"))

    rdf_path = Path("rdf.out")
    if not rdf_path.is_file():
        print("Error: rdf.out not found.", file=sys.stderr)
        sys.exit(1)

    r, g_all, _ = load_gpumd_rdf(rdf_path)

    if col < 2 or col > g_all.shape[1] + 1:
        print(f"Error: column {col} invalid (valid: 2 to {g_all.shape[1]+1})", file=sys.stderr)
        sys.exit(1)

    g = g_all[:, col - 2]
    pmf = compute_pmf(g, T)
    info = analyze_peaks_min(r, g, pmf)

    plot_rdf_pmf(r, g, pmf, info, T, save_flag)


if __name__ == "__main__":
    main()
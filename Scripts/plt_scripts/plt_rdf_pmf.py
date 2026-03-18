"""
RDF + PMF visualization for GPUMD rdf.out

Author: Qilin Guo (guoqilin@buaa.edu.cn)
Modified by Zihan Yan (yanzihan@westlake.edu.cn)
Last modified: 2026-03-18

Usage:
    python plt_rdf_pmf.py [temperature] [column_index] [save]
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
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(7.5, 6.5), sharex=True, gridspec_kw={"hspace": 0.18}
    )

    # ─── RDF panel ───────────────────────────────────────
    ax1.plot(r, g, lw=2.0, color="#1f77b4")

    # Light gray dashed lines
    for pos, style in [
        (info["r_peak1"], (0, (5, 2))),
        (info["r_min"],   (0, (3, 1.5))),
        (info["r_peak2"], (0, (5, 2))),
    ]:
        ax1.axvline(pos, color="0.68", ls=style, lw=1.0, alpha=0.8)

    ax1.axhline(1.0, color="gray", ls="--", lw=0.9, alpha=0.5)

    # Key positions in top-right corner (multi-line, no box)
    pos_text = (
        f"1st peak: {info['r_peak1']:.2f} Å\n"
        f"1st min : {info['r_min']:.2f} Å\n"
        f"2nd peak: {info['r_peak2']:.2f} Å"
    )
    ax1.text(0.98, 0.96, pos_text,
             transform=ax1.transAxes,
             ha="right", va="top",
             fontsize=10.0,
             color="0.2")

    ax1.set_ylabel("g(r)", fontsize=11)
    ax1.set_xlabel("Radial distance r (Å)", fontsize=11)
    ax1.grid(True, alpha=0.16, which="both")
    ax1.minorticks_on()
    ax1.tick_params(which="both", direction="in", top=True, right=True, labelsize=9.5)

    # ─── PMF panel ───────────────────────────────────────
    ax2.plot(r, pmf, lw=2.0, color="#d62728")

    # Same dashed lines
    for pos, style in [
        (info["r_peak1"], (0, (5, 2))),
        (info["r_min"],   (0, (3, 1.5))),
        (info["r_peak2"], (0, (5, 2))),
    ]:
        ax2.axvline(pos, color="0.68", ls=style, lw=1.0, alpha=0.8)

    ax2.axhline(0, color="gray", ls="--", lw=0.9, alpha=0.5)

    # ΔG bracket
    ax2.vlines(info["r_min"], info["pmf_peak1"], info["pmf_min"],
               color="0.35", lw=1.4)
    ax2.hlines([info["pmf_peak1"], info["pmf_min"]],
               info["r_min"] - 0.20, info["r_min"] + 0.20,
               color="0.35", lw=1.4)

    # ΔG text in top-right (no box)
    ax2.text(0.98, 0.96,
             f"ΔG = {info['delta_G']:.2f} kJ/mol\n(T = {temperature:.0f} K)",
             transform=ax2.transAxes,
             ha="right", va="top",
             fontsize=10.5,
             color="#d62728")

    ax2.set_ylabel("PMF (kJ/mol)", fontsize=11)
    ax2.set_xlabel("Radial distance r (Å)", fontsize=11)
    ax2.grid(True, alpha=0.16, which="both")
    ax2.minorticks_on()
    ax2.tick_params(which="both", direction="in", top=True, right=True, labelsize=9.5)

    # Reasonable view range for PMF
    view_mask = r < info["r_peak2"] + 5
    y_vals = pmf[view_mask]
    y_min, y_max = np.min(y_vals), np.max(y_vals)
    margin = (y_max - y_min) * 0.12
    ax2.set_ylim(y_min - margin, y_max + margin * 1.5)

    # Adjust overall margins (this controls whitespace around the figure)
    plt.subplots_adjust(
        left=0.10,    # left margin
        right=0.96,   # right margin
        top=0.96,     # top margin
        bottom=0.10,  # bottom margin
        hspace=0.18   # height space between subplots
    )

    # plt.tight_layout(rect=[0, 0, 1, 0.98])  # slight extra safety

    if save:
        plt.savefig("rdf_pmf.png", dpi=400, bbox_inches="tight")
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
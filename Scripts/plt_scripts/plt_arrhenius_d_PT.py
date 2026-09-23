"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_arrhenius_d_PT.py
Category:   Plot Scripts
Purpose:    Compute diffusivity from temperature folders and fit separate
            Arrhenius lines below and above a user-specified transition
            temperature.
Usage:      gpumdkit.sh -plt D_PT <transition_temperature> [save]
            python3 plt_arrhenius_d_PT.py <transition_temperature> [save]
Arguments:
  transition_temperature  Transition temperature in K; a suffix 'K' is allowed
  save                    Save the plot as 'Arrhenius_D_PT.png'
Output:
  Arrhenius_D_PT.png      (if save is used, or if backend is non-interactive)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-08-28
=============================================================================
"""

import os
import re
import sys


args = sys.argv[1:]
if not args or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -plt D_PT <transition_temperature> [save]")
    print("    or: python3 plt_arrhenius_d_PT.py <transition_temperature> [save]")
    print("")
    print(" Arguments:")
    print("   transition_temperature  Transition temperature in K; '380K' is also accepted")
    print("   save                    Save the plot as 'Arrhenius_D_PT.png'")
    print("")
    print(" Example: gpumdkit.sh -plt D_PT 380 save")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

if len(args) > 2 or (len(args) == 2 and args[1] != "save"):
    print(" Error: expected <transition_temperature> followed by optional 'save'.")
    sys.exit(1)

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib import get_backend
from scipy.stats import linregress


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "font.size": 11,
    "axes.labelsize": 11.5,
    "axes.titlesize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.figsize": (4.3, 3.8),
    "axes.spines.top": True,
    "axes.spines.right": True,
    "axes.linewidth": 0.8,
})


BOLTZMANN_EV = 8.617333262145e-5
TARGET_FRACTION_START = 0.4
TARGET_FRACTION_END = 0.8


def parse_transition_temperature(value):
    """Parse a positive transition temperature expressed in K."""
    normalized = value.strip()
    if normalized.lower().endswith("k"):
        normalized = normalized[:-1]
    try:
        transition_temperature = float(normalized)
    except ValueError:
        print(f" Error: invalid transition temperature '{value}'.")
        sys.exit(1)
    if not np.isfinite(transition_temperature) or transition_temperature <= 0:
        print(" Error: transition temperature must be a positive finite value in K.")
        sys.exit(1)
    return transition_temperature


def fit_arrhenius(temperatures, values):
    """Fit log10(value) against 1000/T and return fit parameters."""
    valid = np.isfinite(temperatures) & np.isfinite(values) & (values > 0)
    if np.count_nonzero(valid) < 2:
        return None

    x_values = 1000.0 / temperatures[valid]
    y_values = np.log10(values[valid])
    result = linregress(x_values, y_values)
    activation_energy = -result.slope * BOLTZMANN_EV * 1000.0 * np.log(10.0)
    return {
        "activation_energy": activation_energy,
        "intercept": result.intercept,
        "slope": result.slope,
        "r_squared": result.rvalue ** 2,
        "count": int(np.count_nonzero(valid)),
    }


def read_msd_file(msd_file):
    """Read time and total MSD from a GPUMDkit four-column MSD file."""
    data = np.atleast_2d(np.loadtxt(msd_file))
    if data.ndim != 2 or data.shape[1] < 4:
        raise ValueError("msd.out must contain at least four columns")
    times = data[:, 0]
    total_msd = np.sum(data[:, 1:4], axis=1)
    return times, total_msd


def calculate_diffusivity(msd_file):
    """Calculate D using the established middle 40%-80% MSD interval."""
    times, total_msd = read_msd_file(msd_file)
    start = int(len(times) * TARGET_FRACTION_START)
    end = int(len(times) * TARGET_FRACTION_END)
    if end - start < 2:
        raise ValueError("msd.out does not contain enough rows for the 40%-80% fit")

    valid = np.isfinite(times[start:end]) & np.isfinite(total_msd[start:end])
    if np.count_nonzero(valid) < 2:
        raise ValueError("the 40%-80% MSD interval contains fewer than two finite rows")
    slope, _ = np.polyfit(times[start:end][valid], total_msd[start:end][valid], 1)
    return slope * 1e-4 / 6.0


def collect_diffusivities(base_dir):
    """Collect finite-temperature diffusivities from integer K directories."""
    temperatures = []
    diffusivities = []
    folders = sorted(
        (folder for folder in os.listdir(base_dir)
         if re.fullmatch(r"\d+K", folder)
         and os.path.isdir(os.path.join(base_dir, folder))),
        key=lambda folder: int(folder[:-1]),
    )

    for folder in folders:
        msd_file = os.path.join(base_dir, folder, "msd.out")
        if not os.path.isfile(msd_file):
            print(f" [Warning] Missing msd.out in {folder}, skipping")
            continue
        try:
            diffusivity = calculate_diffusivity(msd_file)
        except (OSError, ValueError) as error:
            print(f" [Warning] Failed to process {folder}: {error}")
            continue
        temperatures.append(int(folder[:-1]))
        diffusivities.append(diffusivity)
        if not np.isfinite(diffusivity) or diffusivity <= 0:
            print(f" [Warning] {folder}: invalid diffusivity {diffusivity:.3e}; excluded from fitting")

    if not temperatures:
        raise ValueError("no usable integer temperature folders with msd.out were found")
    order = np.argsort(temperatures)
    return np.asarray(temperatures, dtype=float)[order], np.asarray(diffusivities)[order]


def format_temperature(value):
    """Format a temperature without an unnecessary decimal suffix."""
    return f"{value:.0f}" if value.is_integer() else f"{value:g}"


def plot_branch(ax, temperatures, values, mask, label, color):
    """Plot one temperature branch and its fitted Arrhenius line."""
    valid = mask & np.isfinite(temperatures) & np.isfinite(values) & (values > 0)
    fit = fit_arrhenius(temperatures[mask], values[mask])
    if fit is None:
        return None

    x_data = 1000.0 / temperatures[valid]
    y_data = np.log10(values[valid])
    x_fit = np.linspace(np.min(x_data), np.max(x_data), 100)
    y_fit = fit["slope"] * x_fit + fit["intercept"]
    legend_label = f"{label} ({fit['activation_energy']:.3f} eV)"
    ax.plot(
        x_data, y_data, "o", color=color, markersize=8,
        markerfacecolor="none", markeredgewidth=1.5, label=legend_label,
    )
    ax.plot(x_fit, y_fit, "--", color=color, linewidth=1.5)
    return fit


def main():
    """Load data, fit both branches, print results, and draw the plot."""
    transition_temperature = parse_transition_temperature(args[0])
    base_dir = os.getcwd()
    try:
        temperatures, diffusivities = collect_diffusivities(base_dir)
    except (OSError, ValueError) as error:
        print(f" Error: {error}.")
        sys.exit(1)

    low_mask = temperatures <= transition_temperature
    high_mask = temperatures >= transition_temperature
    low_fit = fit_arrhenius(temperatures[low_mask], diffusivities[low_mask])
    high_fit = fit_arrhenius(temperatures[high_mask], diffusivities[high_mask])
    if low_fit is None or high_fit is None:
        print(
            f" Error: at least two positive finite temperatures are required in both "
            f"LowT (T <= {format_temperature(transition_temperature)} K) and "
            f"HighT (T >= {format_temperature(transition_temperature)} K) branches."
        )
        sys.exit(1)

    transition_label = format_temperature(transition_temperature)
    print(f"\n Transition temperature: {transition_label} K")
    print(
        f" LowT (T <= {transition_label} K): Ea = "
        f"{low_fit['activation_energy']:.3f} eV, R^2 = {low_fit['r_squared']:.4f}, "
        f"N = {low_fit['count']}"
    )
    print(
        f" HighT (T >= {transition_label} K): Ea = "
        f"{high_fit['activation_energy']:.3f} eV, R^2 = {high_fit['r_squared']:.4f}, "
        f"N = {high_fit['count']}"
    )

    width_temperature, width_diffusivity = 10, 15
    line = f"+{'-' * (width_temperature + 2)}-{'-' * (width_diffusivity + 2)}+"
    print(f" {line}")
    print(f" | {'Diffusivity (unit: cm2/s)':^{len(line) - 4}} |")
    print(f" {line}")
    print(f" | {'T (K)':^{width_temperature}} | {'D_total':^{width_diffusivity}} |")
    print(f" {line}")
    for temperature, diffusivity in zip(temperatures, diffusivities):
        print(f" | {temperature:^{width_temperature}.0f} | {diffusivity:^{width_diffusivity}.3e} |")
    print(f" {line}")

    fig, ax = plt.subplots(figsize=(4.3, 3.8), dpi=150)
    plot_branch(ax, temperatures, diffusivities, high_mask, "HighT", "#D62828")
    plot_branch(ax, temperatures, diffusivities, low_mask, "LowT", "#D62828")
    ax.set_xlabel("1000/T (1/K)", labelpad=7)
    ax.set_ylabel(r"log10(D) (cm$^2$/s)")
    ax.legend(loc="lower left", frameon=False, fontsize=11)

    ax_top = ax.secondary_xaxis("top")
    ax_top.set_xlabel("Temperature (K)", labelpad=7)
    xticks = ax.get_xticks()
    ax_top.set_xticks(xticks)
    ax_top.set_xticklabels([f"{int(1000 / tick)}" if tick > 0 else "" for tick in xticks])
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    plt.tight_layout()

    output_file = "Arrhenius_D_PT.png"
    if len(args) > 1 and args[-1] == "save":
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
    elif get_backend().lower() in ["agg", "cairo", "pdf", "ps", "svg"]:
        print(" Unable to display the plot due to the non-interactive backend.")
        print(f" The plot has been automatically saved as '{output_file}'.")
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    main()

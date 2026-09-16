"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_arrhenius_sigma_PT.py
Category:   Plot Scripts
Purpose:    Calculate ionic conductivity from temperature folders and fit
            separate Arrhenius lines below and above a user-specified
            transition temperature.
Usage:      gpumdkit.sh -plt sigma_PT <transition_temperature> [save]
            python3 plt_arrhenius_sigma_PT.py <transition_temperature> [save]
Arguments:
  transition_temperature  Transition temperature in K; a suffix 'K' is allowed
  save                    Save the plot as 'Arrhenius_sigma_PT.png'
Output:
  Arrhenius_sigma_PT.png  (if save is used, or if backend is non-interactive)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-08-28
=============================================================================
"""

import os
import re
import sys


args = sys.argv[1:]
if not args or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -plt sigma_PT <transition_temperature> [save]")
    print("    or: python3 plt_arrhenius_sigma_PT.py <transition_temperature> [save]")
    print("")
    print(" Arguments:")
    print("   transition_temperature  Transition temperature in K; '380K' is also accepted")
    print("   save                    Save the plot as 'Arrhenius_sigma_PT.png'")
    print("")
    print(" Example: gpumdkit.sh -plt sigma_PT 380 save")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

if len(args) > 2 or (len(args) == 2 and args[1] != "save"):
    print(" Error: expected <transition_temperature> followed by optional 'save'.")
    sys.exit(1)

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import scipy.constants as consts
from ase.io import read
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
ROOM_TEMPERATURE = 300.0


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


def fit_arrhenius(temperatures, sigma_times):
    """Fit log10(sigma*T) against 1000/T and return fit parameters."""
    valid = np.isfinite(temperatures) & np.isfinite(sigma_times) & (sigma_times > 0)
    if np.count_nonzero(valid) < 2:
        return None

    x_values = 1000.0 / temperatures[valid]
    y_values = np.log10(sigma_times[valid])
    result = linregress(x_values, y_values)
    activation_energy = -result.slope * BOLTZMANN_EV * 1000.0 * np.log(10.0)
    return {
        "activation_energy": activation_energy,
        "intercept": result.intercept,
        "slope": result.slope,
        "r_squared": result.rvalue ** 2,
        "count": int(np.count_nonzero(valid)),
    }


def read_volume(thermo_file):
    """Read and average the volume from a 12- or 18-column thermo.out."""
    data = np.atleast_2d(np.loadtxt(thermo_file))
    if data.ndim != 2:
        raise ValueError("thermo.out is not a two-dimensional numeric table")
    if data.shape[1] == 12:
        volumes = data[:, 9] * data[:, 10] * data[:, 11]
    elif data.shape[1] == 18:
        vector_a = data[:, 9:12]
        vector_b = data[:, 12:15]
        vector_c = data[:, 15:18]
        volumes = np.abs(np.einsum("ij,ij->i", vector_a, np.cross(vector_b, vector_c)))
    else:
        raise ValueError("thermo.out must contain 12 or 18 columns")
    if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0):
        raise ValueError("thermo.out contains invalid volume values")
    return float(np.mean(volumes))


def read_msd_file(msd_file):
    """Read time and total MSD from a GPUMDkit four-column MSD file."""
    data = np.atleast_2d(np.loadtxt(msd_file))
    if data.ndim != 2 or data.shape[1] < 4:
        raise ValueError("msd.out must contain at least four columns")
    return data[:, 0], np.sum(data[:, 1:4], axis=1)


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


def read_ion_count(first_folder):
    """Count Li/Na ions and apply the first folder's GPUMD replicate factor."""
    xyz_file = os.path.join(first_folder, "model.xyz")
    if not os.path.isfile(xyz_file):
        raise FileNotFoundError(f"missing {xyz_file}")
    atoms = read(xyz_file, format="extxyz")
    num_ions = sum(1 for atom in atoms if atom.symbol in ["Li", "Na"])
    if num_ions <= 0:
        raise ValueError("model.xyz contains no Li or Na ions")

    run_in_path = os.path.join(first_folder, "run.in")
    if not os.path.isfile(run_in_path):
        print(" [Note] No run.in file found, assuming no replication")
        return num_ions

    with open(run_in_path, "r", encoding="utf-8") as run_file:
        for line in run_file:
            if line.strip().startswith("replicate"):
                parts = line.split()
                if len(parts) < 4:
                    print(" [Warning] Invalid replicate format, assuming no replication")
                    return num_ions
                try:
                    replicate_x, replicate_y, replicate_z = (int(parts[index]) for index in (1, 2, 3))
                except ValueError:
                    print(" [Warning] Invalid replicate format, assuming no replication")
                    return num_ions
                factor = replicate_x * replicate_y * replicate_z
                if factor <= 0:
                    print(" [Warning] Invalid replicate factor, assuming no replication")
                    return num_ions
                print(
                    f" [Info] Detected replicate {replicate_x} {replicate_y} {replicate_z}, "
                    f"multiplying ion count by {factor}"
                )
                return num_ions * factor
    return num_ions


def collect_conductivities(base_dir, num_ions):
    """Collect conductivity values from integer K directories."""
    temperatures = []
    conductivities = []
    sigma_times = []
    folders = sorted(
        (folder for folder in os.listdir(base_dir)
         if re.fullmatch(r"\d+K", folder)
         and os.path.isdir(os.path.join(base_dir, folder))),
        key=lambda folder: int(folder[:-1]),
    )

    for folder in folders:
        folder_path = os.path.join(base_dir, folder)
        thermo_file = os.path.join(folder_path, "thermo.out")
        msd_file = os.path.join(folder_path, "msd.out")
        if not os.path.isfile(thermo_file) or not os.path.isfile(msd_file):
            print(f" [Warning] Missing thermo.out or msd.out in {folder}, skipping")
            continue
        try:
            temperature = float(folder[:-1])
            volume = read_volume(thermo_file)
            diffusivity = calculate_diffusivity(msd_file)
            volume_cm3 = volume * 1e-24
            conversion_factor = (
                num_ions / (volume_cm3 * consts.N_A)
                * (consts.N_A * consts.e) ** 2
                / (consts.R * temperature)
            )
            conductivity = conversion_factor * diffusivity
            sigma_time = conductivity * temperature
            if not np.isfinite(sigma_time) or sigma_time <= 0:
                print(f" [Warning] {folder}: invalid sigma*T = {sigma_time}; skipping")
                continue
            temperatures.append(temperature)
            conductivities.append(conductivity)
            sigma_times.append(sigma_time)
        except (OSError, ValueError, IndexError) as error:
            print(f" [Warning] Failed to process {folder}: {error}; skipping")

    if not temperatures:
        raise ValueError("no usable temperature folders with thermo.out and msd.out were found")
    order = np.argsort(temperatures)
    return (
        np.asarray(temperatures, dtype=float)[order],
        np.asarray(conductivities, dtype=float)[order],
        np.asarray(sigma_times, dtype=float)[order],
    )


def format_temperature(value):
    """Format a temperature without an unnecessary decimal suffix."""
    return f"{value:.0f}" if value.is_integer() else f"{value:g}"


def plot_branch(ax, temperatures, sigma_times, mask, label, color):
    """Plot one conductivity branch and its fitted Arrhenius line."""
    valid = mask & np.isfinite(temperatures) & np.isfinite(sigma_times) & (sigma_times > 0)
    fit = fit_arrhenius(temperatures[mask], sigma_times[mask])
    if fit is None:
        return None

    x_data = 1000.0 / temperatures[valid]
    y_data = np.log(sigma_times[valid])
    x_fit = np.linspace(np.min(x_data), np.max(x_data), 100)
    y_fit = (fit["slope"] * x_fit + fit["intercept"]) * np.log(10.0)
    legend_label = f"{label} ({fit['activation_energy']:.3f} eV)"
    ax.plot(
        x_data, y_data, "o", color=color, markersize=8,
        markerfacecolor="none", markeredgewidth=1.5, label=legend_label,
    )
    ax.plot(x_fit, y_fit, "--", color=color, linewidth=1.5)
    return fit


def main():
    """Load data, fit both branches, extrapolate sigma at 300 K, and plot."""
    transition_temperature = parse_transition_temperature(args[0])
    base_dir = os.getcwd()
    folders = sorted(
        (folder for folder in os.listdir(base_dir)
         if re.fullmatch(r"\d+K", folder)
         and os.path.isdir(os.path.join(base_dir, folder))),
        key=lambda folder: int(folder[:-1]),
    )
    if not folders:
        print(" Error: no valid temperature folders (for example, 400K or 500K) found.")
        sys.exit(1)

    try:
        num_ions = read_ion_count(os.path.join(base_dir, folders[0]))
        temperatures, conductivities, sigma_times = collect_conductivities(base_dir, num_ions)
    except (OSError, ValueError) as error:
        print(f" Error: {error}.")
        sys.exit(1)

    low_mask = temperatures <= transition_temperature
    high_mask = temperatures >= transition_temperature
    low_fit = fit_arrhenius(temperatures[low_mask], sigma_times[low_mask])
    high_fit = fit_arrhenius(temperatures[high_mask], sigma_times[high_mask])
    transition_label = format_temperature(transition_temperature)
    if low_fit is None or high_fit is None:
        print(
            f" Error: at least two positive finite temperatures are required in both "
            f"LowT (T <= {transition_label} K) and HighT (T >= {transition_label} K) branches."
        )
        sys.exit(1)

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

    width_temperature, width_sigma, width_sigma_time = 10, 16, 18
    line = (
        f"+{'-' * (width_temperature + 2)}-"
        f"{'-' * (width_sigma + 2)}-{'-' * (width_sigma_time + 2)}+"
    )
    print(f" {line}")
    print(f" | {'Conductivity (unit: S/cm, Sigma*T: K S/cm)':^{len(line) - 4}} |")
    print(f" {line}")
    print(
        f" | {'T (K)':^{width_temperature}} | {'Sigma':^{width_sigma}} | "
        f"{'Sigma*T':^{width_sigma_time}} |"
    )
    print(f" {line}")
    for temperature, conductivity, sigma_time in zip(
        temperatures, conductivities, sigma_times
    ):
        print(
            f" | {temperature:^{width_temperature}.0f} | "
            f"{conductivity:^{width_sigma}.3e} | "
            f"{sigma_time:^{width_sigma_time}.3e} |"
        )
    print(f" {line}")

    fig, ax = plt.subplots(figsize=(4.3, 3.8), dpi=150)
    plot_branch(ax, temperatures, sigma_times, high_mask, "HighT", "#457B9D")
    plot_branch(ax, temperatures, sigma_times, low_mask, "LowT", "#457B9D")
    ax.set_xlabel("1000/T (1/K)", labelpad=7)
    ax.set_ylabel(r"ln($\sigma$T) (S$\cdot$K/cm)")
    ax.legend(loc="lower left", frameon=False, fontsize=11)

    ax_top = ax.secondary_xaxis("top")
    ax_top.set_xlabel("Temperature (K)", labelpad=7)
    xticks = ax.get_xticks()
    ax_top.set_xticks(xticks)
    ax_top.set_xticklabels([f"{int(1000 / tick)}" if tick > 0 else "" for tick in xticks])
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    plt.tight_layout()

    selected_label = "LowT" if ROOM_TEMPERATURE <= transition_temperature else "HighT"
    selected_fit = low_fit if selected_label == "LowT" else high_fit
    sigma_300 = (
        10 ** (selected_fit["slope"] * (1000.0 / ROOM_TEMPERATURE) + selected_fit["intercept"])
        / ROOM_TEMPERATURE
    )
    if ROOM_TEMPERATURE == transition_temperature:
        print(" [Info] 300 K equals the transition temperature; using the LowT branch by default.")
    print(
        f"\n Room-temperature extrapolation uses {selected_label}: "
        f"Sigma(300 K) = {sigma_300:.3e} S/cm = {sigma_300 * 1000:.3f} mS/cm"
    )

    cell = os.path.basename(base_dir)
    temp_array = temperatures.astype(int).tolist()
    conductivity_array = [float(f"{value:.3e}") for value in conductivities]
    print("\n Exportable arrays for Python:")
    print(f" temperatures = {temp_array}")
    print(f" conductivity_values = {conductivity_array}")

    output_file = "Arrhenius_sigma_PT.png"
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

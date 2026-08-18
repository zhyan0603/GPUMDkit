"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_emd2.py
Category:   Plot Scripts
Purpose:    Plot EMD HAC and thermal conductivity convergence in all three
            directions.
Usage:      gpumdkit.sh -plt emd2 [save]
            python plt_emd2.py [save]
Arguments:
  save       Save the plot as 'emd2.png' instead of displaying it
Output:
  emd2.png   (if save is used)
Author:     Xin Wu (xinwuchn97@gmial.com)
Last-modified: 2026-08-18
=============================================================================
"""

import os
import sys


args = sys.argv[1:]


def print_usage():
    """Print usage instructions."""
    print("Usage: gpumdkit.sh -plt emd2 [save]")
    print("   or: python3 plt_emd2.py [save]")
    print("")
    print("Arguments:")
    print("  save     Optional, save the plot as 'emd2.png'")
    print("")
    print("Example: gpumdkit.sh -plt emd2 save")


if args and args[0] in ("-h", "--help", "help"):
    print_usage()
    sys.exit(0)

if any(argument != "save" for argument in args):
    print_usage()
    sys.exit(1)

# Matplotlib and NumPy are imported after help handling so that -h remains
# available when plotting dependencies are not installed.
import matplotlib.pyplot as plt
import numpy as np


aw, fs = 1.2, 12
plt.rc("font", size=fs)
plt.rc("font", family="sans-serif")
plt.rc("font", **{"sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"]})
plt.rc("axes", linewidth=aw)


def set_fig_properties(ax):
    """Apply the standard GPUMDkit transport-plot axis styling."""
    ax.tick_params(which="both", length=4, width=1.2, direction="in",
                   right=True, top=True)
    ax.tick_params(which="minor", length=4)


class EMD2Processor:
    """Read HAC output and plot total thermal conductivity for x, y, and z."""

    COLUMNS = [
        "time", "jx_in", "jx_out", "jy_in", "jy_out", "jz_tot",
        "kx_in", "kx_out", "ky_in", "ky_out", "kz_tot",
    ]
    DIRECTIONS = ("x", "y", "z")

    def __init__(self, directory):
        """Initialize the processor for a GPUMD output directory."""
        self.directory = directory
        self.run_path = os.path.join(directory, "run.in")
        self.hac_path = os.path.join(directory, "hac.out")

    def _read_run_parameters(self):
        """Read the time step and HAC sampling parameters from run.in."""
        time_step = None
        output_interval = None
        run_steps = None
        sample_interval = None

        with open(self.run_path, "r", encoding="utf-8") as file:
            for line in file:
                tokens = line.split()
                if not tokens or tokens[0].startswith("#"):
                    continue
                if tokens[0] == "time_step" and len(tokens) >= 2:
                    time_step = float(tokens[1])
                elif tokens[0] == "compute_hac" and len(tokens) >= 4:
                    output_interval = int(tokens[1])
                    run_steps = int(tokens[2])
                    sample_interval = int(tokens[3])

        if None in (time_step, output_interval, run_steps, sample_interval):
            raise ValueError(
                "run.in must contain time_step and compute_hac commands."
            )
        if run_steps <= 0 or sample_interval <= 0:
            raise ValueError("run.in contains invalid compute_hac intervals.")
        if run_steps % sample_interval != 0:
            raise ValueError(
                "The compute_hac run length must be divisible by its sample interval."
            )

        return time_step, output_interval, run_steps, run_steps // sample_interval

    def _load_data(self):
        """Load and reshape HAC output into one column per repeat."""
        if not os.path.isfile(self.run_path):
            raise FileNotFoundError(f"Input file '{self.run_path}' does not exist.")
        if not os.path.isfile(self.hac_path):
            raise FileNotFoundError(f"Input file '{self.hac_path}' does not exist.")

        time_step, output_interval, run_steps, n_hac_data = self._read_run_parameters()
        raw_data = np.loadtxt(self.hac_path, ndmin=2)
        if raw_data.ndim != 2 or raw_data.shape[1] != len(self.COLUMNS):
            raise ValueError(
                f"'{self.hac_path}' must contain {len(self.COLUMNS)} columns."
            )
        if not np.isfinite(raw_data).all():
            raise ValueError(f"'{self.hac_path}' contains NaN or Inf values.")
        if len(raw_data) == 0 or len(raw_data) % n_hac_data != 0:
            raise ValueError(
                "The MD calculation seems to be incomplete; "
                "check that hac.out contains complete repeats."
            )

        n_repeat = len(raw_data) // n_hac_data
        data = {}
        for index, column in enumerate(self.COLUMNS):
            repeats = raw_data[:, index].reshape(n_hac_data, n_repeat, order="F")
            data[column] = np.column_stack((repeats, repeats.mean(axis=1)))

        data["kx_tot"] = data["kx_in"] + data["kx_out"]
        data["ky_tot"] = data["ky_in"] + data["ky_out"]
        data["Results"] = self._compute_results(data, n_repeat)
        time_upper = output_interval * run_steps * time_step * 1e-6
        return data, time_upper, n_repeat

    @staticmethod
    def _compute_results(data, n_repeat):
        """Calculate the mean and standard error using the original EMD rule."""
        results = {}
        for direction in EMD2Processor.DIRECTIONS:
            key = f"k{direction}_tot"
            values = data[key][len(data[key]) // 2:, -1]
            results[f"{key}_ave"] = values.mean()
            results[f"{key}_std"] = values.std() / np.sqrt(n_repeat)
        return results

    @staticmethod
    def _print_results(results):
        """Print all directional thermal conductivity results."""
        print("\n" + "=" * 60)
        print("EMD Thermal Conductivity Results")
        print("=" * 60)
        for direction in EMD2Processor.DIRECTIONS:
            key = f"k{direction}_tot"
            print(
                f"Direction {direction}: "
                f"kappa = {results[f'{key}_ave']:.4f} +/- "
                f"{results[f'{key}_std']:.4f} W/mK"
            )
        print("=" * 60 + "\n")

    @staticmethod
    def _plot_results(data, time_upper, n_repeat, save_plot):
        """Plot HAC and total thermal conductivity for all directions."""
        time_data = data["time"] * 1e-3  # ns
        results = data["Results"]
        colors = {"x": "C0", "y": "C1", "z": "C2"}
        hac_columns = {
            "x": (("jx_in", "HAC in", "C0"), ("jx_out", "HAC out", "C1")),
            "y": (("jy_in", "HAC in", "C0"), ("jy_out", "HAC out", "C1")),
            "z": (("jz_tot", "HAC total", "C2"),),
        }
        figure, axes = plt.subplots(2, 3, figsize=(12, 8))

        for column, direction in enumerate(EMD2Processor.DIRECTIONS):
            hac_axis, conductivity_axis = axes[:, column]
            set_fig_properties(hac_axis)
            for key, label, color in hac_columns[direction]:
                for repeat in range(n_repeat):
                    values = data[key][:, repeat]
                    scale = values.max()
                    if scale == 0:
                        raise ValueError(f"HAC column '{key}' is identically zero.")
                    hac_axis.loglog(
                        time_data[:, repeat], values / scale,
                        color="k", alpha=0.25, linewidth=0.8,
                    )
                mean_values = data[key][:, -1]
                scale = mean_values.max()
                if scale == 0:
                    raise ValueError(f"HAC column '{key}' is identically zero.")
                hac_axis.loglog(
                    time_data[:, -1], mean_values / scale,
                    color=color, linewidth=2.0, label=label,
                )
            hac_axis.set_xlim(1e-5, time_upper)
            hac_axis.set_ylabel("Normalized HAC")
            hac_axis.set_title(f"{direction}: heat-current correlation")
            hac_axis.legend(frameon=False, fontsize=9)

            axis = conductivity_axis
            set_fig_properties(axis)
            key = f"k{direction}_tot"
            color = colors[direction]
            for repeat in range(n_repeat):
                axis.plot(time_data[:, repeat], data[key][:, repeat],
                          color="k", alpha=0.25, linewidth=0.8)
            axis.plot(time_data[:, -1], data[key][:, -1], color=color, linewidth=2.5,
                      label="Running average")
            average = results[f"{key}_ave"]
            standard_error = results[f"{key}_std"]
            axis.axhline(y=average, color="0.25", linestyle="--", linewidth=1.0,
                         label="Half-window average")
            axis.fill_between(
                time_data[:, -1], average - standard_error,
                average + standard_error, color=color, alpha=0.18,
                label="Standard error",
            )
            axis.set_xlim(0, time_upper)
            axis.set_xlabel("Correlation Time (ns)")
            axis.set_ylabel(r"$\kappa_{%s}$ (W/mK)" % direction)
            axis.set_title(
                rf"{direction}: $\kappa$ = {average:.2f} +/- {standard_error:.2f} W/mK"
            )
            axis.legend(frameon=False, fontsize=9)

        figure.tight_layout()
        if save_plot:
            figure.savefig("emd2.png", dpi=300, bbox_inches="tight")
            print("Saved plot to emd2.png")
        else:
            plt.show()

    def process(self, save_plot=False):
        """Process EMD output, print results, and display or save the plot."""
        data, time_upper, n_repeat = self._load_data()
        self._print_results(data["Results"])
        self._plot_results(data, time_upper, n_repeat, save_plot)


if __name__ == "__main__":
    processor = EMD2Processor(os.getcwd())
    try:
        processor.process(save_plot="save" in args)
    except (FileNotFoundError, OSError, ValueError) as error:
        print(f" Error: {error}")
        sys.exit(1)

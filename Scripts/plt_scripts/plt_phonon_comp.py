#!/usr/bin/env python3
"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_phonon_comp.py
Category:   Plot Scripts
Purpose:    Compare two or more phonon band-structure data files.
Usage:      gpumdkit.sh -plt phonon_comp <file1> <file2> [file3 ...] [save]
            python plt_phonon_comp.py <file1> <file2> [file3 ...] [save]
            ... [--qpoints FILE]
Arguments:
  file1...   Two or more phonon data files, normally named phonon_<label>.dat
  save       Save the figure as 'phonon_comp.png' instead of displaying it
Options:
  --qpoints FILE  Use FILE instead of the default 'QPOINTS' path definition.
Output:
  phonon_comp.png  Comparison figure when saving is requested or the selected
                   Matplotlib backend cannot display figures.
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-08-01
=============================================================================
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import List, Optional, Tuple

from plt_phonon import (
    _tick_positions,
    parse_phonon_data,
    parse_qpoints,
    set_qpoint_ticks,
)


def print_help() -> None:
    """Print command-line usage without importing plotting dependencies."""
    print(" Usage: gpumdkit.sh -plt phonon_comp <file1> <file2> [file3 ...] [save]")
    print("    or: python plt_phonon_comp.py <file1> <file2> [file3 ...] [save]")
    print("")
    print(" file1...       Two or more phonon data files")
    print(" save           Save the figure as 'phonon_comp.png'")
    print(" --qpoints FILE Use FILE instead of the default 'QPOINTS'")
    print("")
    print(" Labels are read from filenames: phonon_NEP.dat -> NEP")


def parse_arguments() -> Optional[Tuple[List[Path], Path, bool]]:
    """Parse two or more phonon files and the optional QPOINTS/save values."""
    arguments = sys.argv[1:]
    if arguments and arguments[0] in ("-h", "--help"):
        print_help()
        return None

    save = False
    if arguments and arguments[-1].lower() == "save":
        save = True
        arguments = arguments[:-1]

    files: List[Path] = []
    qpoints_file = Path("QPOINTS")
    index = 0
    while index < len(arguments):
        argument = arguments[index]
        if argument in ("--qpoints", "-q"):
            if index + 1 >= len(arguments):
                raise ValueError("--qpoints requires a file path")
            qpoints_file = Path(arguments[index + 1])
            index += 2
            continue
        if argument.startswith("-"):
            raise ValueError(f"unknown option '{argument}'")
        files.append(Path(argument))
        index += 1

    if len(files) < 2:
        raise ValueError("phonon_comp requires at least two phonon data files")
    return files, qpoints_file, save


def label_from_filename(filename: Path) -> str:
    """Derive a legend label from a conventional phonon_<label>.dat name."""
    stem = filename.stem
    if stem.lower().startswith("phonon_"):
        stem = stem[len("phonon_"):]
    return stem or filename.name


def plot_comparison(
    q_lengths,
    datasets,
    npoints: int,
    tick_positions,
    labels: List[str],
    save: bool,
) -> None:
    """Display or save a multi-model phonon comparison plot."""
    import matplotlib.pyplot as plt
    from matplotlib import get_backend

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "font.size": 12,
        "axes.labelsize": 14,
    })
    figure, axis = plt.subplots(figsize=(5.5, 4.0), dpi=150)
    palette = ["cornflowerblue", "firebrick", "goldenrod", "0.6"]
    linestyles = ["-", "--", "-.", ":"]
    legend_handles = []
    legend_labels = []
    two_file_comparison = len(datasets) == 2
    for model_index, (filename, frequencies) in enumerate(datasets):
        label = label_from_filename(filename)
        upper_label = label.upper()
        if two_file_comparison and upper_label == "DFT":
            color, linewidth, linestyle = "0.6", 1.8, "-"
        elif two_file_comparison and upper_label == "NEP":
            color, linewidth, linestyle = "firebrick", 1.5, "--"
        else:
            color = palette[model_index % len(palette)]
            linewidth = 1.4 if two_file_comparison else 1.2
            linestyle = linestyles[model_index % len(linestyles)]
        if two_file_comparison and upper_label not in {"DFT", "NEP"}:
            color = palette[model_index]
        representative_line = None
        for segment_index in range(len(q_lengths) // npoints):
            start_index = segment_index * npoints
            stop_index = (segment_index + 1) * npoints
            lines = axis.plot(
                q_lengths[start_index:stop_index],
                frequencies[start_index:stop_index],
                color=color,
                linewidth=linewidth,
                linestyle=linestyle,
                label="_nolegend_",
                alpha=0.8,
            )
            if representative_line is None:
                representative_line = lines[0]
        legend_handles.append(representative_line)
        legend_labels.append(label)
        model_index += 1

    for position in tick_positions[1:-1]:
        axis.axvline(position, color="0.8", linewidth=0.8, zorder=0)
    axis.axhline(0.0, color="gray", linewidth=1.0, linestyle=":", zorder=0)
    axis.set_xlim(q_lengths[0], q_lengths[-1])
    all_frequencies = [frequencies for _, frequencies in datasets]
    minimum = min(float(frequencies.min()) for frequencies in all_frequencies)
    maximum = max(float(frequencies.max()) for frequencies in all_frequencies)
    margin = max((maximum - minimum) * 0.05, 0.1)
    axis.set_ylim(minimum - margin, maximum + margin)
    set_qpoint_ticks(axis, tick_positions, labels)
    axis.set_xlabel("Wave vector path")
    axis.set_ylabel("Frequency (THz)")
    axis.legend(
        legend_handles,
        legend_labels,
        loc="upper right",
        ncol=min(3, len(datasets)),
        frameon=False,
        fontsize="small",
    )

    conversion = 33.35641
    right_axis = axis.twinx()
    lower, upper = axis.get_ylim()
    right_axis.set_ylim(lower * conversion, upper * conversion)
    right_axis.set_ylabel(r"Frequency (cm$^{-1}$)")
    axis.tick_params(axis="x", top=True)
    axis.tick_params(axis="y", right=False)
    right_axis.tick_params(top=True, left=False)
    figure.tight_layout()

    if save or get_backend().lower() in {"agg", "cairo", "pdf", "ps", "svg"}:
        figure.savefig("phonon_comp.png", dpi=300, bbox_inches="tight")
        print(" Saved phonon comparison plot to 'phonon_comp.png'.")
    else:
        plt.show()
    plt.close(figure)


def main() -> int:
    """Read all requested files and plot their phonon bands together."""
    try:
        parsed = parse_arguments()
        if parsed is None:
            return 0
        files, qpoints_file, save = parsed
        npoints, endpoint_labels = parse_qpoints(qpoints_file)
        datasets = []
        reference_q_lengths = None
        reference_rows = None
        for filename in files:
            q_lengths, frequencies = parse_phonon_data(filename)
            if reference_q_lengths is None:
                reference_q_lengths = q_lengths
                reference_rows = len(q_lengths)
            elif len(q_lengths) != reference_rows:
                raise ValueError(
                    f"'{filename}' has {len(q_lengths)} q-points; all comparison files "
                    f"must have {reference_rows}"
                )
            datasets.append((filename, frequencies))

        plot_q_lengths, tick_positions, path_labels = _tick_positions(
            reference_q_lengths, npoints, endpoint_labels
        )
        plot_comparison(
            plot_q_lengths,
            datasets,
            npoints,
            tick_positions,
            path_labels,
            save,
        )
    except (OSError, ValueError) as error:
        print(f" Error: {error}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

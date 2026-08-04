#!/usr/bin/env python3
"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_phonon.py
Category:   Plot Scripts
Purpose:    Plot a phonon band structure from calc_phonon.py output.
Usage:      gpumdkit.sh -plt phonon [phonon_file] [QPOINTS] [save]
            python plt_phonon.py [phonon_file] [QPOINTS] [save]
Arguments:
  phonon_file  Phonon data file; defaults to 'phonon_NEP.dat'
  QPOINTS      Line-mode q-point path; defaults to 'QPOINTS'
  save         Save the figure as 'phonon.png' instead of displaying it
Output:
  phonon.png   Phonon band-structure figure when saving is requested or the
               selected Matplotlib backend cannot display figures.
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-08-01
=============================================================================
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple


def print_help() -> None:
    """Print command-line usage without importing plotting dependencies."""
    print(" Usage: gpumdkit.sh -plt phonon [phonon_file] [QPOINTS] [save]")
    print("    or: python plt_phonon.py [phonon_file] [QPOINTS] [save]")
    print("")
    print(" phonon_file  Phonon data file; defaults to 'phonon_NEP.dat'")
    print(" QPOINTS      Line-mode q-point path; defaults to 'QPOINTS'")
    print(" save         Save the figure as 'phonon.png'")


def parse_arguments() -> Optional[Tuple[Path, Path, bool]]:
    """Parse the optional phonon data, QPOINTS, and save arguments."""
    arguments = sys.argv[1:]
    if arguments and arguments[0] in ("-h", "--help"):
        print_help()
        return None

    if len(arguments) > 3:
        raise ValueError("expected at most a phonon file, a QPOINTS file, and 'save'")

    save = False
    if arguments and arguments[-1].lower() == "save":
        save = True
        arguments = arguments[:-1]
    if arguments and arguments[0].lower() == "save":
        raise ValueError("'save' must be the last argument")

    phonon_file = Path(arguments[0]) if arguments else Path("phonon_NEP.dat")
    qpoints_file = Path(arguments[1]) if len(arguments) == 2 else Path("QPOINTS")
    return phonon_file, qpoints_file, save


def _is_integer_token(token: str) -> bool:
    """Return whether a token is an integer q-point index."""
    return re.fullmatch(r"[+-]?\d+", token) is not None


def format_qpoint_label(label: str) -> str:
    """Format Gamma and underscore labels for compact scientific tick labels."""
    if label.upper() == "GAMMA":
        return r"$\Gamma$"
    if "_" in label:
        base, subscript = label.split("_", 1)
        return f"{base}$_{{{subscript}}}$"
    return label


def parse_phonon_data(filename: Path):
    """Read q-path distances and rectangular frequency data from a phonon file."""
    import numpy as np

    if not filename.is_file():
        raise OSError(f"phonon data file not found: '{filename}'")
    try:
        lines = filename.read_text(encoding="utf-8").splitlines()
    except UnicodeError as error:
        raise ValueError(f"phonon data file '{filename}' is not valid UTF-8 text") from error

    records: List[Tuple[float, List[float]]] = []
    current: Optional[Tuple[float, List[float]]] = None
    for line in lines:
        fields = line.split()
        if not fields:
            continue
        if len(fields) >= 6 and _is_integer_token(fields[0]):
            if current is not None:
                records.append(current)
            try:
                current = (float(fields[4]), [float(value) for value in fields[5:]])
            except ValueError as error:
                raise ValueError(f"invalid numeric data in '{filename}'") from error
        elif current is not None:
            try:
                current[1].extend(float(value) for value in fields)
            except ValueError as error:
                raise ValueError(f"invalid frequency data in '{filename}'") from error
    if current is not None:
        records.append(current)

    if not records:
        raise ValueError(f"phonon data file '{filename}' contains no q-point rows")
    band_counts = {len(frequencies) for _, frequencies in records}
    if len(band_counts) != 1 or not next(iter(band_counts)):
        raise ValueError(f"phonon data file '{filename}' has inconsistent band counts")

    q_lengths = np.asarray([distance for distance, _ in records], dtype=float)
    frequencies = np.asarray([frequencies for _, frequencies in records], dtype=float)
    if not np.all(np.isfinite(q_lengths)) or not np.all(np.isfinite(frequencies)):
        raise ValueError(f"phonon data file '{filename}' contains NaN or Inf values")
    return q_lengths, frequencies


def parse_qpoints(filename: Path) -> Tuple[int, List[str]]:
    """Read the points per segment and endpoint labels from QPOINTS."""
    if not filename.is_file():
        raise OSError(f"QPOINTS file not found: '{filename}'")
    try:
        lines = filename.read_text(encoding="utf-8").splitlines()
    except UnicodeError as error:
        raise ValueError(f"QPOINTS file '{filename}' is not valid UTF-8 text") from error
    if len(lines) < 5:
        raise ValueError(f"QPOINTS file '{filename}' is too short")
    try:
        npoints = int(lines[1].split()[0])
    except (IndexError, ValueError) as error:
        raise ValueError("the second QPOINTS line must contain the points per segment") from error
    if npoints <= 0:
        raise ValueError("QPOINTS points per segment must be positive")

    labels_raw: List[str] = []
    for line in lines[4:]:
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            [float(fields[0]), float(fields[1]), float(fields[2])]
        except ValueError:
            continue
        labels_raw.append(fields[-1])
    if not labels_raw or len(labels_raw) % 2:
        raise ValueError("QPOINTS line-mode endpoints must occur in start/end pairs")

    labels = [format_qpoint_label(label) for label in labels_raw]
    return npoints, labels


def _tick_positions(q_lengths, npoints: int, endpoint_labels: List[str]):
    """Normalize disconnected segments and return q-path ticks and labels."""
    import numpy as np

    if len(endpoint_labels) % 2:
        raise ValueError("QPOINTS endpoint labels must occur in start/end pairs")
    segment_count = len(endpoint_labels) // 2
    expected_rows = npoints * segment_count
    if len(q_lengths) != expected_rows:
        raise ValueError(
            f"phonon data contains {len(q_lengths)} q-points, but QPOINTS defines "
            f"{expected_rows} points"
        )
    normalized_lengths = np.asarray(q_lengths, dtype=float).copy()
    indices = []
    labels = []
    previous_end = None
    for segment_index in range(segment_count):
        start_index = segment_index * npoints
        end_index = (segment_index + 1) * npoints - 1
        start_label = endpoint_labels[2 * segment_index]
        end_label = endpoint_labels[2 * segment_index + 1]
        if segment_index == 0 or start_label != previous_end:
            if segment_index > 0:
                normalized_lengths[start_index:end_index + 1] += (
                    normalized_lengths[start_index - 1] - normalized_lengths[start_index]
                )
                labels[-1] = f"{labels[-1]}|{start_label}"
            else:
                indices.append(start_index)
                labels.append(start_label)
        indices.append(end_index)
        labels.append(end_label)
        previous_end = end_label
    return normalized_lengths, normalized_lengths[indices], labels


def set_qpoint_ticks(axis, tick_positions, labels: List[str]) -> None:
    """Set same-height q-point labels and separate crowded neighbors laterally."""
    axis.set_xticks(tick_positions)
    tick_texts = axis.set_xticklabels(labels)
    # Keep all labels at one height. When two neighboring labels collide,
    # anchor them outwards from their ticks instead of moving either label to
    # another row or shrinking the font.
    axis.figure.canvas.draw()
    renderer = axis.figure.canvas.get_renderer()
    for index, (left, right) in enumerate(zip(tick_texts, tick_texts[1:])):
        left_box = left.get_window_extent(renderer=renderer)
        right_box = right.get_window_extent(renderer=renderer)
        if left_box.x1 > right_box.x0:
            left.set_ha("right")
            right.set_ha("left")


def plot_phonon(
    q_lengths,
    frequencies,
    npoints: int,
    tick_positions,
    labels: List[str],
    save: bool,
) -> None:
    """Display or save one phonon band-structure figure."""
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
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
    })
    figure, axis = plt.subplots(figsize=(5.0, 4.0), dpi=150)
    for segment_index in range(len(q_lengths) // npoints):
        start_index = segment_index * npoints
        stop_index = (segment_index + 1) * npoints
        axis.plot(
            q_lengths[start_index:stop_index],
            frequencies[start_index:stop_index],
            color="cornflowerblue",
            linewidth=1.2,
            alpha=0.8,
        )
    for position in tick_positions[1:-1]:
        axis.axvline(position, color="0.8", linewidth=0.8)
    axis.axhline(0.0, color="gray", linewidth=1.0, linestyle="--")
    axis.set_xlim(q_lengths[0], q_lengths[-1])
    set_qpoint_ticks(axis, tick_positions, labels)
    axis.set_xlabel("Wave vector path")
    axis.set_ylabel("Frequency (THz)")
    axis.tick_params(axis="x", top=True)
    axis.tick_params(axis="y", right=True)
    figure.tight_layout()

    if save or get_backend().lower() in {"agg", "cairo", "pdf", "ps", "svg"}:
        figure.savefig("phonon.png", dpi=300, bbox_inches="tight")
        print(" Saved phonon plot to 'phonon.png'.")
    else:
        plt.show()
    plt.close(figure)


def main() -> int:
    """Read inputs, validate the QPOINTS/data relationship, and plot it."""
    try:
        parsed = parse_arguments()
        if parsed is None:
            return 0
        phonon_file, qpoints_file, save = parsed
        npoints, endpoint_labels = parse_qpoints(qpoints_file)
        q_lengths, frequencies = parse_phonon_data(phonon_file)
        plot_q_lengths, ticks, labels = _tick_positions(q_lengths, npoints, endpoint_labels)
        plot_phonon(plot_q_lengths, frequencies, npoints, ticks, labels, save)
    except (OSError, ValueError) as error:
        print(f" Error: {error}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     calc_phonon.py
Category:   Calculator Scripts
Purpose:    Calculate a phonon band structure from a primitive cell and a
            Calorine CPUNEP calculator.
Usage:      gpumdkit.sh -> 4) Calculators -> 414) Phonon band structure
            python calc_phonon.py
Arguments:  None. All calculation parameters are entered interactively.
Output:
  <output>   Single-line phonon data with q-point coordinates, path lengths,
             and frequencies in THz.
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-08-01
=============================================================================
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import List, Tuple


def print_help() -> None:
    """Print the interactive calculator usage without importing dependencies."""
    print(" Usage: gpumdkit.sh -> 4) Calculators -> 414) Phonon band structure")
    print("    or: python calc_phonon.py")
    print("")
    print(" Interactive inputs:")
    print("   primitive_cell  Primitive-cell structure file (default: POSCAR)")
    print("   nep_model       Calorine-compatible NEP model (default: nep.txt)")
    print("   qpoints         VASP/VASPKIT line-mode QPOINTS file (default: QPOINTS)")
    print("   supercell       Three positive integers, for example: 1 1 1")
    print("   displacement    Displacement distance in Angstrom (default: 0.015)")
    print("   output          Output data file (default: phonon_NEP.dat)")
    print("")
    print(" The number of points per segment and the q-point path are read from QPOINTS.")


args = sys.argv[1:]
if args and args[0] in ("-h", "--help"):
    print_help()
    sys.exit(0)
if args:
    print(" Error: calc_phonon.py is interactive and does not accept positional arguments.")
    print_help()
    sys.exit(1)


try:
    import numpy as np
    from ase.io import read
    from calorine.calculators import CPUNEP
    from calorine.tools import get_force_constants
except ImportError as error:
    print(f" Error: phonon calculation dependencies are unavailable: {error}")
    print("        Use an environment with ASE, Calorine, phonopy, and NumPy installed.")
    sys.exit(1)


def prompt_value(label: str, default: str = "") -> str:
    """Read one interactive value, applying a displayed default when blank."""
    suffix = f" [{default}]" if default else ""
    try:
        value = input(f" {label}{suffix}: ").strip()
    except EOFError:
        print(" Error: input closed while reading phonon parameters.")
        sys.exit(1)
    return value or default


def parse_supercell(value: str) -> Tuple[int, int, int]:
    """Parse and validate the three diagonal supercell factors."""
    fields = value.split()
    if len(fields) != 3:
        raise ValueError("supercell must contain exactly three positive integers")
    try:
        factors = tuple(int(field) for field in fields)
    except ValueError as error:
        raise ValueError("supercell must contain exactly three positive integers") from error
    if any(factor <= 0 for factor in factors):
        raise ValueError("supercell factors must be positive")
    return factors


def parse_qpoints(filename: Path) -> Tuple[int, List[Tuple[np.ndarray, np.ndarray]], List[str], np.ndarray]:
    """Read a line-mode QPOINTS file and construct its complete q-point path."""
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

    endpoints: List[Tuple[np.ndarray, str]] = []
    for line in lines[4:]:
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            coordinate = np.array([float(fields[0]), float(fields[1]), float(fields[2])])
        except ValueError:
            continue
        label = fields[-1]
        endpoints.append((coordinate, label))

    if not endpoints:
        raise ValueError(f"QPOINTS file '{filename}' contains no q-point endpoints")
    if len(endpoints) % 2:
        raise ValueError("QPOINTS line-mode endpoints must occur in start/end pairs")

    segments: List[Tuple[np.ndarray, np.ndarray]] = []
    labels = [endpoints[0][1]]
    for index in range(0, len(endpoints), 2):
        start, _start_label = endpoints[index]
        end, end_label = endpoints[index + 1]
        segments.append((start, end))
        labels.append(end_label)

    kpath = np.vstack(
        [np.linspace(start, end, npoints, endpoint=True) for start, end in segments]
    )
    if not np.all(np.isfinite(kpath)):
        raise ValueError("QPOINTS contains non-finite q-point coordinates")
    return npoints, segments, labels, kpath


def save_phonon_data(
    q_lengths: np.ndarray,
    kpath: np.ndarray,
    frequencies: np.ndarray,
    filename: Path,
) -> None:
    """Write one complete q-point and all phonon frequencies per output line."""
    if q_lengths.ndim != 1 or kpath.ndim != 2 or kpath.shape[1] != 3:
        raise ValueError("phonon q-point arrays have unexpected dimensions")
    if len(q_lengths) != len(kpath) or frequencies.ndim != 2:
        raise ValueError("phonon q-point and frequency arrays have inconsistent dimensions")
    if len(frequencies) != len(q_lengths):
        raise ValueError("phonon q-point and frequency arrays have different lengths")
    if not np.all(np.isfinite(q_lengths)) or not np.all(np.isfinite(frequencies)):
        raise ValueError("phonon output contains non-finite values")
    if filename.parent != Path(".") and not filename.parent.is_dir():
        raise OSError(f"output directory does not exist: '{filename.parent}'")

    with filename.open("w", encoding="utf-8") as handle:
        handle.write(
            " q-point No.           q_x           q_y           q_z  "
            "q-path length  frequencies [THz]\n"
        )
        for index, (qpoint, distance, band_values) in enumerate(
            zip(kpath, q_lengths, frequencies), start=1
        ):
            line = (
                f"{index:12d}{qpoint[0]:14.6f}{qpoint[1]:14.6f}"
                f"{qpoint[2]:14.6f}{distance:15.6f}"
            )
            line += "".join(f"{frequency:14.6f}" for frequency in band_values)
            handle.write(line + "\n")


def main() -> int:
    """Collect inputs, calculate the band structure, and save the result."""
    primitive_cell = Path(prompt_value("Primitive-cell structure file", "POSCAR"))
    nep_model = Path(prompt_value("NEP model file", "nep.txt"))
    qpoints_file = Path(prompt_value("QPOINTS file", "QPOINTS"))
    supercell_text = prompt_value("Supercell factors (a b c)", "1 1 1")
    displacement_text = prompt_value("Displacement distance in Angstrom", "0.015")
    output_file = Path(prompt_value("Output file", "phonon_NEP.dat"))

    for input_file, description in (
        (primitive_cell, "primitive-cell structure"),
        (nep_model, "NEP model"),
        (qpoints_file, "QPOINTS"),
    ):
        if not input_file.is_file():
            print(f" Error: {description} file not found: '{input_file}'")
            return 1

    try:
        supercell = parse_supercell(supercell_text)
        displacement = float(displacement_text)
    except ValueError as error:
        print(f" Error: invalid phonon parameter: {error}")
        return 1
    if displacement <= 0 or not np.isfinite(displacement):
        print(" Error: displacement distance must be a finite positive number.")
        return 1

    try:
        npoints, segments, labels, kpath = parse_qpoints(qpoints_file)
        structure = read(str(primitive_cell))
        calculator = CPUNEP(str(nep_model))
        phonon = get_force_constants(
            structure,
            calculator,
            np.diag(supercell),
            # Preserve the user-supplied primitive/input cell. This keeps the
            # band count and QPOINTS path consistent with the supplied data.
            kwargs_phonopy={"primitive_matrix": "P"},
            kwargs_generate_displacements={"distance": displacement},
        )
        phonon.run_band_structure([kpath])
        band = phonon.band_structure
        q_lengths = np.asarray(band.distances[0], dtype=float)
        frequencies = np.asarray(band.frequencies[0], dtype=float)
        save_phonon_data(q_lengths, kpath, frequencies, output_file)
    except Exception as error:
        print(f" Error: phonon calculation failed: {error}")
        return 1

    print(f" Read {len(segments)} q-point segments from '{qpoints_file}'.")
    print(f" Path: {' - '.join(labels)} ({npoints} points per segment)")
    print(f" Supercell: {supercell}; displacement: {displacement:g} Angstrom")
    print(f" Saved phonon data to '{output_file}'.")
    print(f" {len(q_lengths)} q-points; {frequencies.shape[1]} bands")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

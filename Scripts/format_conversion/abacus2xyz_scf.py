"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     abacus2xyz_scf.py
Category:   Format Conversion Scripts
Purpose:    Convert ABACUS SCF logs or abacus.json files to extended XYZ format.
Usage:      python abacus2xyz_scf.py <dir> <extxyz>
            (invoked interactively via gpumdkit.sh menu 104)
Arguments:
  dir      Directory containing ABACUS SCF outputs
  extxyz   Output extxyz file
Output:
  <extxyz>  (converted structures in extxyz format)
Author:     Benrui Tang (tang070205@proton.me)
Last-modified: 2026-11-17
=============================================================================
"""

import json
import os
import re
import sys
from pathlib import Path


KBAR_PER_EV_A3 = 1602.1766208


def read_input(folder):
    path = folder / "INPUT"
    if not path.is_file():
        path = folder.parent / "INPUT"
    if not path.is_file():
        return {}
    return dict(line.split()[:2] for line in path.read_text().splitlines()
                if len(line.split()) >= 2 and not line.lstrip().startswith("#"))


def parse_force_rows(force_section, natoms):
    """Extract natoms force rows from an ABACUS TOTAL-FORCE section."""
    rows = []
    for line in force_section.splitlines():
        parts = line.split()
        if len(parts) == 4 and re.fullmatch(r"[A-Za-z]+\d*", parts[0]):
            try:
                [float(val) for val in parts[1:4]]
            except ValueError:
                if rows:
                    break
                continue
            rows.append(parts)
            if len(rows) == natoms:
                break
        elif rows:
            break
    return rows


def read_scf_log(folder, scf_nmax):
    import numpy as np
    from ase.units import Bohr

    log = folder / "running_scf.log"
    contents = log.read_text(errors="replace")
    lower = contents.lower()
    energies = re.findall(r"FINAL_ETOT_IS\s+([^\s]+)", contents)
    if not energies or "convergence has not been achieved" in lower:
        return None
    if ("charge density convergence is achieved" not in lower
            and scf_nmax is not None and contents.count("ALGORITHM") >= scf_nmax):
        return None

    natoms = int(re.findall(r"TOTAL ATOM NUMBER\s*=\s*(\d+)", contents)[-1])
    if "TOTAL-FORCE (eV/Angstrom)" not in contents:
        raise ValueError("no forces found in running_scf.log")
    force_section = contents.rsplit("TOTAL-FORCE (eV/Angstrom)", 1)[-1]
    force_rows = parse_force_rows(force_section, natoms)
    if len(force_rows) != natoms:
        raise ValueError("no forces found in running_scf.log")
    symbols = [row[0].rstrip("0123456789") for row in force_rows]
    forces = [row[1:4] for row in force_rows]

    # Use the cell and positions from the same Cartesian basis as the logged forces and stress.
    lattice_constant = re.findall(r"lattice constant \(Bohr\)\s*=\s*(\S+)", contents)
    lattice = re.search(r"(?im)^\s*Lattice vectors:[^\n]*\n((?:[^\n]+\n){3})", contents)
    if not lattice_constant or lattice is None:
        raise ValueError("no complete lattice found in running_scf.log")
    scale = float(lattice_constant[-1]) * Bohr
    cell = np.array([line.split() for line in lattice[1].splitlines()], dtype=float) * scale
    coords = re.findall(r"(?m)^\s*tau([cd])_\S+[^\S\n]+(\S+)[^\S\n]+(\S+)[^\S\n]+(\S+)", contents)
    if len(coords) < natoms:
        raise ValueError("no complete atomic positions found in running_scf.log")
    positions = np.array([row[1:] for row in coords[:natoms]], dtype=float)
    positions = positions @ cell if coords[0][0] == "d" else positions * scale

    virial = None
    if "TOTAL-STRESS" in contents:
        stress_section = contents.rsplit("TOTAL-STRESS", 1)[-1]
        stress_lines = [line.split() for line in stress_section.splitlines() if len(line.split()) == 3][:3]
        volume = abs(np.linalg.det(cell))
        virial = np.array(stress_lines, dtype=float) * (volume / KBAR_PER_EV_A3)
    return symbols, positions, forces, cell, float(energies[-1]), virial


def read_scf_json(path, scf_nmax):
    import numpy as np

    data = json.loads(path.read_text())
    result = data["output"][0]
    if len(result["scf"]) == scf_nmax:
        return None
    cell = np.asarray(result["cell"])
    virial = None
    if result.get("stress") is not None:
        virial = np.asarray(result["stress"]) * (abs(np.linalg.det(cell)) / KBAR_PER_EV_A3)
    return data["init"]["label"], result["coordinate"], result["force"], cell, result["energy"], virial


def format_frame(symbols, positions, forces, cell, energy, virial, configuration):
    import numpy as np

    if len(symbols) != len(positions) or np.shape(forces) != (len(symbols), 3):
        raise ValueError("atom, position, and force counts do not match")
    lattice = " ".join(f"{value:.10f}" for value in np.asarray(cell).ravel())
    config = configuration.replace("\\", "\\\\").replace('"', '\\"')
    header = f'Energy={energy} Lattice="{lattice}" Config_type="{config}" Weight=1.0'
    if virial is not None:
        values = " ".join(f"{value:.10f}" for value in np.asarray(virial).ravel())
        header += f' Virial="{values}"'
    header += " Properties=species:S:1:pos:R:3:forces:R:3"
    rows = [f"{symbol} " + " ".join(f"{float(value):.10f}" for value in (*position, *force))
            for symbol, position, force in zip(symbols, positions, forces)]
    return "\n".join([str(len(symbols)), header, *rows]) + "\n"


def main():
    args = sys.argv[1:]
    if len(args) != 2 or args[0] in ("-h", "--help"):
        print(" Usage: python abacus2xyz_scf.py <dir> <extxyz>")
        print(" Converts converged ABACUS SCF results under <dir> to an extxyz dataset.")
        return 0 if args and args[0] in ("-h", "--help") else 1

    root, output = map(Path, args)
    if not root.is_dir():
        print(f" Error: input directory '{root}' does not exist.")
        return 1

    converted = 0
    for directory, _, files in os.walk(root, followlinks=True):
        folder = Path(directory)
        if "running_scf.log" not in files and "abacus.json" not in files:
            continue
        try:
            scf_nmax = read_input(folder).get("scf_nmax")
            scf_nmax = int(scf_nmax) if scf_nmax is not None else None
            if "running_scf.log" in files:
                result = read_scf_log(folder, scf_nmax)
            else:
                result = read_scf_json(folder / "abacus.json", scf_nmax)
            if result is None:
                print(f" Skipping incomplete or unconverged SCF in {folder}.")
                continue
            frame = format_frame(*result, "/".join(folder.parts[-3:]))
            output.parent.mkdir(parents=True, exist_ok=True)
            with output.open("a" if converted else "w") as stream:
                stream.write(frame)
        except (OSError, ValueError, KeyError, IndexError, ImportError) as exc:
            print(f" Error: {folder}: {exc}")
            return 1
        converted += 1

    if not converted:
        print(f" Error: no converged ABACUS SCF structures found under '{root}'.")
        return 1
    print(f" Converted {converted} SCF structures to {output}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

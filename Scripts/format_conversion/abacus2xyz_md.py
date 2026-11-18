"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     abacus2xyz_md.py
Category:   Format Conversion Scripts
Purpose:    Convert an ABACUS MD trajectory to extxyz.
Usage:      python abacus2xyz_md.py <dir> <extxyz>
            (invoked interactively via gpumdkit.sh menu 104)
Arguments:
  dir      Directory containing running_md.log and MD_dump
  extxyz   Output extxyz file
Output:
  <extxyz>  Converted MD structures with energies and forces
Author:     Benrui Tang (tang070205@proton.me)
Last-modified: 2026-11-17
=============================================================================
"""

import os
import re
import sys
from pathlib import Path

from abacus2xyz_scf import KBAR_PER_EV_A3, format_frame


def convert_md(folder, output):
    import numpy as np
    from ase.units import Bohr

    log = (folder / "running_md.log").read_text(errors="replace")
    energies = []
    for line in log.splitlines():
        if "final etot is" in line or "#TOTAL ENERGY#" in line:
            energies.append(float(line.split()[-2]))
        elif "!! convergence has not been achieved" in line:
            energies.append(None)

    first_step = re.search(r"STEP OF MOLECULAR DYNAMICS\s*:\s*(\d+)", log)
    start = int(first_step[1]) if first_step else 0
    dump = (folder / "MD_dump").read_text(errors="replace")
    blocks = re.split(r"(?m)^MDSTEP:\s*(\d+)[^\S\n]*\n", dump)
    steps, frames = blocks[1::2], blocks[2::2]
    natoms = int(re.findall(r"TOTAL ATOM NUMBER\s*=\s*(\d+)", log)[-1])

    converted = 0
    for step, frame in zip(steps[1:], frames[1:]):  # Preserve the menu's initial-frame skip.
        energy_index = int(step) - start
        if not 0 <= energy_index < len(energies):
            raise ValueError(f"no matching energy in running_md.log for MD step {step}")
        energy = energies[energy_index]
        if energy is None:
            print(f" Skipping unconverged MD step {step}.")
            continue

        lines = frame.splitlines()
        lattice_constant = lines[0].split()
        new_format = "angstrom" in lines[0].lower()
        scale = float(lattice_constant[1]) * (1 if new_format else Bohr)
        cell = np.array([line.split() for line in lines[2:5]], dtype=float) * scale

        header = next(i for i, line in enumerate(lines) if line.startswith("INDEX"))
        if "FORCE" not in lines[header]:
            raise ValueError("force columns are missing from MD_dump")
        rows = [line.split() for line in lines[header + 1:header + 1 + natoms]]
        if len(rows) != natoms:
            raise ValueError(f"incomplete MD_dump frame at step {step}")
        symbols = [row[1] for row in rows]
        positions = np.array([row[2:5] for row in rows], dtype=float)
        if not new_format:
            positions *= scale
        forces = [row[5:8] for row in rows]

        virial = None
        if lines[5].upper().startswith("VIRIAL"):
            stress = np.array([line.split() for line in lines[6:9]], dtype=float)
            virial = stress * (abs(np.linalg.det(cell)) / KBAR_PER_EV_A3)
        config = "/".join(folder.parts[-3:]) + f"-{step}"
        record = format_frame(symbols, positions, forces, cell, energy, virial, config)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("a" if converted else "w") as stream:
            stream.write(record)
        converted += 1

    return converted


def main():
    args = sys.argv[1:]
    if len(args) != 2 or args[0] in ("-h", "--help"):
        print(" Usage: python abacus2xyz_md.py <dir> <extxyz>")
        print(" Converts the first ABACUS MD run found under <dir>.")
        return 0 if args and args[0] in ("-h", "--help") else 1

    root, output = map(Path, args)
    if not root.is_dir():
        print(f" Error: input directory '{root}' does not exist.")
        return 1

    for directory, _, files in os.walk(root, followlinks=True):
        if "running_md.log" in files and "MD_dump" in files:
            try:
                converted = convert_md(Path(directory), output)
            except (OSError, ValueError, IndexError, ImportError, StopIteration) as exc:
                print(f" Error: failed to convert ABACUS MD in '{directory}': {exc}")
                return 1
            if not converted:
                print(f" Error: no converged MD frames to convert in '{directory}'.")
                return 1
            print(f" Converted {converted} MD frames to {output}.")
            return 0

    print(f" Error: no directory with both running_md.log and MD_dump found under '{root}'.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

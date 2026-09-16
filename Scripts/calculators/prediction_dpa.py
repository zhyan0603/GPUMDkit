#!/usr/bin/env python3
"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     prediction_dpa.py
Category:   Calculator Scripts
Purpose:    Evaluate a DeepMD DPA model on a labeled training set and write
            NEP-compatible prediction output files.
Usage:      gpumdkit.sh -prediction_dpa <input.xyz> <dpa_model>
            python prediction_dpa.py <input.xyz> <dpa_model>
Arguments:
  input.xyz  Labeled extended-XYZ training set
  dpa_model  DeepMD PyTorch checkpoint (.pt) or frozen model (.pt2)
Output:
  energy_train.out  Predicted and target energy, eV/atom
  force_train.out   Predicted and target forces, eV/Angstrom
  virial_train.out  Predicted and target virial, eV/atom
  stress_train.out  Predicted and target stress, GPa
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-09-12
=============================================================================

The four output files follow the NEP/GPUMD parity-file convention. The model
is evaluated with the public ``deepmd.infer.DeepPot`` interface, which can
load a PyTorch checkpoint (``.pt``) and, when available, a frozen
PyTorch-exportable model (``.pt2``).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np


EV_PER_A3_TO_GPA = 160.21766208
FLOAT_FORMAT = "%.10e"


def _parse_comment(comment: str) -> dict[str, str]:
    """Parse quoted and unquoted key=value fields from an XYZ comment line."""
    fields: dict[str, str] = {}
    pattern = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=(?:\"([^\"]*)\"|(\S+))")
    for match in pattern.finditer(comment):
        fields[match.group(1).lower()] = (
            match.group(2) if match.group(2) is not None else match.group(3)
        )
    return fields


def _parse_properties(specification: str) -> dict[str, tuple[int, int]]:
    """Return column offsets and widths from an extended-XYZ Properties field."""
    properties: dict[str, tuple[int, int]] = {}
    offset = 0
    tokens = specification.split(":")
    if len(tokens) % 3 != 0:
        raise ValueError(f"Malformed Properties field: {specification!r}")
    for name, _kind, width_text in zip(tokens[0::3], tokens[1::3], tokens[2::3]):
        width = int(width_text)
        properties[name.lower()] = (offset, width)
        offset += width
    return properties


def _read_frames(path: Path):
    """Yield frames from the labeled extended-XYZ file used by this project."""
    with path.open() as handle:
        frame_index = 0
        while True:
            line = handle.readline()
            if not line:
                return
            if not line.strip():
                continue
            try:
                natoms = int(line.strip())
            except ValueError as exc:
                raise ValueError(
                    f"Expected atom count before frame {frame_index}, got {line!r}"
                ) from exc
            comment = handle.readline().rstrip("\n")
            tags = _parse_comment(comment)
            if "lattice" not in tags:
                raise ValueError(f"Frame {frame_index} has no Lattice field")
            cell = np.fromstring(tags["lattice"], sep=" ", dtype=np.float64)
            if cell.size != 9:
                raise ValueError(f"Frame {frame_index} Lattice must have 9 values")
            cell = cell.reshape(3, 3)

            properties = _parse_properties(tags.get("properties", ""))
            required = {"species", "pos"}
            missing = required.difference(properties)
            if missing:
                raise ValueError(f"Frame {frame_index} is missing Properties {missing}")
            force_key = "forces" if "forces" in properties else "force"
            if force_key not in properties:
                raise ValueError(f"Frame {frame_index} has no force field")

            species: list[str] = []
            positions = np.empty((natoms, 3), dtype=np.float64)
            forces = np.empty((natoms, 3), dtype=np.float64)
            for atom_index in range(natoms):
                atom_line = handle.readline()
                if not atom_line:
                    raise ValueError(f"Unexpected EOF in frame {frame_index}")
                values = atom_line.split()
                species_offset, species_width = properties["species"]
                if species_width != 1:
                    raise ValueError("Only scalar species fields are supported")
                species.append(values[species_offset])
                pos_offset, pos_width = properties["pos"]
                force_offset, force_width = properties[force_key]
                if pos_width != 3 or force_width != 3:
                    raise ValueError("Position and force fields must have width 3")
                positions[atom_index] = np.asarray(
                    values[pos_offset : pos_offset + 3], dtype=np.float64
                )
                forces[atom_index] = np.asarray(
                    values[force_offset : force_offset + 3], dtype=np.float64
                )

            energy_text = tags.get("energy")
            energy = None if energy_text is None else float(energy_text)
            virial_text = tags.get("virial")
            virial = None
            if virial_text is not None:
                virial = np.fromstring(virial_text, sep=" ", dtype=np.float64)
                if virial.size != 9:
                    raise ValueError(f"Frame {frame_index} virial must have 9 values")

            yield {
                "index": frame_index,
                "natoms": natoms,
                "species": species,
                "positions": positions,
                "forces": forces,
                "cell": cell,
                "energy": energy,
                "virial": None if virial is None else virial.reshape(3, 3),
            }
            frame_index += 1


def _voigt(tensor: np.ndarray) -> np.ndarray:
    """Convert a 3x3 tensor to GPUMD order xx, yy, zz, xy, yz, zx."""
    return np.asarray(
        [
            tensor[0, 0],
            tensor[1, 1],
            tensor[2, 2],
            tensor[0, 1],
            tensor[1, 2],
            tensor[2, 0],
        ],
        dtype=np.float64,
    )


def _write(path: Path, rows: list[np.ndarray]) -> None:
    if not rows:
        path.write_text("")
        return
    np.savetxt(path, np.vstack(rows), fmt=FLOAT_FORMAT)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Predict a DeepMD DPA model and write GPUMD-style parity files."
    )
    parser.add_argument(
        "xyz",
        type=Path,
        help="labeled extended-XYZ training set",
    )
    parser.add_argument(
        "model",
        type=Path,
        help="DeepMD PyTorch checkpoint (.pt) or frozen model (.pt2)",
    )
    args = parser.parse_args()

    if not args.model.is_file():
        print(f" Error: model file '{args.model}' does not exist.")
        sys.exit(1)
    if not args.xyz.is_file():
        print(f" Error: XYZ file '{args.xyz}' does not exist.")
        sys.exit(1)

    try:
        from deepmd.infer import DeepPot
    except ImportError as error:
        print(" Error: this command requires the 'deepmd-kit' package.")
        print(f" Details: {error}")
        sys.exit(1)

    model = DeepPot(str(args.model), auto_batch_size=True)
    type_map = list(model.get_type_map())
    type_index = {name: index for index, name in enumerate(type_map)}
    print(f" model: {args.model}")
    print(f" type_map: {type_map}")
    print(f" rcut: {model.get_rcut():g} A")

    energy_rows: list[np.ndarray] = []
    force_rows: list[np.ndarray] = []
    virial_rows: list[np.ndarray] = []
    stress_rows: list[np.ndarray] = []
    energy_errors: list[float] = []
    force_errors: list[np.ndarray] = []
    virial_errors: list[np.ndarray] = []
    nframes = 0
    natoms_total = 0

    for frame in _read_frames(args.xyz):
        unknown = sorted(set(frame["species"]).difference(type_index))
        if unknown:
            raise ValueError(
                f"Frame {frame['index']} contains elements absent from model: {unknown}"
            )
        atom_types = np.asarray([type_index[x] for x in frame["species"]], dtype=np.int32)
        pred_energy, pred_force, pred_virial = model.eval(
            frame["positions"][None, :, :],
            frame["cell"][None, :, :],
            atom_types,
        )
        pred_energy = float(np.asarray(pred_energy).reshape(-1)[0])
        pred_force = np.asarray(pred_force, dtype=np.float64).reshape(frame["natoms"], 3)
        pred_virial = np.asarray(pred_virial, dtype=np.float64).reshape(3, 3)

        if frame["energy"] is None:
            raise ValueError(f"Frame {frame['index']} has no target energy")
        target_energy = float(frame["energy"])
        energy_rows.append(np.asarray([pred_energy / frame["natoms"], target_energy / frame["natoms"]]))
        force_rows.append(np.column_stack((pred_force, frame["forces"])))

        target_virial = frame["virial"]
        if target_virial is None:
            target_virial = np.full((3, 3), -1.0e6, dtype=np.float64)
        pred_virial_voigt = _voigt(pred_virial)
        target_virial_voigt = _voigt(target_virial)
        virial_rows.append(
            np.concatenate(
                (pred_virial_voigt / frame["natoms"], target_virial_voigt / frame["natoms"])
            )
        )

        volume = abs(float(np.linalg.det(frame["cell"])))
        if volume <= 0.0:
            raise ValueError(f"Frame {frame['index']} has non-positive cell volume")
        stress_rows.append(
            np.concatenate(
                (
                    -pred_virial_voigt * EV_PER_A3_TO_GPA / volume,
                    -target_virial_voigt * EV_PER_A3_TO_GPA / volume,
                )
            )
        )

        energy_errors.append(pred_energy / frame["natoms"] - target_energy / frame["natoms"])
        force_errors.append((pred_force - frame["forces"]).reshape(-1))
        virial_errors.append(pred_virial_voigt / frame["natoms"] - target_virial_voigt / frame["natoms"])
        nframes += 1
        natoms_total += frame["natoms"]
        if nframes % 100 == 0:
            print(f" predicted {nframes} frames ({natoms_total} atoms)")

    output_dir = Path.cwd()
    _write(output_dir / "energy_train.out", energy_rows)
    _write(output_dir / "force_train.out", force_rows)
    _write(output_dir / "virial_train.out", virial_rows)
    _write(output_dir / "stress_train.out", stress_rows)

    energy_rmse = float(np.sqrt(np.mean(np.square(energy_errors)))) if energy_errors else float("nan")
    force_rmse = (
        float(np.sqrt(np.mean(np.square(np.concatenate(force_errors)))))
        if force_errors
        else float("nan")
    )
    virial_rmse = (
        float(np.sqrt(np.mean(np.square(np.vstack(virial_errors)))))
        if virial_errors
        else float("nan")
    )
    print(f" wrote {nframes} frames to {output_dir}")
    print(f" RMSE energy/atom = {energy_rmse:.6e} eV")
    print(f" RMSE force       = {force_rmse:.6e} eV/A")
    print(f" RMSE virial/atom = {virial_rmse:.6e} eV")


if __name__ == "__main__":
    main()

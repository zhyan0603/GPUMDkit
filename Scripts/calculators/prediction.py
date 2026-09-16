"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     prediction.py
Category:   Calculator Scripts
Purpose:    Evaluate every frame in an extended XYZ file with a Calorine
            CPUNEP calculator and write NEP-compatible prediction outputs.
Usage:      gpumdkit.sh -prediction <input.xyz> <nep.txt> [workers]
            python prediction.py <input.xyz> <nep.txt> [workers]
Arguments:
  input.xyz  Extended XYZ file containing energy, force, and optional
             virial/stress targets
  nep.txt    Calorine-compatible NEP model file
  workers    Number of CPU workers (default: 1)
Output:
  energy_<input-stem>.out  Predicted and target energy, eV/atom
  force_<input-stem>.out   Predicted and target forces, eV/Angstrom
  stress_<input-stem>.out  Predicted and target stress, GPa
  virial_<input-stem>.out  Predicted and target virial, eV/atom
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-08-28
=============================================================================
"""

import os
import sys


args = sys.argv[1:]


def print_usage():
    """Print command-line usage and output-file conventions."""
    print(" Usage: gpumdkit.sh -prediction <input.xyz> <nep.txt> [workers]")
    print("    or: python prediction.py <input.xyz> <nep.txt> [workers]")
    print("")
    print(" Arguments:")
    print("   input.xyz  Extended XYZ file with energy, force, and optional")
    print("              virial/stress targets")
    print("   nep.txt    Calorine-compatible NEP model file")
    print("   workers    Number of CPU workers (default: 1)")
    print("")
    print(" Output:")
    print("   energy_<input-stem>.out, force_<input-stem>.out")
    print("   stress_<input-stem>.out, virial_<input-stem>.out")
    print("")
    print(" Example: gpumdkit.sh -prediction train.xyz nep.txt 8")
    print("")


if __name__ == "__main__" and (
    len(args) < 2 or args[0] in ("-h", "--help")
):
    print_usage()
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

if __name__ == "__main__" and len(args) > 3:
    print(" Error: expected <input.xyz> <nep.txt> and optional [workers].")
    print_usage()
    sys.exit(1)

try:
    requested_workers = int(args[2]) if len(args) == 3 else 1
except ValueError:
    print(" Error: workers must be a positive integer.")
    sys.exit(1)

if requested_workers < 1:
    print(" Error: workers must be a positive integer.")
    sys.exit(1)

# One native OpenMP thread per process makes the worker count the total CPU
# parallelism and keeps the default mode genuinely single-core.
os.environ["OMP_NUM_THREADS"] = "1"

try:
    from concurrent.futures import ProcessPoolExecutor, as_completed

    import numpy as np
    from ase.io import read
    from ase.units import GPa
    from calorine.calculators import CPUNEP
    from tqdm import tqdm
except ImportError as error:
    print(" Error: this command requires ASE, calorine, and tqdm.")
    print(f" Details: {error}")
    print(" Install them with: pip install ase calorine tqdm")
    sys.exit(1)


SENTINEL = -1.0e6
NEP_REDUCED6_INDEX = np.array([0, 1, 2, 5, 3, 4])


def print_dependency_notice():
    """Print the dependency and citation notice for Calorine."""
    print(" This function requires the calorine package.")
    print(" If you use this function, we recommend citing:")
    print(" Lindgren et al., J. Open Source Softw. 9, 6264 (2024).")
    print(" https://doi.org/10.21105/joss.06264")


def _calculator_result(atoms, name):
    """Return a calculator result or extxyz info value for one property."""
    calculator = atoms.calc
    if calculator is not None:
        results = getattr(calculator, "results", {})
        for key in (name, name.capitalize(), name.upper()):
            result = results.get(key)
            if result is not None:
                return result
    for key in (name, name.capitalize(), name.upper()):
        result = atoms.info.get(key)
        if result is not None:
            return result
    return None


def _as_scalar(value, property_name, frame_index):
    """Convert one scalar target and reject malformed or non-finite values."""
    values = np.asarray(value, dtype=float).reshape(-1)
    if values.size != 1:
        raise ValueError(
            f"frame {frame_index}: {property_name} must contain one value."
        )
    if not np.isfinite(values[0]):
        raise ValueError(f"frame {frame_index}: {property_name} contains NaN or Inf.")
    return float(values[0])


def _as_nep_reduced6(value, property_name, frame_index):
    """Convert a symmetric tensor to NEP's ``xx yy zz xy yz xz`` order."""
    if isinstance(value, str):
        values = np.fromstring(value, sep=" ", dtype=float)
    else:
        values = np.asarray(value, dtype=float)

    if values.size == 9:
        full_tensor = values.reshape(3, 3)
        reduced = np.array(
            [
                full_tensor[0, 0],
                full_tensor[1, 1],
                full_tensor[2, 2],
                full_tensor[0, 1],
                full_tensor[1, 2],
                full_tensor[0, 2],
            ],
            dtype=float,
        )
    elif values.size == 6:
        # ASE stores stress in xx, yy, zz, yz, xz, xy order. This also makes
        # six-component virial inputs consistent with the ASE representation.
        reduced = values.reshape(6)[NEP_REDUCED6_INDEX]
    else:
        raise ValueError(
            f"frame {frame_index}: {property_name} must contain 6 or 9 values."
        )

    if not np.all(np.isfinite(reduced)):
        raise ValueError(f"frame {frame_index}: {property_name} contains NaN or Inf.")
    return reduced


def _target_stress_virial(atoms, frame_index):
    """Return per-atom virial and stress targets, deriving a missing form."""
    natoms = len(atoms)
    volume = atoms.get_volume()
    raw_virial = _calculator_result(atoms, "virial")
    raw_stress = _calculator_result(atoms, "stress")

    target_virial = None
    target_stress = None

    if raw_virial is not None:
        virial_total = _as_nep_reduced6(raw_virial, "virial", frame_index)
        target_virial = virial_total / natoms
        if raw_stress is None:
            target_stress = virial_total / volume / GPa

    if raw_stress is not None:
        stress_density = _as_nep_reduced6(raw_stress, "stress", frame_index)
        target_stress = stress_density / GPa
        if raw_virial is None:
            target_virial = stress_density * volume / natoms

    if target_virial is None:
        target_virial = np.full(6, SENTINEL, dtype=float)
    if target_stress is None:
        target_stress = np.full(6, SENTINEL, dtype=float)

    return target_stress, target_virial


def _prepare_frame(atoms, frame_index):
    """Validate one input frame and detach its reference labels for workers."""
    if len(atoms) == 0:
        raise ValueError(f"frame {frame_index}: an empty structure is not supported.")
    if atoms.cell.rank != 3 or atoms.get_volume() <= 0:
        raise ValueError(
            f"frame {frame_index}: a valid three-dimensional cell is required."
        )

    energy = _calculator_result(atoms, "energy")
    if energy is None:
        raise ValueError(f"frame {frame_index}: energy target is missing.")
    target_energy = _as_scalar(energy, "energy", frame_index) / len(atoms)

    forces = _calculator_result(atoms, "forces")
    if forces is None:
        forces = atoms.arrays.get("forces", atoms.arrays.get("force"))
    if forces is None:
        raise ValueError(f"frame {frame_index}: force target is missing.")
    target_forces = np.asarray(forces, dtype=float)
    if target_forces.shape != (len(atoms), 3):
        raise ValueError(
            f"frame {frame_index}: force target must have shape ({len(atoms)}, 3)."
        )
    if not np.all(np.isfinite(target_forces)):
        raise ValueError(f"frame {frame_index}: force target contains NaN or Inf.")

    target_stress, target_virial = _target_stress_virial(atoms, frame_index)

    structure = atoms.copy()
    structure.calc = None
    targets = {
        "energy": target_energy,
        "forces": target_forces,
        "stress": target_stress,
        "virial": target_virial,
    }
    return frame_index, structure, targets


def _predict_chunk(payload, show_progress=False):
    """Predict a contiguous frame chunk with one independent CPUNEP object."""
    model_path, frames = payload
    calculator = CPUNEP(model_path)
    predictions = []

    frame_iterator = frames
    if show_progress:
        frame_iterator = tqdm(frames, total=len(frames), desc="Predicting", unit="frame")

    for frame_index, atoms, targets in frame_iterator:
        atoms.calc = calculator
        natoms = len(atoms)
        energy = float(atoms.get_potential_energy()) / natoms
        forces = np.asarray(atoms.get_forces(), dtype=float)
        stress_ase = np.asarray(atoms.get_stress(), dtype=float)
        if stress_ase.shape != (6,) or not np.all(np.isfinite(stress_ase)):
            raise ValueError(
                f"frame {frame_index}: Calorine returned an invalid stress tensor."
            )

        # Calorine/ASE uses the opposite stress sign from NEP's stress/virial
        # output convention. NEP stress and virial therefore have the same
        # sign, with stress expressed as an energy density.
        stress_density_nep = -stress_ase[NEP_REDUCED6_INDEX]
        stress = stress_density_nep / GPa
        virial = stress_density_nep * atoms.get_volume() / natoms
        predictions.append(
            (
                frame_index,
                np.array([energy, targets["energy"]], dtype=float),
                np.column_stack((forces, targets["forces"])),
                np.concatenate((stress, targets["stress"])),
                np.concatenate((virial, targets["virial"])),
            )
        )

    return predictions


def _run_predictions(frames, model_path, workers):
    """Run predictions serially or in ordered CPU worker processes."""
    worker_count = min(workers, len(frames))
    chunks = []
    chunk_count = min(len(frames), worker_count * 4)
    chunk_size = (len(frames) + chunk_count - 1) // chunk_count
    for start in range(0, len(frames), chunk_size):
        chunks.append(frames[start:start + chunk_size])

    if worker_count == 1:
        return _predict_chunk((model_path, frames), show_progress=True), worker_count

    predictions = []
    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        payloads = [(model_path, chunk) for chunk in chunks]
        futures = [executor.submit(_predict_chunk, payload) for payload in payloads]
        with tqdm(total=len(frames), desc="Predicting", unit="frame") as progress:
            for future in as_completed(futures):
                chunk_predictions = future.result()
                predictions.extend(chunk_predictions)
                progress.update(len(chunk_predictions))
    return predictions, worker_count


def _write_outputs(predictions, output_dir, input_stem):
    """Write the four NEP-compatible predicted/target output files."""
    ordered = sorted(predictions, key=lambda item: item[0])
    energy = np.vstack([item[1] for item in ordered])
    force = np.vstack([item[2] for item in ordered])
    stress = np.vstack([item[3] for item in ordered])
    virial = np.vstack([item[4] for item in ordered])

    output_data = {
        "energy": energy,
        "force": force,
        "stress": stress,
        "virial": virial,
    }
    output_paths = {}
    for property_name, data in output_data.items():
        output_path = os.path.join(
            output_dir, f"{property_name}_{input_stem}.out"
        )
        np.savetxt(output_path, data, fmt="%.16e")
        output_paths[property_name] = output_path
    return output_paths


def main():
    """Validate inputs, calculate all frames, and write NEP output files."""
    input_file = os.path.abspath(args[0])
    model_file = os.path.abspath(args[1])
    if not os.path.isfile(input_file):
        print(f" Error: input file '{args[0]}' does not exist.")
        sys.exit(1)
    if not os.path.isfile(model_file):
        print(f" Error: NEP model file '{args[1]}' does not exist.")
        sys.exit(1)

    print_dependency_notice()
    try:
        structures = read(input_file, index=":")
        if not isinstance(structures, list):
            structures = [structures]
        if not structures:
            print(" Error: input file contains no structures.")
            sys.exit(1)

        frames = [
            _prepare_frame(atoms, frame_index)
            for frame_index, atoms in enumerate(structures)
        ]
        print(
            f" Processing {len(frames)} frame(s) with {requested_workers} "
            "CPU worker(s)."
        )
        predictions, worker_count = _run_predictions(
            frames, model_file, requested_workers
        )
        input_stem = os.path.splitext(os.path.basename(input_file))[0] or "input"
        _write_outputs(predictions, os.getcwd(), input_stem)
    except Exception as error:
        print(f" Error: prediction failed: {error}")
        sys.exit(1)

    print(
        f" Prediction completed for {len(predictions)} frame(s) using "
        f"{worker_count} CPU worker(s)."
    )


if __name__ == "__main__":
    main()

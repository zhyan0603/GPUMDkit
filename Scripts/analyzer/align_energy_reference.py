"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     align_energy_reference.py
Category:   Analyzer Scripts
Purpose:    Interactively shift atomic energy references in an extxyz file.
            Three alignment modes are available: reference-group alignment,
            zero atomic baseline alignment, and DFT-to-NEP alignment.
Usage:      gpumdkit.sh -shift_energy
            python3 align_energy_reference.py
Arguments:  None. Run without arguments to open the interactive interface.
Output:     A new extxyz file with updated energy= values.
Author:     Zherui Chen (chenzherui0124@foxmail.com)
Last-modified: 2026-07-29
=============================================================================
"""

import argparse
import os
import re
import sys
from collections import Counter
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np


REF_GROUP_ALIGNMENT = "REF_GROUP_ALIGNMENT"
ZERO_BASELINE_ALIGNMENT = "ZERO_BASELINE_ALIGNMENT"
DFT_TO_NEP_ALIGNMENT = "DFT_TO_NEP_ALIGNMENT"

DEFAULT_INPUT_FILE = "train.xyz"
DEFAULT_OUTPUT_FILE = "output_aligned.xyz"
DEFAULT_REFERENCE_GROUP = "cp2k2xyz"
DEFAULT_NEP_MODEL_FILE = "nep.txt"
DEFAULT_MAX_GENERATIONS = 10000
DEFAULT_POPULATION_SIZE = 40
DEFAULT_CONVERGENCE_TOL = 1e-8
DEFAULT_RANDOM_SEED = 42
DEFAULT_NEP_BATCH_SIZE = 32

ENERGY_PATTERN = re.compile(
    r"(\benergy\s*=\s*)[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",
    re.IGNORECASE,
)
CONFIG_TYPE_PATTERN = re.compile(
    r"\bconfig_type\s*=\s*(?:\"([^\"]+)\"|'([^']+)'|([^\s]+))",
    re.IGNORECASE,
)


class AlignmentError(Exception):
    """Expected input or processing error shown to the user."""


class AlignmentCancelled(Exception):
    """User cancelled the interactive operation."""


@dataclass
class FrameMetadata:
    """Metadata needed for optimization and energy rewriting."""

    frame_index: int
    n_atoms: int
    primary_energy: float
    config_type: str
    element_counts: Counter


@dataclass
class AlignmentSettings:
    """Interactive settings for one alignment operation."""

    input_file: str
    output_file: str
    mode: str
    reference_group: str = DEFAULT_REFERENCE_GROUP
    shift_groups: Optional[List[str]] = None
    nep_model_file: str = DEFAULT_NEP_MODEL_FILE
    max_generations: int = DEFAULT_MAX_GENERATIONS
    population_size: int = DEFAULT_POPULATION_SIZE
    convergence_tol: float = DEFAULT_CONVERGENCE_TOL
    random_seed: int = DEFAULT_RANDOM_SEED
    nep_batch_size: int = DEFAULT_NEP_BATCH_SIZE


def build_parser() -> argparse.ArgumentParser:
    """Build the intentionally small command-line interface."""
    return argparse.ArgumentParser(
        description="Interactively shift atomic energy references in an extxyz file.",
        epilog=(
            "Run without arguments to open the interactive interface. "
            "The same interface is available from menu 5 -> 509."
        ),
    )


def split_line_ending(line: str) -> Tuple[str, str]:
    """Return line content and its original line-ending sequence."""
    if line.endswith("\r\n"):
        return line[:-2], "\r\n"
    if line.endswith("\n") or line.endswith("\r"):
        return line[:-1], line[-1:]
    return line, ""


def parse_config_type(header: str) -> str:
    """Extract config_type while accepting quoted and unquoted values."""
    match = CONFIG_TYPE_PATTERN.search(header)
    if match is None:
        return "default_group"
    return next(value for value in match.groups() if value is not None)


def parse_xyz_metadata(input_file: str) -> List[FrameMetadata]:
    """Read frame metadata without retaining the atom lines in memory."""
    frames: List[FrameMetadata] = []

    try:
        with open(input_file, "r", encoding="utf-8", newline="") as fin:
            frame_index = 0
            while True:
                n_atoms_line = fin.readline()
                if not n_atoms_line:
                    break
                if not n_atoms_line.strip():
                    continue

                try:
                    n_atoms = int(n_atoms_line.strip())
                except ValueError as exc:
                    raise AlignmentError(
                        f"invalid atom count at frame {frame_index}: {n_atoms_line.strip()}"
                    ) from exc
                if n_atoms <= 0:
                    raise AlignmentError(
                        f"frame {frame_index} has an invalid atom count: {n_atoms}"
                    )

                header_line = fin.readline()
                if not header_line:
                    raise AlignmentError(f"missing header for frame {frame_index}")
                header, _ = split_line_ending(header_line)

                energy_match = ENERGY_PATTERN.search(header)
                if energy_match is None:
                    raise AlignmentError(
                        f"missing energy= in the header of frame {frame_index}"
                    )
                try:
                    primary_energy = float(
                        re.search(
                            r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",
                            energy_match.group(0).split("=", 1)[1],
                        ).group(0)
                    )
                except (AttributeError, ValueError) as exc:
                    raise AlignmentError(
                        f"invalid energy= value in frame {frame_index}"
                    ) from exc

                element_counts: Counter = Counter()
                for atom_index in range(n_atoms):
                    atom_line = fin.readline()
                    if not atom_line:
                        raise AlignmentError(
                            f"frame {frame_index} ended before atom {atom_index + 1}"
                        )
                    atom_fields = atom_line.split()
                    if not atom_fields:
                        raise AlignmentError(
                            f"empty atom line at frame {frame_index}, atom {atom_index + 1}"
                        )
                    element_counts[atom_fields[0]] += 1

                frames.append(
                    FrameMetadata(
                        frame_index=frame_index,
                        n_atoms=n_atoms,
                        primary_energy=primary_energy,
                        config_type=parse_config_type(header),
                        element_counts=element_counts,
                    )
                )
                frame_index += 1
    except OSError as exc:
        raise AlignmentError(f"failed to read '{input_file}': {exc}") from exc

    if not frames:
        raise AlignmentError(f"no XYZ frames found in '{input_file}'")
    return frames


def calculate_nep_batch(
    batch_start: int,
    batch_end: int,
    nep_model_file: str,
    input_file: str,
) -> List[float]:
    """Calculate NEP energies for one contiguous frame batch."""
    from ase.io import read
    from calorine.calculators import CPUNEP

    atoms_batch = read(input_file, index=f"{batch_start}:{batch_end}")
    if not isinstance(atoms_batch, list):
        atoms_batch = [atoms_batch]

    calculator = CPUNEP(nep_model_file)
    energies = []
    for atoms in atoms_batch:
        atoms.calc = calculator
        energies.append(float(atoms.get_potential_energy()))
    return energies


def calculate_nep_energies(
    input_file: str,
    nep_model_file: str,
    n_frames: int,
    batch_size: int,
) -> np.ndarray:
    """Calculate NEP energies in ordered batches using worker processes."""
    from multiprocessing import Pool, cpu_count

    jobs = [
        (
            start,
            min(start + batch_size, n_frames),
            nep_model_file,
            input_file,
        )
        for start in range(0, n_frames, batch_size)
    ]
    process_count = min(cpu_count(), len(jobs))
    print(
        f" Calculating NEP energies in {len(jobs)} batch(es) "
        f"using {process_count} worker process(es)..."
    )

    try:
        with Pool(processes=process_count) as pool:
            batch_results = pool.starmap(calculate_nep_batch, jobs)
    except Exception as exc:
        raise AlignmentError(f"NEP energy calculation failed: {exc}") from exc

    energies = [energy for batch in batch_results for energy in batch]
    if len(energies) != n_frames:
        raise AlignmentError(
            f"NEP returned {len(energies)} energies for {n_frames} frames"
        )
    return np.asarray(energies, dtype=float)


def atomic_baseline_cost(
    population: np.ndarray,
    source_energies: np.ndarray,
    element_counts: np.ndarray,
    target_energies: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Evaluate the MSE used by the original NES optimizer."""
    shifted_energies = source_energies[None, :] - np.dot(
        population, element_counts.T
    )
    if target_energies is None:
        residual = shifted_energies
    else:
        residual = shifted_energies - target_energies[None, :]
    return np.mean(residual**2, axis=1).reshape(-1, 1)


def optimize_atomic_baseline(
    num_variables: int,
    max_generations: int,
    source_energies: np.ndarray,
    element_counts: np.ndarray,
    target_energies: Optional[np.ndarray],
    population_size: int,
    convergence_tol: float,
    random_seed: int,
    print_every: int = 100,
) -> np.ndarray:
    """Optimize atomic baselines using the original reproducible NES method."""
    rng = np.random.RandomState(random_seed)
    best_fitness_history = np.ones(max_generations)
    elite_solutions = np.zeros((max_generations, num_variables))
    mean = -1 * rng.rand(1, num_variables)
    stddev = 0.1 * np.ones((1, num_variables))
    learning_rate_mean = 1.0
    learning_rate_stddev = (
        (3 + np.log(num_variables)) / (5 * np.sqrt(num_variables)) / 2
    )
    selection_weights = np.maximum(
        0,
        np.log(population_size / 2 + 1)
        - np.log(np.arange(1, population_size + 1)),
    )
    selection_weights = (
        selection_weights / np.sum(selection_weights) - 1 / population_size
    )

    for generation in range(max_generations):
        z_samples = rng.randn(population_size, num_variables)
        population = mean + stddev * z_samples
        fitness = atomic_baseline_cost(
            population, source_energies, element_counts, target_energies
        )
        sorted_indices = np.argsort(fitness.flatten())
        fitness = fitness[sorted_indices]
        z_samples = z_samples[sorted_indices, :]
        population = population[sorted_indices, :]

        best_fitness_history[generation] = fitness[0, 0]
        elite_solutions[generation, :] = population[0, :]
        mean += learning_rate_mean * stddev * (
            selection_weights @ z_samples
        )
        stddev *= np.exp(
            learning_rate_stddev
            * (selection_weights @ (z_samples**2 - 1))
        )

        if generation % print_every == 0:
            print(
                f" Generation = {generation}, "
                f"best fitness = {fitness[0, 0]:.8f}"
            )
        if (
            generation > 0
            and abs(
                best_fitness_history[generation]
                - best_fitness_history[generation - 1]
            )
            < convergence_tol
        ):
            print(f" Converged at generation {generation}.")
            return elite_solutions[generation, :]

    return elite_solutions[-1, :]


def prompt_text(prompt: str, default: Optional[str] = None) -> str:
    """Read one line of text, handling EOF and optional defaults."""
    suffix = f" [{default}]" if default is not None else ""
    try:
        value = input(f" {prompt}{suffix}: ").strip()
    except EOFError as exc:
        raise AlignmentCancelled("input closed") from exc
    return value if value else (default or "")


def prompt_yes_no(prompt: str, default: bool = True) -> bool:
    """Read a yes/no response with a default."""
    default_text = "Y/n" if default else "y/N"
    while True:
        answer = prompt_text(f"{prompt} ({default_text})")
        if not answer:
            return default
        if answer.lower() in ("y", "yes"):
            return True
        if answer.lower() in ("n", "no"):
            return False
        print(" Please enter y or n.")


def prompt_positive_int(prompt: str, default: int) -> int:
    """Read a positive integer."""
    while True:
        value = prompt_text(prompt, str(default))
        try:
            parsed = int(value)
        except ValueError:
            parsed = 0
        if parsed > 0:
            return parsed
        print(" Please enter a positive integer.")


def prompt_positive_float(prompt: str, default: float) -> float:
    """Read a positive floating-point value."""
    while True:
        value = prompt_text(prompt, str(default))
        try:
            parsed = float(value)
        except ValueError:
            parsed = 0.0
        if parsed > 0:
            return parsed
        print(" Please enter a positive number.")


def choose_mode() -> str:
    """Select one of the three supported alignment modes."""
    print(" Alignment modes:")
    print("   1) Reference-group alignment")
    print("   2) Zero atomic baseline alignment")
    print("   3) DFT-to-NEP alignment")
    while True:
        choice = prompt_text("Select mode", "1")
        if choice == "1":
            return REF_GROUP_ALIGNMENT
        if choice == "2":
            return ZERO_BASELINE_ALIGNMENT
        if choice == "3":
            return DFT_TO_NEP_ALIGNMENT
        print(" Please select 1, 2, or 3.")


def interactive_settings() -> Tuple[AlignmentSettings, List[FrameMetadata]]:
    """Collect all settings through the Python interactive interface."""
    print(" +------------------------------------------------------+")
    print(" |             GPUMDkit Energy Reference Shifter        |")
    print(" +------------------------------------------------------+")
    print(" This tool changes only energy= values in the output XYZ file.")

    input_file = prompt_text("Input XYZ file", DEFAULT_INPUT_FILE)
    while not os.path.isfile(input_file):
        print(f" Error: input file '{input_file}' does not exist.")
        input_file = prompt_text("Input XYZ file")

    frames = parse_xyz_metadata(input_file)
    elements = sorted({element for frame in frames for element in frame.element_counts})
    groups = sorted({frame.config_type for frame in frames})
    print(f" Detected {len(frames)} frame(s).")
    print(f" Elements: {' '.join(elements)}")
    print(f" config_type groups: {', '.join(groups)}")

    output_file = prompt_text("Output XYZ file", DEFAULT_OUTPUT_FILE)
    if os.path.abspath(input_file) == os.path.abspath(output_file):
        raise AlignmentError("input and output files must be different")
    output_directory = os.path.dirname(os.path.abspath(output_file))
    if not os.path.isdir(output_directory):
        raise AlignmentError(f"output directory '{output_directory}' does not exist")
    if os.path.exists(output_file) and not prompt_yes_no(
        f"Output file '{output_file}' exists. Overwrite", default=False
    ):
        raise AlignmentCancelled("output overwrite declined")

    mode = choose_mode()
    settings = AlignmentSettings(
        input_file=input_file,
        output_file=output_file,
        mode=mode,
        shift_groups=[],
    )

    if mode == REF_GROUP_ALIGNMENT:
        settings.reference_group = prompt_text(
            "Reference group", DEFAULT_REFERENCE_GROUP
        )
        while settings.reference_group not in groups:
            print(
                f" Error: reference group '{settings.reference_group}' "
                "was not found in the input file."
            )
            settings.reference_group = prompt_text("Reference group")
        raw_shift_groups = prompt_text(
            "Shift groups, separated by spaces (blank = all non-reference groups)"
        )
        settings.shift_groups = raw_shift_groups.split() if raw_shift_groups else []
    elif mode == DFT_TO_NEP_ALIGNMENT:
        settings.nep_model_file = prompt_text(
            "NEP model file", DEFAULT_NEP_MODEL_FILE
        )
        if not os.path.isfile(settings.nep_model_file):
            raise AlignmentError(
                f"NEP model file '{settings.nep_model_file}' does not exist"
            )

    if prompt_yes_no("Change advanced NES settings", default=False):
        settings.max_generations = prompt_positive_int(
            "Maximum generations", DEFAULT_MAX_GENERATIONS
        )
        settings.population_size = prompt_positive_int(
            "Population size", DEFAULT_POPULATION_SIZE
        )
        settings.convergence_tol = prompt_positive_float(
            "Convergence tolerance", DEFAULT_CONVERGENCE_TOL
        )
        settings.random_seed = prompt_positive_int(
            "Random seed", DEFAULT_RANDOM_SEED
        )
        if mode == DFT_TO_NEP_ALIGNMENT:
            settings.nep_batch_size = prompt_positive_int(
                "NEP batch size", DEFAULT_NEP_BATCH_SIZE
            )

    print(" Using NES settings:")
    print(f"   max generations: {settings.max_generations}")
    print(f"   population size: {settings.population_size}")
    print(f"   convergence tolerance: {settings.convergence_tol}")
    print(f"   random seed: {settings.random_seed}")
    if mode == DFT_TO_NEP_ALIGNMENT:
        print(f"   NEP batch size: {settings.nep_batch_size}")
    return settings, frames


def optimize_groups(
    frames: List[FrameMetadata],
    settings: AlignmentSettings,
) -> Dict[str, np.ndarray]:
    """Optimize one atomic baseline for each selected config_type group."""
    all_elements = sorted(
        {element for frame in frames for element in frame.element_counts}
    )
    all_groups = sorted({frame.config_type for frame in frames})
    nep_energies: Optional[np.ndarray] = None

    if settings.mode == REF_GROUP_ALIGNMENT:
        reference_frames = [
            frame for frame in frames if frame.config_type == settings.reference_group
        ]
        if not reference_frames:
            raise AlignmentError(
                f"reference group '{settings.reference_group}' was not found"
            )
        groups_to_process = (
            [
                group
                for group in all_groups
                if group != settings.reference_group
            ]
            if not settings.shift_groups
            else [
                group
                for group in settings.shift_groups
                if group != settings.reference_group
            ]
        )
        reference_mean = float(
            np.mean([frame.primary_energy for frame in reference_frames])
        )
        print(
            f" Reference group '{settings.reference_group}' mean energy: "
            f"{reference_mean:.8f}"
        )
    elif settings.mode == ZERO_BASELINE_ALIGNMENT:
        groups_to_process = all_groups
        reference_mean = None
    elif settings.mode == DFT_TO_NEP_ALIGNMENT:
        groups_to_process = all_groups
        reference_mean = None
        nep_energies = calculate_nep_energies(
            settings.input_file,
            settings.nep_model_file,
            len(frames),
            settings.nep_batch_size,
        )
    else:
        raise AlignmentError(f"unsupported alignment mode '{settings.mode}'")

    optimized_baselines: Dict[str, np.ndarray] = {}
    for group_name in groups_to_process:
        group_frames = [frame for frame in frames if frame.config_type == group_name]
        if not group_frames:
            print(
                f" Warning: no frames found for config_type '{group_name}'. "
                "Skipping."
            )
            continue

        source_energies = np.asarray(
            [frame.primary_energy for frame in group_frames], dtype=float
        )
        element_counts = np.asarray(
            [
                [frame.element_counts.get(element, 0) for element in all_elements]
                for frame in group_frames
            ],
            dtype=float,
        )

        if settings.mode == REF_GROUP_ALIGNMENT:
            target_energies = np.full_like(source_energies, reference_mean)
            print(
                f" Optimizing '{group_name}' against "
                f"'{settings.reference_group}'..."
            )
        elif settings.mode == ZERO_BASELINE_ALIGNMENT:
            target_energies = None
            print(f" Optimizing '{group_name}' against zero...")
        else:
            target_energies = np.asarray(
                [nep_energies[frame.frame_index] for frame in group_frames],
                dtype=float,
            )
            print(f" Optimizing '{group_name}' against NEP energies...")

        optimized_baselines[group_name] = optimize_atomic_baseline(
            len(all_elements),
            settings.max_generations,
            source_energies,
            element_counts,
            target_energies,
            settings.population_size,
            settings.convergence_tol,
            settings.random_seed,
        )
        print(f" Optimized atomic baselines for '{group_name}':")
        for element, value in zip(all_elements, optimized_baselines[group_name]):
            print(f"   {element}: {value:.8f}")

    return optimized_baselines


def rewrite_aligned_xyz(
    settings: AlignmentSettings,
    frames: List[FrameMetadata],
    optimized_baselines: Dict[str, np.ndarray],
) -> None:
    """Write a new XYZ file while preserving atom lines and other headers."""
    all_elements = sorted(
        {element for frame in frames for element in frame.element_counts}
    )
    frame_index = 0

    try:
        with open(settings.input_file, "r", encoding="utf-8", newline="") as fin, open(
            settings.output_file, "w", encoding="utf-8", newline=""
        ) as fout:
            while frame_index < len(frames):
                n_atoms_line = fin.readline()
                if not n_atoms_line:
                    raise AlignmentError(
                        f"input ended while writing frame {frame_index}"
                    )
                if not n_atoms_line.strip():
                    fout.write(n_atoms_line)
                    continue

                header_line = fin.readline()
                if not header_line:
                    raise AlignmentError(
                        f"input ended before the header of frame {frame_index}"
                    )
                header, line_ending = split_line_ending(header_line)

                frame = frames[frame_index]
                shifted_energy = frame.primary_energy
                baseline = optimized_baselines.get(frame.config_type)
                if baseline is not None:
                    counts = np.asarray(
                        [frame.element_counts.get(element, 0) for element in all_elements],
                        dtype=float,
                    )
                    shifted_energy -= float(np.dot(counts, baseline))

                new_header = ENERGY_PATTERN.sub(
                    lambda match: f"{match.group(1)}{shifted_energy:.8f}",
                    header,
                    count=1,
                )
                fout.write(n_atoms_line)
                fout.write(new_header + line_ending)
                for _ in range(frame.n_atoms):
                    atom_line = fin.readline()
                    if not atom_line:
                        raise AlignmentError(
                            f"input ended while writing frame {frame_index}"
                        )
                    fout.write(atom_line)
                frame_index += 1

            for remaining_line in fin:
                fout.write(remaining_line)
    except OSError as exc:
        raise AlignmentError(f"failed to write '{settings.output_file}': {exc}") from exc


def run_interactive() -> None:
    """Run one complete interactive alignment operation."""
    settings, frames = interactive_settings()
    optimized_baselines = optimize_groups(frames, settings)
    print(f" Writing aligned energies to '{settings.output_file}'...")
    rewrite_aligned_xyz(settings, frames, optimized_baselines)
    print(" Energy alignment completed.")
    print(f" Output written to '{settings.output_file}'.")


def main() -> int:
    """Parse the help option and launch the interactive interface."""
    parser = build_parser()
    parser.parse_args()
    try:
        run_interactive()
    except AlignmentCancelled as exc:
        print(f" Operation canceled: {exc}.")
        return 0
    except AlignmentError as exc:
        print(f" Error: {exc}.")
        return 1
    except KeyboardInterrupt:
        print(" Operation canceled by user.")
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())

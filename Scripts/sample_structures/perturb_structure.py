"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     perturb_structure.py
Category:   Sample Structure Scripts
Purpose:    Generate perturbed structures from a POSCAR/CONTCAR file using
            dpdata. Supports cell and atom perturbations with different
            perturbation styles (normal, uniform, const) and writes a concise
            statistical report of cell changes and atomic displacements.
Usage:      gpumdkit.sh
            choose 204) Perturb structure
            python perturb_structure.py <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>
Arguments:
  input.vasp           Input POSCAR/CONTCAR file
  pert_num             Number of perturbed structures to generate
  cell_pert_fraction   Fraction of cell perturbation
  atom_pert_distance   Distance of atom perturbation (Angstrom)
  atom_pert_style      Style: normal, uniform, or const
Output:
  POSCAR_*.vasp  (perturbed structures in VASP POSCAR format)
  perturb.txt    (statistical perturbation report)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-08-12
=============================================================================
"""

import argparse
import os
import sys

import numpy as np


REPORT_FILE = "perturb.txt"


def print_dependency_notice():
    print(" This function requires the dpdata package.")
    print(" If you use this function, we recommend citing dpdata according to its official documentation.")


def print_usage():
    print(" Usage: gpumdkit.sh")
    print("        choose 204) Perturb structure")
    print("    or: python perturb_structure.py <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>")
    print("")
    print(" Arguments:")
    print("   input.vasp           Input POSCAR/CONTCAR file")
    print("   pert_num             Number of perturbed structures to generate")
    print("   cell_pert_fraction   Fraction of cell perturbation")
    print("   atom_pert_distance   Distance of atom perturbation (Angstrom)")
    print("   atom_pert_style      Style: normal, uniform, or const")
    print("")
    print(" Output:")
    print("   POSCAR_*.vasp")
    print(f"   {REPORT_FILE}")
    print("")
    print(" Example: in interactive mode, enter: POSCAR 20 0.03 0.2 uniform")
    print("          python perturb_structure.py POSCAR 20 0.03 0.2 uniform")
    print("")


def parse_args():
    parser = argparse.ArgumentParser(
             formatter_class=argparse.RawDescriptionHelpFormatter,
             description="""Generate perturbed structures from a POSCAR/CONTCAR file.

Usage:
  gpumdkit.sh
  choose 204) Perturb structure
  python perturb_structure.py <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>

Example:
  in interactive mode, enter: POSCAR 20 0.03 0.2 uniform
  python perturb_structure.py POSCAR 20 0.03 0.2 uniform
""")
    parser.add_argument('input_file', help='The path to POSCAR/CONTCAR file')
    parser.add_argument('pert_num', type=int, default=20, help='The perturbation number')
    parser.add_argument('cell_pert_fraction', type=float, default=0.03, help='The fraction of cell perturbation')
    parser.add_argument('atom_pert_distance', type=float, default=0.2, help='The distance of atom perturbation')
    parser.add_argument('atom_pert_style', type=str, default='uniform', choices=['normal', 'uniform', 'const'], help='The style for atom perturbation')
    return parser.parse_args()


def vector_angle(vector_1, vector_2):
    """Return the angle between two vectors in degrees."""
    denominator = np.linalg.norm(vector_1) * np.linalg.norm(vector_2)
    cosine = np.dot(vector_1, vector_2) / denominator
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def cell_metrics(cell):
    """Return cell lengths, angles, and volume."""
    lengths = np.linalg.norm(cell, axis=1)
    angles = np.array([
        vector_angle(cell[1], cell[2]),
        vector_angle(cell[0], cell[2]),
        vector_angle(cell[0], cell[1]),
    ])
    volume = abs(np.linalg.det(cell))
    return lengths, angles, volume


def cell_strain_magnitude(reference_cell, perturbed_cell):
    """Return the rotation-invariant Green-Lagrange strain magnitude."""
    deformation = np.linalg.solve(reference_cell, perturbed_cell)
    strain = 0.5 * (np.matmul(deformation, deformation.T) - np.eye(3))
    return np.linalg.norm(strain)


def minimum_distance(coords, cell, distance_function):
    """Return the shortest periodic distance and its atom pair."""
    if distance_function is None or len(coords) < 2:
        return None, None
    shortest_distance = np.inf
    shortest_pair = None
    for atom_i in range(len(coords) - 1):
        _, distances = distance_function(
            coords[atom_i], coords[atom_i + 1:], cell=cell, pbc=True
        )
        local_atom = int(np.argmin(distances[0]))
        local_distance = float(distances[0, local_atom])
        if local_distance < shortest_distance:
            shortest_distance = local_distance
            shortest_pair = (atom_i, atom_i + 1 + local_atom)
    return shortest_distance, shortest_pair


def summary(values):
    """Return minimum, mean, and maximum values."""
    values = np.asarray(values, dtype=float)
    return np.min(values), np.mean(values), np.max(values)


def histogram_lines(values, atom_pert_distance, atom_pert_style, bins=10):
    """Build an ASCII histogram for atomic displacement magnitudes."""
    observed_max = float(np.max(values))
    if observed_max == 0.0:
        return ["0.0000          | {:7d} 100.00% ########################################".format(len(values))]

    if atom_pert_style in ("uniform", "const") and atom_pert_distance > 0:
        upper = max(atom_pert_distance, observed_max)
    else:
        upper = observed_max

    counts, edges = np.histogram(values, bins=np.linspace(0.0, upper, bins + 1))
    max_count = max(int(np.max(counts)), 1)
    total = len(values)
    lines = []
    for index, count in enumerate(counts):
        bar_length = int(round(40 * count / max_count)) if count else 0
        percentage = 100.0 * count / total
        lines.append(
            f"{edges[index]:7.4f} - {edges[index + 1]:7.4f} | "
            f"{count:7d} {percentage:6.2f}% {'#' * bar_length}"
        )
    return lines


def collect_statistics(system, perturbed_systems, output_files, distance_function):
    """Collect aggregate cell, displacement, and distance statistics."""
    reference_cell = np.asarray(system["cells"][0], dtype=float)
    reference_coords = np.asarray(system["coords"][0], dtype=float)
    reference_fractional = np.linalg.solve(reference_cell.T, reference_coords.T).T
    reference_lengths, reference_angles, reference_volume = cell_metrics(reference_cell)

    atom_names = system["atom_names"]
    atom_types = np.asarray(system["atom_types"], dtype=int)
    symbols = [atom_names[atom_type] for atom_type in atom_types]

    length_changes = []
    angle_changes = []
    volume_changes = []
    strain_magnitudes = []
    displacement_frames = []
    minimum_distances = []
    minimum_pairs = []

    distance_error = None
    if distance_function is None:
        reference_minimum, reference_pair = None, None
        distance_error = "ASE is not installed."
    elif len(reference_coords) < 2:
        reference_minimum, reference_pair = None, None
        distance_error = "the structure contains fewer than two atoms."
    else:
        try:
            reference_minimum, reference_pair = minimum_distance(
                reference_coords, reference_cell, distance_function
            )
        except Exception as exc:
            reference_minimum, reference_pair = None, None
            distance_error = str(exc)

    for frame_index in range(len(perturbed_systems)):
        cell = np.asarray(perturbed_systems["cells"][frame_index], dtype=float)
        coords = np.asarray(perturbed_systems["coords"][frame_index], dtype=float)
        lengths, angles, volume = cell_metrics(cell)

        length_changes.append(100.0 * (lengths / reference_lengths - 1.0))
        angle_changes.append(angles - reference_angles)
        volume_changes.append(100.0 * (volume / reference_volume - 1.0))
        strain_magnitudes.append(cell_strain_magnitude(reference_cell, cell))

        affine_coords = np.matmul(reference_fractional, cell)
        displacement_frames.append(coords - affine_coords)

        if distance_function is not None and distance_error is None:
            try:
                distance, pair = minimum_distance(coords, cell, distance_function)
            except Exception as exc:
                distance_error = str(exc)
                minimum_distances = []
                minimum_pairs = []
            else:
                minimum_distances.append(distance)
                minimum_pairs.append(pair)

    displacement_frames = np.asarray(displacement_frames)
    displacement_magnitudes = np.linalg.norm(displacement_frames, axis=2)
    max_frame, max_atom = np.unravel_index(
        np.argmax(displacement_magnitudes), displacement_magnitudes.shape
    )

    return {
        "reference_lengths": reference_lengths,
        "reference_angles": reference_angles,
        "reference_volume": reference_volume,
        "reference_minimum": reference_minimum,
        "reference_pair": reference_pair,
        "symbols": symbols,
        "length_changes": np.asarray(length_changes),
        "angle_changes": np.asarray(angle_changes),
        "volume_changes": np.asarray(volume_changes),
        "strain_magnitudes": np.asarray(strain_magnitudes),
        "displacement_magnitudes": displacement_magnitudes,
        "max_frame": max_frame,
        "max_atom": max_atom,
        "minimum_distances": np.asarray(minimum_distances),
        "minimum_pairs": minimum_pairs,
        "distance_error": distance_error,
        "output_files": output_files,
    }


def format_change_row(label, values):
    """Format one minimum/mean/maximum report row."""
    minimum, mean, maximum = summary(values)
    return f"{label:<27}{minimum:>11.4f}{mean:>12.4f}{maximum:>12.4f}"


def write_report(args, dpdata_version, system, statistics):
    """Write the concise perturbation report."""
    length_changes = statistics["length_changes"]
    angle_changes = statistics["angle_changes"]
    volume_changes = statistics["volume_changes"]
    strain_magnitudes = statistics["strain_magnitudes"]
    displacement_magnitudes = statistics["displacement_magnitudes"]
    flat_displacements = displacement_magnitudes.reshape(-1)
    output_files = statistics["output_files"]

    contraction_index = int(np.argmin(volume_changes))
    expansion_index = int(np.argmax(volume_changes))
    strain_index = int(np.argmax(strain_magnitudes))
    max_frame = statistics["max_frame"]
    max_atom = statistics["max_atom"]

    lines = [
        "GPUMDkit Perturbation Report",
        "============================",
        "",
        f"Input structure : {args.input_file}",
        f"Structures      : {args.pert_num}",
        f"Atoms/structure : {system.get_natoms()}",
        f"Elements        : {' '.join(system['atom_names'])}",
        "",
        "Perturbation settings",
        "---------------------",
        f"Cell fraction   : {args.cell_pert_fraction:.6f}",
        f"Atom distance   : {args.atom_pert_distance:.6f} Angstrom",
        f"Atom style      : {args.atom_pert_style}",
        f"dpdata version  : {dpdata_version}",
        "",
        "Reference cell",
        "--------------",
        f"Lengths (A)     : {statistics['reference_lengths'][0]:.6f} "
        f"{statistics['reference_lengths'][1]:.6f} {statistics['reference_lengths'][2]:.6f}",
        f"Angles (degree) : {statistics['reference_angles'][0]:.6f} "
        f"{statistics['reference_angles'][1]:.6f} {statistics['reference_angles'][2]:.6f}",
        f"Volume (A^3)    : {statistics['reference_volume']:.6f}",
        "",
        "Cell changes",
        "------------",
        f"{'Quantity':<27}{'Minimum':>11}{'Mean':>12}{'Maximum':>12}",
        format_change_row("a change (%)", length_changes[:, 0]),
        format_change_row("b change (%)", length_changes[:, 1]),
        format_change_row("c change (%)", length_changes[:, 2]),
        format_change_row("volume change (%)", volume_changes),
        format_change_row("alpha change (degree)", angle_changes[:, 0]),
        format_change_row("beta change (degree)", angle_changes[:, 1]),
        format_change_row("gamma change (degree)", angle_changes[:, 2]),
        "",
        f"Largest contraction : {output_files[contraction_index]} "
        f"({volume_changes[contraction_index]:+.4f}%)",
        f"Largest expansion   : {output_files[expansion_index]} "
        f"({volume_changes[expansion_index]:+.4f}%)",
        f"Largest cell strain : {output_files[strain_index]}",
        "",
        "Atomic displacements",
        "--------------------",
        f"Mean displacement   : {np.mean(flat_displacements):.6f} Angstrom",
        f"Median displacement : {np.median(flat_displacements):.6f} Angstrom",
        f"RMS displacement    : {np.sqrt(np.mean(flat_displacements ** 2)):.6f} Angstrom",
        f"95th percentile     : {np.percentile(flat_displacements, 95):.6f} Angstrom",
        f"Maximum displacement: {np.max(flat_displacements):.6f} Angstrom",
        f"Maximum location    : {output_files[max_frame]}, "
        f"{statistics['symbols'][max_atom]} atom {max_atom + 1}",
        "",
        "Distribution (Angstrom)",
    ]
    lines.extend(
        histogram_lines(
            flat_displacements,
            args.atom_pert_distance,
            args.atom_pert_style,
        )
    )

    lines.extend(["", "Minimum distances", "-----------------"])
    minimum_distances = statistics["minimum_distances"]
    if statistics["reference_minimum"] is None or minimum_distances.size == 0:
        lines.append(f"Not calculated: {statistics['distance_error']}")
    else:
        shortest_index = int(np.argmin(minimum_distances))
        shortest_pair = statistics["minimum_pairs"][shortest_index]
        atom_i, atom_j = shortest_pair
        lines.extend([
            f"Original structure : {statistics['reference_minimum']:.6f} Angstrom",
            f"Generated minimum  : {np.min(minimum_distances):.6f} Angstrom",
            f"Generated mean     : {np.mean(minimum_distances):.6f} Angstrom",
            f"Generated maximum  : {np.max(minimum_distances):.6f} Angstrom",
            f"Shortest structure : {output_files[shortest_index]}",
            f"Shortest pair      : {statistics['symbols'][atom_i]} atom {atom_i + 1} - "
            f"{statistics['symbols'][atom_j]} atom {atom_j + 1}",
        ])

    lines.extend([
        "",
        "Notes",
        "-----",
        "Atomic displacements exclude affine motion caused by cell deformation.",
        "No structural acceptance threshold was applied.",
        "",
    ])

    with open(REPORT_FILE, "w", encoding="utf-8") as report_file:
        report_file.write("\n".join(lines))

def main():
    args_in = sys.argv[1:]
    if len(args_in) != 5 and not (args_in and args_in[0] in ("-h", "--help")):
        print_usage()
        sys.exit(1)

    try:
        args = parse_args()
    except SystemExit as exc:
        if exc.code == 0:
            sys.exit(0)
        print(f"Default values: pert_num=20, cell_pert_fraction=0.03, atom_pert_distance=0.2, atom_pert_style=uniform")
        print("atom_pert_style options: 'normal', 'uniform', 'const'")
        print("dpdata documentation: https://docs.deepmodeling.com/projects/dpdata/en/master/index.html \n")
        sys.exit(exc.code)

    print_dependency_notice()

    if not os.path.isfile(args.input_file):
        print(f" Error: file '{args.input_file}' does not exist.")
        sys.exit(1)
    if args.pert_num <= 0:
        print(" Error: pert_num must be a positive integer.")
        sys.exit(1)
    if not np.isfinite(args.cell_pert_fraction) or args.cell_pert_fraction < 0:
        print(" Error: cell_pert_fraction must be a finite non-negative number.")
        sys.exit(1)
    if not np.isfinite(args.atom_pert_distance) or args.atom_pert_distance < 0:
        print(" Error: atom_pert_distance must be a finite non-negative number.")
        sys.exit(1)

    try:
        import dpdata
    except ImportError:
        print(" Error: dpdata is not installed or cannot be imported.")
        print(" Please install dpdata before using this function.")
        sys.exit(1)

    try:
        from ase.geometry import get_distances
    except ImportError:
        get_distances = None

    # Read the POSCAR file and perform perturbation
    try:
        system = dpdata.System(args.input_file, fmt='vasp/poscar')
        perturbed_systems = system.perturb(pert_num=args.pert_num,
                                           cell_pert_fraction=args.cell_pert_fraction,
                                           atom_pert_distance=args.atom_pert_distance,
                                           atom_pert_style=args.atom_pert_style,)
    except Exception as exc:
        print(f" Error: failed to generate perturbed structures: {exc}")
        sys.exit(1)

    # Save the perturbed structures
    pert_num = args.pert_num
    width = len(str(pert_num))
    output_files = []
    for i in range(pert_num):
        output_file = f"POSCAR_{str(i + 1).zfill(width)}.vasp"
        try:
            perturbed_systems.sub_system(i).to('vasp/poscar', output_file)
        except Exception as exc:
            print(f" Error: failed to write '{output_file}': {exc}")
            sys.exit(1)
        output_files.append(output_file)

    try:
        statistics = collect_statistics(
            system, perturbed_systems, output_files, get_distances
        )
        write_report(
            args,
            getattr(dpdata, "__version__", "unknown"),
            system,
            statistics,
        )
    except Exception as exc:
        print(f" Error: failed to write perturbation report: {exc}")
        sys.exit(1)

    print(f" Generated {pert_num} perturbed structures.")
    print(f" Saved perturbation report to {REPORT_FILE}.")

if __name__ == "__main__":
     main()

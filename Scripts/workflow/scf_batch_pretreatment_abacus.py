"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     scf_batch_pretreatment_abacus.py
Category:   Workflow Scripts
Purpose:    Batch pretreatment of extxyz structures for ABACUS SCF
            calculations. Writes one STRU file per frame, matches each
            element to a pseudopotential and orbital file by
            filename prefix, and links the shared fp/INPUT and fp/KPT.
Usage:      python scf_batch_pretreatment_abacus.py <extxyz_file> <prefix> <basis_dir>
            (invoked interactively via gpumdkit.sh menu 301 -> 3)
Arguments:
  extxyz_file  Input extxyz file with lattice information
  prefix       Prefix for output directories (<prefix>_<index>)
  basis_dir    Directory with exactly one .upf and one .orb per element
Output:
  <prefix>_<index>/STRU          ABACUS structure file
  <prefix>_<index>/INPUT, KPT    Links to ../fp/INPUT and ../fp/KPT
  fp/                            Directory for the shared INPUT and KPT
  presub.sh                      Example loop for running ABACUS
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-11-17
=============================================================================
"""

import os
import sys

args = sys.argv[1:]
if len(args) != 3 or args[0] in ("-h", "--help"):
    print(" Usage: python scf_batch_pretreatment_abacus.py <extxyz_file> <prefix> <basis_dir>")
    print("  (invoked interactively via gpumdkit.sh menu 301 -> 3)")
    print("")
    print(" Arguments:")
    print("   extxyz_file  Input extxyz file with lattice information")
    print("   prefix       Prefix for output directories (<prefix>_<index>)")
    print("   basis_dir    Directory with exactly one .upf and one .orb per element")
    print("")
    print(" Example: python scf_batch_pretreatment_abacus.py train.xyz LiF_iter01 ./basis")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

from ase.data import atomic_masses, atomic_numbers
from ase.io import read

# ABACUS reads LATTICE_CONSTANT in Bohr; this value makes the lattice vectors
# and Cartesian positions below Angstrom.
ANGSTROM_TO_BOHR = 1.889726125457828


def match_element_file(directory, element, extensions):
    """Return candidate files whose names start with the element symbol."""
    candidates = []
    for name in sorted(os.listdir(directory)):
        if not os.path.isfile(os.path.join(directory, name)):
            continue
        if not name.lower().endswith(extensions):
            continue
        if not name.startswith(element):
            continue
        # Separate "S" from "Si", "B" from "Bi": the next character must not
        # continue the element symbol.
        rest = name[len(element):]
        if rest and rest[0].islower():
            continue
        candidates.append(name)
    return candidates


def build_species_map(elements, directory, extensions, label):
    """Map every element to exactly one file, or exit with a clear message."""
    mapping = {}
    problems = []
    for element in elements:
        candidates = match_element_file(directory, element, extensions)
        if len(candidates) == 1:
            mapping[element] = candidates[0]
        elif not candidates:
            problems.append(f"   {element}: no {label} file found")
        else:
            problems.append(f"   {element}: several {label} files found: {', '.join(candidates)}")
    if problems:
        print(f" Error: could not assign one {label} file per element in '{directory}'.")
        for line in problems:
            print(line)
        print(f" Keep exactly one {label} file per element whose name starts with the")
        print(" element symbol (for example Li_ONCV_PBE-1.0.upf or Li_gga_8au_100Ry_4s1p.orb).")
        sys.exit(1)
    return mapping


def species_in_order(atoms):
    """Return element symbols in order of first appearance."""
    order = []
    for symbol in atoms.get_chemical_symbols():
        if symbol not in order:
            order.append(symbol)
    return order


def write_stru(path, atoms, pp_map, orb_map):
    """Write one ABACUS STRU file with Angstrom lattice vectors and positions."""
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    species = species_in_order(atoms)
    lines = ["ATOMIC_SPECIES"]
    for element in species:
        mass = atomic_masses[atomic_numbers[element]]
        lines.append(f"{element} {mass:.6f} {pp_map[element]}")
    lines += ["", "NUMERICAL_ORBITAL"]
    lines += [orb_map[element] for element in species]
    lines += ["", "LATTICE_CONSTANT", f"{ANGSTROM_TO_BOHR:.15f}", "", "LATTICE_VECTORS"]
    for vector in atoms.get_cell():
        lines.append(" ".join(f"{value:18.10f}" for value in vector))
    lines += ["", "ATOMIC_POSITIONS", "Cartesian    # Angstrom (LATTICE_CONSTANT = 1 Angstrom)"]
    for element in species:
        indices = [i for i, symbol in enumerate(symbols) if symbol == element]
        lines += [element, "0.0", str(len(indices))]
        for i in indices:
            coord = " ".join(f"{value:18.10f}" for value in positions[i])
            lines.append(f"{coord} 1 1 1")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def replace_link(target, link_path):
    """Create or refresh a relative symbolic link."""
    if os.path.islink(link_path) or os.path.exists(link_path):
        os.remove(link_path)
    os.symlink(target, link_path)


def main():
    """Validate inputs, write STRU directories, and print next steps."""
    extxyz_file, prefix, basis_dir = args

    if not os.path.isfile(extxyz_file):
        print(f" Error: file '{extxyz_file}' does not exist.")
        sys.exit(1)
    if not os.path.isdir(basis_dir):
        print(f" Error: basis directory '{basis_dir}' does not exist.")
        sys.exit(1)

    try:
        frames = read(extxyz_file, index=":")
    except Exception as error:
        print(f" Error: failed to read {extxyz_file}: {error}")
        sys.exit(1)
    if not isinstance(frames, list):
        frames = [frames]
    if not frames:
        print(f" Error: no structures found in {extxyz_file}.")
        sys.exit(1)

    elements = []
    for index, atoms in enumerate(frames, start=1):
        if atoms.cell.rank != 3 or atoms.get_volume() <= 0:
            print(f" Error: structure {index} has no valid three-dimensional cell.")
            sys.exit(1)
        for element in species_in_order(atoms):
            if element not in elements:
                elements.append(element)

    pp_map = build_species_map(elements, basis_dir, (".upf",), "pseudopotential (.upf)")
    orb_map = build_species_map(elements, basis_dir, (".orb",), "orbital (.orb)")

    print(f" Found {len(frames)} structures with elements: {' '.join(elements)}")
    for element in elements:
        print(f"   {element:<3} {pp_map[element]}  {orb_map[element]}")

    os.makedirs("fp", exist_ok=True)
    for index, atoms in enumerate(frames, start=1):
        dir_name = f"{prefix}_{index}"
        os.makedirs(dir_name, exist_ok=True)
        write_stru(os.path.join(dir_name, "STRU"), atoms, pp_map, orb_map)
        replace_link("../fp/INPUT", os.path.join(dir_name, "INPUT"))
        replace_link("../fp/KPT", os.path.join(dir_name, "KPT"))

    with open("presub.sh", "w") as handle:
        handle.write(
            "#!/bin/bash\n\n"
            "# You can cat it to your submit script.\n\n"
            f"for dir in {prefix}_*; do\n"
            "    cd $dir\n"
            "    echo \"Running ABACUS SCF in $dir...\"\n"
            "    mpirun -n X abacus > log\n"
            "    cd ..\n"
            "done\n"
        )
    os.chmod("presub.sh", 0o755)

    print(f"\n Prepared {len(frames)} ABACUS calculation directories with prefix '{prefix}'.")
    print(" Before running ABACUS, provide your own fp/INPUT and fp/KPT (unless using kspacing).")
    print(f" Set pseudo_dir and orbital_dir in fp/INPUT to {os.path.abspath(basis_dir)}")


if __name__ == "__main__":
    main()

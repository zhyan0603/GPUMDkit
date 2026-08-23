#!/usr/bin/env python3
"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     xyz2dp.py
Category:   Format Conversion Scripts
Purpose:    Convert a labeled extended XYZ file to a DeepMD-kit npy dataset.
            The user supplies the element order used by the DeepMD type map.
Usage:      gpumdkit.sh -xyz2dp
            python3 xyz2dp.py <input.xyz> <type1> <type2> ...
Arguments:
  input.xyz     Labeled extended XYZ input containing energy and force data
  type1 ...      Element order for the DeepMD type_map.raw file
Output:
  deepmd_data    Fixed output directory; grouped by composition when needed
Dependencies:
  dpdata       DeepMD-kit data conversion library
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-08-23
=============================================================================
"""

import os
import re
import sys
import tempfile

OUTPUT_DIR = "deepmd_data"


def print_usage():
    """Print the standard GPUMDkit and direct-Python usage text."""
    print(" Usage: gpumdkit.sh -xyz2dp")
    print("    or: python3 xyz2dp.py <input.xyz> <type1> <type2> ...")
    print("")
    print(" Arguments:")
    print("   input.xyz  Labeled extended XYZ input containing energy and force data")
    print("   type1 ...  Element order for the DeepMD type_map.raw file")
    print("")
    print(" Example: python3 xyz2dp.py train.xyz Li P S")
    print("")


def normalize_extxyz_header(line):
    """Normalize common GPUMDkit extxyz header spellings for dpdata."""
    key_names = {
        "Energy": "energy",
        "energy": "energy",
        "Virial": "virial",
        "virial": "virial",
        "Lattice": "Lattice",
        "lattice": "Lattice",
        "Properties": "Properties",
        "properties": "Properties",
    }
    for source, target in key_names.items():
        line = re.sub(
            rf"(?<![\w]){source}\s*=",
            f"{target}=",
            line,
        )

    properties_match = re.search(
        r"(?<![\w])Properties=(\"[^\"]*\"|'[^']*'|\S+)", line
    )
    if properties_match:
        properties_value = properties_match.group(1)
        quote = properties_value[0] if properties_value[:1] in ("'", '\"') else ""
        if quote:
            properties_value = properties_value[1:-1]
        properties_value = re.sub(
            r"(?<=:)forces(?=:|$)",
            "force",
            properties_value,
        )
        if quote:
            properties_value = f"{quote}{properties_value}{quote}"
        line = (
            line[: properties_match.start(1)]
            + properties_value
            + line[properties_match.end(1) :]
        )

    return line


def write_normalized_xyz(input_file):
    """Write a temporary dpdata-compatible copy and return its path."""
    temporary_file = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        suffix=".xyz",
        prefix="gpumdkit_xyz2dp_",
        delete=False,
    )
    temporary_path = temporary_file.name
    try:
        with temporary_file, open(input_file, encoding="utf-8") as source:
            for line in source:
                temporary_file.write(normalize_extxyz_header(line))
    except Exception:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)
        raise
    return temporary_path


def load_dpdata():
    """Import dpdata after command-line help and basic checks are complete."""
    try:
        from dpdata.system import MultiSystems
    except ImportError:
        print(" Error: the 'dpdata' package is not installed.")
        print(" Please install dpdata before using this conversion.")
        sys.exit(1)
    return MultiSystems


def validate_type_map(type_map):
    """Validate the user-provided DeepMD element order."""
    if not type_map or any(not item.strip() for item in type_map):
        print(" Error: type-map must contain at least one element symbol.")
        sys.exit(1)
    if len(type_map) != len(set(type_map)):
        print(" Error: type-map contains duplicate element symbols.")
        sys.exit(1)


def print_system_summary(multi_systems, type_map):
    """Print the parsed systems and the exact type-map used for output."""
    total_frames = sum(system.get_nframes() for system in multi_systems.systems.values())
    print(f" Type map: {' '.join(type_map)}")
    print(f" Found {len(multi_systems.systems)} composition system(s).")
    for name, system in multi_systems.systems.items():
        species = [str(item) for item in system.get_atom_names()]
        print(
            f"  {name}: frames={system.get_nframes()}, "
            f"atoms={system.get_natoms()}, species={species}"
        )
    print(f" Total frames: {total_frames}")


def main(args):
    """Convert the input extxyz file to DeepMD-kit npy format."""
    if not os.path.isfile(args[0]):
        print(f" Error: input file '{args[0]}' does not exist.")
        return 1

    input_file = args[0]
    type_map = args[1:]
    validate_type_map(type_map)
    print(" This function requires the 'dpdata' package.")

    MultiSystems = load_dpdata()
    temporary_path = None
    try:
        try:
            temporary_path = write_normalized_xyz(input_file)
        except Exception as error:
            print(f" Error: failed to prepare extxyz file '{input_file}': {error}")
            return 1
        try:
            multi_systems = MultiSystems.from_file(temporary_path, fmt="extxyz")
        except Exception as error:
            print(f" Error: failed to parse extxyz file '{input_file}': {error}")
            return 1

        if not multi_systems.systems:
            print(f" Error: no structures were found in '{input_file}'.")
            return 1

        input_species = {
            species
            for system in multi_systems.systems.values()
            for species in system.get_atom_names()
        }
        missing_species = sorted(input_species.difference(type_map))
        if missing_species:
            print(
                " Error: type-map is missing species found in the input: "
                + ", ".join(missing_species)
                + "."
            )
            return 1

        multi_systems = MultiSystems(
            *multi_systems.systems.values(), type_map=type_map
        )
        print_system_summary(multi_systems, type_map)
        try:
            multi_systems.to_deepmd_npy(OUTPUT_DIR)
        except Exception as error:
            print(f" Error: failed to write DeepMD npy data to '{OUTPUT_DIR}': {error}")
            return 1
    finally:
        if temporary_path and os.path.exists(temporary_path):
            os.unlink(temporary_path)

    print(f" Saved DeepMD npy dataset(s) to '{OUTPUT_DIR}'.")
    print(
        " Use the same type-map order in DeepMD input.json and any downstream "
        "pair-table configuration."
    )
    return 0


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) < 2 or args[0] in ("-h", "--help"):
        print_usage()
        sys.exit(0 if args and args[0] in ("-h", "--help") else 1)
    sys.exit(main(args))

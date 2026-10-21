"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     exyz2pos_scf.py
Category:   Workflow Scripts
Purpose:    Convert extxyz frames to validated POSCAR files for VASP 301,
            grouping frames by their element sets when needed.
Usage:      python3 exyz2pos_scf.py <input.xyz> <struct_fp_dir>
Arguments:
  input.xyz      Input extended XYZ trajectory
  struct_fp_dir  Directory in which POSCAR files will be created
Output:
  One POSCAR_N.vasp per frame, flat for one element set or grouped by set
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-09-26
=============================================================================
"""

import os
import sys
import tempfile
from collections import Counter
from pathlib import Path


def print_help():
    print(" Usage: python3 exyz2pos_scf.py <input.xyz> <struct_fp_dir>")
    print(" Input:")
    print("   input.xyz       Input extended XYZ trajectory")
    print("   struct_fp_dir  Directory for generated POSCAR files")
    print(" Output:")
    print("   One POSCAR_N.vasp per frame; multiple element sets are grouped")
    print(" Notes:")
    print("   Elements follow their first-seen order in the input trajectory.")


def unique_symbol_order(atoms):
    """Return the distinct symbols in the order they occur in an Atoms object."""
    return list(dict.fromkeys(atoms.get_chemical_symbols()))


def main():
    args = sys.argv[1:]
    if args and args[0] in ("-h", "--help"):
        print_help()
        return 0
    if len(args) != 2:
        print_help()
        print(" Error: provide an input extxyz file and an output directory.")
        return 1

    input_file = Path(args[0])
    output_dir = Path(args[1])
    if not input_file.is_file():
        print(f" Error: input file not found: {input_file}")
        return 1
    if output_dir.name in ("", ".", ".."):
        print(f" Error: invalid output directory: {output_dir}")
        return 1
    if os.path.lexists(output_dir):
        if output_dir.is_symlink() or not output_dir.is_dir():
            print(f" Error: output path already exists and is not a directory: {output_dir}")
            return 1
        if any(output_dir.iterdir()):
            print(f" Error: output directory is not empty; refusing to overwrite: {output_dir}")
            return 1

    # Reuse the public converter's atom ordering instead of defining a second
    # element-order rule for the 301 workflow.
    converter_dir = Path(__file__).resolve().parents[1] / "format_conversion"
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(converter_dir))
    try:
        from exyz2pos import first_seen_element_order, reorder_atoms
    except ImportError as exc:
        print(f" Error: failed to load the shared exyz2pos ordering helpers: {exc}")
        return 1

    # Import ASE after handling help and argument errors.
    try:
        from ase.io import read, write
    except ImportError as exc:
        print(f" Error: ASE is required for VASP 301 extxyz conversion: {exc}")
        return 1

    try:
        frames = read(str(input_file), index=":", format="extxyz")
    except Exception as exc:
        print(f" Error: failed to read {input_file}: {exc}")
        return 1

    if not frames:
        print(f" Error: no frames found in {input_file}")
        return 1

    element_order = first_seen_element_order(frames)
    if not element_order:
        print(f" Error: no atoms found in {input_file}")
        return 1

    frame_signatures = []
    grouped_frames = {}
    for frame_index, atoms in enumerate(frames, start=1):
        frame_symbols = set(atoms.get_chemical_symbols())
        # POTCAR selection depends on the species set; atom-count ratios do not
        # split structures into separate groups.
        signature = tuple(symbol for symbol in element_order if symbol in frame_symbols)
        if not signature:
            print(f" Error: frame {frame_index} contains no atoms.")
            return 1
        frame_signatures.append(signature)
        grouped_frames.setdefault(signature, []).append(frame_index)

    use_group_directories = len(grouped_frames) > 1
    try:
        output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{output_dir.name}.stage-", dir=str(output_dir.parent)
        ) as temporary_directory:
            staged_output = Path(temporary_directory) / output_dir.name
            staged_output.mkdir()

            for frame_index, (atoms, signature) in enumerate(
                zip(frames, frame_signatures), start=1
            ):
                frame_dir = staged_output
                if use_group_directories:
                    frame_dir = staged_output / "_".join(signature)
                    frame_dir.mkdir(exist_ok=True)
                output_file = frame_dir / f"POSCAR_{frame_index}.vasp"
                reordered = reorder_atoms(atoms, element_order)
                try:
                    write(str(output_file), reordered, format="vasp")
                    written_atoms = read(str(output_file), format="vasp")
                except Exception as exc:
                    print(f" Error: failed to write or read back {output_file}: {exc}")
                    return 1

                source_counts = Counter(atoms.get_chemical_symbols())
                output_counts = Counter(written_atoms.get_chemical_symbols())
                written_order = unique_symbol_order(written_atoms)
                if output_counts != source_counts or written_order != list(signature):
                    print(f" Error: POSCAR validation failed for frame {frame_index}.")
                    print(f" Expected order: {' '.join(signature)}")
                    print(f" Read-back order: {' '.join(written_order)}")
                    print(f" Expected counts: {dict(source_counts)}")
                    print(f" Read-back counts: {dict(output_counts)}")
                    return 1

            # Keep the final output all-or-nothing and never replace existing files.
            if os.path.lexists(output_dir):
                if output_dir.is_symlink() or not output_dir.is_dir() or any(output_dir.iterdir()):
                    print(f" Error: output directory changed during conversion: {output_dir}")
                    return 1
                output_dir.rmdir()
            os.replace(staged_output, output_dir)
    except Exception as exc:
        print(f" Error: failed to create POSCAR output in {output_dir}: {exc}")
        return 1

    print(f" Converted {len(frames)} frame(s) using element order: {' '.join(element_order)}")
    if use_group_directories:
        print(" Multiple element sets detected. Prepare these POTCAR files in fp/:")
        for signature, frame_indices in grouped_frames.items():
            key = "_".join(signature)
            print(
                f"   fp/POTCAR_{key}  (POSCAR order: {' '.join(signature)}; "
                f"{len(frame_indices)} frame(s))"
            )
    else:
        signature = next(iter(grouped_frames))
        print(f" One element set detected: {' '.join(signature)}")
        print(" Prepare the shared fp/POTCAR with the same element order.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

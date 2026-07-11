"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     split_train_test.py
Category:   Sample Structure Scripts
Purpose:    Inspect an extxyz data set and interactively split it into training
            and test sets using uniform, random, or NEP-descriptor FPS selection.
Usage:      gpumdkit.sh
            choose 206) Split training and test sets
            python split_train_test.py <input.xyz>
Arguments:
  input.xyz  Input extxyz data set
Output:
  <input>_train.xyz  Frames not selected for the test set
  <input>_test.xyz   Frames selected for the test set
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-07-11
=============================================================================
"""

import os
import sys
from pathlib import Path


def print_usage():
    print(" Usage: gpumdkit.sh")
    print("        choose 206) Split training and test sets")
    print("    or: python split_train_test.py <input.xyz>")
    print("")
    print(" Arguments:")
    print("   input.xyz  Input extxyz data set")
    print("")
    print(" Output:")
    print("   <input>_train.xyz and <input>_test.xyz")
    print("")


def read_prompt(message):
    try:
        return input(message).strip()
    except (EOFError, KeyboardInterrupt):
        print("\n Input closed. Exiting.")
        sys.exit(1)


def print_section(title):
    print(" +------------------------------------------------------+")
    print(f" | {title:<52} |")
    print(" +------------------------------------------------------+")


def parse_test_size(value, total_frames):
    try:
        if any(character in value.lower() for character in (".", "e")):
            ratio = float(value)
            if not 0 < ratio < 1:
                raise ValueError
            number = max(1, int(total_frames * ratio + 0.5))
        else:
            number = int(value)
            if number <= 0:
                raise ValueError
    except ValueError:
        print(" Error: enter a ratio between 0 and 1, or a positive integer frame count.")
        sys.exit(1)

    if number >= total_frames:
        print(f" Error: the test set size ({number}) must be smaller than the total frame count ({total_frames}).")
        sys.exit(1)
    return number


def fps_indices(frames, number, model_file, np):
    if not os.path.isfile(model_file):
        print(f" Error: NEP model file '{model_file}' does not exist.")
        sys.exit(1)

    try:
        from NepTrain.core.nep import Nep3Calculator
        from scipy.spatial.distance import cdist
    except ImportError:
        print(" Error: FPS requires NepTrain and scipy, but they cannot be imported.")
        print(" Please install NepTrain and the scientific Python stack before using FPS.")
        sys.exit(1)

    try:
        calculator = Nep3Calculator(model_file)
        descriptors = []
        for index, frame in enumerate(frames):
            per_atom = calculator.get_descriptors(frame)
            descriptors.append(per_atom.mean(axis=0))
            print(f"\r Calculating NEP descriptors: {index + 1}/{len(frames)}", end="", flush=True)
        print()
        descriptors = np.asarray(descriptors)
    except Exception as error:
        print(f"\n Error: failed to calculate NEP descriptors: {error}")
        sys.exit(1)

    selected = [0]
    distances = cdist(descriptors, descriptors[[0]], metric="euclidean")[:, 0]
    distances[0] = -1.0
    while len(selected) < number:
        next_index = int(np.argmax(distances))
        selected.append(next_index)
        new_distances = cdist(descriptors, descriptors[[next_index]], metric="euclidean")[:, 0]
        distances = np.minimum(distances, new_distances)
        distances[selected] = -1.0
    return selected


def main():
    args = sys.argv[1:]
    if len(args) != 1 or args[0] in ("-h", "--help"):
        print_usage()
        sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

    input_file = args[0]
    if not os.path.isfile(input_file):
        print(f" Error: input file '{input_file}' does not exist.")
        sys.exit(1)

    try:
        import numpy as np
        from ase.io import read, write
    except ImportError:
        print(" Error: this function requires numpy and ASE.")
        sys.exit(1)

    try:
        frames = read(input_file, index=":")
    except Exception as error:
        print(f" Error: failed to read '{input_file}': {error}")
        sys.exit(1)

    if len(frames) < 2:
        print(" Error: the input data set must contain at least two frames.")
        sys.exit(1)

    elements = sorted({symbol for frame in frames for symbol in frame.get_chemical_symbols()})
    atom_counts = [len(frame) for frame in frames]
    print("")
    print_section("DATASET SUMMARY")
    print(f"   File: {input_file}")
    print(f"   Elements: {' '.join(elements)}")
    print(f"   Frames: {len(frames)}")
    print(f"   Atoms per frame: {min(atom_counts)}-{max(atom_counts)}")
    print("")
    print_section("STEP 1/2: TEST-SET SIZE")
    print("   0  < value <  1 : fraction of all frames")
    print("                     Example: 0.1 selects 10% for testing")
    print("   Positive integer: exact number of test frames")
    print("                     Example: 100 selects 100 frames")
    print(" ------------>>")
    size_value = read_prompt(" Enter the test-set fraction or frame count:\n ")
    test_count = parse_test_size(size_value, len(frames))
    print(f" Test-set size: {test_count} frames; training-set size: {len(frames) - test_count} frames")

    print("")
    print_section("STEP 2/2: SELECTION METHOD")
    print("   1) Uniform")
    print("   2) Random")
    print("   3) FPS using NepTrain and a NEP model")
    print(" ------------>>")
    method = read_prompt(" Enter the selection method:\n ").lower()

    if method in ("1", "uniform"):
        test_indices = np.linspace(0, len(frames) - 1, test_count, dtype=int).tolist()
        method_name = "uniform"
    elif method in ("2", "random"):
        print(" A fixed integer seed makes the random split reproducible.")
        seed_value = read_prompt(" Enter a random seed, or press Enter to use no fixed seed:\n ")
        try:
            seed = None if seed_value == "" else int(seed_value)
        except ValueError:
            print(" Error: random seed must be an integer or empty.")
            sys.exit(1)
        generator = np.random.default_rng(seed)
        test_indices = generator.choice(len(frames), test_count, replace=False).tolist()
        method_name = "random"
    elif method in ("3", "fps"):
        print(" FPS uses mean atomic NEP descriptors, following function 203.")
        model_file = read_prompt(" Enter the NEP model file (e.g., nep.txt):\n ")
        test_indices = fps_indices(frames, test_count, model_file, np)
        method_name = "FPS"
    else:
        print(" Error: selection method must be 1/uniform, 2/random, or 3/FPS.")
        sys.exit(1)

    test_index_set = set(test_indices)
    train_frames = [frame for index, frame in enumerate(frames) if index not in test_index_set]
    test_frames = [frame for index, frame in enumerate(frames) if index in test_index_set]

    input_path = Path(input_file)
    output_stem = input_path.with_suffix("")
    train_file = f"{output_stem}_train.xyz"
    test_file = f"{output_stem}_test.xyz"
    try:
        write(train_file, train_frames, format="extxyz")
        write(test_file, test_frames, format="extxyz")
    except Exception as error:
        print(f" Error: failed to write split data sets: {error}")
        sys.exit(1)

    print("")
    print_section("SPLIT COMPLETED")
    print(f"   Method: {method_name}")
    print(f"   Training frames: {len(train_frames)} -> {train_file}")
    print(f"   Test frames: {len(test_frames)} -> {test_file}")


if __name__ == "__main__":
    main()

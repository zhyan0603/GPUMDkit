"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     dp2xyz.py
Category:   Format Conversion Scripts
Purpose:    Convert DeepMD npy datasets under a directory to extended XYZ
            format using dpdata. Recursively finds type.raw + set.000
            directories and converts all frames into a single extxyz file.
Usage:      gpumdkit.sh -dp2xyz <input_dir/> [output.xyz]
            python3 dp2xyz.py <input_dir/> [output.xyz]
Arguments:
  input_dir/    Directory to scan recursively for DeepMD npy datasets
  output.xyz    Output extxyz file (default: train.xyz)
Output:
  <output.xyz>  (converted structures in extxyz format)
Author:     Denan LI (lidenan@westlake.edu.cn)
Last-modified: 2026-07-01
=============================================================================
"""

import os
import sys

args = sys.argv[1:]
if len(args) < 1 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -dp2xyz <input_dir/> [output.xyz]")
    print("    or: python3 dp2xyz.py <input_dir/> [output.xyz]")
    print("")
    print(" Arguments:")
    print("   input_dir/    Directory to scan recursively for DeepMD npy datasets")
    print("   output.xyz    Output extxyz file (default: train.xyz)")
    print("")
    print(" Example: gpumdkit.sh -dp2xyz database train.xyz")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

import dpdata
from ase.io import write


def is_deepmd_npy_dataset(path):
    return (
        os.path.isfile(os.path.join(path, "type.raw"))
        and os.path.isfile(os.path.join(path, "type_map.raw"))
        and os.path.isdir(os.path.join(path, "set.000"))
    )


def find_deepmd_npy_datasets(input_dir):
    dataset_dirs = []
    for root, dirs, _ in os.walk(input_dir):
        dirs.sort()
        if is_deepmd_npy_dataset(root):
            dataset_dirs.append(root)
    return dataset_dirs


if __name__ == "__main__":
    input_dir = args[0]
    output_file = args[1] if len(args) == 2 else "train.xyz"

    if not os.path.isdir(input_dir):
        print(f" Error: input directory '{input_dir}' does not exist.")
        sys.exit(1)

    dataset_dirs = find_deepmd_npy_datasets(input_dir)
    if not dataset_dirs:
        print(f" Error: no DeepMD npy datasets found under {input_dir}")
        sys.exit(1)

    all_frames = []
    print(f" Found {len(dataset_dirs)} DeepMD npy dataset(s) under {input_dir}")

    for dataset_dir in dataset_dirs:
        try:
            system = dpdata.LabeledSystem(dataset_dir, fmt="deepmd/npy")
        except Exception as e:
            print(f" Error: failed to load {dataset_dir}: {e}")
            sys.exit(1)

        frames = system.to_ase_structure()
        all_frames.extend(frames)
        print(f" Loaded {len(frames)} frames from {dataset_dir}")

    try:
        write(output_file, all_frames, format="extxyz", append=False)
    except Exception as e:
        print(f" Error: failed to write {output_file}: {e}")
        sys.exit(1)

    print(f" Converted {len(all_frames)} frames from {len(dataset_dirs)} dataset(s)")
    print(f" Saved extxyz file to {output_file}")

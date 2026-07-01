"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     dp2xyz.py
Category:   Format Conversion Scripts
Purpose:    Convert DeepMD npy datasets under a directory to extended XYZ
            format using dpdata.
Usage:      python3 dp2xyz.py <input/> [output.xyz]
Arguments:
  input/      Directory to scan recursively for DeepMD npy datasets
  output.xyz  Output extxyz file, default: train.xyz
Output:
  <output.xyz>  (converted structures in extxyz format)
Author:     GPUMDkit Contributors
Last-modified: 2026-07-01
=============================================================================
"""

import os
import sys

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
    if len(sys.argv) not in [2, 3]:
        print(" Usage: python3 dp2xyz.py <input/> [output.xyz]")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) == 3 else "train.xyz"

    if not os.path.isdir(input_dir):
        print(f" Error: input directory '{input_dir}' does not exist.")
        sys.exit(1)

    dataset_dirs = find_deepmd_npy_datasets(input_dir)
    if not dataset_dirs:
        print(f" No DeepMD npy datasets found under {input_dir}")
        sys.exit(1)

    all_frames = []
    print(f" Found {len(dataset_dirs)} DeepMD npy datasets under {input_dir}")

    for dataset_dir in dataset_dirs:
        try:
            system = dpdata.LabeledSystem(dataset_dir, fmt="deepmd/npy")
        except Exception as e:
            print(f" Error loading {dataset_dir}: {e}")
            sys.exit(1)

        frames = system.to_ase_structure()
        all_frames.extend(frames)
        print(f" Loaded {len(frames)} frames from {dataset_dir}")

    try:
        write(output_file, all_frames, format="extxyz", append=False)
    except Exception as e:
        print(f" Error writing extxyz file {output_file}: {e}")
        sys.exit(1)

    print(f" Converted {len(all_frames)} frames from {len(dataset_dirs)} datasets")
    print(f" Saved extxyz file to {output_file}")

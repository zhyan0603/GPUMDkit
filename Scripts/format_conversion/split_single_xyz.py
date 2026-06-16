"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     split_single_xyz.py
Category:   Format Conversion Scripts
Purpose:    Split a multi-frame XYZ file into individual model_${i}.xyz files
            (one frame per file).
Usage:      python split_single_xyz.py <input.xyz>
Arguments:
  input.xyz   Input multi-frame XYZ file
Output:
  model_1.xyz, model_2.xyz, ...  (individual frame files)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
from ase.io import read, write

def split_xyz_to_individual_models(input_xyz_filename):
    """
    Split an XYZ file into individual model files.

    Parameters:
    - input_xyz_filename (str): Name of the input XYZ file to read.

    Reads the input XYZ file, splits each frame into individual files named model_${i}.xyz.
    """

    # Read all frames from input XYZ file
    frames = read(input_xyz_filename, format='extxyz', index=':')

    # Iterate over each frame and write to separate model files
    for i, frame in enumerate(frames):
        model_filename = f'model_{i + 1}.xyz'  
        write(model_filename, frame, format='extxyz')

    print(f' All frames from "{input_xyz_filename}" have been split into individual model files.')

if __name__ == '__main__':
    # Check if the number of arguments is correct
    if len(sys.argv) < 2:
        print("Usage: python split_single_xyz.py <input.xyz>")
        sys.exit(1)

    # Get input XYZ filename from command line argument
    input_xyz_filename = sys.argv[1]

    # Call function to split XYZ file into individual models
    split_xyz_to_individual_models(input_xyz_filename)

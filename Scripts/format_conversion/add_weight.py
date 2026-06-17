"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     add_weight.py
Category:   Format Conversion Scripts
Purpose:    Add a weight value to the 'Weight' attribute of each structure
            in an input file and save the modified structures.
Usage:      python add_weight.py <input_file> <output_file> <new_weight>
Arguments:
  input_file   Input structure file
  output_file  Output structure file with added weight
  new_weight   Weight value to assign
Output:
  <output_file>  (structures with updated Weight attribute)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import os
import sys
from ase.io import read, write

def main():
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_file> <output_file> <new_weight>")
        sys.exit(1)

    # Parse command line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    try:
        new_weight = float(sys.argv[3])
    except ValueError:
        print("Error: The third argument (new weight) must be a number.")
        sys.exit(1)

    # Read all structures from the input file
    try:
        structures = read(input_file, index=':')
    except Exception as e:
        print(f"Error: Could not read the file {input_file}. Reason: {e}")
        sys.exit(1)

    # Modify the 'Weight' in the info dictionary for each structure
    for structure in structures:
        structure.info['Weight'] = new_weight

    # Write the modified structures to the output file
    try:
        write(output_file, structures)
    except Exception as e:
        print(f"Error: Could not write to the file {output_file}. Reason: {e}")
        sys.exit(1)

    # Print success message
    print(f"Weights in {input_file} have been updated and saved to {output_file}")

if __name__ == "__main__":
    main()


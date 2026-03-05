import argparse
import sys
from ase.io import read, write

# 1. Setup Argument Parser for positional arguments
parser = argparse.ArgumentParser(
    description=" Convert VASP XDATCAR to extxyz format."
)
parser.add_argument("input", help=" Path to the input XDATCAR file")
parser.add_argument("output", help=" Path to the output .xyz file")

args = parser.parse_args()

# 2. Execute conversion logic
try:
    print(f" Reading XDATCAR: {args.input} ...")
    
    # Read all frames from the trajectory (index=":")
    # ASE automatically parses scaling, lattice, and species
    trajectory = read(args.input, index=":", format="vasp-xdatcar")
    
    num_frames = len(trajectory)
    print(f" Detected {num_frames} frames.")
    
    print(f" Writing to extxyz: {args.output} ...")
    
    # Writing in extxyz format preserves the Lattice matrix in the header
    write(args.output, trajectory, format="extxyz")
    
    print(" Conversion complete.")

except FileNotFoundError:
    print(f" Error: Could not find '{args.input}'. Please check the path.")
    sys.exit(1)
except Exception as e:
    print(f" An unexpected error occurred: {e}")
    sys.exit(1)
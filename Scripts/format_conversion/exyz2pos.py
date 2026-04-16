#!/usr/bin/env python3
"""
Convert extended XYZ file to POSCAR format.
Developer: Zihan YAN

Usage: exyz2pos.py <input.xyz>

Example:
    python exyz2pos.py train.xyz
    python exyz2pos.py dump.xyz

Output: Creates POSCAR_1.vasp, POSCAR_2.vasp, etc. for each frame.
"""

import sys
import os
from ase.io import read, write


def print_error(message, usage=None, example=None):
    """Print formatted error message to stderr."""
    print(f"Error: {message}", file=sys.stderr)
    
    if usage:
        print(f"Usage: {usage}", file=sys.stderr)
    
    if example:
        print(f"Example: {example}", file=sys.stderr)


def check_file_exists(filepath, description="File"):
    """Check if a file exists."""
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"{description} '{filepath}' does not exist")
    return filepath


def print_progress_bar(iteration, total, length=50):
    """Print progress bar to stdout."""
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = chr(9608) * filled_length + '-' * (length - filled_length)
    print(f'\r Progress: |{bar}| {percent}% Complete', end='\r')
    if iteration == total:
        print()


def validate_arguments():
    """Validate command line arguments."""
    if len(sys.argv) < 2:
        print_error(
            "Missing required argument: input file",
            usage="exyz2pos.py <input.xyz>",
            example="python exyz2pos.py train.xyz"
        )
        sys.exit(1)
    
    if sys.argv[1] in ['-h', '--help']:
        print(__doc__)
        sys.exit(0)
    
    input_file = sys.argv[1]
    
    try:
        check_file_exists(input_file, "Input file")
    except FileNotFoundError as e:
        print_error(
            str(e),
            usage="exyz2pos.py <input.xyz>",
            example="python exyz2pos.py train.xyz"
        )
        sys.exit(1)
    
    return input_file


def main():
    """Main function to convert XYZ to POSCAR."""
    # Validate arguments
    input_file = validate_arguments()
    
    # Read all frames
    try:
        frames = read(input_file, index=':')
    except Exception as e:
        print_error(
            f"Failed to read file '{input_file}': {str(e)}",
            usage="exyz2pos.py <input.xyz>"
        )
        sys.exit(1)
    
    if len(frames) == 0:
        print_error(
            f"No frames found in '{input_file}'",
            usage="exyz2pos.py <input.xyz>"
        )
        sys.exit(1)
    
    total_frames = len(frames)
    
    # Save to POSCAR
    for i, frame in enumerate(frames):
        poscar_filename = f'POSCAR_{i + 1}.vasp'
        write(poscar_filename, frame)
        print_progress_bar(i + 1, total_frames)
    
    print(f' All frames have been converted.')


if __name__ == "__main__":
    main()

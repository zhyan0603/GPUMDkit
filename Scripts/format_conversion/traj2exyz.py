import sys
import os
from ase.io import read, write

def convert_traj_to_extxyz(input_file, output_file):
    """
    Converts an ASE .traj file to .extxyz format.
    """
    if not os.path.exists(input_file):
        print(f" Error: The file '{input_file}' does not exist.")
        return

    try:
        print(f" Reading frames from: {input_file} ...")
        # index=':' ensures all images in the trajectory are loaded
        images = read(input_file, index=':')
        
        print(f" Writing {len(images)} frames to: {output_file} ...")
        write(output_file, images, format='extxyz')
        
        print(" Conversion completed successfully.")
        
    except Exception as e:
        print(f" An unexpected error occurred: {e}")

if __name__ == "__main__":
    # Check if correct number of arguments are provided
    if len(sys.argv) != 3:
        print(" Usage: python traj2exyz.py <input.traj> <output.xyz>")
        sys.stdout.flush()
    else:
        input_path = sys.argv[1]
        output_path = sys.argv[2]
        convert_traj_to_extxyz(input_path, output_path)
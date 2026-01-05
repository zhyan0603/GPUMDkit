import sys
import glob
from ase.io import read, write

def convert_poscar_to_extxyz(poscar_filenames, extxyz_filename):

    # Create an empty list to store the frames
    all_frames = []

    # Read the structure from each POSCAR file and append to the list
    for poscar_filename in poscar_filenames:
        frames = read(poscar_filename, format='vasp')
        
        # If 'frames' is a list (multiple frames), we should iterate over it
        if isinstance(frames, list):
            all_frames.extend(frames)  # Add all frames from this file
        else:
            all_frames.append(frames)  # Add the single frame

    # Write all frames to extxyz file (using append=False ensures overwrite)
    write(extxyz_filename, all_frames, format='extxyz', append=False)

if __name__ == '__main__':
    # Check if the number of arguments is correct
    if len(sys.argv) < 3:
        print("Usage: python pos2exyz.py <POSCAR> <extxyz_filename>")
        print("Example 1: python pos2exyz.py POSCAR model.xyz")
        # print("Example 2: python pos2exyz.py 'POSCAR*' train.xyz")
        sys.exit(1)

    poscar_filename_pattern = sys.argv[1]
    extxyz_filename = sys.argv[2]

    # Use glob to handle file patterns (with wildcards)
    poscar_filenames = glob.glob(poscar_filename_pattern)

    # If no files matched the pattern, print a message and exit
    if not poscar_filenames:
        print(f"No files matched the pattern: {poscar_filename_pattern}")
        sys.exit(1)

    # Call function to convert POSCAR to extxyz
    convert_poscar_to_extxyz(poscar_filenames, extxyz_filename)

    print(f"Conversion completed, structures saved to {extxyz_filename}")

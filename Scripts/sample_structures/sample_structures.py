import argparse
import numpy as np
import textwrap
from ase.io import read, write
from pathlib import Path


class MyFormatter(argparse.RawDescriptionHelpFormatter,
                  argparse.ArgumentDefaultsHelpFormatter):
    pass

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=MyFormatter,
                                     description=textwrap.dedent("""
    Purpose:
        This script samples structures from an extxyz file using either 'uniform' or 'random' sampling methods. 
        The sampled structures are then written to the 'sampled_structures.xyz' file.

    Author:
        Zihan YAN <yanzihan@westlake.edu.cn>
    Modified by:
        Dr. Huan Wang <huan.wang@whut.edu.cn>
    Usage:
        python sample_structures.py <extxyz_file> <sampling_method> <num_samples> [skip_initial])
    """))
    parser.add_argument('extxyz_file', 
                        type=Path, 
                        metavar='extXYZ_File',
                        help='Path to the extxyz file for sampling.'
                        )
    parser.add_argument('sampling_method', 
                        type=str, 
                        metavar='Method',
                        choices=['uniform', 'random'], 
                        help='Sampling method to use. Please choose either "uniform" or "random".'
                        )
    parser.add_argument('num_samples', 
                        type=int, 
                        metavar='Num_Samples',
                        help='Number of frames to sample.'
                        )
    parser.add_argument('skip_initial', 
                        type=int, 
                        metavar='NUM',
                        nargs='?',
                        default=0,
                        help='Number of initial frames to skip.'
                        )
    return parser.parse_args() 

def main():
    args = parse_args()
    file_path = args.extxyz_file
    sampling_method = args.sampling_method
    num_samples = args.num_samples
    skip_initial = args.skip_initial

    # Read all frames
    frames = read(file_path, index=':')

    # Total number of frames
    num_frames = len(frames)

    # Adjust the frames list if skip_initial is specified
    if skip_initial:
        frames = frames[skip_initial:]
        num_frames = len(frames)  # Update the number of frames after skipping

    # Generate evenly spaced indices
    if sampling_method == 'uniform':
        sampled_indices = np.linspace(0, num_frames-1, num_samples, dtype=int)
    # Generate random indices
    else:
        sampled_indices = np.random.choice(num_frames, num_samples, replace=False)

    # Initialize an empty list to store the sampled frames
    sampled_frames = []

    # Collect the frame for the sample.xyz file
    for i, idx in enumerate(sampled_indices):
        sampled_frames.append(frames[idx])

    # Write the sampled frames to sample.xyz
    write('sampled_structures.xyz', sampled_frames)

    print('All sampled frames written to sampled_structures.xyz')


if __name__ == "__main__":
    main()
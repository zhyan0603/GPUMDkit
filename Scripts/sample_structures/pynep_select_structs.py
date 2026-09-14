"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     pynep_select_structs.py
Category:   Sample Structure Scripts
Purpose:    Select diverse structures from a sampledata set using
            farthest-point sampling on pynep descriptors, relative to
            an existing training set.
Usage:      python pynep_select_structs.py <sampledata_file> <traindata_file> <nep_model_file>
Arguments:
  sampledata_file  Extxyz file with candidate structures
  traindata_file   Extxyz file with existing training structures
  nep_model_file   NEP model file (e.g., nep.txt)
Output:
  selected.xyz     Selected diverse structures
  select.png       PCA visualization of descriptor space
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import os
import sys

args = sys.argv[1:]
if len(args) < 3 or args[0] in ("-h", "--help"):
    print(" Usage: python pynep_select_structs.py <sampledata_file> <traindata_file> <nep_model_file>")
    print("        (superseded by parallel_pynep_select_structs.py / NepTrain sampling)")
    print("")
    print(" Example: python pynep_select_structs.py dump.xyz train.xyz nep.txt")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

import numpy as np
from ase.io import read, write
from pynep.calculate import NEP
from pynep.select import FarthestPointSample


def print_dependency_notice():
    print(" This function requires the pynep package.")
    print(" This PyNEP sampling entry is deprecated. We recommend using NepTrain sampling instead.")


def read_prompt(message):
    """Read one line from stdin with EOF-safe handling."""
    try:
        return input(message).strip()
    except (EOFError, KeyboardInterrupt):
        print("\n Input closed. Exiting.")
        sys.exit(1)

def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

# Calculate descriptors with progress bar
def calculate_descriptors():
    total_sample = len(sampledata)
    total_train = len(traindata)

    des_sample = []
    for i in range(total_sample):
        des_sample.append(np.mean(calc.get_property('descriptor', sampledata[i]), axis=0))
        print_progress_bar(i + 1, total_sample, prefix=' Processing sampledata:', suffix='Complete', length=50)
    #des_sample = np.load('des_sample.npy')
    des_sample = np.array(des_sample)
    #np.save('des_sample.npy', des_sample)

    des_train = []
    for i in range(total_train):
        des_train.append(np.mean(calc.get_property('descriptor', traindata[i]), axis=0))
        print_progress_bar(i + 1, total_train, prefix=' Processing traindata: ', suffix='Complete', length=50)
    #des_train = np.load('des_train.npy')
    des_train = np.array(des_train)
    #np.save('des_train.npy', des_train)
    
    return des_sample, des_train

print_dependency_notice()

for input_file in (args[0], args[1], args[2]):
    if not os.path.isfile(input_file):
        print(f" Error: file '{input_file}' does not exist.")
        sys.exit(1)

# Load data
sampledata = read(args[0], ':')
traindata = read(args[1], ':')

# Initialize NEP calculator
calc = NEP(args[2])
print(calc)

# Interactive selection method
print(" Choose selection method:")
print(" 1) Select structures based on minimum distance")
print(" 2) Select structures based on number of structures")
choice = read_prompt(" ------------>>\n ")

sampler = FarthestPointSample()
if choice == '1':
    min_dist = float(read_prompt(" Enter min_dist (e.g., 0.01): "))
    des_sample, des_train = calculate_descriptors()
    selected = sampler.select(des_sample, des_train, min_distance=min_dist, max_select=None)
elif choice == '2':
    try:
        min_max_input = read_prompt(" Enter min_select and max_select (e.g., '50 100'): ")
        min_select, max_select = map(int, min_max_input.split())
        if min_select < 1 or max_select < min_select:
            print(" Error: min_select must be >= 1 and max_select must be >= min_select.")
            sys.exit(1)
        des_sample, des_train = calculate_descriptors()
        selected = sampler.select(des_sample, des_train, min_select=min_select, max_select=max_select)
    except ValueError:
        print(" Error: Please enter two integers separated by a space (e.g., '50 100').")
        sys.exit(1)
else:
    print(" Invalid choice. Exiting.")
    sys.exit(1)

write('selected.xyz', [sampledata[i] for i in selected])

# Check if seaborn is installed
try:
    import seaborn as sns
    sns_installed = True
except ImportError:
    sns_installed = False

# PCA for dimensionality reduction and visualization
reducer = PCA(n_components=2)
reducer.fit(des_sample)
proj_sample = reducer.transform(des_sample)
proj_train = reducer.transform(des_train)
proj_selected = reducer.transform(np.array([des_sample[i] for i in selected]))

# Create the figure
plt.figure(figsize=(5, 5), dpi=200)

# Add the main scatter plot
main_ax = plt.gca()
main_ax.scatter(proj_sample[:, 0], proj_sample[:, 1], color='C0', label=sys.argv[1], alpha=0.4)
main_ax.scatter(proj_train[:, 0], proj_train[:, 1], color='C1', label=sys.argv[2], alpha=0.4)
main_ax.scatter(proj_selected[:, 0], proj_selected[:, 1], color='C2', label='selected.xyz', alpha=0.4)
main_ax.set_xlabel('PC1')
main_ax.set_ylabel('PC2')
main_ax.set_xticks([])
main_ax.set_yticks([])
main_ax.legend()

# Add projections if seaborn is available
if sns_installed:
    # Add density plots
    top_kde = main_ax.inset_axes([0, 1.05, 1, 0.2], transform=main_ax.transAxes)
    sns.kdeplot(x=proj_sample[:, 0], color='C0', ax=top_kde, fill=True, alpha=0.4)
    sns.kdeplot(x=proj_train[:, 0], color='C1', ax=top_kde, fill=True, alpha=0.4)
    sns.kdeplot(x=proj_selected[:, 0], color='C2', ax=top_kde, fill=True, alpha=0.4)
    top_kde.set_xticks([])
    top_kde.set_yticks([])

    side_kde = main_ax.inset_axes([1.05, 0, 0.2, 1], transform=main_ax.transAxes)
    sns.kdeplot(y=proj_sample[:, 1], color='C0', ax=side_kde, fill=True, alpha=0.4)
    sns.kdeplot(y=proj_train[:, 1], color='C1', ax=side_kde, fill=True, alpha=0.4)
    sns.kdeplot(y=proj_selected[:, 1], color='C2', ax=side_kde, fill=True, alpha=0.4)
    side_kde.set_xticks([])
    side_kde.set_yticks([])

plt.tight_layout()
#plt.show()
plt.savefig('select.png')

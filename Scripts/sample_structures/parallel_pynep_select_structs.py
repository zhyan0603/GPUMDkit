"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    parallel_pynep_select_structs.py
"""

import sys
import numpy as np
from ase.io import read, write
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pynep.calculate import NEP
from pynep.select import FarthestPointSample
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import argparse
import time
import os

def get_num_frames(file_path):
    count = 0
    with open(file_path, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.strip().isdigit():
                count += 1
                num_atoms = int(line.strip())
                for _ in range(num_atoms + 1):  # comment line + atoms
                    f.readline()
    return count

def read_chunk(file_path, start, end):
    return read(file_path, index=f'{start}:{end}')

def write_chunk(images, temp_file):
    write(temp_file, images, format='extxyz')

def parallel_read(file_path, max_workers, chunk_size=100):
    num_frames = get_num_frames(file_path)
    chunks = [(i, min(i + chunk_size, num_frames)) for i in range(0, num_frames, chunk_size)]
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(read_chunk, file_path, start, end) for start, end in chunks]
        results = []
        for future in as_completed(futures):
            results.extend(future.result())
    return results

def parallel_write(images, output_file, max_workers, chunk_size=100):
    chunks = [images[i:i + chunk_size] for i in range(0, len(images), chunk_size)]
    temp_files = [f'temp_{i}.xyz' for i in range(len(chunks))]
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(write_chunk, chunk, temp_file) for chunk, temp_file in zip(chunks, temp_files)]
        for future in as_completed(futures):
            future.result()
    
    # Merge temp files
    with open(output_file, 'w') as outfile:
        for temp_file in temp_files:
            with open(temp_file, 'r') as infile:
                outfile.write(infile.read())
            os.remove(temp_file)

# Custom progress bar function
def print_progress(completed, total, prefix='', length=50, fill='█'):
    percent = (100 * (completed / float(total)))
    filled_length = int(length * completed // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    sys.stdout.write(f'\r{prefix} |{bar}| {percent:.1f}% Complete')
    sys.stdout.flush()
    if completed == total:
        print()

# Calculate descriptors with progress bar
def compute_descriptor(atom):
    return np.mean(calc.get_property('descriptor', atom), axis=0)

def calculate_descriptors(data, desc_type, max_workers=None):
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()
    
    total = len(data)
    completed = 0
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(compute_descriptor, atom) for atom in data]
        results = []
        
        print_progress(completed, total, prefix=f'Processing {desc_type}:', length=50)
        
        for future in as_completed(futures):
            results.append(future.result())
            completed += 1
            print_progress(completed, total, prefix=f'Processing {desc_type}:', length=50)
    
    return np.array(results)

# Check command line arguments
if len(sys.argv) < 4:
    print(" Usage: python pynep_select_structs.py <sampledata_file> <traindata_file> <nep_model_file>")
    print(" Examp: python pynep_select_structs.py dump.xyz train.xyz nep.txt")
    sys.exit(1)

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Select structures using NEP descriptors.')
parser.add_argument('sampledata_file', type=str, help='Path to sampledata file')
parser.add_argument('traindata_file', type=str, help='Path to traindata file')
parser.add_argument('nep_model_file', type=str, help='Path to NEP model file')
parser.add_argument('threads', type=int, default=multiprocessing.cpu_count(), help='Number of threads for parallel processing')
args = parser.parse_args()

# Function to load data from file
def load_data(file_path):
    return read(file_path, ':')

# Parallel loading of data files
with ProcessPoolExecutor(max_workers=2) as executor:  # Only 2 files, so max_workers=2
    future_sample = executor.submit(load_data, args.sampledata_file)
    future_train = executor.submit(load_data, args.traindata_file)
    
    sampledata = future_sample.result()
    traindata = future_train.result()

# Initialize NEP calculator
calc = NEP(args.nep_model_file)
print(calc)

# Interactive selection method
print(" Choose selection method:")
print(" 1) Select structures based on minimum distance")
print(" 2) Select structures based on number of structures")
choice = input(" ------------>>\n ").strip()

sampler = FarthestPointSample()
if choice == '1':
    min_dist = float(input(" Enter min_dist (e.g., 0.01): ").strip())
    des_sample = calculate_descriptors(sampledata, 'sampledata', max_workers=args.threads)
    des_train = calculate_descriptors(traindata, 'traindata', max_workers=args.threads)
    selected = sampler.select(des_sample, des_train, min_distance=min_dist, max_select=None)
elif choice == '2':
    try:
        min_max_input = input(" Enter min_select and max_select (e.g., '50 100'): ").strip()
        min_select, max_select = map(int, min_max_input.split())
        if min_select < 1 or max_select < min_select:
            print(" Error: min_select must be >= 1 and max_select must be >= min_select.")
            sys.exit(1)
        des_sample = calculate_descriptors(sampledata, 'sampledata', max_workers=args.threads)
        des_train = calculate_descriptors(traindata, 'traindata', max_workers=args.threads)
        selected = sampler.select(des_sample, des_train, min_select=min_select, max_select=max_select)
    except ValueError:
        print(" Error: Please enter two integers separated by a space (e.g., '50 100').")
        sys.exit(1)
else:
    print(" Invalid choice. Exiting.")
    sys.exit(1)

parallel_write([sampledata[i] for i in selected], 'selected.xyz', args.threads)

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

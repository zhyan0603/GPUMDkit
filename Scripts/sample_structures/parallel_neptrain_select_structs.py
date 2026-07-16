"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     parallel_neptrain_select_structs.py
Category:   Sample Structure Scripts
Purpose:    Select diverse structures using farthest-point sampling on
            NepTrain descriptors, with parallel descriptor calculation.
Usage:      gpumdkit.sh
            choose 203) FPS sampling by NepTrain
            python parallel_neptrain_select_structs.py <sample.xyz> <train.xyz> <nep.txt> [threads]
Arguments:
  sample.xyz  Extxyz file with candidate structures
  train.xyz   Extxyz file with existing training structures
  nep.txt     NEP model file used to calculate descriptors
  threads     Number of parallel descriptor workers (default: 1)
Output:
  selected.xyz     Selected diverse structures
  select.png       PCA visualization of descriptor space
  pca_sample.txt / pca_train.txt / pca_selected.txt  PCA projections
Author:     Benrui TANG (tang070205@proton.me)
Last-modified: 2026-07-14
=============================================================================
"""

import importlib.util
import os
import sys
from concurrent.futures import ProcessPoolExecutor


_worker_calculator = None


def print_dependency_notice():
    """Print the NepTrain dependency and citation notice."""
    print(" This function requires the NepTrain package.")
    print(" If you use this function, we recommend citing:")
    print(" Chen et al., Comput. Phys. Commun. 317, 109859 (2025).")
    print(" https://doi.org/10.1016/j.cpc.2025.109859")


def print_usage():
    """Print command usage without importing optional dependencies."""
    print(" Usage: gpumdkit.sh")
    print("        choose 203) FPS sampling by NepTrain")
    print("    or: python parallel_neptrain_select_structs.py <sample.xyz> <train.xyz> <nep.txt> [threads]")
    print("")
    print(" Arguments:")
    print("   sample.xyz  Extxyz file with candidate structures")
    print("   train.xyz   Extxyz file with existing training structures")
    print("   nep.txt     NEP model file used to calculate descriptors")
    print("   threads     Number of parallel descriptor workers (default: 1)")
    print("")
    print(" Output:")
    print("   selected.xyz")
    print("   select.png")
    print("   pca_sample.txt, pca_train.txt, pca_selected.txt")
    print("")
    print(" Example: in interactive mode, enter: dump.xyz train.xyz nep.txt 4")
    print("          python parallel_neptrain_select_structs.py dump.xyz train.xyz nep.txt 4")
    print("")


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='#'):
    """Print a terminal progress bar."""
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end='\r')
    if iteration == total:
        print()


def initialize_worker(nep_model_file):
    """Load one independent NepTrain calculator in each worker process."""
    global _worker_calculator
    os.environ["OMP_NUM_THREADS"] = "1"
    from NepTrain.core.nep import Nep3Calculator
    _worker_calculator = Nep3Calculator(nep_model_file)


def calculate_mean_descriptor(atoms):
    """Calculate the mean atomic descriptor for one structure in a worker."""
    descriptors_per_atom = _worker_calculator.get_descriptors(atoms)
    return descriptors_per_atom.mean(axis=0)


def calculate_descriptor_set(data, data_name, calculator, executor, threads, np):
    """Calculate descriptors in input order and report progress."""
    total = len(data)
    if executor is None:
        descriptors = []
        for index, atoms in enumerate(data, start=1):
            descriptors_per_atom = calculator.get_descriptors(atoms)
            descriptors.append(descriptors_per_atom.mean(axis=0))
            print_progress_bar(
                index,
                total,
                prefix=f' Processing {data_name}:',
                suffix='Complete',
                length=50,
            )
        return np.asarray(descriptors)

    chunk_size = max(1, total // (threads * 4))
    descriptors = []
    for index, descriptor in enumerate(
        executor.map(calculate_mean_descriptor, data, chunksize=chunk_size),
        start=1,
    ):
        descriptors.append(descriptor)
        print_progress_bar(
            index,
            total,
            prefix=f' Processing {data_name}:',
            suffix='Complete',
            length=50,
        )
    return np.asarray(descriptors)


def farthest_point_sample(
    new_data,
    now_data=None,
    min_distance=None,
    min_select=1,
    max_select=None,
    metric='euclidean',
):
    """Select new descriptor vectors using farthest-point sampling."""
    max_select = max_select or len(new_data)
    to_add = []
    if len(new_data) == 0:
        return to_add
    if now_data is None:
        now_data = []
    if len(now_data) == 0:
        to_add.append(0)
        now_data = [new_data[0]]
    distances = np.min(cdist(new_data, now_data, metric=metric), axis=1)
    while np.max(distances) > min_distance or len(to_add) < min_select:
        index = np.argmax(distances)
        to_add.append(index)
        if len(to_add) >= max_select:
            break
        distances = np.minimum(
            distances,
            cdist([new_data[index]], new_data, metric=metric)[0],
        )
    return to_add


def main():
    """Run parallel NepTrain descriptor calculation and FPS selection."""
    args = sys.argv[1:]
    if len(args) not in (3, 4) or args[0] in ("-h", "--help"):
        print_usage()
        sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

    sample_file, train_file, model_file = args[:3]
    try:
        threads = int(args[3]) if len(args) == 4 else 1
    except ValueError:
        print(" Error: threads must be a positive integer.")
        sys.exit(1)
    if threads < 1:
        print(" Error: threads must be a positive integer.")
        sys.exit(1)

    # NepTrain's native backend can use OpenMP. Keep one native thread in each
    # worker so the user-selected worker count is the total parallelism.
    os.environ["OMP_NUM_THREADS"] = "1"

    for file_path in (sample_file, train_file, model_file):
        if not os.path.isfile(file_path):
            print(f" Error: file '{file_path}' does not exist.")
            sys.exit(1)

    print_dependency_notice()
    if len(args) == 3:
        print(" Descriptor threads: 1 (default).")
    else:
        print(f" Descriptor threads: {threads}.")
    print(" Set the optional fourth argument to change the thread count.")

    if importlib.util.find_spec("NepTrain") is None:
        print(" Error: required Python dependencies are not installed or cannot be imported.")
        print(" Please install NepTrain and the scientific Python stack before using this function.")
        sys.exit(1)

    try:
        global np, cdist
        import numpy as np
        from ase.io import read, write
        import matplotlib.pyplot as plt
        from sklearn.decomposition import PCA
        from scipy.spatial.distance import cdist
    except ImportError:
        print(" Error: required Python dependencies are not installed or cannot be imported.")
        print(" Please install NepTrain and the scientific Python stack before using this function.")
        sys.exit(1)

    sampledata = read(sample_file, ':')
    traindata = read(train_file, ':')
    calculator = None
    if threads == 1:
        try:
            initialize_worker(model_file)
            calculator = _worker_calculator
        except ImportError:
            print(" Error: required Python dependencies are not installed or cannot be imported.")
            print(" Please install NepTrain and the scientific Python stack before using this function.")
            sys.exit(1)
    if calculator is not None:
        print(calculator)

    print(" Choose selection method:")
    print(" 1) Select structures based on minimum distance")
    print(" 2) Select structures based on number of structures")
    choice = input(" ------------>>\n ").strip()

    min_dist = None
    min_select = 1
    max_select = None
    if choice == '1':
        try:
            min_dist = float(input(" Enter min_dist (e.g., 0.01): ").strip())
        except ValueError:
            print(" Error: min_dist must be a number.")
            sys.exit(1)
    elif choice == '2':
        try:
            min_max_input = input(" Enter min_select and max_select (e.g., '50 100'): ").strip()
            min_select, max_select = map(int, min_max_input.split())
        except ValueError:
            print(" Error: Please enter two integers separated by a space (e.g., '50 100').")
            sys.exit(1)
        if min_select < 1 or max_select < min_select:
            print(" Error: min_select must be >= 1 and max_select must be >= min_select.")
            sys.exit(1)
    else:
        print(" Error: selection method must be 1 or 2.")
        sys.exit(1)

    try:
        if threads == 1:
            des_sample = calculate_descriptor_set(sampledata, 'sampledata', calculator, None, threads, np)
            des_train = calculate_descriptor_set(traindata, 'traindata', calculator, None, threads, np)
        else:
            with ProcessPoolExecutor(
                max_workers=threads,
                initializer=initialize_worker,
                initargs=(model_file,),
            ) as executor:
                des_sample = calculate_descriptor_set(sampledata, 'sampledata', None, executor, threads, np)
                des_train = calculate_descriptor_set(traindata, 'traindata', None, executor, threads, np)
    except Exception as error:
        print(f" Error: descriptor calculation failed: {error}")
        sys.exit(1)

    if choice == '1':
        selected = farthest_point_sample(
            des_sample,
            des_train,
            min_distance=min_dist,
            max_select=None,
        )
    elif choice == '2':
        selected = farthest_point_sample(
            des_sample,
            des_train,
            min_distance=-1,
            min_select=min_select,
            max_select=max_select,
        )
    write('selected.xyz', [sampledata[index] for index in selected])

    try:
        import seaborn as sns
        sns_installed = True
    except ImportError:
        sns_installed = False

    reducer = PCA(n_components=2)
    reducer.fit(des_sample)
    proj_sample = reducer.transform(des_sample)
    proj_train = reducer.transform(des_train)
    proj_selected = reducer.transform(np.asarray([des_sample[index] for index in selected]))

    plt.figure(figsize=(5, 5), dpi=200)
    main_ax = plt.gca()
    main_ax.scatter(proj_sample[:, 0], proj_sample[:, 1], color='C0', label=sample_file, alpha=0.4)
    main_ax.scatter(proj_train[:, 0], proj_train[:, 1], color='C1', label=train_file, alpha=0.4)
    main_ax.scatter(proj_selected[:, 0], proj_selected[:, 1], color='C2', label='selected.xyz', alpha=0.4)
    main_ax.set_xlabel('PC1')
    main_ax.set_ylabel('PC2')
    main_ax.set_xticks([])
    main_ax.set_yticks([])
    main_ax.legend()

    if sns_installed:
        top_kde = main_ax.inset_axes([0, 1.05, 1, 0.2], transform=main_ax.transAxes)
        sns.kdeplot(x=proj_sample[:, 0], color='C0', ax=top_kde, fill=True, alpha=0.4, warn_singular=False)
        sns.kdeplot(x=proj_train[:, 0], color='C1', ax=top_kde, fill=True, alpha=0.4, warn_singular=False)
        sns.kdeplot(x=proj_selected[:, 0], color='C2', ax=top_kde, fill=True, alpha=0.4, warn_singular=False)
        top_kde.set_xticks([])
        top_kde.set_yticks([])

        side_kde = main_ax.inset_axes([1.05, 0, 0.2, 1], transform=main_ax.transAxes)
        sns.kdeplot(y=proj_sample[:, 1], color='C0', ax=side_kde, fill=True, alpha=0.4, warn_singular=False)
        sns.kdeplot(y=proj_train[:, 1], color='C1', ax=side_kde, fill=True, alpha=0.4, warn_singular=False)
        sns.kdeplot(y=proj_selected[:, 1], color='C2', ax=side_kde, fill=True, alpha=0.4, warn_singular=False)
        side_kde.set_xticks([])
        side_kde.set_yticks([])

    np.savetxt('pca_sample.txt', proj_sample, fmt='%.8f', header='sample_x sample_y', comments='')
    np.savetxt('pca_train.txt', proj_train, fmt='%.8f', header='train_x train_y', comments='')
    np.savetxt('pca_selected.txt', proj_selected, fmt='%.8f', header='selected_x selected_y', comments='')

    plt.tight_layout()
    plt.savefig('select.png')


if __name__ == "__main__":
    main()

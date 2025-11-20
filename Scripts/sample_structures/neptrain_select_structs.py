import sys
import numpy as np
from ase.io import read, write
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from NepTrain.core.nep import *
from scipy.spatial.distance import cdist


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
        descriptors_per_atom = calc.get_descriptors(sampledata[i])    
        mean_descriptor = descriptors_per_atom.mean(axis=0)        
        des_sample.append(mean_descriptor)
        print_progress_bar(i + 1, total_sample, prefix=' Processing sampledata:', suffix='Complete', length=50)
    #des_sample = np.load('des_sample.npy')
    des_sample = np.array(des_sample)
    #np.save('des_sample.npy', des_sample)

    des_train = []
    for i in range(total_train):
        descriptors_per_atom = calc.get_descriptors(traindata[i])
        mean_descriptor = descriptors_per_atom.mean(axis=0)
        des_train.append(mean_descriptor)
        print_progress_bar(i + 1, total_train, prefix=' Processing traindata: ', suffix='Complete', length=50)
    #des_train = np.load('des_train.npy')
    des_train = np.array(des_train)
    #np.save('des_train.npy', des_train)
    
    return des_sample, des_train

def FarthestPointSample(new_data, now_data=[], min_distance=None, min_select=1, max_select=None, metric='euclidean'):
    max_select = max_select or len(new_data)
    to_add = []
    if len(new_data) == 0:
        return to_add
    if len(now_data) == 0:
        to_add.append(0)
        now_data.append(new_data[0])
    distances = np.min(cdist(new_data, now_data, metric=metric), axis=1)
    while np.max(distances) > min_distance or len(to_add) < min_select:
        i = np.argmax(distances)
        to_add.append(i)
        if len(to_add) >= max_select:
            break
        distances = np.minimum(distances, cdist([new_data[i]], new_data, metric=metric)[0])
    return to_add


# Check command line arguments
if len(sys.argv) < 4:
    print(" Usage: python pynep_select_structs.py <sampledata_file> <traindata_file> <nep_model_file>")
    print(" Examp: python pynep_select_structs.py dump.xyz train.xyz nep.txt")
    sys.exit(1)

# Load data
sampledata = read(sys.argv[1], ':')
traindata = read(sys.argv[2], ':')

# Initialize NEP calculator
calc = Nep3Calculator(sys.argv[3])
print(calc)

# Interactive selection method
print(" Choose selection method:")
print(" 1) Select structures based on minimum distance")
print(" 2) Select structures based on number of structures")
choice = input(" ------------>>\n ").strip()

if choice == '1':
    min_dist = float(input(" Enter min_dist (e.g., 0.01): ").strip())
    des_sample, des_train = calculate_descriptors()
    selected = FarthestPointSample(des_sample, des_train, min_distance=min_dist, max_select=None)
elif choice == '2':
    try:
        min_max_input = input(" Enter min_select and max_select (e.g., '50 100'): ").strip()
        min_select, max_select = map(int, min_max_input.split())
        if min_select < 1 or max_select < min_select:
            print(" Error: min_select must be >= 1 and max_select must be >= min_select.")
            sys.exit(1)
        des_sample, des_train = calculate_descriptors()
        selected = FarthestPointSample(des_sample, des_train, min_distance=-1, min_select=min_select, max_select=max_select)
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

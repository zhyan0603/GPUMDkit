from pynep.calculate import NEP
from pynep.select import FarthestPointSample
from ase.io import read, write
import numpy as np
import matplotlib.pyplot as plt
import sys
from sklearn.decomposition import PCA

def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

# Check command line arguments
if len(sys.argv) < 5:
    print(" Usage: python pynep_select_structs.py <sampledata_file> <traindata_file> <nep_model_file> <min_distance>")
    print(" Examp: python pynep_select_structs.py dump.xyz train.xyz ./nep.txt 0.01")
    sys.exit(1)

# Load data
sampledata = read(sys.argv[1], ':')
traindata = read(sys.argv[2], ':')

# Initialize NEP calculator
calc = NEP(sys.argv[3])
print(calc)

# Calculate descriptors with progress bar
total_sample = len(sampledata)
total_train = len(traindata)

des_sample = []
for i in range(total_sample):
    des_sample.append(np.mean(calc.get_property('descriptor', sampledata[i]), axis=0))
    print_progress_bar(i + 1, total_sample, prefix=' Processing sampledata:', suffix='Complete', length=50)
#des_train = np.load('des_sample.npy')
des_sample = np.array(des_sample)
np.save('des_sample.npy', des_sample)

des_train = []
for i in range(total_train):
    des_train.append(np.mean(calc.get_property('descriptor', traindata[i]), axis=0))
    print_progress_bar(i + 1, total_train, prefix=' Processing traindata:', suffix='Complete', length=50)
#des_train = np.load('des_train.npy')
des_train = np.array(des_train)
np.save('des_train.npy', des_train)

# Farthest Point Sampling
min_dist = float(sys.argv[4])
sampler = FarthestPointSample(min_distance=min_dist)
selected_i = sampler.select(des_sample, des_train)
write('selected.xyz', [sampledata[i] for i in selected_i])

# PCA for dimensionality reduction and visualization
reducer = PCA(n_components=2)
reducer.fit(des_sample)
proj = reducer.transform(des_sample)
proj_old = reducer.transform(des_train)
plt.figure(figsize=(6,4), dpi=200)
plt.xlabel('PCA_1')
plt.ylabel('PCA_2')
plt.scatter(proj[:,0], proj[:,1], label='dump.xyz')
plt.scatter(proj_old[:,0], proj_old[:,1], label='train.xyz')
selected_proj = reducer.transform(np.array([des_sample[i] for i in selected_i]))
plt.scatter(selected_proj[:,0], selected_proj[:,1], label='selected.xyz')
plt.legend()
plt.axis('off')
plt.tight_layout()
#plt.show()
plt.savefig('select.png')

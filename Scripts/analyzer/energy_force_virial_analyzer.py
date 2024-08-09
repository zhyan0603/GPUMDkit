import sys
import numpy as np
from ase.io import read
import matplotlib.pyplot as plt

def calculate_range(frames, property_name):
    values = []
    for frame in frames:
        if property_name in ["energy", "Energy"]:
            values.append(frame.get_potential_energy())
        elif property_name in ["force", "forces"]:
            forces = frame.get_forces()
            values.extend(np.linalg.norm(forces, axis=1))
        elif property_name in ["virial", "Virial"]:
            if 'virial' in frame.info:
                virial = frame.info['virial']
                values.append(np.linalg.norm(virial))
            else:
                raise ValueError("Virial information not found in frame.")
        else:
            raise ValueError("Invalid property. Choose from 'energy', 'force', or 'virial'.")
    return np.min(values), np.max(values), values

def plot_histogram(values, property_name):
    plt.hist(values, bins=30, edgecolor='black')
    plt.title(f'{property_name.capitalize()} Histogram')
    plt.xlabel(f'{property_name.capitalize()}')
    plt.ylabel('Frequency')
    plt.show()

if __name__ == "__main__":
    # Check if the required arguments are provided
    if len(sys.argv) < 3:
        print("Usage: python script.py <filename> <property> [hist]")
        sys.exit(1)
    
    filename = sys.argv[1]
    property_name = sys.argv[2]
    plot_hist = len(sys.argv) > 3 and sys.argv[3] == 'hist'
    
    # Read the extxyz file
    frames = read(filename, index=":")
    
    # Calculate the range of the specified property
    min_val, max_val, values = calculate_range(frames, property_name)
    
    # Print the range
    print(f"{property_name.capitalize()} range: {min_val} to {max_val}")
    
    # Plot histogram if requested
    if plot_hist:
        plot_histogram(values, property_name)


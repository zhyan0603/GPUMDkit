"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    get_volume.py
"""

import os
import numpy as np

# Function to calculate volume
def calculate_volume(a, b, c):
    volume = np.einsum('ij,ij->i', a, np.cross(b, c))
    return np.abs(volume)

# Function to extract volume information from the thermo.out file
def extract_volume_from_thermo(file_path):
    """
    Parse the thermo.out file and extract volume information.
    """
    data = np.loadtxt(file_path)
    num_columns = data.shape[1]

    if num_columns == 12:
        # For orthogonal box, directly calculate the volume
        box_length_x = data[:, 9]
        box_length_y = data[:, 10]
        box_length_z = data[:, 11]
        volume = box_length_x * box_length_y * box_length_z
    elif num_columns == 18:
        # For non-orthogonal box, calculate the volume using three vectors
        ax, ay, az = data[:, 9], data[:, 10], data[:, 11]
        bx, by, bz = data[:, 12], data[:, 13], data[:, 14]
        cx, cy, cz = data[:, 15], data[:, 16], data[:, 17]

        a_vectors = np.column_stack((ax, ay, az))
        b_vectors = np.column_stack((bx, by, bz))
        c_vectors = np.column_stack((cx, cy, cz))

        volume = calculate_volume(a_vectors, b_vectors, c_vectors)
    else:
        raise ValueError("Unsupported number of columns in thermo.out. Expected 12 or 18.")

    return volume

# Function to calculate the average volume for each temperature in *K folders
def calculate_average_volumes_for_temperature(folder):
    temperature_volume_dict = {}

    for subfolder in os.listdir(folder):
        if subfolder.endswith('K') and os.path.isdir(os.path.join(folder, subfolder)):
            thermo_file_path = os.path.join(folder, subfolder, 'thermo.out')
            if os.path.exists(thermo_file_path):
                # Extract volume data from thermo.out file
                volume = extract_volume_from_thermo(thermo_file_path)

                if len(volume) > 0:
                    try:
                        # Extract temperature from the subfolder name
                        temperature = int(subfolder.replace('K', ''))

                        # Calculate the average volume for this temperature
                        avg_volume = np.mean(volume)

                        # Store the result in the dictionary
                        temperature_volume_dict[temperature] = avg_volume
                    except ValueError:
                        continue  # Skip if the subfolder name does not represent a valid temperature
    return temperature_volume_dict

# Main function
def main():
    current_folder = os.getcwd()  # Get the current folder path
    temperature_volume_dict = calculate_average_volumes_for_temperature(current_folder)

    # Output the results in one line, with temperature and volume separated by a colon, and each pair separated by a comma
    output = ', '.join([f"{temp}: {volume:.3f}" for temp, volume in sorted(temperature_volume_dict.items())])
    print(output)

if __name__ == "__main__":
    main()

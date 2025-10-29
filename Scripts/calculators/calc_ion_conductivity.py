"""
This script calculates diffusivity and ionic conductivity from MSD data (msd.out)
obtained via GPUMD. It computes directional (x, y, z) and total diffusivities
and converts them to ionic conductivities by Nernst-Einstein equation.

Usage:
1. Run the script with:
   python script.py <element> <charge>
   - <element>: Chemical species (e.g., Li).
   - <charge>: Ion charge (e.g., 1 for Li⁺).
2. If required files (`thermo.out` and `model.xyz`) are not found, the script
   will prompt for manual input of structure volume, temperature, and number of ions.

Author: Zihan YAN
Modified by: Shengjie Tang
Date: Sep 1, 2025
"""

import os
import sys
import numpy as np
import scipy.constants as consts


# Function to calculate the volume of a triclinic box
def calculate_volume(a, b, c):
    return np.abs(np.einsum('ij,ij->i', a, np.cross(b, c)))


def get_conversion_factor(structure_volume, species_charge, num_ions, temperature):
    """
    Compute the conversion factor to relate diffusivity (cm^2/s) to ionic conductivity (mS/cm).

    Parameters:
        structure_volume (float): Volume of the structure in Å^3.
        species_charge (int): Charge of the diffusing ion.
        num_ions (int): Number of ions in the structure.
        temperature (float): Temperature in Kelvin.

    Returns:
        float: Conversion factor to compute ionic conductivity from diffusivity.
    """
    z = species_charge
    n = num_ions

    vol = structure_volume * 1e-24  # units cm^3
    return (
            1000
            * n
            / (vol * consts.N_A)
            * z ** 2
            * (consts.N_A * consts.e) ** 2
            / (consts.R * temperature)
    )


def read_msd_file(msd_file):
    """
    Read MSD data from a file.

    Parameters:
        msd_file (str): Path to the msd.out file.

    Returns:
        tuple: Time steps (dts) and MSD values for x, y, z, and total.
    """
    data = np.loadtxt(msd_file)
    dts = data[:, 0]  # Time steps
    msd_x = data[:, 1]  # MSD in x-direction
    msd_y = data[:, 2]  # MSD in y-direction
    msd_z = data[:, 3]  # MSD in z-direction
    msd_total = msd_x + msd_y + msd_z  # Total MSD
    return dts, msd_x, msd_y, msd_z, msd_total


def calculate_diffusivity_and_conductivity(filepath, structure_volume, species_charge, num_ions, temperature):
    """
    Calculate the diffusivity and ionic conductivity from MSD data.

    Parameters:
        filepath (str): Path to the msd.out file.
        structure_volume (float): Volume of the structure in Å^3.
        species_charge (int): Charge of the diffusing ion.
        num_ions (int): Number of ions in the structure.
        temperature (float): Temperature in Kelvin.

    Returns:
        tuple: Diffusivity (cm^2/s) and ionic conductivity (mS/cm).
    """
    dts, msd_x, msd_y, msd_z, msd_total = read_msd_file(filepath)

    # Use a portion of the MSD data for fitting
    start = int(len(dts) * 0.1)  # Start at 10% of the data
    end = int(len(dts) * 0.4)  # End at 40% of the data

    # Linear fit to obtain slopes for x, y, z, and total MSD
    k_x, _ = np.polyfit(dts[start:end], msd_x[start:end], 1)
    k_y, _ = np.polyfit(dts[start:end], msd_y[start:end], 1)
    k_z, _ = np.polyfit(dts[start:end], msd_z[start:end], 1)
    k_total, _ = np.polyfit(dts[start:end], msd_total[start:end], 1)

    # Compute diffusivity D (cm^2/s) for each direction and total
    D_x = k_x * 1e-4 / 2  # maybe?
    D_y = k_y * 1e-4 / 2
    D_z = k_z * 1e-4 / 2
    D_total = k_total * 1e-4 / 6

    # Compute conversion factor
    conversion_factor = get_conversion_factor(structure_volume, species_charge, num_ions, temperature)

    # Compute ionic conductivity (mS/cm) for each direction and total
    conductivity_x = D_x * conversion_factor
    conductivity_y = D_y * conversion_factor
    conductivity_z = D_z * conversion_factor
    conductivity_total = D_total * conversion_factor

    return (D_x, D_y, D_z, D_total), (conductivity_x, conductivity_y, conductivity_z, conductivity_total)


# Function to extract volume and temperature from thermo.out
def extract_thermo_data():
    if not os.path.exists('./thermo.out'):
        raise FileNotFoundError(" The file 'thermo.out' does not exist.")

    data = np.loadtxt('./thermo.out')
    num_columns = data.shape[1]

    # Extract thermodynamic quantities
    temperature = data[:, 0]

    # Extract box dimensions and volume
    if num_columns == 12:
        box_length_x = data[:, 9]
        box_length_y = data[:, 10]
        box_length_z = data[:, 11]
        volume = box_length_x * box_length_y * box_length_z
    elif num_columns == 18:
        ax, ay, az = data[:, 9], data[:, 10], data[:, 11]
        bx, by, bz = data[:, 12], data[:, 13], data[:, 14]
        cx, cy, cz = data[:, 15], data[:, 16], data[:, 17]

        a_vectors = np.column_stack((ax, ay, az))
        b_vectors = np.column_stack((bx, by, bz))
        c_vectors = np.column_stack((cx, cy, cz))

        volume = calculate_volume(a_vectors, b_vectors, c_vectors)
    else:
        raise ValueError(" Unsupported number of columns in thermo.out. Expected 12 or 18.")

    # Calculate averages after 50% of simulation time
    start_index = int(len(temperature) * 0.5)  # Can adjust this threshold
    avg_temperature = np.mean(temperature[start_index:])
    avg_volume = np.mean(volume[start_index:])

    return avg_temperature, avg_volume


# Function to count ions from model.xyz
def count_ions(atom_type="Li"):
    # Check if model.xyz exists
    if not os.path.exists("model.xyz"):
        raise FileNotFoundError("The file 'model.xyz' does not exist.")

    # Count ions from model.xyz
    num_ions = 0
    with open("model.xyz", 'r') as f:
        lines = f.readlines()
        num_atoms = int(lines[0].strip())  # First line is number of atoms
        for line in lines[2:]:  # Skip header lines
            parts = line.split()
            if parts[0] == atom_type:
                num_ions += 1

    # Check for run.in file and replication
    replication_factor = 1
    try:
        with open("run.in", 'r') as f:
            for line in f:
                if line.strip().startswith("replicate"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            x, y, z = map(int, parts[1:4])
                            replication_factor = x * y * z
                            break
                        except ValueError:
                            print(" Warning: Invalid replicate parameters in run.in")
    except FileNotFoundError:
        print(" +---------------------------------------------------------+")
        print(" | Warning: run.in file not found, assuming no replication |")
        print(" +---------------------------------------------------------+")

    # Return total ions considering replication
    return num_ions * replication_factor


# Main function to orchestrate the process
def main():
    # Check if thermo.out and model.xyz exist to extract or prompt for values
    if os.path.exists("./thermo.out") and os.path.exists("model.xyz"):
        avg_temperature, avg_volume = extract_thermo_data()
        atom_type = sys.argv[1] if len(sys.argv) > 1 else "Li"  # Default to Li if not passed as argument
        num_ions = count_ions(atom_type)
    else:
        # If files don't exist, prompt the user to enter values manually
        print(" Files 'thermo.out' and 'model.xyz' are not found.")
        print(" Please provide the following values:")
        print(" --------------------------->")
        avg_temperature = float(input(" Enter average temperature (in K): "))
        avg_volume = float(input(" Enter system volume (in Å^3): "))
        atom_type = sys.argv[1]
        num_ions = int(input(" Enter number of ions: "))

    # Output the calculated or inputted values
    print(f" Number of ions: {num_ions}")
    print(f" Average Volume: {avg_volume:.3f} Å^3")
    print(f" Average Temperature: {avg_temperature:.3f} K")

    # Calculate diffusivity and conductivity
    msd_file = "./msd.out"  # Path to MSD file
    structure_volume = avg_volume  # Use calculated or inputted volume
    species_charge = int(sys.argv[2])  # Replace with actual ion charge (e.g., Li+ is 1)

    diffusivities, conductivities = calculate_diffusivity_and_conductivity(
        msd_file, structure_volume, species_charge, num_ions, avg_temperature
    )

    # Print results
    D_x, D_y, D_z, D_total = diffusivities
    conductivity_x, conductivity_y, conductivity_z, conductivity_total = conductivities

    print(" ------------------------------")
    print(f" Diffusivity (D):")
    print(f"   D_x: {D_x:.3e} cm^2/s")
    print(f"   D_y: {D_y:.3e} cm^2/s")
    print(f"   D_z: {D_z:.3e} cm^2/s")
    print(f"   D_total: {D_total:.3e} cm^2/s")
    print(" ------------------------------")
    print(f" Ionic Conductivity:")
    print(f"   Sigma_x: {conductivity_x:.3e} mS/cm")
    print(f"   Sigma_y: {conductivity_y:.3e} mS/cm")
    print(f"   Sigma_z: {conductivity_z:.3e} mS/cm")
    print(f"   Sigma_total: {conductivity_total:.3e} mS/cm")
    print(" ------------------------------")

    if D_total < 1e-7:
        print(" WARNING: ")
        print(" The D_total is below the threshold (1e-7) cm^2/s. ")
        print(" Please check the MSD and trajectory to confirm whether diffusion has occurred. ")
        # print("+-----------------------------------------------+")
        # print("|               GPUMDkit WARNING                |")
        # print("|-----------------------------------------------|")
        # print("|         D_total is below 1e-7 cm^2/s          |")
        # print("|        Has diffusion really occurred?         |")
        # print("| Please check MSD and trajectory to confirm it.|")
        # print("+-----------------------------------------------+")

# Run the main function
if __name__ == "__main__":
    main()



"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Find outlier structures in dataset

Usage:
    python find_outliers.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

# Script to find outliers in training data based on RMSE thresholds for energy, force, and stress.
# It is still during development and may not handle all cases perfectly.
# Author: Zihan YAN (yanzihan@westlake.edu.cn)

import sys
import matplotlib.pyplot as plt
import numpy as np

# Function to parse xyz file and get natoms list and frames
def parse_xyz(filename):
    """
    Reads an XYZ file and returns a list of atom counts per structure and a list of frames.
    Each frame contains the number of atoms, comment line, and atom coordinates.
    """
    natoms = []
    frames = []
    with open(filename, 'r') as f:
        while True:
            line = f.readline().strip()
            if not line:
                break
            nat = int(line)
            natoms.append(nat)
            comment = f.readline().strip()
            atoms = []
            for _ in range(nat):
                atoms.append(f.readline().strip())
            frames.append((nat, comment, atoms))
    return natoms, frames

# Function to write xyz file
def write_xyz(filename, selected_frames):
    """
    Writes a list of frames to an XYZ file.
    Each frame includes the number of atoms, comment line, and atom coordinates.
    """
    with open(filename, 'w') as f:
        for nat, comment, atoms in selected_frames:
            f.write(f"{nat}\n")
            f.write(f"{comment}\n")
            for atom in atoms:
                f.write(f"{atom}\n")

# Function to calculate RMSE
def calculate_rmse(pred, actual):
    """
    Calculates the RMSE between predicted and actual values.
    Returns infinity if the input arrays are empty to avoid premature termination.
    """
    if len(pred) == 0:
        return 0.0
    return np.sqrt(np.mean((pred - actual) ** 2))

# Generalized select outliers
def select_outliers(rmse_func, err_func, thresh, n_struct):
    """
    Iteratively removes structures with the highest errors until the RMSE is below the threshold.
    Args:
        rmse_func: Function to compute RMSE for the current set of structures.
        err_func: Function to compute per-structure errors for ranking.
        thresh: RMSE threshold (in meV/atom for energy, meV/Å for force, GPa for stress).
        n_struct: Total number of structures.
    Returns:
        A set of indices of structures to be selected (high-error structures).
    """
    indices = list(range(n_struct))
    selected = set()
    current_rmse = rmse_func(indices)
    if current_rmse <= thresh:
        return selected
    while current_rmse > thresh and len(indices) > 0:
        errors = err_func(indices)
        max_err_idx = np.argmax(errors)
        global_idx = indices.pop(max_err_idx)
        selected.add(global_idx)
        current_rmse = rmse_func(indices)
    return selected

# Load data
energy_train = np.loadtxt('energy_train.out')
force_train = np.loadtxt('force_train.out')
stress_train = np.loadtxt('stress_train.out')
natoms_list, frames = parse_xyz('train.xyz')

# Filter out rows with invalid stress data
valid_struct = ~np.any(np.abs(stress_train[:, :12]) > 1e6, axis=1)

# Filter energy and stress
energy_train = energy_train[valid_struct]
stress_train = stress_train[valid_struct]

# Filter force
cum_old = np.cumsum([0] + natoms_list)
force_keep = []
start = 0
for i, valid in enumerate(valid_struct):
    end = start + natoms_list[i]
    if valid:
        force_keep.extend(range(start, end))
    start = end
force_train = force_train[force_keep]

# Filter natoms and frames
natoms_list = [natoms_list[i] for i in range(len(valid_struct)) if valid_struct[i]]
frames = [frames[i] for i in range(len(valid_struct)) if valid_struct[i]]

n_struct = len(natoms_list)
cum = np.cumsum([0] + natoms_list)

# Input RMSE thresholds (in meV for energy and force, GPa for stress)
if len(sys.argv) >= 4:
    # Use command-line arguments if provided (in order: energy, force, stress)
    energy_thresh = float(sys.argv[1]) / 1000  # Convert meV/atom to eV/atom
    force_thresh = float(sys.argv[2]) / 1000   # Convert meV/Å to eV/Å
    stress_thresh = float(sys.argv[3])         # GPa, no conversion needed
else:
    # Prompt for manual input if command-line arguments are not provided
    energy_thresh = float(input(" Enter energy RMSE threshold (meV/atom): ")) / 1000  # Convert to eV/atom
    force_thresh = float(input(" Enter force RMSE threshold (meV/Å): ")) / 1000      # Convert to eV/Å
    stress_thresh = float(input(" Enter stress RMSE threshold (GPa): "))             # No conversion needed

# Define funcs for energy
energy_rmse_func = lambda inds: calculate_rmse(energy_train[inds, 0], energy_train[inds, 1])
energy_err_func = lambda inds: np.abs(energy_train[inds, 0] - energy_train[inds, 1])

# For stress
stress_rmse_func = lambda inds: np.mean([calculate_rmse(stress_train[inds, i], stress_train[inds, i+6]) for i in range(6)]) if inds else 0.0
stress_err_func = lambda inds: np.array([np.sqrt(np.mean((stress_train[j, 0:6] - stress_train[j, 6:12])**2)) for j in inds])

# For force
force_rmse_func = lambda inds: (
    np.mean([calculate_rmse(fr[:, i], fr[:, i+3]) for i in range(3)]) if (fr := np.vstack([force_train[cum[j]:cum[j+1]] for j in inds] if inds else [])) is not None and len(fr) > 0 else 0.0
)
force_err_func = lambda inds: np.array([np.sqrt(np.mean((force_train[cum[j]:cum[j+1], 0:3] - force_train[cum[j]:cum[j+1], 3:6])**2)) for j in inds])

# Select
selected_energy = select_outliers(energy_rmse_func, energy_err_func, energy_thresh, n_struct)
selected_stress = select_outliers(stress_rmse_func, stress_err_func, stress_thresh, n_struct)
selected_force = select_outliers(force_rmse_func, force_err_func, force_thresh, n_struct)

# Union
selected_all = selected_energy | selected_stress | selected_force
left_set = set(range(n_struct)) - selected_all
left_indices = sorted(left_set)
selected_indices = sorted(selected_all)

# Save files
selected_frames = [frames[i] for i in selected_indices]
left_frames = [frames[i] for i in left_indices]
write_xyz('selected.xyz', selected_frames)
write_xyz('remained.xyz', left_frames)

# Colors
Remained_color = '#237B9F'
selected_color = '#EC817E' 

# Function to calculate dynamic axis limits
def calculate_limits(data1, data2, padding=0.08):
    all_data = np.concatenate((data1, data2))
    if len(all_data) == 0:
        return -1, 1
    data_min = np.min(all_data)
    data_max = np.max(all_data)
    data_range = data_max - data_min
    return data_min - padding * data_range, data_max + padding * data_range

# Create subplot
fig, axs = plt.subplots(1, 3, figsize=(12, 3.3), dpi=100)

# Energy plot
energy_left = energy_train[left_indices]
energy_select = energy_train[selected_indices]
xmin_energy, xmax_energy = calculate_limits(energy_train[:, 1], energy_train[:, 0])
axs[0].set_xlim(xmin_energy, xmax_energy)
axs[0].set_ylim(xmin_energy, xmax_energy)
axs[0].plot(energy_left[:, 1], energy_left[:, 0], '.', markersize=10, label='Remained', color=Remained_color)
axs[0].plot(energy_select[:, 1], energy_select[:, 0], '.', markersize=10, label='Selected', color=selected_color)
axs[0].plot([xmin_energy, xmax_energy], [xmin_energy, xmax_energy], linewidth=1, color='red')
axs[0].set_xlabel('DFT energy (eV/atom)', fontsize=10)
axs[0].set_ylabel('NEP energy (eV/atom)', fontsize=10)
axs[0].legend(frameon=False)
axs[0].tick_params(axis='both', labelsize=10)

# RMSE for energy
energy_rmse = calculate_rmse(energy_left[:, 0], energy_left[:, 1]) * 1000
axs[0].text(0.3, 0.1, f'RMSE (Rem): {energy_rmse:.2f} meV/atom', transform=axs[0].transAxes, fontsize=10, verticalalignment='center')

# Force plot
force_left = np.vstack([force_train[cum[i]:cum[i+1]] for i in left_indices]) if left_indices else np.empty((0,6))
force_select = np.vstack([force_train[cum[i]:cum[i+1]] for i in selected_indices]) if selected_indices else np.empty((0,6))
xmin_force, xmax_force = calculate_limits(force_train[:, 3:6].reshape(-1), force_train[:, 0:3].reshape(-1))
axs[1].set_xlim(xmin_force, xmax_force)
axs[1].set_ylim(xmin_force, xmax_force)
axs[1].plot(force_left[:, 3:6].reshape(-1), force_left[:, 0:3].reshape(-1), '.', markersize=10, label='Remained', color=Remained_color)
axs[1].plot(force_select[:, 3:6].reshape(-1), force_select[:, 0:3].reshape(-1), '.', markersize=10, label='Selected', color=selected_color)
axs[1].plot([xmin_force, xmax_force], [xmin_force, xmax_force], linewidth=1, color='red')
axs[1].set_xlabel(r'DFT force (eV/Å)', fontsize=10)
axs[1].set_ylabel(r'NEP force (eV/Å)', fontsize=10)
axs[1].tick_params(axis='both', labelsize=10)
axs[1].legend(frameon=False)

# RMSE for force
if len(force_left) > 0:
    force_rmse_comps = [calculate_rmse(force_left[:, i], force_left[:, i + 3]) for i in range(3)]
    mean_force_rmse = np.mean(force_rmse_comps) * 1000
else:
    mean_force_rmse = 0.0
axs[1].text(0.35, 0.1, rf'RMSE (Rem): {mean_force_rmse:.2f} meV/Å', transform=axs[1].transAxes, fontsize=10, verticalalignment='center')

# Stress plot
stress_left = stress_train[left_indices]
stress_select = stress_train[selected_indices]
xmin_stress, xmax_stress = calculate_limits(stress_train[:, 6:12].reshape(-1), stress_train[:, 0:6].reshape(-1))
axs[2].set_xlim(xmin_stress, xmax_stress)
axs[2].set_ylim(xmin_stress, xmax_stress)
axs[2].plot(stress_left[:, 6:12].reshape(-1), stress_left[:, 0:6].reshape(-1), '.', markersize=10, label='Remained', color=Remained_color)
axs[2].plot(stress_select[:, 6:12].reshape(-1), stress_select[:, 0:6].reshape(-1), '.', markersize=10, label='Selected', color=selected_color)
axs[2].plot([xmin_stress, xmax_stress], [xmin_stress, xmax_stress], linewidth=1, color='red')
axs[2].set_xlabel('DFT stress (GPa)', fontsize=10)
axs[2].set_ylabel('NEP stress (GPa)', fontsize=10)
axs[2].tick_params(axis='both', labelsize=10)
axs[2].legend(frameon=False)

# RMSE for stress
if len(stress_left) > 0:
    stress_rmse_comps = [calculate_rmse(stress_left[:, i], stress_left[:, i + 6]) for i in range(6)]
    mean_stress_rmse = np.mean(stress_rmse_comps)
else:
    mean_stress_rmse = 0.0
axs[2].text(0.4, 0.1, f'RMSE (Rem): {mean_stress_rmse:.3f} GPa', transform=axs[2].transAxes, fontsize=10, verticalalignment='center')

# Adjust layout
plt.tight_layout()
fig.subplots_adjust(top=0.968, bottom=0.16, left=0.086, right=0.983, hspace=0.2, wspace=0.25)
plt.savefig('slected_remained.png', dpi=300)

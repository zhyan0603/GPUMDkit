import re
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import get_backend
from matplotlib.lines import Line2D
import scipy.constants as consts
from ase.io import read

# Set plot font size
plt.rcParams.update({'font.size': 10})

# Get current directory and temperature folders (e.g., 400K, 500K)
base_dir = os.getcwd()
temp_folders = [f for f in os.listdir(base_dir) if re.match(r'\d+K$', f) and os.path.isdir(f)]
temp_folders = sorted(temp_folders, key=lambda x: int(x[:-1]))  # Sort by temperature

# Optional: Uncomment to manually specify temperatures
# temp_list = [400, 500, 600, 700]
# temp_folders = [f"{t}K" for t in temp_list if os.path.isdir(os.path.join(base_dir, f"{t}K"))]

if not temp_folders:
    raise ValueError("No valid temperature folders (e.g., 400K, 500K) found in current directory")

# Define color for plotting
color = '#FF5D5D'

# Get structure label from current directory name
cell = os.path.basename(base_dir)

# Read model.xyz from the first temperature folder
xyz_file = os.path.join(base_dir, temp_folders[0], 'model.xyz')
if not os.path.exists(xyz_file):
    raise FileNotFoundError(f"Missing {temp_folders[0]}/model.xyz file")

# Count Li or Na ions
atoms = read(xyz_file, format='extxyz')
num_ions = sum(1 for atom in atoms if atom.symbol in ['Li', 'Na'])

# Check for replicate parameter in run.in
run_in_path = os.path.join(base_dir, 'run.in')
rep_factor = 1
if os.path.exists(run_in_path):
    with open(run_in_path, 'r') as f:
        for line in f:
            if line.strip().startswith('replicate'):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        rep_x = int(parts[1])
                        rep_y = int(parts[2])
                        rep_z = int(parts[3])
                        rep_factor = rep_x * rep_y * rep_z
                        print(f"[Info] Detected replicate {rep_x} {rep_y} {rep_z}, multiplying ion count by {rep_factor}")
                        num_ions *= rep_factor
                        break
                    except ValueError:
                        print("[Warning] Invalid replicate format, assuming no replication")
else:
    print("[Note] No run.in file found, assuming no replication")

# Initialize lists for plotting and data storage
group_Ts = []
group_sigmaTs = []
group_sigmas = []
fig, ax1 = plt.subplots(figsize=(5, 3.4), dpi=150)

# Print header for conductivity data
print("\nConductivity Data:")
print(f"{'T (K)':<6} {'Sigma (S/cm)':<15} {'Sigma·T (K·S/cm)':<15}")
print("-" * 36)

# Process each temperature folder
for temp_folder in temp_folders:
    temp = int(temp_folder[:-1])  # Extract temperature value
    folder_path = os.path.join(base_dir, temp_folder)
    thermo_file = os.path.join(folder_path, 'thermo.out')
    msd_file = os.path.join(folder_path, 'msd.out')

    if not os.path.exists(thermo_file) or not os.path.exists(msd_file):
        print(f"[Warning] Missing thermo.out or msd.out in {temp_folder}, skipping")
        continue

    try:
        # Extract volume from thermo.out
        data = np.loadtxt(thermo_file)
        num_columns = data.shape[1]
        if num_columns == 12:
            vol = data[:, 9] * data[:, 10] * data[:, 11]
        elif num_columns == 18:
            a = np.column_stack((data[:, 9], data[:, 10], data[:, 11]))
            b = np.column_stack((data[:, 12], data[:, 13], data[:, 14]))
            c = np.column_stack((data[:, 15], data[:, 16], data[:, 17]))
            vol = np.abs(np.einsum('ij,ij->i', a, np.cross(b, c)))
        else:
            raise ValueError("Unsupported number of columns in thermo.out")
        volume = np.mean(vol)

        # Read MSD data
        data = np.loadtxt(msd_file)
        dts = data[:, 0]
        msds = np.sum(data[:, 1:4], axis=1)
        start = int(len(dts) * 0.4)
        end = int(len(dts) * 0.8)

        # Calculate diffusivity
        k, _ = np.polyfit(dts[start:end], msds[start:end], 1)
        D = k * 1e-4 / 6  # Convert A²/ps to cm²/s

        # Calculate conductivity
        z = 1  # Species charge
        vol_cm3 = volume * 1e-24  # Convert volume to cm³
        conversion_factor = (num_ions / (vol_cm3 * consts.N_A) * z**2 * (consts.N_A * consts.e)**2 / (consts.R * temp))
        conductivity = conversion_factor * D
        sigma_T = conductivity * temp

        if sigma_T <= 0 or np.isnan(sigma_T):
            print(f"[Skip] {temp_folder}: Invalid sigma_T = {sigma_T}")
            continue

        # Store data and print conductivity
        group_Ts.append(temp)
        group_sigmaTs.append(sigma_T)
        group_sigmas.append(conductivity)
        print(f"{temp:<6} {conductivity:<15.3e} {sigma_T:<15.3e}")

        ax1.scatter(1000 / temp, np.log(sigma_T), color=color, s=20, alpha=0.8)

    except Exception as e:
        print(f"[Error] Processing {temp_folder}: {e}")
        continue
print("-" * 36)


if len(group_Ts) < 2:
    raise ValueError("Insufficient data points for fitting, need at least two temperatures")

# Perform linear fitting
invT = 1000 / np.array(group_Ts)
ln_sigmaT = np.log(np.array(group_sigmaTs))
k, b = np.polyfit(invT, ln_sigmaT, 1)
Ea = -k * 1000 * consts.k / consts.e
print(f"\n{cell}, Ea: {Ea:.3f} eV")

# Plot fitted line
ax1.plot(invT, k * invT + b, linestyle='-', linewidth=2, color=color, label=cell, alpha=0.8)

# Set axis labels
ax1.set_xlabel("1000/T (1/K)", fontsize=10)
ax1.set_ylabel(r"ln($\sigma \cdot T$) (K S/cm)", fontsize=10)

# Add legend
legend_element = Line2D([0], [0], marker='o', linestyle='-', color=color, label=cell,
                        markersize=5, linewidth=2, alpha=0.8)
ax1.legend(handles=[legend_element], loc='lower left', bbox_to_anchor=(0.015, 0.02), fontsize=8, frameon=True)

# Extrapolate from lowest temperature to 300K
specified_temperature = 300
x_ext = 1000 / np.arange(min(group_Ts), specified_temperature - 1, -10)  # Start from lowest T to 300K
y_ext = k * x_ext + b
ax1.plot(x_ext, y_ext, linestyle='--', color=color, alpha=0.5)

ln_sigmaT_300 = k * 1000 / specified_temperature + b
sigma_300K = np.exp(ln_sigmaT_300) / specified_temperature
print(f"at {specified_temperature}K, {cell}: Sigma = {sigma_300K:.3e} S/cm")

def tick_function(X):
    V = 1000 / X  
    return ["%.0f" % z for z in V]

# Add top temperature axis
ax2 = ax1.twiny()
x_min, x_max = ax1.get_xlim()
ax2.set_xticks(ax1.get_xticks())
ax2.set_xbound(ax1.get_xbound())
ax2.set_xticklabels(tick_function(ax1.get_xticks()))
ax2.set_xlabel("T (K)", fontsize=10)

# Save and show plot
fig.tight_layout()
if len(sys.argv) > 1 and sys.argv[1] == 'save':
    fig.savefig('Arrhenius_sigma.png', dpi=300)
else:
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'Arrhenius_sigma.png'.")
        fig.savefig('Arrhenius_sigma.png', dpi=300)
    else:
        plt.show()
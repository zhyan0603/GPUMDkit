import re
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import get_backend
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
group_sigmaTs_all = []
group_sigmaTs_x = []
group_sigmaTs_y = []
group_sigmaTs_z = []
group_sigmas_all = []
group_sigmas_x = []
group_sigmas_y = []
group_sigmas_z = []
fig, ax1 = plt.subplots(figsize=(5, 3.4), dpi=150)

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
        msd_x = data[:, 1]
        msd_y = data[:, 2]
        msd_z = data[:, 3]
        msd_all = msd_x + msd_y + msd_z
        start = int(len(dts) * 0.4)
        end = int(len(dts) * 0.8)

        # Calculate diffusivity
        k_all, _ = np.polyfit(dts[start:end], msd_all[start:end], 1)
        k_x, _ = np.polyfit(dts[start:end], msd_x[start:end], 1)
        k_y, _ = np.polyfit(dts[start:end], msd_y[start:end], 1)
        k_z, _ = np.polyfit(dts[start:end], msd_z[start:end], 1)

        D_all = k_all * 1e-4 / 6  # Convert A/ps to cm/s
        D_x = k_x * 1e-4 / 2
        D_y = k_y * 1e-4 / 2
        D_z = k_z * 1e-4 / 2

        # Calculate conductivity
        z = 1  # Species charge
        vol_cm3 = volume * 1e-24  # Convert volume to cm^3
        conversion_factor = (num_ions / (vol_cm3 * consts.N_A) * z**2 * (consts.N_A * consts.e)**2 / (consts.R * temp))
        sigma_all = conversion_factor * D_all
        sigma_x = conversion_factor * D_x
        sigma_y = conversion_factor * D_y
        sigma_z = conversion_factor * D_z
        sigmaT_all = sigma_all * temp
        sigmaT_x = sigma_x * temp
        sigmaT_y = sigma_y * temp
        sigmaT_z = sigma_z * temp

        if any(v <= 0 or np.isnan(v) for v in [sigmaT_all, sigmaT_x, sigmaT_y, sigmaT_z]):
            print(f"[Skip] {temp_folder}: Invalid sigma_T values")
            continue

        group_Ts.append(temp)
        group_sigmaTs_all.append(sigmaT_all)
        group_sigmaTs_x.append(sigmaT_x)
        group_sigmaTs_y.append(sigmaT_y)
        group_sigmaTs_z.append(sigmaT_z)
        group_sigmas_all.append(sigma_all)
        group_sigmas_x.append(sigma_x)
        group_sigmas_y.append(sigma_y)
        group_sigmas_z.append(sigma_z)

    except Exception as e:
        print(f"[Error] Processing {temp_folder}: {e}")
        continue

if len(group_Ts) < 2:
    raise ValueError("Insufficient data points for fitting, need at least two temperatures")

# Sort by temperature
sorted_indices = np.argsort(group_Ts)
group_Ts = np.array(group_Ts)[sorted_indices]
group_sigmaTs_all = np.array(group_sigmaTs_all)[sorted_indices]
group_sigmaTs_x = np.array(group_sigmaTs_x)[sorted_indices]
group_sigmaTs_y = np.array(group_sigmaTs_y)[sorted_indices]
group_sigmaTs_z = np.array(group_sigmaTs_z)[sorted_indices]
group_sigmas_all = np.array(group_sigmas_all)[sorted_indices]
group_sigmas_x = np.array(group_sigmas_x)[sorted_indices]
group_sigmas_y = np.array(group_sigmas_y)[sorted_indices]
group_sigmas_z = np.array(group_sigmas_z)[sorted_indices]

# Print conductivity tables
w_t, w_all, w_xyz = 10, 16, 12
line = f"+{'-'*(w_t+2)}-{'-'*(w_all+2)}-{'-'*(w_xyz+2)}-{'-'*(w_xyz+2)}-{'-'*(w_xyz+2)}+"

print(line)
print(f"| {'Conductivity (unit: S/cm)':^{len(line)-4}} |")
print(line)
print(f"| {'T (K)':^{w_t}} | {'Sigma_all':^{w_all}} | {'Sigma_x':^{w_xyz}} | {'Sigma_y':^{w_xyz}} | {'Sigma_z':^{w_xyz}} |")
print(line)
for temp, s_all, s_x, s_y, s_z in zip(group_Ts, group_sigmas_all, group_sigmas_x, group_sigmas_y, group_sigmas_z):
    print(
        f"| {temp:^{w_t}d} | {s_all:^{w_all}.3e} | {s_x:^{w_xyz}.3e} | {s_y:^{w_xyz}.3e} | {s_z:^{w_xyz}.3e} |"
    )
print(line)

print(line)
print(f"| {'Sigma*T (unit: K S/cm)':^{len(line)-4}} |")
print(line)
print(f"| {'T (K)':^{w_t}} | {'SigmaT_all':^{w_all}} | {'SigmaT_x':^{w_xyz}} | {'SigmaT_y':^{w_xyz}} | {'SigmaT_z':^{w_xyz}} |")
print(line)
for temp, s_all, s_x, s_y, s_z in zip(group_Ts, group_sigmaTs_all, group_sigmaTs_x, group_sigmaTs_y, group_sigmaTs_z):
    print(
        f"| {temp:^{w_t}d} | {s_all:^{w_all}.3e} | {s_x:^{w_xyz}.3e} | {s_y:^{w_xyz}.3e} | {s_z:^{w_xyz}.3e} |"
    )
print(line)

# Plot scatter points
ax1.scatter(1000 / group_Ts, np.log(group_sigmaTs_all), color='C4', s=30, label=r'$\sigma T$ total')
ax1.scatter(1000 / group_Ts, np.log(group_sigmaTs_x), color='C0', s=25, label=r'$\sigma T$ x')
ax1.scatter(1000 / group_Ts, np.log(group_sigmaTs_y), color='C2', s=25, label=r'$\sigma T$ y')
ax1.scatter(1000 / group_Ts, np.log(group_sigmaTs_z), color='C3', s=25, label=r'$\sigma T$ z')

# Perform linear fitting
invT = 1000 / np.array(group_Ts)
ln_sigmaT_all = np.log(np.array(group_sigmaTs_all))
ln_sigmaT_x = np.log(np.array(group_sigmaTs_x))
ln_sigmaT_y = np.log(np.array(group_sigmaTs_y))
ln_sigmaT_z = np.log(np.array(group_sigmaTs_z))

k_all, b_all = np.polyfit(invT, ln_sigmaT_all, 1)
k_x, b_x = np.polyfit(invT, ln_sigmaT_x, 1)
k_y, b_y = np.polyfit(invT, ln_sigmaT_y, 1)
k_z, b_z = np.polyfit(invT, ln_sigmaT_z, 1)

Ea_all = -k_all * 1000 * consts.k / consts.e
Ea_x = -k_x * 1000 * consts.k / consts.e
Ea_y = -k_y * 1000 * consts.k / consts.e
Ea_z = -k_z * 1000 * consts.k / consts.e

print(f"\n{cell}, Ea_all: {Ea_all:.3f} eV")
print(f"{cell}, Ea_x: {Ea_x:.3f} eV")
print(f"{cell}, Ea_y: {Ea_y:.3f} eV")
print(f"{cell}, Ea_z: {Ea_z:.3f} eV")

# Plot fitted lines
ax1.plot(invT, k_all * invT + b_all, linestyle='-', linewidth=2, color='C4')
ax1.plot(invT, k_x * invT + b_x, linestyle='--', linewidth=1.5, color='C0')
ax1.plot(invT, k_y * invT + b_y, linestyle='--', linewidth=1.5, color='C2')
ax1.plot(invT, k_z * invT + b_z, linestyle='--', linewidth=1.5, color='C3')

# Extrapolate from lowest temperature to 300K
specified_temperature = 300
x_ext = 1000 / np.arange(min(group_Ts), specified_temperature - 1, -10)  # Start from lowest T to 300K

ax1.plot(x_ext, k_all * x_ext + b_all, linestyle='--', color='C4', alpha=0.5)
ax1.plot(x_ext, k_x * x_ext + b_x, linestyle=':', color='C0', alpha=0.5)
ax1.plot(x_ext, k_y * x_ext + b_y, linestyle=':', color='C2', alpha=0.5)
ax1.plot(x_ext, k_z * x_ext + b_z, linestyle=':', color='C3', alpha=0.5)

ln_sigmaT_300_all = k_all * 1000 / specified_temperature + b_all
ln_sigmaT_300_x = k_x * 1000 / specified_temperature + b_x
ln_sigmaT_300_y = k_y * 1000 / specified_temperature + b_y
ln_sigmaT_300_z = k_z * 1000 / specified_temperature + b_z

sigma_300K_all = np.exp(ln_sigmaT_300_all) / specified_temperature
sigma_300K_x = np.exp(ln_sigmaT_300_x) / specified_temperature
sigma_300K_y = np.exp(ln_sigmaT_300_y) / specified_temperature
sigma_300K_z = np.exp(ln_sigmaT_300_z) / specified_temperature

print(f"at {specified_temperature}K, {cell}: Sigma_all = {sigma_300K_all:.3e} S/cm")
print(f"at {specified_temperature}K, {cell}: Sigma_x = {sigma_300K_x:.3e} S/cm")
print(f"at {specified_temperature}K, {cell}: Sigma_y = {sigma_300K_y:.3e} S/cm")
print(f"at {specified_temperature}K, {cell}: Sigma_z = {sigma_300K_z:.3e} S/cm")

# Set axis labels
ax1.set_xlabel("1000/T (1/K)", fontsize=10)
ax1.set_ylabel(r"ln($\sigma \cdot T$) (K S/cm)", fontsize=10)
ax1.legend(loc='lower left', fontsize=8, frameon=True)

# Add top temperature axis
def tick_function(X):
    V = 1000 / X
    return ["%.0f" % z for z in V]

ax2 = ax1.twiny()
ax2.set_xticks(ax1.get_xticks())
ax2.set_xbound(ax1.get_xbound())
ax2.set_xticklabels(tick_function(ax1.get_xticks()))
ax2.set_xlabel("T (K)", fontsize=10)

# Save and show plot
fig.tight_layout()
if len(sys.argv) > 1 and sys.argv[1] == 'save':
    fig.savefig('Arrhenius_sigma_xyz.png', dpi=300)
else:
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'Arrhenius_sigma_xyz.png'.")
        fig.savefig('Arrhenius_sigma_xyz.png', dpi=300)
    else:
        plt.show()

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as consts

def read_msd_file(msd_file):
    data = np.loadtxt(msd_file)
    dts = data[:, 0]
    msd_x = data[:, 1]
    msd_y = data[:, 2]
    msd_z = data[:, 3]
    msd_all = msd_x + msd_y + msd_z
    return dts, msd_all, msd_x, msd_y, msd_z

# Plot setup
plt.rcParams.update({'font.size': 10})
plt.figure(figsize=(5, 3.4), dpi=200)

# Automatically find *K folders
temps = []
D_alls = []
D_xs = []
D_ys = []
D_zs = []

# You can specify which temperature folders if you don't want to process all
#specified_folders = ['600K', '700K', '800K', '900K']

current_dir = os.getcwd()
for folder in os.listdir(current_dir):
    if folder.endswith('K') and os.path.isdir(folder):
#    if folder.endswith('K') and os.path.isdir(folder) and folder in specified_folders:
        try:
            temp = int(folder[:-1])  # Remove 'K' and convert to int
            msd_file = os.path.join(folder, 'msd.out')
            if os.path.exists(msd_file):
                dts, msd_all, msd_x, msd_y, msd_z = read_msd_file(msd_file)

                start = int(len(dts) * 0.15)
                end = int(len(dts) * 0.3)
                k_all, _ = np.polyfit(dts[start:end], msd_all[start:end], 1)
                k_x, _ = np.polyfit(dts[start:end], msd_x[start:end], 1)
                k_y, _ = np.polyfit(dts[start:end], msd_y[start:end], 1)
                k_z, _ = np.polyfit(dts[start:end], msd_z[start:end], 1)

                D_all = k_all * 1e-4 / 6
                D_x = k_x * 1e-4 / 2
                D_y = k_y * 1e-4 / 2
                D_z = k_z * 1e-4 / 2

                temps.append(temp)
                D_alls.append(D_all)
                D_xs.append(D_x)
                D_ys.append(D_y)
                D_zs.append(D_z)
        except ValueError:
            continue  # Skip if folder name doesn't convert to int properly

# Sort by temperature
sorted_indices = np.argsort(temps)
temps = np.array(temps)[sorted_indices]
D_alls = np.array(D_alls)[sorted_indices]
D_xs = np.array(D_xs)[sorted_indices]
D_ys = np.array(D_ys)[sorted_indices]
D_zs = np.array(D_zs)[sorted_indices]

# Print diffusion coefficients in temperature order
w_t, w_all, w_xyz = 10, 15, 12
line = f"+{'-'*(w_t+2)}-{'-'*(w_all+2)}-{'-'*(w_xyz+2)}-{'-'*(w_xyz+2)}-{'-'*(w_xyz+2)}+"

print(line)
print(f"| {'Diffusion coefficients (unit: cm2/s)':^{len(line)-4}} |")
print(line)
print(f"| {'T (K)':^{w_t}} | {'D_all':^{w_all}} | {'D_x':^{w_xyz}} | {'D_y':^{w_xyz}} | {'D_z':^{w_xyz}} |")
print(line)

for temp, D_all, D_x, D_y, D_z in zip(temps, D_alls, D_xs, D_ys, D_zs):
    print(
        f"| {temp:^{w_t}d} | {D_all:^{w_all}.3e} | {D_x:^{w_xyz}.3e} | {D_y:^{w_xyz}.3e} | {D_z:^{w_xyz}.3e} |"
    )
print(line)

# Plot scatter points
plt.scatter(1000 / temps, np.log10(D_alls), color='C4', s=40, label=r'D$_{all}$')
plt.scatter(1000 / temps, np.log10(D_xs), color='C0', s=30, label=r'D$_x$')
plt.scatter(1000 / temps, np.log10(D_ys), color='C2', s=30, label=r'D$_y$')
plt.scatter(1000 / temps, np.log10(D_zs), color='C3', s=30, label=r'D$_z$')

# Fit and calculate activation energy
if len(temps) > 1:
    k_all, b_all = np.polyfit(1000 / temps, np.log10(D_alls), 1)
    k_x, b_x = np.polyfit(1000 / temps, np.log10(D_xs), 1)
    k_y, b_y = np.polyfit(1000 / temps, np.log10(D_ys), 1)
    k_z, b_z = np.polyfit(1000 / temps, np.log10(D_zs), 1)

    Ea_all = -k_all * 1000 * (consts.k / consts.e) * np.log(10)
    Ea_x = -k_x * 1000 * (consts.k / consts.e) * np.log(10)
    Ea_y = -k_y * 1000 * (consts.k / consts.e) * np.log(10)
    Ea_z = -k_z * 1000 * (consts.k / consts.e) * np.log(10)

    print(f"Ea_all: {Ea_all:.3f} eV")
    print(f"Ea_x: {Ea_x:.3f} eV")
    print(f"Ea_y: {Ea_y:.3f} eV")
    print(f"Ea_z: {Ea_z:.3f} eV")

    # Plot fit lines
    plt.plot(1000 / temps, k_all * (1000 / temps) + b_all,
             linestyle='-', linewidth=2, color='C4')
    plt.plot(1000 / temps, k_x * (1000 / temps) + b_x,
             linestyle='--', linewidth=1.5, color='C0')
    plt.plot(1000 / temps, k_y * (1000 / temps) + b_y,
             linestyle='--', linewidth=1.5, color='C2')
    plt.plot(1000 / temps, k_z * (1000 / temps) + b_z,
             linestyle='--', linewidth=1.5, color='C3')

# Axes setup
plt.xlabel('1000/T (1/K)', fontsize=11)
plt.ylabel(r'log10(D) (cm$^2$/s)')
plt.legend(loc='lower left', fontsize=9)

# Secondary x-axis
ax1 = plt.gca()
ax2 = ax1.twiny()
ax2.set_xticks(ax1.get_xticks())
ax2.set_xbound(ax1.get_xbound())
ax2.set_xticklabels([f"{1000/x:.0f}" for x in ax1.get_xticks()])
ax2.set_xlabel("T (K)", fontsize=11)

plt.tight_layout()
if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('Arrhenius_D_xyz.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'Arrhenius_D_xyz.png'.")
        plt.savefig('Arrhenius_D_xyz.png', dpi=300)
    else:
        plt.show()

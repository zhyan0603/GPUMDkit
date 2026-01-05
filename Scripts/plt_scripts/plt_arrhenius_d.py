import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as consts

def read_msd_file(msd_file):
    data = np.loadtxt(msd_file)
    dts = data[:, 0]
    msds_all = np.sum(data[:, 1:4], axis=1)
    return dts, msds_all

# Plot setup
plt.rcParams.update({'font.size': 10})
plt.figure(figsize=(5, 3.4), dpi=200)

# Automatically find *K folders
temps = []
D_alls = []

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
                dts, msds_all = read_msd_file(msd_file)
                
                start = int(len(dts) * 0.15)
                end = int(len(dts) * 0.3)
                k_all, _ = np.polyfit(dts[start:end], msds_all[start:end], 1)
                D_all = k_all * 1e-4 / 6
                temps.append(temp)
                D_alls.append(D_all)
                print(f"T: {temp}K, D_total: {D_all:.3e} cm2/s")
                plt.scatter(1000 / temp, np.log10(D_all), color='C4', s=40)
        except ValueError:
            continue  # Skip if folder name doesn't convert to int properly

# Sort by temperature
sorted_indices = np.argsort(temps)
temps = np.array(temps)[sorted_indices]
D_alls = np.array(D_alls)[sorted_indices]

# Fit and calculate activation energy
if len(temps) > 1:
    k_all, b_all = np.polyfit(1000 / temps, np.log10(D_alls), 1)
    Ea_all = -k_all * 1000 * (consts.k / consts.e) * np.log(10)
    print(f"Ea: {Ea_all:.3f} eV")
    
    # Plot fit line
    plt.plot(1000 / temps, k_all * (1000 / temps) + b_all,
             linestyle='-', linewidth=2, color='C4', label=r'D$_{total}$')

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
    plt.savefig('Arrhenius_D.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'Arrhenius_D.png'.")
        plt.savefig('Arrhenius_D.png', dpi=300)
    else:
        plt.show()
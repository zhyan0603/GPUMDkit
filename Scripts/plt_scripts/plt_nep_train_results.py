import sys
import matplotlib.pyplot as plt
import numpy as np
# from cycler import cycler

# custom_colors = ['#9BBBE1', '#EAB883', '#A9CA70', '#DD7C4F', '#F09BA0', '#B58C9A']
# plt.rcParams['axes.prop_cycle'] = cycler(color=custom_colors)

# Load data
loss = np.loadtxt('loss.out')
energy_data = np.loadtxt('energy_train.out')
force_data = np.loadtxt('force_train.out')
#virial_data = np.loadtxt('virial_train.out')
stress_data = np.loadtxt('stress_train.out')

# Function to calculate RMSE
def calculate_rmse(pred, actual):
    return np.sqrt(np.mean((pred - actual) ** 2))

# Function to calculate dynamic axis limits
def calculate_limits(train_data, padding=0.08):
    data_min = np.min(train_data)
    data_max = np.max(train_data)
    data_range = data_max - data_min
    return data_min - padding * data_range, data_max + padding * data_range

# Create a subplot with 2 row and 2 columns
fig, axs = plt.subplots(2, 2, figsize=(9, 7), dpi=100)

# Plotting the loss figure
axs[0, 0].loglog(loss[:, 1:7], '-', linewidth=2)
xlabel = 'Generation/100' if loss[1, 0] - loss[0, 0] == 100 else 'Epoch'
axs[0, 0].set_xlabel(xlabel, fontsize=10)
axs[0, 0].set_ylabel('Loss functions', fontsize=10)
axs[0, 0].tick_params(axis='both', labelsize=10)
axs[0, 0].legend(['Total', 'L1-Reg', 'L2-Reg', 'Energy-train', 'Force-train', 'Virial-train'], prop = {'size':8})
axs[0, 0].axis('tight')
#axs[0, 0].text(-0.07, 1.03, "(a)", transform=axs[0, 0].transAxes, fontsize=12, va='top', ha='right')

# Plotting the energy_data figure
xmin_energy, xmax_energy = calculate_limits(energy_data[:, 1])
axs[0, 1].set_xlim(xmin_energy, xmax_energy)
axs[0, 1].set_ylim(xmin_energy, xmax_energy)
axs[0, 1].plot(energy_data[:, 1], energy_data[:, 0], '.', markersize=10)
axs[0, 1].plot([xmin_energy, xmax_energy], [xmin_energy, xmax_energy], linewidth=2, color='grey', linestyle='--')
axs[0, 1].set_xlabel('DFT energy (eV/atom)', fontsize=10)
axs[0, 1].set_ylabel('NEP energy (eV/atom)', fontsize=10)
axs[0, 1].tick_params(axis='both', labelsize=10)
axs[0, 1].legend(['energy'])
axs[0, 1].axis('tight')

# Calculate and display RMSE for energy
energy_rmse = calculate_rmse(energy_data[:, 0], energy_data[:, 1]) * 1000
axs[0, 1].text(0.5, 0.08, f'RMSE: {energy_rmse:.2f} meV/atom', transform=axs[0, 1].transAxes, fontsize=10, verticalalignment='center')
#axs[0, 1].text(-0.07, 1.03, "(b)", transform=axs[0, 1].transAxes, fontsize=12, va='top', ha='right')

# Plotting the force_data figure
xmin_force, xmax_force = calculate_limits(force_data[:, 3:6].reshape(-1))
axs[1, 0].set_xlim(xmin_force, xmax_force)
axs[1, 0].set_ylim(xmin_force, xmax_force)
axs[1, 0].plot(force_data[:, 3:6], force_data[:, 0:3], '.', markersize=10)
axs[1, 0].plot([xmin_force, xmax_force], [xmin_force, xmax_force], linewidth=2, color='grey', linestyle='--')
axs[1, 0].set_xlabel(r'DFT force (eV/$\mathrm{\AA}$)', fontsize=10)
axs[1, 0].set_ylabel(r'NEP force (eV/$\mathrm{\AA}$)', fontsize=10)
axs[1, 0].tick_params(axis='both', labelsize=10)
axs[1, 0].legend(['fx', 'fy', 'fz'])
axs[1, 0].axis('tight')

# Calculate and display RMSE for forces
force_rmse = [calculate_rmse(force_data[:, i], force_data[:, i + 3]) for i in range(3)]
mean_force_rmse = np.mean(force_rmse) * 1000
axs[1, 0].text(0.5, 0.08, rf'RMSE: {mean_force_rmse:.2f} meV/$\mathrm{{\AA}}$', transform=axs[1, 0].transAxes, fontsize=10, verticalalignment='center')
#axs[1, 0].text(-0.07, 1.03, "(c)", transform=axs[1, 0].transAxes, fontsize=12, va='top', ha='right')

# Plotting the stress figure
xmin_stress, xmax_stress = calculate_limits(stress_data[:, 6:12].reshape(-1))
axs[1, 1].set_xlim(xmin_stress, xmax_stress)
axs[1, 1].set_ylim(xmin_stress, xmax_stress)
axs[1, 1].plot(stress_data[:, 6:12], stress_data[:, 0:6], '.', markersize=10)
axs[1, 1].plot([xmin_stress, xmax_stress], [xmin_stress, xmax_stress], linewidth=2, color='grey', linestyle='--')
axs[1, 1].set_xlabel('DFT stress (GPa)', fontsize=10)
axs[1, 1].set_ylabel('NEP stress (GPa)', fontsize=10)
axs[1, 1].tick_params(axis='both', labelsize=10)
axs[1, 1].legend(['xx', 'yy', 'zz', 'xy', 'yz', 'zx'])
axs[1, 1].axis('tight')

# Calculate and display RMSE for stresses
stress_rmse = [calculate_rmse(stress_data[:, i], stress_data[:, i + 6]) for i in range(6)]
mean_stress_rmse = np.mean(stress_rmse) 
axs[1, 1].text(0.5, 0.08, f'RMSE: {mean_stress_rmse:.4f} GPa', transform=axs[1, 1].transAxes, fontsize=10, verticalalignment='center')
#axs[1, 1].text(-0.07, 1.03, "(d)", transform=axs[1, 1].transAxes, fontsize=12, va='top', ha='right')

# Adjust layout for better spacing
plt.tight_layout()
#fig.subplots_adjust(top=0.968,bottom=0.088,left=0.086,right=0.983,hspace=0.22,wspace=0.24)

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('train.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'train.png'.")
        plt.savefig('train.png', dpi=300)
    else:
        plt.show()

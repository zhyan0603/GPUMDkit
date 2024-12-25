import sys
import matplotlib.pyplot as plt
#import scienceplots
import numpy as np

# Set the plotting style
#plt.style.use(['retro'])

# Load data
energy_train = np.loadtxt('energy_train.out')
force_train = np.loadtxt('force_train.out')
stress_train = np.loadtxt('stress_train.out')
energy_test = np.loadtxt('energy_test.out')
force_test = np.loadtxt('force_test.out')
stress_test = np.loadtxt('stress_test.out')

# Function to calculate RMSE
def calculate_rmse(pred, actual):
    return np.sqrt(np.mean((pred - actual) ** 2))

# Function to calculate dynamic axis limits
def calculate_limits(train_data, test_data, padding=0.08):
    data_min = np.min(np.concatenate((train_data, test_data)))
    data_max = np.max(np.concatenate((train_data, test_data)))
    data_range = data_max - data_min
    return data_min - padding * data_range, data_max + padding * data_range

# Create a subplot with 1 row and 3 columns
fig, axs = plt.subplots(1, 3, figsize=(12, 3.3), dpi=100)

# Energy plot
xmin_energy, xmax_energy = calculate_limits(energy_train[:, 1], energy_test[:, 1])
axs[0].set_xlim(xmin_energy, xmax_energy)
axs[0].set_ylim(xmin_energy, xmax_energy)
axs[0].plot(energy_train[:, 1], energy_train[:, 0], '.', markersize=10, label='Train', color='#237B9F') #377EB9
axs[0].plot(energy_test[:, 1], energy_test[:, 0], '.', markersize=10, label='Test', color='#EC817E')  #F48892
axs[0].plot([xmin_energy, xmax_energy], [xmin_energy, xmax_energy], linewidth=1, color='red')
axs[0].set_xlabel('DFT energy (eV/atom)', fontsize=10)
axs[0].set_ylabel('NEP energy (eV/atom)', fontsize=10)
axs[0].legend(frameon=False)
axs[0].tick_params(axis='both', labelsize=10)

# Calculate and display RMSE for energy
energy_train_rmse = calculate_rmse(energy_train[:, 0], energy_train[:, 1]) * 1000
energy_test_rmse = calculate_rmse(energy_test[:, 0], energy_test[:, 1]) * 1000
axs[0].text(0.3, 0.2, f'RMSE (Train): {energy_train_rmse:.2f} meV/atom', transform=axs[0].transAxes, fontsize=10, verticalalignment='center')
axs[0].text(0.3, 0.1, f'RMSE (Test): {energy_test_rmse:.2f} meV/atom', transform=axs[0].transAxes, fontsize=10, verticalalignment='center')
#axs[0].text(-0.1, 1.03, "(a)", transform=axs[0].transAxes, fontsize=13, va='top', ha='right')

# Force plot
xmin_force, xmax_force = calculate_limits(force_train[:, 3:6].reshape(-1), force_test[:, 3:6].reshape(-1))
axs[1].set_xlim(xmin_force, xmax_force)
axs[1].set_ylim(xmin_force, xmax_force)
axs[1].plot(force_train[:, 3:6].reshape(-1), force_train[:, 0:3].reshape(-1), '.', markersize=10, label='Train', color='#237B9F')
axs[1].plot(force_test[:, 3:6].reshape(-1), force_test[:, 0:3].reshape(-1), '.', markersize=10, label='Test', color='#EC817E')
axs[1].plot([xmin_force, xmax_force], [xmin_force, xmax_force], linewidth=1, color='red')
axs[1].set_xlabel(r'DFT force (eV/Å)', fontsize=10)
axs[1].set_ylabel(r'NEP force (eV/Å)', fontsize=10)
axs[1].tick_params(axis='both', labelsize=10)
axs[1].legend(frameon=False)

# Calculate and display RMSE for forces
force_train_rmse = [calculate_rmse(force_train[:, i], force_train[:, i + 3]) for i in range(3)]
force_test_rmse = [calculate_rmse(force_test[:, i], force_test[:, i + 3]) for i in range(3)]
mean_force_train_rmse = np.mean(force_train_rmse) * 1000
mean_force_test_rmse = np.mean(force_test_rmse) * 1000
axs[1].text(0.35, 0.2, rf'RMSE (Train): {mean_force_train_rmse:.2f} meV/Å', transform=axs[1].transAxes, fontsize=10, verticalalignment='center')
axs[1].text(0.35, 0.1, rf'RMSE (Test): {mean_force_test_rmse:.2f} meV/Å', transform=axs[1].transAxes, fontsize=10, verticalalignment='center')
#axs[1].text(-0.1, 1.03, "(b)", transform=axs[1].transAxes, fontsize=13, va='top', ha='right')

# Stress plot
xmin_stress, xmax_stress = calculate_limits(stress_train[:, 6:12].reshape(-1), stress_test[:, 6:12].reshape(-1))
axs[2].set_xlim(xmin_stress, xmax_stress)
axs[2].set_ylim(xmin_stress, xmax_stress)
axs[2].plot(stress_train[:, 6:12].reshape(-1), stress_train[:, 0:6].reshape(-1), '.', markersize=10, label='Train', color='#237B9F')
axs[2].plot(stress_test[:, 6:12].reshape(-1), stress_test[:, 0:6].reshape(-1), '.', markersize=10, label='Test', color='#EC817E')
axs[2].plot([xmin_stress, xmax_stress], [xmin_stress, xmax_stress], linewidth=1, color='red')
axs[2].set_xlabel('DFT stress (GPa)', fontsize=10)
axs[2].set_ylabel('NEP stress (GPa)', fontsize=10)
axs[2].tick_params(axis='both', labelsize=10)
axs[2].legend(frameon=False)

# Calculate and display RMSE for stresses
stress_train_rmse = [calculate_rmse(stress_train[:, i], stress_train[:, i + 6]) for i in range(6)]
stress_test_rmse = [calculate_rmse(stress_test[:, i], stress_test[:, i + 6]) for i in range(6)]
mean_stress_train_rmse = np.mean(stress_train_rmse)
mean_stress_test_rmse = np.mean(stress_test_rmse)
axs[2].text(0.4, 0.2, f'RMSE (Train): {mean_stress_train_rmse:.3f} GPa', transform=axs[2].transAxes, fontsize=10, verticalalignment='center')
axs[2].text(0.4, 0.1, f'RMSE (Test): {mean_stress_test_rmse:.3f} GPa', transform=axs[2].transAxes, fontsize=10, verticalalignment='center')
#axs[2].text(-0.1, 1.03, "(c)", transform=axs[2].transAxes, fontsize=13, va='top', ha='right')

# Adjust layout for better spacing
plt.tight_layout()
fig.subplots_adjust(top=0.968, bottom=0.16, left=0.086, right=0.983, hspace=0.2, wspace=0.25)

# Show plot
if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('train_test.png', dpi=300)
else:
    plt.show()

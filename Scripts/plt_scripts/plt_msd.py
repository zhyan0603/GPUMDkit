import sys
import numpy as np
import matplotlib.pyplot as plt

# Function to read data from the given file
def read_data(file_name):
    data = np.loadtxt(file_name)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# Function to calculate the slope of the first 50% of the data
def calculate_slope(x, y):
    half_idx = len(x) // 2
    x_half, y_half = x[:half_idx], y[:half_idx]
    coeffs = np.polyfit(x_half, y_half, 1)  # Linear fit
    return coeffs[0]  # Slope of the linear fit

# Determine the input file
#input_file = sys.argv[1] if len(sys.argv) > 1 else './msd.out'
input_file = './msd.out'

# Read the data from the file
time, msd_x, msd_y, msd_z = read_data(input_file)

slope_x = calculate_slope(time, msd_x)
slope_y = calculate_slope(time, msd_y)
slope_z = calculate_slope(time, msd_z)

# Calculate the total MSD and mean MSD
#msd_total = msd_x + msd_y + msd_z
#msd_mean = (msd_x + msd_y + msd_z)/3

# Create a plot for the MSD data
plt.figure(figsize=(4.5, 3.5), dpi=150)
plt.plot(time, msd_x, label='x')
plt.plot(time, msd_y, label='y')
plt.plot(time, msd_z, label='z')
#plt.plot(time, msd_mean, label='mean', color='C4')
#plt.plot(time, msd_total, label='total', color='C5')

plt.text(0.5, 0.3, f'Slope (x): {slope_x:.3f} Å$^2$/ps', fontsize=10, color='C0', verticalalignment='top', transform=plt.gca().transAxes)
plt.text(0.5, 0.21, f'Slope (y): {slope_y:.3f} Å$^2$/ps', fontsize=10, color='C1', verticalalignment='top', transform=plt.gca().transAxes)
plt.text(0.5, 0.12, f'Slope (z): {slope_z:.3f} Å$^2$/ps', fontsize=10, color='C2', verticalalignment='top', transform=plt.gca().transAxes)

# Add titles and labels
#plt.title('MSD vs dt')
plt.xlabel('dt (ps)')
plt.ylabel(r'MSD ($\AA^2$)')
plt.legend()
plt.tight_layout()
#plt.grid(True)

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('msd.png', dpi=300)
else:
    plt.show()
"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    plt_sdc.py
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# Function to read data from the given file
def read_data(file_name):
    data = np.loadtxt(file_name)
    return data[:, 0], data[:, 4], data[:, 5], data[:, 6]

input_file = './msd.out'

# Read the data from the file
time, sdc_x, sdc_y, sdc_z = read_data(input_file)

# Calculate the total MSD and mean MSD
#sdc_total = sdc_x + sdc_y + sdc_z
#sdc_mean = (sdc_x + sdc_y + sdc_z)/3

# Create a plot for the MSD data
plt.figure(figsize=(4.5, 3.5), dpi=150)
plt.plot(time, sdc_x, label='x')
plt.plot(time, sdc_y, label='y')
plt.plot(time, sdc_z, label='z')
#plt.plot(time, sdc_mean, label='mean', color='C4')
#plt.plot(time, sdc_total, label='total', color='C5')

#plt.title('SDC vs dt')
plt.xlabel('dt (ps)')
plt.ylabel(r'SDC ($\AA^2$/ps)')
plt.legend()
plt.tight_layout()
#plt.grid(True)

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('sdc.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'sdc.png'.")
        plt.savefig('sdc.png', dpi=300)
    else:
        plt.show()
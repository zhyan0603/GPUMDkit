import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ase.io import read
from collections import defaultdict

def load_xyz_and_charges(xyz_file, charge_file):
    """
    Load atomic charges from charge.out and match with atom symbols from an extxyz file.
    
    Args:
        xyz_file (str): Path to the extxyz file (e.g., train.xyz).
        charge_file (str): Path to the charge file (e.g., charge.out).
    
    Returns:
        dict: Dictionary with element symbols as keys and lists of corresponding charges as values.
    """
    # Read extxyz file
    structures = read(xyz_file, index=':', format='extxyz')
    
    # Read charge file
    charges = np.loadtxt(charge_file, ndmin=1)
    
    # Verify that the number of charges matches the total number of atoms
    total_atoms = sum(len(struct) for struct in structures)
    if len(charges) != total_atoms:
        raise ValueError(f"Charge data length ({len(charges)}) does not match total atoms in XYZ ({total_atoms})")
    
    # Collect charges by element
    charges_by_element = defaultdict(list)
    charge_index = 0
    for struct in structures:
        symbols = struct.get_chemical_symbols()  # Get element symbols for the structure
        for symbol in symbols:
            charges_by_element[symbol].append(charges[charge_index])
            charge_index += 1
    
    return charges_by_element

# File paths
xyz_file = "train.xyz"
charge_file = "charge.out"

# Load data and plot
charges_by_element = load_xyz_and_charges(xyz_file, charge_file)

plt.figure(figsize=(5, 3.5), dpi=150)
# Plot charge distribution for each element
for element, charge_values in charges_by_element.items():
    sns.kdeplot(charge_values, fill=True, label=element, alpha=0.6)

plt.xlabel("Charge (e)")
plt.ylabel("Density")
plt.legend()
#plt.grid()
plt.tight_layout()

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('charge.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'charge.png'.")
        plt.savefig('charge.png', dpi=300)
    else:
        plt.show()
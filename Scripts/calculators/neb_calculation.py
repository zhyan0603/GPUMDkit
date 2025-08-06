"""
Performing NEB calculations with NEP model.

Usage:
    python neb_calculation.py <initial_structure> <final_structure> <n_image> <nep_model>

Arguments:
    initial_structure : Path to the initial structure XYZ file
    final_structure   : Path to the final structure XYZ file
    n_image           : Number of intermediate images
    nep_model         : Path to the NEP model file (nep.txt)

Author:
    Zhoulin LIU <>

Modified by Zihan YAN <yanzihan@westlake.edu.cn>
"""
import os
import sys
import numpy as np
from ase.io import read
from ase.mep import NEB
from ase.optimize.fire import FIRE as QuasiNewton
from calorine.calculators import CPUNEP
import matplotlib.pyplot as plt
from ase.constraints import FixAtoms
from ase.io import Trajectory, write

# Check command-line arguments
if len(sys.argv) != 5:
    print("python neb_calculation.py <initial_structure> <final_structure> <n_image> <nep_model>")
    sys.exit(1)

initial_file = sys.argv[1]
final_file = sys.argv[2]
n_images = int(sys.argv[3])
nep_model = sys.argv[4]

# Read initial and final structures
initial = read(initial_file)
final = read(final_file)

# Create intermediate images
images = [initial] + [initial.copy() for _ in range(n_images)] + [final]
print(f"Created {len(images)} images (1 initial + {n_images} intermediate + 1 final)")

# Get fixing method interactively (only once, applied to all images)
fix_method = None
fixed_indices = None
fixed_element = None
fixed_range = None

if len(images) > 0:
    # Prompt for fixing method
    while fix_method not in ['none', 'index', 'element', 'position']:
        print("\nAvailable atom fixing methods:")
        print("  none: No atoms will be fixed")
        print("  index: Fix atoms by their 0-based indices")
        print("  element: Fix all atoms of a specific element")
        print("  position: Fix atoms within a specified x, y, z range")
        fix_method = input("Enter fixing method (none/index/element/position): ").strip().lower()
        if fix_method not in ['none', 'index', 'element', 'position']:
            print("Invalid choice. Please choose: none, index, element, or position")

    # Handle additional input based on fixing method
    if fix_method == 'index':
        while True:
            try:
                indices_input = input("Enter atom indices to fix (0-based, space-separated, e.g., '0 1 2'): ").strip()
                fixed_indices = [int(idx) for idx in indices_input.split()]
                # Validate indices
                for idx in fixed_indices:
                    if idx < 0 or idx >= len(images[0]):
                        raise ValueError(f"Index {idx} is out of range (0 to {len(images[0])-1})")
                break
            except ValueError as e:
                print(f"Error: {e}. Please enter valid integer indices.")
    elif fix_method == 'element':
        while True:
            fixed_element = input("Enter element symbol to fix (e.g., Mg): ").strip()
            # Check if element exists in the structure
            if any(atom.symbol == fixed_element for atom in images[0]):
                break
            print(f"Error: Element {fixed_element} not found in structure. Try again.")
    elif fix_method == 'position':
        # Display coordinate ranges of the initial structure
        positions = images[0].get_positions()
        x_coords, y_coords, z_coords = positions[:, 0], positions[:, 1], positions[:, 2]
        print(f"\nCoordinate ranges of the initial structure (in Ångströms):")
        print(f"  X: {min(x_coords):.3f} to {max(x_coords):.3f}")
        print(f"  Y: {min(y_coords):.3f} to {max(y_coords):.3f}")
        print(f"  Z: {min(z_coords):.3f} to {max(z_coords):.3f}")
        while True:
            try:
                range_input = input("Enter x_min x_max y_min y_max z_min z_max (space-separated floats): ").strip()
                fixed_range = [float(x) for x in range_input.split()]
                if len(fixed_range) != 6:
                    raise ValueError("Exactly 6 values are required for x_min, x_max, y_min, y_max, z_min, z_max")
                # Validate ranges
                x_min, x_max, y_min, y_max, z_min, z_max = fixed_range
                if x_min > x_max or y_min > y_max or z_min > z_max:
                    raise ValueError("Min values must be less than or equal to max values")
                break
            except ValueError as e:
                print(f"Error: {e}. Please enter 6 valid floats.")

# Set up NEP calculator and constraints for each image
print("\nSetting up NEP calculator...")
# Prepare to save fixed atom indices
os.makedirs('images', exist_ok=True)
fixed_atoms_file = 'images/fixed_atoms.txt'
with open(fixed_atoms_file, 'w') as f:
    f.write("Fixed Atoms for Each Image\n")
    f.write("-" * 30 + "\n")

for i, image in enumerate(images):
    image.calc = CPUNEP(nep_model)

    # Apply constraints based on fix_method
    if fix_method == 'none':
        print(f"Image {i}: No atoms fixed")
        with open(fixed_atoms_file, 'a') as f:
            f.write(f"Image {i}: No atoms fixed\n")
    elif fix_method == 'index':
        fixed_mask = [False] * len(image)
        for idx in fixed_indices:
            fixed_mask[idx] = True
        image.set_constraint(FixAtoms(mask=fixed_mask))
        print(f"Image {i}: Fixed {sum(fixed_mask)} atoms by index: {fixed_indices}")
        with open(fixed_atoms_file, 'a') as f:
            f.write(f"Image {i}: Fixed atoms by index: {fixed_indices}\n")
    elif fix_method == 'element':
        fixed_mask = [atom.symbol == fixed_element for atom in image]
        fixed_indices = [idx for idx, fixed in enumerate(fixed_mask) if fixed]
        image.set_constraint(FixAtoms(mask=fixed_mask))
        print(f"Image {i}: Fixed {sum(fixed_mask)} atoms of element {fixed_element}: {fixed_indices}")
        with open(fixed_atoms_file, 'a') as f:
            f.write(f"Image {i}: Fixed atoms of element {fixed_element}: {fixed_indices}\n")
    elif fix_method == 'position':
        x_min, x_max, y_min, y_max, z_min, z_max = fixed_range
        positions = image.get_positions()
        fixed_mask = [(x_min <= pos[0] <= x_max) and (y_min <= pos[1] <= y_max) and (z_min <= pos[2] <= z_max) 
                     for pos in positions]
        fixed_indices = [idx for idx, fixed in enumerate(fixed_mask) if fixed]
        image.set_constraint(FixAtoms(mask=fixed_mask))
        print(f"Image {i}: Fixed {sum(fixed_mask)} atoms within position range: {fixed_indices}")
        with open(fixed_atoms_file, 'a') as f:
            f.write(f"Image {i}: Fixed atoms within position range (x: {x_min} to {x_max}, y: {y_min} to {y_max}, z: {z_min} to {z_max}): {fixed_indices}\n")

    print(f"Image {i} initialization completed, number of atoms: {len(image)}")

# Initialize NEB
print(f"\nInitializing NEB path with {len(images)} images...")
neb = NEB(
    images,
    climb=True,
    k=0.5,
    parallel=True
)

# Interpolate path with constraints applied
neb.interpolate(method='idpp', apply_constraint=True)
print("Path interpolation completed, starting NEB optimization...")

# Run NEB optimization
optimizer = QuasiNewton(neb, trajectory='neb.traj', logfile='neb.log')
optimizer.run(fmax=0.01)

# Restore standard output
sys.stdout = open('/dev/stdout', 'w')

print("NEB optimization completed!")

# Save NEB trajectory to a single XYZ file
neb_traj = Trajectory('neb.traj')
write('neb.xyz', neb_traj)

# Analyze and plot NEB energy profile
e0 = initial.get_potential_energy()
print("\nEnergy Profile:")
print(f"{'Image':<5} {'Energy (eV)':<15} {'Delta E (eV)':<10}")

barrier_height = 0
transition_state_idx = 0
energies = []

for i, img in enumerate(images):
    energy = img.get_potential_energy()
    delta_energy = energy - e0
    energies.append(delta_energy)

    if delta_energy > barrier_height:
        barrier_height = delta_energy
        transition_state_idx = i

    print(f"{i:<5} {energy:<15.6f} {delta_energy:<10.6f}")

# Output key results
print("\nKey Results:")
print("-" * 40)
print(f"{'IS energy':<20}: {e0:>10.6f} eV")
print(f"{'TS energy':<20}: {images[transition_state_idx].get_potential_energy():>10.6f} eV")
print(f"{'FS energy':<20}: {images[-1].get_potential_energy():>10.6f} eV")
print(f"{'Reaction energy':<20}: {energies[-1]:>10.6f} eV")
print(f"{'Energy barrier':<20}: {barrier_height:>10.6f} eV")
print(f"{'Transition state':<20}: image_{transition_state_idx}.xyz")
print("-" * 40)

# Save structures
os.makedirs('images', exist_ok=True)
transition_state = images[transition_state_idx]
transition_state.write('images/transition_state.xyz')
print("Transition state saved as images/transition_state.xyz")

for i, img in enumerate(images):
    img.write(f'images/image_{i}.xyz')
print("All NEB images saved as images/image_*.xyz")

# Plot: Image index vs relative energy
plt.figure(figsize=(6, 4))
image_indices = list(range(len(images)))
plt.plot(image_indices, energies, 'o-', c='C0', linewidth=2, alpha=0.8)

plt.xlabel('Image Index')
plt.ylabel('Delta Energy (eV)')
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig('nep_barrier.png', dpi=300)
print("Energy profile plot saved as nep_barrier.png")
print(f"Fixed atoms indices saved to {fixed_atoms_file}")
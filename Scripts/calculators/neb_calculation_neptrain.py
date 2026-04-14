"""
Performing NEB calculations with NEP model (Using NepTrainKit).

Usage:
    python neb_calculation.py <initial_structure> <final_structure> <n_image> <nep_model>
"""
import os
import sys
import numpy as np
from ase.io import read
from ase.mep import NEB
from ase.optimize.fire import FIRE as QuasiNewton
from NepTrainKit.core.calculator import NepAseCalculator
import matplotlib.pyplot as plt
from ase.constraints import FixAtoms
from ase.io import Trajectory, write
from loguru import logger

logger.remove()
logger.add(sys.stderr, level="INFO")
# ----------------------------------------------

if len(sys.argv) != 5:
    print("python neb_calculation.py <initial_structure> <final_structure> <n_image> <nep_model>")
    sys.exit(1)

initial_file = sys.argv[1]
final_file = sys.argv[2]
n_images = int(sys.argv[3])
nep_model = sys.argv[4]

initial = read(initial_file)
final = read(final_file)

images = [initial] + [initial.copy() for _ in range(n_images)] + [final]
print(f"Created {len(images)} images (1 initial + {n_images} intermediate + 1 final)")

fix_method = None
fixed_indices = None
fixed_element = None
fixed_range = None

if len(images) > 0:
    while fix_method not in ['none', 'index', 'element', 'position']:
        print("\nAvailable atom fixing methods:")
        print("  none: No atoms will be fixed")
        print("  index: Fix atoms by their 0-based indices")
        print("  element: Fix all atoms of a specific element")
        print("  position: Fix atoms within a specified x, y, z range")
        fix_method = input("Enter fixing method (none/index/element/position): ").strip().lower()
        if fix_method not in ['none', 'index', 'element', 'position']:
            print("Invalid choice. Please choose: none, index, element, or position")

    if fix_method == 'index':
        while True:
            try:
                indices_input = input("Enter atom indices to fix (0-based, space-separated): ").strip()
                fixed_indices = [int(idx) for idx in indices_input.split()]
                for idx in fixed_indices:
                    if idx < 0 or idx >= len(images[0]):
                        raise ValueError(f"Index {idx} is out of range")
                break
            except ValueError as e:
                print(f"Error: {e}. Please try again.")
    elif fix_method == 'element':
        while True:
            fixed_element = input("Enter element symbol to fix (e.g., Mg): ").strip()
            if any(atom.symbol == fixed_element for atom in images[0]):
                break
            print(f"Error: Element {fixed_element} not found. Try again.")
    elif fix_method == 'position':
        positions = images[0].get_positions()
        print(f"\nCoordinate ranges (Å):")
        print(f"  X: {min(positions[:,0]):.3f} to {max(positions[:,0]):.3f}")
        print(f"  Y: {min(positions[:,1]):.3f} to {max(positions[:,1]):.3f}")
        print(f"  Z: {min(positions[:,2]):.3f} to {max(positions[:,2]):.3f}")
        while True:
            try:
                range_input = input("Enter x_min x_max y_min y_max z_min z_max: ").strip()
                fixed_range = [float(x) for x in range_input.split()]
                if len(fixed_range) != 6:
                    raise ValueError("Need 6 values")
                break
            except ValueError as e:
                print(f"Error: {e}. Please try again.")

print("\nSetting up calculator...")
os.makedirs('images', exist_ok=True)
fixed_atoms_file = 'images/fixed_atoms.txt'
with open(fixed_atoms_file, 'w') as f:
    f.write("Fixed Atoms for Each Image\n")

for i, image in enumerate(images):
    image.calc = NepAseCalculator(model_file=nep_model)
    
    if fix_method == 'index':
        fixed_mask = [False] * len(image)
        for idx in fixed_indices:
            fixed_mask[idx] = True
        image.set_constraint(FixAtoms(mask=fixed_mask))
        print(f"Image {i}: Fixed {sum(fixed_mask)} atoms by index")
    elif fix_method == 'element':
        fixed_mask = [atom.symbol == fixed_element for atom in image]
        image.set_constraint(FixAtoms(mask=fixed_mask))
        print(f"Image {i}: Fixed {sum(fixed_mask)} atoms of {fixed_element}")
    elif fix_method == 'position':
        x_min, x_max, y_min, y_max, z_min, z_max = fixed_range
        positions = image.get_positions()
        fixed_mask = [(x_min <= pos[0] <= x_max) and 
                     (y_min <= pos[1] <= y_max) and 
                     (z_min <= pos[2] <= z_max) for pos in positions]
        image.set_constraint(FixAtoms(mask=fixed_mask))
        print(f"Image {i}: Fixed {sum(fixed_mask)} atoms by position")
    else:
        print(f"Image {i}: No atoms fixed")

print(f"\nInitializing NEB path...")
neb = NEB(images, climb=True, k=0.5, parallel=False)
neb.interpolate(method='idpp', apply_constraint=True)

print("Checking initial energies...")
for i, img in enumerate(images):
    e = img.get_potential_energy()
    print(f"  Image {i}: {e:.4f} eV")

print("\n" + "="*60)
print("Starting NEB optimization...")
print("="*60 + "\n")

optimizer = QuasiNewton(neb, trajectory='neb.traj', logfile='-')
optimizer.run(fmax=0.01)

print("\n" + "="*60)
print("NEB optimization completed!")
print("="*60 + "\n")

neb_traj = Trajectory('neb.traj')
write('neb.xyz', neb_traj)

e0 = initial.get_potential_energy()
barrier_height = 0
transition_state_idx = 0
energies = []

print("\nFinal Energy Profile:")
print(f"{'Image':<6} {'Energy (eV)':<15} {'Delta E (eV)':<12}")
print("-" * 35)

for i, img in enumerate(images):
    energy = img.get_potential_energy()
    delta_energy = energy - e0
    energies.append(delta_energy)
    if delta_energy > barrier_height:
        barrier_height = delta_energy
        transition_state_idx = i
    print(f"{i:<6} {energy:<15.6f} {delta_energy:<12.6f}")

print("-" * 35)
print(f"\n{'='*60}")
print("Key Results:")
print(f"{'='*60}")
print(f"  Initial State (IS) Energy : {e0:>12.6f} eV")
print(f"  Transition State (TS)     : {images[transition_state_idx].get_potential_energy():>12.6f} eV")
print(f"  Final State (FS) Energy   : {images[-1].get_potential_energy():>12.6f} eV")
print(f"  {'-'*56}")
print(f"  Reaction Energy           : {energies[-1]:>12.6f} eV")
print(f"  Energy Barrier            : {barrier_height:>12.6f} eV")
print(f"  Transition State Image    : {transition_state_idx}")
print(f"{'='*60}\n")

transition_state = images[transition_state_idx]
transition_state.write('images/transition_state.xyz')
print("Transition state saved: images/transition_state.xyz")

for i, img in enumerate(images):
    img.write(f'images/image_{i}.xyz')
print("All images saved: images/image_*.xyz")

plt.figure(figsize=(6, 4))
plt.plot(range(len(images)), energies, 'o-', linewidth=2, alpha=0.8)
plt.xlabel('Image Index')
plt.ylabel('Delta Energy (eV)')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('nep_barrier.png', dpi=300)
print("Energy profile saved: nep_barrier.png")

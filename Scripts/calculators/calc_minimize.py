"""
This script performs structure relaxation using calorine.
It reads a single structure from POSCAR or extxyz file, performs energy minimization,
and outputs the optimization trajectory to minimize.xyz in extxyz format.
The force convergence threshold and maximum steps can be specified via command line arguments.
"""

import sys
from tqdm import tqdm
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from calorine.calculators import CPUNEP
import os

class TrajectoryCallback:
    """Callback class to save trajectory during optimization."""
    def __init__(self, atoms):
        self.atoms = atoms
        self.trajectory = [atoms.copy()]
    
    def __call__(self):
        self.trajectory.append(self.atoms.copy())

def main():
    # Parse command line arguments
    if len(sys.argv) < 3:
        print(" Usage: python calc_minimize.py <structure_file> <nep.txt> [fmax=0.01] [max_steps=1000]")
        print(" Supported structure formats: POSCAR, .xyz")
        print(" Optional arguments: ")
        print("  fmax: Force convergence threshold in eV/Ang (default: 0.01)")
        print("  max_steps: Maximum number of optimization steps (default: 1000)")
        sys.exit(1)

    structure_file = sys.argv[1]  # Structure file (POSCAR or extxyz)
    model_path = sys.argv[2]  # Path to the model for CPUNEP

    # Check if structure file exists
    if not os.path.exists(structure_file):
        print(f" Error: Structure file '{structure_file}' does not exist.")
        sys.exit(1)

    # Check if model file exists
    if not os.path.exists(model_path):
        print(f" Error: Model file '{model_path}' does not exist.")
        sys.exit(1)

    # Set default values and parse optional arguments
    fmax = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01  # Default max force threshold
    max_steps = int(sys.argv[4]) if len(sys.argv) > 4 else 1000  # Default max steps

    print(f" Force convergence threshold: {fmax} eV/Ang")
    print(f" Maximum steps: {max_steps}")
    print(f" Structure file: {structure_file}")
    print(f" Model path: {model_path}")

    # Determine file format and read structure
    file_ext = os.path.splitext(structure_file)[1].lower()
    if file_ext in ['.xyz']:
        atoms = read(structure_file)
    elif 'POSCAR' in structure_file or 'CONTCAR' in structure_file:
        atoms = read(structure_file, format='vasp')
    else:
        print(f" Warning: Unknown file extension '{file_ext}'. Attempting to read as extxyz.")
        atoms = read(structure_file)

    print(f" Read structure with {len(atoms)} atoms")

    # Create calculator and assign to atoms
    calc = CPUNEP(model_path)
    atoms.calc = calc

    # Initialize trajectory callback
    traj_callback = TrajectoryCallback(atoms)

    # Perform structure optimization with trajectory tracking
    opt = BFGS(atoms, logfile='minimize.log')
    opt.attach(traj_callback)  # Attach callback to save trajectory
    opt.run(fmax=fmax, steps=max_steps)  # Convergence criterion based on input

    # Save the complete trajectory to minimize.xyz
    write('minimize.xyz', traj_callback.trajectory)
    print(f" Optimization completed after {len(traj_callback.trajectory)} steps.")
    print(f" Final max force: {np.max(np.abs(atoms.get_forces())):.6f} eV/Ang")
    print(f" Optimization trajectory saved to minimize.xyz")

if __name__ == "__main__":
    main()

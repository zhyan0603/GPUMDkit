"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     calc_minimize_fixed.py
Category:   Calculator Scripts
Purpose:    Structure relaxation with a specified group of atoms fixed.
            Uses GPUMD-style grouping (group_method + group_id).
Usage:      python calc_minimize_fixed.py <structure> <nep.txt> <group_method> <group_id> [fmax] [max_steps]
Arguments:
  structure     Input structure file (POSCAR or .xyz)
  nep.txt       Path to the NEP model file
  group_method  GPUMD group method: 1=by group, 2=by type, 3=group+type
  group_id      1-based group ID to fix (following GPUMD convention)
  fmax          Force convergence threshold (eV/A, default: 0.01)
  max_steps     Max optimization steps (default: 1000)
Output:
  minimized_fixed.xyz  (optimization trajectory in extxyz format)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-07-18
"""

import os
import sys


class TrajectoryCallback:
    """Save optimization trajectory at each step."""
    def __init__(self, atoms):
        self.atoms = atoms
        self.trajectory = [atoms.copy()]

    def __call__(self):
        self.trajectory.append(self.atoms.copy())


def build_fix_mask(atoms, group_method, group_id):
    """Build boolean mask for atoms to fix, following GPUMD grouping convention.

    GPUMD group_method:
      1 = group by 'group' property in extxyz
      2 = group by atom 'type' (0-based type index)
      3 = group by both 'group' and 'type' (combined)

    group_id is 1-based (GPUMD convention):
      group_method=1: group_id maps to group property value + 1
      group_method=2: group_id maps to type value + 1
      group_method=3: group_id = group * n_types + type + 1
    """
    import numpy as np

    if group_method == 1:
        if 'group' not in atoms.arrays:
            raise ValueError(
                "Structure has no 'group' field. "
                "Use 'gpumdkit.sh -addgroup' to add group labels."
            )
        groups = atoms.arrays['group']
        target = group_id - 1
        mask = groups == target
        if not mask.any():
            available = sorted(set(groups.tolist()))
            raise ValueError(
                f"No atoms in group_id={group_id} (group value {target}). "
                f"Available group values: {available}"
            )

    elif group_method == 2:
        if 'type' in atoms.arrays:
            types = atoms.arrays['type']
        else:
            unique_numbers = sorted(set(atoms.numbers.tolist()))
            type_map = {num: i for i, num in enumerate(unique_numbers)}
            types = np.array([type_map[num] for num in atoms.numbers])
        target = group_id - 1
        mask = types == target
        if not mask.any():
            available = sorted(set(types.tolist()))
            raise ValueError(
                f"No atoms with type_id={group_id} (type value {target}). "
                f"Available type values: {available}"
            )

    elif group_method == 3:
        if 'group' not in atoms.arrays:
            raise ValueError(
                "Structure has no 'group' field. "
                "Use 'gpumdkit.sh -addgroup' to add group labels."
            )
        if 'type' in atoms.arrays:
            types = atoms.arrays['type']
        else:
            unique_numbers = sorted(set(atoms.numbers.tolist()))
            type_map = {num: i for i, num in enumerate(unique_numbers)}
            types = np.array([type_map[num] for num in atoms.numbers])
        groups = atoms.arrays['group']
        n_types = len(set(types.tolist()))
        target_group = (group_id - 1) // n_types
        target_type = (group_id - 1) % n_types
        mask = (groups == target_group) & (types == target_type)
        if not mask.any():
            raise ValueError(
                f"No atoms with group_id={group_id} "
                f"(group={target_group}, type={target_type})."
            )

    else:
        raise ValueError(
            f"Invalid group_method={group_method}. Must be 1, 2, or 3."
        )

    return mask


def print_usage():
    print("Usage: python calc_minimize_fixed.py <structure> <nep.txt> "
          "<group_method> <group_id> [fmax] [max_steps]")
    print("")
    print("Arguments:")
    print("  structure     Input structure file (POSCAR or .xyz)")
    print("  nep.txt       Path to the NEP model file")
    print("  group_method  GPUMD group method: 1=by group, 2=by type, 3=group+type")
    print("  group_id      1-based group ID to fix (following GPUMD convention)")
    print("  fmax          Force convergence threshold (eV/A, default: 0.01)")
    print("  max_steps     Max optimization steps (default: 1000)")
    print("")
    print("Examples:")
    print("  # Fix atoms in group 1 (by 'group' property in extxyz)")
    print("  python calc_minimize_fixed.py POSCAR nep.txt 1 1")
    print("")
    print("  # Fix atoms of type 1 (1-based)")
    print("  python calc_minimize_fixed.py POSCAR nep.txt 2 1")
    print("")
    print("  # Fix atoms with combined group+type (group_method=3)")
    print("  python calc_minimize_fixed.py POSCAR nep.txt 3 3")
    print("")
    print("Note: calorine requires the gpumdkit conda environment:")
    print("  conda activate gpumdkit")


def main():
    args = sys.argv[1:]

    if not args or args[0] in ('-h', '--help'):
        print_usage()
        sys.exit(0 if args and args[0] in ('-h', '--help') else 1)

    if len(args) < 4:
        print(f"Error: expected at least 4 arguments, got {len(args)}")
        print_usage()
        sys.exit(1)

    structure = args[0]
    nep_path = args[1]

    try:
        group_method = int(args[2])
        group_id = int(args[3])
    except ValueError:
        print(f"Error: group_method and group_id must be integers, "
              f"got '{args[2]}' and '{args[3]}'")
        sys.exit(1)

    fmax = float(args[4]) if len(args) > 4 else 0.01
    max_steps = int(args[5]) if len(args) > 5 else 1000

    # Check files
    if not os.path.exists(structure):
        sys.exit(f"Error: Structure file '{structure}' does not exist.")
    if not os.path.exists(nep_path):
        sys.exit(f"Error: NEP model file '{nep_path}' does not exist.")

    print(f" Structure: {structure}")
    print(f" NEP model: {nep_path}")
    print(f" group_method: {group_method}, group_id: {group_id}")
    print(f" fmax: {fmax} eV/A, max_steps: {max_steps}")

    # Lazy imports: calorine needs gpumdkit conda env
    try:
        import numpy as np
        from ase.constraints import FixAtoms
        from ase.io import read, write
        from ase.optimize import BFGS
        from calorine.calculators import CPUNEP
    except ImportError as e:
        sys.exit(
            f"Error: {e}\n"
            "calorine requires the gpumdkit conda environment:\n"
            "  conda activate gpumdkit"
        )

    # Read structure
    name = os.path.basename(structure)
    if name == 'POSCAR' or name == 'CONTCAR' or 'POSCAR' in structure:
        atoms = read(structure, format='vasp')
    else:
        atoms = read(structure)
    print(f" Read structure with {len(atoms)} atoms")
    print(f" Composition: {atoms.get_chemical_formula()}")

    # Build fix mask
    try:
        mask = build_fix_mask(atoms, group_method, group_id)
    except ValueError as e:
        sys.exit(f"Error: {e}")

    n_fixed = int(mask.sum())
    pct = 100 * n_fixed / len(atoms) if len(atoms) > 0 else 0
    print(f" Fixing {n_fixed}/{len(atoms)} atoms ({pct:.1f}%)")

    atoms.set_constraint(FixAtoms(mask=mask))

    # Setup calculator and run optimization
    calc = CPUNEP(nep_path)
    atoms.calc = calc

    traj_callback = TrajectoryCallback(atoms)
    opt = BFGS(atoms, logfile='minimize_fixed.log')
    opt.attach(traj_callback)
    opt.run(fmax=fmax, steps=max_steps)

    # Save trajectory
    write('minimized_fixed.xyz', traj_callback.trajectory)

    print(f" Optimization completed after {len(traj_callback.trajectory)} steps.")
    print(f" Final max force: {np.max(np.abs(atoms.get_forces())):.6f} eV/A")
    print(f" Trajectory saved to: minimized_fixed.xyz")


if __name__ == "__main__":
    main()

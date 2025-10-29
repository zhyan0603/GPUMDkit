import sys
import os
import math
import numpy as np
from ase.io import read, write
from ase.build import make_supercell
from ase import Atoms

def find_nearest_supercell(atoms, target):
    """
    Find a supercell configuration (a, b, c) based on target atom number such that:
    - The number of atoms is as close as possible to the target
    - For equal atom number differences, prioritize balanced lattice constants
    """
    n_atoms = len(atoms)
    if target < n_atoms:
        print("Warning: Target atom count is less than original; returning 1x1x1 supercell.")
        return (1, 1, 1)
    
    cell = atoms.cell.lengths()
    ratio = (target / n_atoms) ** (1/3)
    max_mult = max(1, int(math.ceil(ratio * 1.5)))
    
    candidates = [(a, b, c) for a in range(1, max_mult + 1)
                  for b in range(1, max_mult + 1)
                  for c in range(1, max_mult + 1)]
    candidates.sort(key=lambda x: abs(n_atoms * x[0] * x[1] * x[2] - target))
    
    best = None
    best_diff = float("inf")
    best_spread = float("inf")
    
    for a, b, c in candidates:
        atoms_num = n_atoms * a * b * c
        diff = abs(atoms_num - target)
        if diff > best_diff:
            break
        new_lengths = (cell[0] * a, cell[1] * b, cell[2] * c)
        spread = np.std(new_lengths)
        if (diff < best_diff) or (diff == best_diff and spread < best_spread):
            best_diff = diff
            best_spread = spread
            best = (a, b, c)
        if diff == 0:
            break
    
    return best

def reorder_atoms_by_input_species_order(atoms, original_atoms):
    """
    Reorder atoms in the supercell to group by species in the order of the original input file.
    """
    # Get the order of species from the original atoms
    original_symbols = original_atoms.get_chemical_symbols()
    # Get unique species in order of appearance
    seen = set()
    species_order = [s for s in original_symbols if not (s in seen or seen.add(s))]
    
    # Group atoms by species in the supercell, maintaining input order
    new_positions = []
    new_symbols = []
    for species in species_order:
        indices = [i for i, sym in enumerate(atoms.get_chemical_symbols()) if sym == species]
        new_positions.extend(atoms.positions[indices])
        new_symbols.extend([species] * len(indices))
    
    # Create a new Atoms object with reordered atoms
    reordered_atoms = Atoms(
        symbols=new_symbols,
        positions=new_positions,
        cell=atoms.cell,
        pbc=atoms.pbc
    )
    return reordered_atoms

def main():
    if len(sys.argv) < 4:
        print(" Usage 1: python replicate.py input.vasp output.vasp a b c")
        print(" Usage 2: python replicate.py input.vasp output.vasp target_num")
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]
    args = sys.argv[3:]

    # Infer input/output format
    in_ext = os.path.splitext(infile)[1][1:].lower()
    out_ext = os.path.splitext(outfile)[1][1:].lower()
    in_format = "extxyz" if in_ext == "xyz" else in_ext or "extxyz"
    out_format = "extxyz" if out_ext == "xyz" else out_ext or "extxyz"
    
    try:
        atoms = read(infile, format=in_format)
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)
    
    n_atoms = len(atoms)
    
    if len(args) == 1:  # Target atom number
        try:
            target_num = int(args[0])
        except ValueError:
            print(" Error: num must be an integer")
            sys.exit(1)
        a, b, c = find_nearest_supercell(atoms, target_num)
        total_atoms = n_atoms * a * b * c
        print(f"Original atoms: {n_atoms}, Target: {target_num}, Selected: {a}x{b}x{c} (total {total_atoms} atoms)")
    elif len(args) == 3:  # a b c
        try:
            a, b, c = map(int, args)
        except ValueError:
            print(" Error: a, b, c must be integers")
            sys.exit(1)
        if any(x <= 0 for x in (a, b, c)):
            print("Error: a, b, c must be positive integers")
            sys.exit(1)
        total_atoms = n_atoms * a * b * c
        print(f"Supercell: {a}x{b}x{c}, Total atoms: {total_atoms}")
    else:
        print(" Invalid number of arguments, use num or a b c")
        sys.exit(1)
    
    # Construct supercell matrix
    P = np.diag([a, b, c])
    new_atoms = make_supercell(atoms, P)
    
    # Reorder atoms to group by species in input order
    new_atoms = reorder_atoms_by_input_species_order(new_atoms, atoms)
    
    # Write output file, preserving order
    try:
        if out_format in ["vasp", "poscar"]:
            write(outfile, new_atoms, format="vasp", vasp5=True, direct=True, sort=False)
        else:
            write(outfile, new_atoms, format=out_format)
        print(f"Supercell created and saved to {outfile}")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
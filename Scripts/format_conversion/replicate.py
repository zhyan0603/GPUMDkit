import sys
import os
import math
import numpy as np
from ase.io import read, write
from ase.build import make_supercell

def find_nearest_supercell(atoms, target):
    """
    Find a supercell configuration (a, b, c) based on target atom number such that:
    - The number of atoms is as close as possible to the target
    - For equal atom number differences, prioritize balanced lattice constants
    """
    n_atoms = len(atoms)
    cell = atoms.cell.lengths()  # Original lattice constants (a0, b0, c0)

    ratio = (target / n_atoms) ** (1/3)  # Ideal average scaling factor
    max_mult = max(1, int(math.ceil(ratio * 2)))  # Search range

    best = None
    best_diff = float("inf")
    best_spread = float("inf")

    for a in range(1, max_mult+1):
        for b in range(1, max_mult+1):
            for c in range(1, max_mult+1):
                atoms_num = n_atoms * a * b * c
                diff = abs(atoms_num - target)

                # Scaled lattice constants
                new_lengths = (cell[0]*a, cell[1]*b, cell[2]*c)
                spread = np.std(new_lengths)

                # Prioritize atom number closeness, then lattice constant balance
                if (diff < best_diff) or (diff == best_diff and spread < best_spread):
                    best_diff = diff
                    best_spread = spread
                    best = (a, b, c)

    return best

def main():
    if len(sys.argv) < 4:
        print(" Usage 1: python supercell.py inputfile outputfile a b c")
        print(" Usage 2: python supercell.py inputfile outputfile target_num")
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]
    args = sys.argv[3:]

    # Infer input/output format
    in_ext = os.path.splitext(infile)[1][1:].lower()
    out_ext = os.path.splitext(outfile)[1][1:].lower()
    
    # Handle .xyz files as extxyz format
    in_format = "extxyz" if in_ext == "xyz" else in_ext or None
    out_format = "extxyz" if out_ext == "xyz" else out_ext or None

    atoms = read(infile, format=in_format)
    n_atoms = len(atoms)

    if len(args) == 1:  # Target atom number
        try:
            target_num = int(args[0])
        except ValueError:
            print(" Error: num must be an integer")
            sys.exit(1)
        a, b, c = find_nearest_supercell(atoms, target_num)
        print(f" Original atom num: {n_atoms}, Target atom num: {target_num}, Selected supercell: {a}x{b}x{c} (total {n_atoms*a*b*c} atoms)")
    elif len(args) == 3:  # a b c
        try:
            a, b, c = map(int, args)
        except ValueError:
            print(" Error: a, b, c must be integers")
            sys.exit(1)
        print(f" Supercell: {a}x{b}x{c}, Atom count: {n_atoms*a*b*c}")
    else:
        print(" Invalid number of arguments, use num or a b c")
        sys.exit(1)

    # Construct supercell matrix
    P = np.diag([a, b, c])
    new_atoms = make_supercell(atoms, P)

    # Write output file, maintaining consistent input/output format
    if out_format in ["vasp", "poscar"]:
        write(outfile, new_atoms, format="vasp", vasp5=True, direct=True, sort=True)
    else:
        write(outfile, new_atoms, format=out_format)

    print(f" Supercell created, saved to {outfile}")

if __name__ == "__main__":
    main()
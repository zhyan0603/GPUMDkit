"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    scf_batch_pretreatment_cp2k.py
"""

import os
import argparse
from ase.io import read, write

def modify_cp2k_template(template_content, cell):
    """
    Modify the &CELL section in the CP2K template with the given cell vectors.
    cell is a 3x3 matrix from ASE Atoms.get_cell().
    """
    lines = template_content.splitlines()
    in_cell_section = False
    modified_lines = []

    for line in lines:
        stripped = line.strip()

        if stripped.startswith('&CELL'):
            in_cell_section = True
            modified_lines.append(line)
            continue

        if in_cell_section:
            if stripped.startswith('&END CELL'):
                in_cell_section = False
                modified_lines.append(line)
                continue
            elif stripped.startswith('A'):
                a_vec = ' '.join(f'{v:14.8f}' for v in cell[0])
                modified_lines.append(f'      A    {a_vec}')
            elif stripped.startswith('B'):
                b_vec = ' '.join(f'{v:14.8f}' for v in cell[1])
                modified_lines.append(f'      B    {b_vec}')
            elif stripped.startswith('C'):
                c_vec = ' '.join(f'{v:14.8f}' for v in cell[2])
                modified_lines.append(f'      C    {c_vec}')
            else:
                # Keep other lines in &CELL, like PERIODIC
                modified_lines.append(line)
        else:
            modified_lines.append(line)

    return '\n'.join(modified_lines) + '\n'

def main():
    parser = argparse.ArgumentParser(description="Process extxyz file for CP2K single-point calculations.")
    parser.add_argument('extxyz_file', type=str, help="Path to the input extxyz file.")
    parser.add_argument('template_inp', type=str, help="Path to the CP2K template inp file.")
    parser.add_argument('prefix', type=str, default='structure', help="Prefix for output directories.")
    args = parser.parse_args()

    # Read the extxyz file into a list of Atoms objects
    atoms_list = read(args.extxyz_file, format='extxyz', index=':')

    # Read the template inp content
    with open(args.template_inp, 'r') as f:
        template_content = f.read()

    num_structures = len(atoms_list)
    print(f" Found {num_structures} structures in {args.extxyz_file}.")

    # Create directories and files for each structure
    for i, atoms in enumerate(atoms_list, start=1):
        dir_name = f"{args.prefix}_{i}"
        os.makedirs(dir_name, exist_ok=True)

        # Get the cell 
        cell = atoms.get_cell()
        if cell is None or (cell == 0).all():
            raise ValueError(f" Structure {i} has no cell defined. Ensure extxyz includes lattice info.")

        # Modify the template with the current cell
        modified_inp = modify_cp2k_template(template_content, cell)

        # Write the modified inp file
        inp_path = os.path.join(dir_name, 'input.inp')  # You can change the filename if needed
        with open(inp_path, 'w') as f:
            f.write(modified_inp)

        # Write the pos.xyz file (standard XYZ format)
        xyz_path = os.path.join(dir_name, 'pos.xyz')
        write(xyz_path, atoms, format='xyz')

if __name__ == "__main__":
    main()

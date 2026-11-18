"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     scf_batch_pretreatment_cp2k.py
Category:   Workflow Scripts
Purpose:    Batch pretreatment of structures for CP2K SCF calculations.
            Reads structures from extxyz/VASP files, copies a CP2K template,
            and populates the CELL section with the correct lattice vectors.
Usage:      python scf_batch_pretreatment_cp2k.py <extxyz_file> <template_inp> <prefix>
Arguments:
  extxyz_file   Input extxyz file with structures
  template_inp  CP2K input template
  prefix        Prefix for output directories
Output:
  <prefix>_<index>/ directories with input.inp and pos.xyz files
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import os
import argparse

def modify_cp2k_template(template_content, cell):
    """
    Modify the &CELL section in the CP2K template with the given cell vectors.
    cell is a 3x3 matrix from ASE Atoms.get_cell().
    """
    modified_lines = []
    depth = 0
    found_cell = False
    geometry_keys = {'A', 'B', 'C', 'ABC', 'ALPHA_BETA_GAMMA',
                     'CELL_FILE_NAME', 'CELL_FILE_FORMAT'}
    for line in template_content.splitlines():
        key = next(iter(line.upper().split()), '')
        if key == '&CELL' and not depth:
            found_cell = True
            depth = 1
            modified_lines.append(line)
            for name, vector in zip('ABC', cell):
                values = ' '.join(f'{value:14.8f}' for value in vector)
                modified_lines.append(f'      {name}    {values}')
            continue
        if depth:
            if key == '&END':
                depth -= 1
            elif key.startswith('&'):
                depth += 1
            elif depth == 1 and key in geometry_keys:
                continue
        modified_lines.append(line)

    if not found_cell or depth:
        raise ValueError("The CP2K template must contain a complete &CELL section.")
    return '\n'.join(modified_lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description="Process extxyz file for CP2K single-point calculations.")
    parser.add_argument('extxyz_file', type=str, help="Path to the input extxyz file.")
    parser.add_argument('template_inp', type=str, help="Path to the CP2K template inp file.")
    parser.add_argument('prefix', type=str, default='structure', help="Prefix for output directories.")
    args = parser.parse_args()

    for path in (args.extxyz_file, args.template_inp):
        if not os.path.isfile(path):
            parser.exit(1, f" Error: file '{path}' does not exist.\n")

    from ase.io import read, write

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
        cell = atoms.get_cell()
        if cell.rank != 3 or cell.volume == 0:
            parser.exit(1, f" Error: structure {i} has no valid 3D cell.\n")
        try:
            modified_inp = modify_cp2k_template(template_content, cell)
        except ValueError as exc:
            parser.exit(1, f" Error: {exc}\n")
        os.makedirs(dir_name, exist_ok=True)

        # Write the modified inp file
        inp_path = os.path.join(dir_name, 'input.inp')  # You can change the filename if needed
        with open(inp_path, 'w') as f:
            f.write(modified_inp)

        # Write the pos.xyz file (standard XYZ format)
        xyz_path = os.path.join(dir_name, 'pos.xyz')
        write(xyz_path, atoms, format='xyz')

if __name__ == "__main__":
    main()

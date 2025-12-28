"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Perturb atomic positions in structures

Usage:
    python perturb_structure.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

import argparse
import sys
import dpdata

def parse_args():
    parser = argparse.ArgumentParser(
             description="Usage: python script.py <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>")
    parser.add_argument('input_file', help='The path to POSCAR/CONTCAR file')
    parser.add_argument('pert_num', type=int, default=20, help='The perturbation number')
    parser.add_argument('cell_pert_fraction', type=float, default=0.03, help='The fraction of cell perturbation')
    parser.add_argument('atom_pert_distance', type=float, default=0.2, help='The distance of atom perturbation')
    parser.add_argument('atom_pert_style', type=str, default='uniform', choices=['normal', 'uniform', 'const'], help='The style for atom perturbation')
    return parser.parse_args()

def main():
    try:
        args = parse_args()
    except SystemExit:
        print(f"Default values: pert_num=20, cell_pert_fraction=0.03, atom_pert_distance=0.2, atom_pert_style=uniform")
        print("atom_pert_style options: 'normal', 'uniform', 'const'")
        print("dpdata documentation: https://docs.deepmodeling.com/projects/dpdata/en/master/index.html \n")
        sys.exit(1)

    # Read the POSCAR file and perform perturbation
    system = dpdata.System(args.input_file, fmt='vasp/poscar')
    perturbed_systems = system.perturb(pert_num=args.pert_num,
                                       cell_pert_fraction=args.cell_pert_fraction,
                                       atom_pert_distance=args.atom_pert_distance,
                                       atom_pert_style=args.atom_pert_style,)

    # Save the perturbed structures
    pert_num = args.pert_num
    width = len(str(pert_num))
    for i in range(pert_num):
        output_file = f"POSCAR_{str(i + 1).zfill(width)}.vasp"
        perturbed_systems.sub_system(i).to('vasp/poscar', output_file)

if __name__ == "__main__":
     main()

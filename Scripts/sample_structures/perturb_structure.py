"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     perturb_structure.py
Category:   Sample Structure Scripts
Purpose:    Generate perturbed structures from a POSCAR/CONTCAR file using
            dpdata. Supports cell and atom perturbations with different
            perturbation styles (normal, uniform, const).
Usage:      gpumdkit.sh
            choose 204) Perturb structure
            python perturb_structure.py <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>
Arguments:
  input.vasp           Input POSCAR/CONTCAR file
  pert_num             Number of perturbed structures to generate
  cell_pert_fraction   Fraction of cell perturbation
  atom_pert_distance   Distance of atom perturbation (Angstrom)
  atom_pert_style      Style: normal, uniform, or const
Output:
  POSCAR_*.vasp  (perturbed structures in VASP POSCAR format)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import argparse
import sys


def print_dependency_notice():
    print(" This function requires the dpdata package.")
    print(" If you use this function, we recommend citing dpdata according to its official documentation.")


def print_usage():
    print(" Usage: gpumdkit.sh")
    print("        choose 204) Perturb structure")
    print("    or: python perturb_structure.py <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>")
    print("")
    print(" Arguments:")
    print("   input.vasp           Input POSCAR/CONTCAR file")
    print("   pert_num             Number of perturbed structures to generate")
    print("   cell_pert_fraction   Fraction of cell perturbation")
    print("   atom_pert_distance   Distance of atom perturbation (Angstrom)")
    print("   atom_pert_style      Style: normal, uniform, or const")
    print("")
    print(" Output:")
    print("   POSCAR_*.vasp")
    print("")
    print(" Example: in interactive mode, enter: POSCAR 20 0.03 0.2 uniform")
    print("          python perturb_structure.py POSCAR 20 0.03 0.2 uniform")
    print("")


def parse_args():
    parser = argparse.ArgumentParser(
             formatter_class=argparse.RawDescriptionHelpFormatter,
             description="""Generate perturbed structures from a POSCAR/CONTCAR file.

Usage:
  gpumdkit.sh
  choose 204) Perturb structure
  python perturb_structure.py <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>

Example:
  in interactive mode, enter: POSCAR 20 0.03 0.2 uniform
  python perturb_structure.py POSCAR 20 0.03 0.2 uniform
""")
    parser.add_argument('input_file', help='The path to POSCAR/CONTCAR file')
    parser.add_argument('pert_num', type=int, default=20, help='The perturbation number')
    parser.add_argument('cell_pert_fraction', type=float, default=0.03, help='The fraction of cell perturbation')
    parser.add_argument('atom_pert_distance', type=float, default=0.2, help='The distance of atom perturbation')
    parser.add_argument('atom_pert_style', type=str, default='uniform', choices=['normal', 'uniform', 'const'], help='The style for atom perturbation')
    return parser.parse_args()

def main():
    args_in = sys.argv[1:]
    if len(args_in) != 5 and not (args_in and args_in[0] in ("-h", "--help")):
        print_usage()
        sys.exit(1)

    try:
        args = parse_args()
    except SystemExit as exc:
        if exc.code == 0:
            sys.exit(0)
        print(f"Default values: pert_num=20, cell_pert_fraction=0.03, atom_pert_distance=0.2, atom_pert_style=uniform")
        print("atom_pert_style options: 'normal', 'uniform', 'const'")
        print("dpdata documentation: https://docs.deepmodeling.com/projects/dpdata/en/master/index.html \n")
        sys.exit(exc.code)

    print_dependency_notice()

    try:
        import dpdata
    except ImportError:
        print(" Error: dpdata is not installed or cannot be imported.")
        print(" Please install dpdata before using this function.")
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

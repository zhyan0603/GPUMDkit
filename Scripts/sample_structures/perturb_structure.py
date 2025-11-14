import argparse
import sys
import dpdata

def parse_args():
    parser = argparse.ArgumentParser(
             description="Usage: python script.py <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>")
    parser.add_argument('input_file', help='输入 VASP 文件路径')
    parser.add_argument('pert_num', type=int, default=20, help='微扰结构数量')
    parser.add_argument('cell_pert_fraction', type=float, default=0.03, help='晶胞微扰比例')
    parser.add_argument('atom_pert_distance', type=float, default=0.2, help='原子微扰距离')
    parser.add_argument('atom_pert_style', default='uniform', choices=['normal', 'uniform', 'const'], help='原子微扰样式')
    return parser.parse_args()

def main():
    try:
        args = parse_args()
    except SystemExit:
        print(f"Default values: pert_num={DEFAULT_PERT_NUM}, cell_pert_fraction={DEFAULT_CELL_PERT_FRACTION}, atom_pert_distance={DEFAULT_ATOM_PERT_DISTANCE}, atom_pert_style={DEFAULT_ATOM_PERT_STYLE}")
        print("atom_pert_style options: 'normal', 'uniform', 'const'")
        print("dpdata documentation: https://docs.deepmodeling.com/projects/dpdata/en/master/index.html")
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

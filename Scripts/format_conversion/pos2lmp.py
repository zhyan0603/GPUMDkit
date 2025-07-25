import sys
from ovito.io import import_file, export_file

def convert_poscar_to_lammps(poscar_path, lammps_data_path):
    pipeline = import_file(poscar_path)
    
    export_file(pipeline, lammps_data_path, "lammps/data", multiple_frames=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pos2lmp.py <poscar_file> <lammps_data_file>")
        sys.exit(1)
    
    poscar_file = sys.argv[1]
    lammps_data_file = sys.argv[2]
    
    convert_poscar_to_lammps(poscar_file, lammps_data_file)
    print(f" Converted {poscar_file} to {lammps_data_file} successfully.")

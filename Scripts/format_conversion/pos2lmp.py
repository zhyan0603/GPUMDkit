"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     pos2lmp.py
Category:   Format Conversion Scripts
Purpose:    Convert VASP POSCAR file to LAMMPS data format using OVITO.
Usage:      gpumdkit.sh -pos2lmp <poscar_file> <lammps_data_file>
            python pos2lmp.py <poscar_file> <lammps_data_file>
Arguments:
  poscar_file      Input VASP POSCAR file
  lammps_data_file Output LAMMPS data file
Output:
  <lammps_data_file>  (LAMMPS data format structure)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys
from ovito.io import import_file, export_file

def convert_poscar_to_lammps(poscar_path, lammps_data_path):
    pipeline = import_file(poscar_path)
    
    export_file(pipeline, lammps_data_path, "lammps/data", multiple_frames=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(" Usage: gpumdkit.sh -pos2lmp <poscar_file> <lammps_data_file>")
        print("    or: python pos2lmp.py <poscar_file> <lammps_data_file>")
        sys.exit(1)
    
    poscar_file = sys.argv[1]
    lammps_data_file = sys.argv[2]
    
    convert_poscar_to_lammps(poscar_file, lammps_data_file)
    print(f" Converted {poscar_file} to {lammps_data_file} successfully.")

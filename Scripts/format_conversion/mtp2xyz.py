"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     mtp2xyz.py
Category:   Format Conversion Scripts
Purpose:    Convert MTP (Machine Learning Interatomic Potential) input file
            format to extended XYZ using a custom cfg parser and ASE.
Usage:      Interactive: gpumdkit.sh -> 1 -> 102
            python mtp2xyz.py <train.cfg> <Symbol1> <Symbol2> ...
Arguments:
  train.cfg   MTP training data file
  SymbolX     Chemical element symbols in order
Output:
  XYZ/mtp2xyz.xyz   Converted structures in extxyz format
Author:     Ke XU (kickhsu@gmail.com)
Last-modified: 2026-09-14
=============================================================================
"""

import os
import sys
import numpy as np
from collections import defaultdict


def load_cfg(filename, type_to_symbol):
    frames = []
    with open(filename) as f:
        line = 'chongchongchong!'
        while line:
            line = f.readline()
            if 'BEGIN_CFG' in line:
                cell = np.zeros((3, 3))
            if 'Size' in line:
                line = f.readline()
                natoms = int(line.split()[0])
                positions = np.zeros((natoms, 3))
                forces = np.zeros((natoms, 3))
                energies = np.zeros(natoms)
                symbols = ['X'] * natoms
            if 'Supercell' in line: 
                for i in range(3):
                    line = f.readline()
                    for j in range(3):
                        cell[i, j] = float(line.split()[j])
            if 'AtomData' in line:
                d = defaultdict(int)
                for (i, x) in enumerate(line.split()[1:]):
                    d[x] = i

                for _ in range(natoms):
                    line = f.readline()
                    fields = line.split()
                    i = int(fields[d['id']]) - 1
                    symbols[i] = type_to_symbol[int(fields[d['type']])]
                    positions[i] = [float(fields[d[attr]]) for attr in ['cartes_x', 'cartes_y' ,'cartes_z']]
                    forces[i] = [float(fields[d[attr]]) for attr in ['fx', 'fy' ,'fz']]
                    energies[i] = float(fields[d['site_en']])

                atoms = Atoms(symbols=symbols, cell=cell, positions=positions, pbc=True)
                if d['fx'] != 0:
                    atoms.info['forces'] = forces
                if d['site_en'] != 0:
                    atoms.info['energies'] = energies

            if 'Energy' in line and 'Weight' not in line:
                line = f.readline()
                atoms.info['energy'] = float(line.split()[0])
            if 'PlusStress' in line:
                line = f.readline()
                plus_stress = np.array(list(map(float, line.split())))
                atoms.info['virial'] = plus_stress
            if 'END_CFG' in line:
                frames.append(atoms)
            if 'EnergyWeight' in line:
                line = f.readline()
                atoms.info['energy_weight'] = float(line.split()[0])
            if 'identification' in line:
                atoms.info['identification'] = int(line.split()[2])

    return frames


def dump_xyz(frames):

    Out_string = ""
    n_frames = len(frames)
    for atoms in frames:
        Out_string += str(len(atoms)) + '\n'
        Out_string += "energy=" + str(atoms.info['energy']) + " "
        Out_string += "config_type=mtp2nep "
        Out_string += "pbc=\"T T T\" "
        Out_string += "virial=\"" + " ".join(list(map(str, atoms.info['virial']))) + "\" "
        Out_string += "Lattice=\"" + " ".join(list(map(str, atoms.get_cell().reshape(-1)))) + "\" "
        Out_string += "Properties=species:S:1:pos:R:3:force:R:3\n"

        s = atoms.get_chemical_symbols()
        p = atoms.get_positions()
        forces = atoms.info['forces']
        for i in range(len(atoms)):
            Out_string += '{:2} {:>15.8e} {:>15.8e} {:>15.8e} {:>15.8e} {:>15.8e} {:>15.8e}\n'.format(s[i], *p[i], *forces[i])

    os.makedirs("XYZ", exist_ok=True)
    fo = open(os.path.join("XYZ", 'mtp2xyz.xyz'), 'w')
    fo.write(Out_string)
    fo.close()


if __name__ == "__main__":
    # Check arguments
    args = sys.argv[1:]
    if len(args) < 2 or args[0] in ("-h", "--help"):
        print(" Usage: gpumdkit.sh -> 1 -> 102 (interactive)")
        print("    or: python mtp2xyz.py <train.cfg> <Symbol1> <Symbol2> ...")
        print("")
        print(" Arguments:")
        print("   train.cfg   MTP training data file")
        print("   SymbolX     Chemical element symbols in order")
        print("")
        print(" Example: python mtp2xyz.py train.cfg Pd Ag")
        print("")
        sys.exit(0 if args and args[0] in ("-h", "--help") else 1)
    if not os.path.isfile(args[0]):
        print(f" Error: file '{args[0]}' does not exist.")
        sys.exit(1)

    from ase.atoms import Atoms

    type_to_symbol = {i: s for i, s in enumerate(args[1:])}
    frames = load_cfg(args[0], type_to_symbol)
    dump_xyz(frames)

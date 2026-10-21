"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     probability_density_analysis.py
Category:   Analyzer Scripts
Purpose:    Compute the 3D probability density of mobile ions from an AIMD
            trajectory using pymatgen, and save the result as a CHGCAR file.
Usage:      gpumdkit.sh -pda <ref_struct> <trajectory_file> <species> <interval>
            python probability_density_analysis.py <ref_struct> <trajectory_file> <species> <interval>
Arguments:
  ref_struct       Reference structure file (e.g., POSCAR)
  trajectory_file  AIMD trajectory in extxyz format
  species          Mobile species symbol (e.g., Li)
  interval         Grid interval for probability density (Angstrom)
Output:
  probability_density_<interval>.vasp (CHGCAR-like file)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import sys

# Check if required arguments are provided before heavy imports
args = sys.argv[1:]
if len(args) < 4 or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -pda <ref_struct> <trajectory_file> <species> <interval>")
    print("    or: python probability_density_analysis.py <ref_struct> <trajectory_file> <species> <interval>")
    print("")
    print(" Arguments:")
    print("   ref_struct       Reference structure file (e.g., POSCAR)")
    print("   trajectory_file  AIMD trajectory in extxyz format")
    print("   species          Mobile species symbol (e.g., Li)")
    print("   interval         Grid interval for probability density (Angstrom)")
    print("")
    print(" Output:")
    print("   probability_density_<interval>.vasp (CHGCAR-like file)")
    print("")
    print(" Example: gpumdkit.sh -pda POSCAR dump.xyz Li 0.1")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

import numpy as np
from pymatgen.analysis.diffusion.aimd.pathway import ProbabilityDensityAnalysis
from pymatgen.core import Structure
from ase.io import read

# Parse command-line arguments
structure_file = sys.argv[1]
trajectory_file = sys.argv[2]
species = (sys.argv[3],)  # Convert species to tuple
interval = float(sys.argv[4])

# Load structure from current directory
structure = Structure.from_file(structure_file)

# Read trajectory from dump.xyz (non-ULM format)
trajectory = read(trajectory_file, index=":")
trajectory_npy = [atom.get_scaled_positions() for atom in trajectory]

# Perform probability density analysis
pda = ProbabilityDensityAnalysis(
    structure=structure,
    trajectories=trajectory_npy,
    interval=interval,
    species=species
)

# Save probability density to CHGCAR file
pda.to_chgcar(filename=f"probability_density_{interval}.vasp")

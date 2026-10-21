"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     rdf_calculator_ovito.py
Category:   Calculator Scripts
Purpose:    Calculate the radial distribution function (RDF) using the
            OVITO library.
Usage:      python rdf_calculator_ovito.py <exyzfile> <cutoff> <bins>
Arguments:
  exyzfile  Path to the input extxyz file
  cutoff    Cutoff distance for coordination analysis
  bins      Number of bins for the RDF histogram
Output:
  rdf.txt   (RDF values exported to text file)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
=============================================================================
"""

import os
import sys

args = sys.argv[1:]
if len(args) != 3 or args[0] in ("-h", "--help"):
    print(" Usage: python rdf_calculator_ovito.py <exyzfile> <cutoff> <bins>")
    print("")
    print(" Arguments:")
    print("   exyzfile   Input extended XYZ trajectory file")
    print("   cutoff     Cutoff distance in Angstrom for the RDF")
    print("   bins       Number of RDF histogram bins")
    print("")
    print(" Output:")
    print("   rdf.txt    Averaged (partial) RDF table from ovito")
    print("")
    print(" Example: python rdf_calculator_ovito.py dump.xyz 6.0 100")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

if not os.path.isfile(args[0]):
    print(f" Error: file '{args[0]}' does not exist.")
    sys.exit(1)

print(" This function requires the ovito package.")
print(" If you use this function, please cite:")
print("   A. Stukowski, Model. Simul. Mater. Sci. Eng. 18, 015012 (2010)")

from ovito.io import import_file, export_file
from ovito.modifiers import CoordinationAnalysisModifier, TimeAveragingModifier

exyzfile = args[0]
cutoff = float(args[1])
bins = int(args[2])

pipeline = import_file(exyzfile)
modifier = CoordinationAnalysisModifier(cutoff=cutoff, number_of_bins=bins, partial=True)
pipeline.modifiers.append(modifier)
pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))
export_file(pipeline,"rdf.txt","txt/table",key="coordination-rdf[average]")

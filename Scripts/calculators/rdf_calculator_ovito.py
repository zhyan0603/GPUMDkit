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

import sys
from ovito.io import import_file, export_file
from ovito.modifiers import CoordinationAnalysisModifier,TimeAveragingModifier

# Check if the number of variables is three
if len(sys.argv) != 4:
    print("Error: Invalid number of arguments.")
    print("Usage: python rdf_calculator_ovito.py exyzfile cutoff bins")
    sys.exit(1)

exyzfile = sys.argv[1]
cutoff = float(sys.argv[2])
bins = int(sys.argv[3])

pipeline = import_file(exyzfile)
modifier = CoordinationAnalysisModifier(cutoff=cutoff, number_of_bins=bins, partial=True)
pipeline.modifiers.append(modifier)
pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))
export_file(pipeline,"rdf.txt","txt/table",key="coordination-rdf[average]")

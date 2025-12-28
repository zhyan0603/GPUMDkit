"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Calculate radial distribution function using OVITO

Usage:
    python rdf_calculator_ovito.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
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
cutoff = sys.argv[2]
bins = sys.argv[3]

pipeline = import_file(exyzfile)
modifier = CoordinationAnalysisModifier(cutoff=cutoff,number_of_bins=bins,partial=True)
pipeline.modifiers.append(modifier)
pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))
export_file(pipeline,"rdf.txt","txt/table",key="coordination-rdf[average]")

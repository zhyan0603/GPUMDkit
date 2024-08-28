"""
This script is used to add a weight value to the 'Weight' attribute of each structure in an input file and save the modified structures to an output file.

Usage:
    python add_weight.py <input_file> <output_file> <new_weight>

Arguments:
    input_file (str): Path to the input file containing the structures.
    output_file (str): Path to the output file where the modified structures will be saved.
    new_weight (str): The new weight value to be assigned to the 'Weight' attribute of each structure.

Example:
    python add_weight.py input.xyz output.xyz 1.5
"""

import os
import sys
from ase.io import read, write

input_file = sys.argv[1]
output_file = sys.argv[2]
new_weight = sys.argv[3]  

structures = read(input_file, index=':')  

for structure in structures:
    structure.info['Weight'] = new_weight

write(output_file, structures)


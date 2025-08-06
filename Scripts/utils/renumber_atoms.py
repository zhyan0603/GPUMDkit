"""
Renumber atom IDs in a LAMMPS dump file.

Usage:
    python renumber_atoms.py <input_file> <output_file>

Author:
    Dian HUANG <huangdian@stu.xjtu.edu.cn>
"""

import sys
from tqdm import tqdm

# Check command-line arguments
if len(sys.argv) != 3:
    print("Usage: python renumber_atoms.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Count total lines for progress bar
print("Counting lines in input file...")
with open(input_file, 'r') as f:
    total_lines = sum(1 for _ in f)

# Process the file
print("Processing file...")
with open(input_file, "r") as in_f, open(output_file, "w") as out_f:
    pbar = tqdm(total=total_lines, desc="Processing", unit="lines")
    in_atoms_section = False

    for line in in_f:
        if "ITEM: ATOMS" in line:
            in_atoms_section = True
            new_id = 1  # Reset ID for each new frame
            out_f.write(line)
        elif in_atoms_section and line.strip() and not line.startswith("ITEM:"):
            parts = line.split()
            parts[0] = str(new_id)
            out_f.write(" ".join(parts) + "\n")
            new_id += 1
        else:
            in_atoms_section = False
            out_f.write(line)
        pbar.update(1)

    pbar.close()

print(f"Processing complete! Output saved to {output_file}")
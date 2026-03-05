import sys
from ase.io import read
from ase.data import chemical_symbols, atomic_numbers

# Load all frames from the provided extxyz file
try:
    atoms_list = read(sys.argv[1], index=':')
except IndexError:
    print("Usage: python analyzer.py <filename.extxyz>")
    sys.exit(1)

# Identify all unique elements present across all frames
element_set = set()
for atoms in atoms_list:
    element_set.update(atoms.get_chemical_symbols())

# Sort elements according to their atomic numbers (Periodic Table order)
sorted_elements = sorted(element_set, key=lambda x: atomic_numbers[x])

# Display results
print(f'Total unique elements: {len(sorted_elements)}')
print('Element list (Periodic Table order):')
print(' '.join(sorted_elements))
"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Analyze atomic composition of structures

Usage:
    python analyze_composition.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

import re
import sys
from ase.io import read, write
from ase.data import atomic_numbers
from collections import Counter

def get_composition(atoms):
    """Get the chemical composition of a single structure, sorted by atomic number"""
    symbols = atoms.get_chemical_symbols()
    # Count the occurrence of each element
    composition = Counter(symbols)
    # Sort elements by atomic number (periodic table order)
    elements = sorted(composition.keys(), key=lambda x: atomic_numbers.get(x, float('inf')))
    # Format the composition string
    comp_str = ''
    for elem in elements:
        count = composition[elem]
        comp_str += f"{elem}{count}"
    # Calculate total number of atoms
    total_atoms = sum(composition.values())
    return comp_str, composition, total_atoms

def parse_export_input(input_str, max_index):
    """Parse user input for export indices (e.g., '1,2', '2-3', 'all')"""
    if input_str.lower() == 'all':
        return list(range(1, max_index + 1))
    indices = set()
    parts = input_str.replace(' ', '').split(',')
    for part in parts:
        if '-' in part:
            try:
                start, end = map(int, part.split('-'))
                if 1 <= start <= end <= max_index:
                    indices.update(range(start, end + 1))
            except ValueError:
                print(f" Invalid range format: {part}, skipping.")
        else:
            try:
                idx = int(part)
                if 1 <= idx <= max_index:
                    indices.add(idx)
            except ValueError:
                print(f" Invalid index: {part}, skipping.")
    return sorted(indices)

def analyze_xyz_compositions(filename):
    """Analyze compositions in the xyz file and allow export"""
    # Read all structures from the extxyz file
    structures = read(filename, index=':', format='extxyz')
    
    # Group structures by composition
    composition_data = {}
    for struct in structures:
        comp_str, comp_dict, total_atoms = get_composition(struct)
        if comp_str not in composition_data:
            composition_data[comp_str] = {
                'comp_dict': comp_dict,
                'total_atoms': total_atoms,
                'structures': []
            }
        composition_data[comp_str]['structures'].append(struct)
    
    # Assign indices to unique compositions, sorted by element atomic numbers
    unique_compositions = sorted(composition_data.keys(), key=lambda x: tuple(atomic_numbers.get(elem[0], float('inf')) for elem in re.findall(r'([A-Z][a-z]?)(\d+)', x)))
    
    # Print results with aligned columns
    print(f" {'Index':<8} {'Compositions':<22} {'N atoms':<12} {'Count':<10}")
    print(" "+"-" * 51)
    comp_to_index = {}
    for idx, comp in enumerate(unique_compositions, 1):
        data = composition_data[comp]
        comp_to_index[comp] = idx
        print(f" {idx:<8} {comp:<22} {data['total_atoms']:<12} {len(data['structures']):<10}")
    print(" "+"-" * 51)
    # Prompt for export
    print(" Enter index to export (e.g., '1,2', '2-3', 'all'), or press Enter to skip:")
    user_input = input(" ").strip()
    if user_input:
        indices = parse_export_input(user_input, len(unique_compositions))
        if not indices:
            print(" No valid indices provided. Exiting without exporting.")
            return
        
        # Export selected compositions
        for comp, idx in comp_to_index.items():
            if idx in indices:
                output_file = f"{comp}.xyz"
                write(output_file, composition_data[comp]['structures'], format='extxyz')
                print(f" Exported {comp} to {output_file} ({len(composition_data[comp]['structures'])} structures)")

# Run the analysis
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(" Usage: python analyze_composition.py <input.xyz>")
        sys.exit(1)
    filename = sys.argv[1]  # your extxyz file
    analyze_xyz_compositions(filename)

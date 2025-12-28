"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Check charge balance in structures

Usage:
    python charge_balance_check.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""

import sys
from ase.io import read, write
from ase.data import atomic_numbers
from collections import Counter
from pymatgen.core.composition import Composition
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial

def get_composition(atoms):
    """Get the reduced chemical composition of a single structure."""
    symbols = atoms.get_chemical_symbols()
    # Count the occurrence of each element
    composition = Counter(symbols)
    # Calculate GCD to reduce the formula
    def gcd(a, b):
        while b:
            a, b = b, a % b
        return a
    counts = list(composition.values())
    if not counts:
        return "", composition, 0
    divisor = counts[0]
    for count in counts[1:]:
        divisor = gcd(divisor, count)
    # Create reduced formula sorted by atomic number
    reduced_comp = {elem: count // divisor for elem, count in composition.items()}
    elements = sorted(reduced_comp.keys(), key=lambda x: atomic_numbers.get(x, float('inf')))
    comp_str = ''.join(f"{elem}{reduced_comp[elem] if reduced_comp[elem] > 1 else ''}" for elem in elements)
    total_atoms = sum(composition.values())
    return comp_str, composition, total_atoms

def check_oxidation_state(comp_str):
    """Check if the composition can have balanced oxidation states."""
    try:
        comp = Composition(comp_str)
        # Check if it's a single element (pure substance)
        if len(comp.elements) == 1:
            return comp_str, True  # Single elements are balanced (oxidation state = 0)
        oxi_states = comp.oxi_state_guesses()
        return comp_str, bool(oxi_states)
    except Exception as e:
        return comp_str, False, str(e)

def process_structure(idx_atoms):
    idx, atoms = idx_atoms
    comp_str, _, _ = get_composition(atoms)
    atoms.info["reduced_formula"] = comp_str
    # Store original_index temporarily for indexing, will not be written to output
    atoms.info["temp_index"] = idx
    return atoms

def main():
    # Check if input file is provided
    if len(sys.argv) != 2:
        print("Usage: python process_structures.py <input.extxyz>", file=sys.stderr)
        sys.exit(1)
    
    input_file = sys.argv[1]
    balanced_file = "balanced.xyz"
    unbalanced_file = "unbalanced.xyz"
    index_file = "indices.txt"
    
    # Read all structures from the extxyz file
    atoms_list = read(input_file, index=":", format="extxyz")
    
    # Compute compositions in parallel
    idx_atoms_pairs = list(enumerate(atoms_list))
    with Pool() as pool:
        structures = []
        for structure in tqdm(pool.imap_unordered(process_structure, idx_atoms_pairs),
                             total=len(idx_atoms_pairs),
                             desc="Computing compositions"):
            structures.append(structure)
    
    # Get unique compositions
    unique_compositions = set(atoms.info["reduced_formula"] for atoms in structures if atoms.info["reduced_formula"])
    
    # Check oxidation states for unique compositions
    with Pool() as pool:
        oxi_results = dict(tqdm(pool.imap_unordered(check_oxidation_state, unique_compositions),
                                total=len(unique_compositions),
                                desc="Checking oxidation states"))
    
    # Assign oxidation state results to structures
    balanced_structures = []
    unbalanced_structures = []
    balanced_indices = []
    unbalanced_indices = []
    
    for atoms in structures:
        idx = atoms.info["temp_index"]
        comp_str = atoms.info["reduced_formula"]
        # Remove temp_index before writing to output
        del atoms.info["temp_index"]
        if not comp_str:
            atoms.info["oxidation_state"] = "unbalanced"
            atoms.info["error"] = "Empty composition"
            unbalanced_structures.append(atoms)
            unbalanced_indices.append(idx)
            continue
        
        result = oxi_results[comp_str]
        if isinstance(result, tuple):
            # Error case
            comp_str, is_balanced, error = result
            atoms.info["oxidation_state"] = "unbalanced"
            atoms.info["error"] = error
            unbalanced_structures.append(atoms)
            unbalanced_indices.append(idx)
        elif result:
            atoms.info["oxidation_state"] = "balanced"
            balanced_structures.append(atoms)
            balanced_indices.append(idx)
        else:
            atoms.info["oxidation_state"] = "unbalanced"
            unbalanced_structures.append(atoms)
            unbalanced_indices.append(idx)
    
    # Write to separate output files
    if balanced_structures:
        write(balanced_file, balanced_structures, format="extxyz")
    if unbalanced_structures:
        write(unbalanced_file, unbalanced_structures, format="extxyz")
    
    # Write indices to a text file
    with open(index_file, "w") as f:
        f.write("Balanced structures (indices): {}\n".format(sorted(balanced_indices)))
        f.write("Unbalanced structures (indices): {}\n".format(sorted(unbalanced_indices)))
        f.write("Balanced structures written to {}\n".format(balanced_file))
        f.write("Unbalanced structures written to {}\n".format(unbalanced_file))

if __name__ == "__main__":
    main()

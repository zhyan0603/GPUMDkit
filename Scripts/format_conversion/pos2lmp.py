import sys
from ase.io import read

def parse_element_order_from_poscar(poscar_file):
    """
    Parse the element order from the POSCAR file.
    The element order is defined on the 6th line of the POSCAR file.

    Parameters:
        poscar_file (str): Path to the POSCAR file.

    Returns:
        list: List of element symbols in the order they appear in the POSCAR.
    """
    with open(poscar_file, 'r') as f:
        lines = f.readlines()
    # Ensure at least 6 lines in POSCAR
    if len(lines) < 6:
        raise ValueError("Invalid POSCAR format: fewer than 6 lines.")
    element_order = lines[5].split()
    return element_order

def convert_poscar_to_lammps_data(poscar_file, lammps_data_file, element_order=None):
    """
    Convert a POSCAR file to LAMMPS data file format, with optional custom element ordering.

    Parameters:
        poscar_file (str): Path to the input POSCAR file.
        lammps_data_file (str): Path to the output LAMMPS data file.
        element_order (list, optional): List of element symbols in the desired order.
                                        If None, elements are ordered as in the POSCAR file.
    """
    # Read the POSCAR file
    atoms = read(poscar_file, format='vasp')

    # Get element symbols and masses
    symbols = atoms.get_chemical_symbols()
    masses = atoms.get_masses()

    # Parse element order from POSCAR if not provided
    if element_order is None:
        element_order = parse_element_order_from_poscar(poscar_file)

    # Ensure element order includes only elements present in the structure
    unique_symbols = list(set(symbols))
    element_order = [symbol for symbol in element_order if symbol in unique_symbols]

    # Create a dictionary to map element symbols to their type numbers
    symbol_to_type = {symbol: i + 1 for i, symbol in enumerate(element_order)}

    # Create LAMMPS data file
    with open(lammps_data_file, 'w') as f:
        f.write("LAMMPS data file\n\n")

        # Write the number of atoms and atom types
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{len(element_order)} atom types\n\n")

        # Write the box dimensions
        cell = atoms.get_cell()
        f.write(f"0.0 {cell[0][0]:.6f} xlo xhi\n")
        f.write(f"0.0 {cell[1][1]:.6f} ylo yhi\n")
        f.write(f"0.0 {cell[2][2]:.6f} zlo zhi\n")
        f.write(f"{cell[1][0]:.6f} {cell[2][0]:.6f} {cell[2][1]:.6f} xy xz yz\n\n")

        # Write the masses
        f.write("Masses\n\n")
        for i, symbol in enumerate(element_order):
            # Use the mass corresponding to the first occurrence of the symbol
            mass = masses[symbols.index(symbol)]
            f.write(f"{i+1} {mass:.6f} # {symbol}\n")
        f.write("\n")

        # Write the atomic coordinates
        f.write("Atoms\n\n")
        for i, atom in enumerate(atoms):
            atom_id = i + 1
            atom_type = symbol_to_type[atom.symbol]
            x, y, z = atom.position
            f.write(f"{atom_id} {atom_type} {x:.6f} {y:.6f} {z:.6f}\n")

    print(f"Conversion complete! {poscar_file} has been converted to {lammps_data_file}")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python pos2lmp.py <poscar_file> <lammps_data_file> [element_order]")
        sys.exit(1)

    poscar_file = sys.argv[1]
    lammps_data_file = sys.argv[2]

    # Parse element order from command line or default to POSCAR order
    element_order = sys.argv[3:] if len(sys.argv) > 3 else None

    convert_poscar_to_lammps_data(poscar_file, lammps_data_file, element_order)
"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     lmp2exyz.py
Category:   Format Conversion Scripts
Purpose:    Convert LAMMPS text dump or data files to extended XYZ format
            with proper element type mapping.
Usage:      gpumdkit.sh -lmp2exyz <input_file> <element1> <element2> ...
            python lmp2exyz.py <input_file> <element1> <element2> ...
Arguments:
  input_file Input LAMMPS text dump or data file (units = metal)
  elementX   Chemical element symbols in order of atomic types
Output:
  dump.xyz   Converted structures in extxyz format, in the working directory
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-10-07
=============================================================================
"""

import re
import sys


def _inspect_input(filename):
    """Recognize file contents and recover the data-file box origin."""
    origin = [0.0, 0.0, 0.0]
    natoms = None
    ntypes = None
    with open(filename, encoding='utf-8-sig') as stream:
        for line in stream:
            text = line.split('#', 1)[0].strip()
            if text.startswith('ITEM: TIMESTEP'):
                return 'lammps-dump-text', None, None, None
            match = re.fullmatch(r'(\d+)\s+atoms', text)
            if match:
                natoms = int(match.group(1))
            match = re.fullmatch(r'(\d+)\s+atom\s+types', text)
            if match:
                ntypes = int(match.group(1))
            fields = text.split()
            if len(fields) == 4 and fields[2:] in (
                    ['xlo', 'xhi'], ['ylo', 'yhi'], ['zlo', 'zhi']):
                origin['xyz'.index(fields[2][0])] = float(fields[0])
            elif len(fields) == 5 and fields[3:] == ['abc', 'origin']:
                origin = [float(value) for value in fields[:3]]
            if natoms is not None and text == 'Atoms':
                return 'lammps-data', natoms, ntypes, origin
    raise ValueError('Input is empty or is not a supported LAMMPS text dump/data file.')


def lmp2exyz(dump_file, elements):
    from ase import units
    from ase.data import atomic_numbers
    from ase.io import read, write

    if not elements:
        raise ValueError('Provide element symbols in order of LAMMPS atom types.')
    for element in elements:
        if element not in atomic_numbers or atomic_numbers[element] == 0:
            raise ValueError(f'Invalid chemical element symbol: {element}')

    input_format, natoms, ntypes, origin = _inspect_input(dump_file)
    if input_format == 'lammps-data':
        if ntypes is not None and ntypes > len(elements):
            raise ValueError(
                f'Input declares {ntypes} atom types, but only {len(elements)} '
                'element symbols were provided.')
        # Supply the mapping explicitly: ASE otherwise guesses elements from
        # Masses, whose atomic numbers are not LAMMPS type IDs.
        frame = read(
            dump_file, format='lammps-data', units='metal',
            Z_of_type={i + 1: atomic_numbers[e] for i, e in enumerate(elements)},
            sort_by_id=True, read_image_flags=False)
        if len(frame) != natoms:
            raise ValueError(f'Expected {natoms} atoms, but read {len(frame)}.')
        # GPUMD uses a cell with zero origin. Translate the whole structure,
        # then wrap it into the cell without applying trajectory image flags.
        frame.positions -= origin
        frame.set_celldisp([0.0, 0.0, 0.0])
        frame.wrap()
        # Preserve ASE arrays and also expose GPUMD's recognized property names.
        # ASE velocities * units.fs converts to GPUMD's Angstrom/fs.
        if frame.has('momenta'):
            frame.set_array('vel', frame.get_velocities() * units.fs)
        if frame.has('masses'):
            frame.set_array('mass', frame.get_masses())
        if frame.has('initial_charges'):
            frame.set_array('charge', frame.get_initial_charges())
        frames = [frame]
    else:
        # Keep the original dump reader and all frames/additional properties.
        frames = read(dump_file, format='lammps-dump-text', index=':')
        type_to_element = {i + 1: e for i, e in enumerate(elements)}
        for frame in frames:
            # Prefer actual type IDs when ASE provides them (e.g. dumps that
            # also contain a mass or element column).
            types = frame.arrays.get('type', frame.get_atomic_numbers())
            if len(types) == 0:
                raise ValueError('A LAMMPS dump frame contains no atoms.')
            if min(types) < 1 or max(types) > len(elements):
                raise ValueError(
                    f'Found atomic type outside 1..{len(elements)}: '
                    f'min={min(types)}, max={max(types)}.')
            frame.set_chemical_symbols([type_to_element[t] for t in types])

    # Do not create/truncate dump.xyz if parsing produced no structures.
    if not frames or any(len(frame) == 0 for frame in frames):
        raise ValueError('No atoms/frames were read; dump.xyz was not written.')
    write('dump.xyz', frames, format='extxyz')


def main():
    args = sys.argv[1:]
    if len(args) < 2 or args[0] in ('-h', '--help'):
        print(' Usage: gpumdkit.sh -lmp2exyz <input_file> <element1> <element2> ...')
        print('    or: python lmp2exyz.py <input_file> <element1> <element2> ...')
        print('\n Arguments:')
        print('   input_file Input LAMMPS text dump or data file (units = metal)')
        print('   elementX   Chemical element symbols in order of atomic types')
        print('\n Output:')
        print('   dump.xyz   Converted structures in extxyz format')
        print('\n Examples:')
        print('   gpumdkit.sh -lmp2exyz dump.lammps Si O')
        print('   python lmp2exyz.py mixSi.data B N O H C Si\n')
        return 0 if args and args[0] in ('-h', '--help') else 1
    try:
        lmp2exyz(args[0], args[1:])
    except (OSError, ValueError, KeyError, IndexError) as error:
        print(f' Error: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

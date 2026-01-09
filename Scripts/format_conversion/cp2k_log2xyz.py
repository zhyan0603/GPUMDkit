#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert validated CP2K calculations into a extended XYZ file.

Supports two processing modes:

**MODE 1 - Manual Selection:**
  Scans for .inp or .xyz files, lets you manually select which to process,
  finds matching .log files, and extracts data.

**MODE 2 - Auto Batch Processing (Recommended):**
  Automatically detects all folders containing .log files.
  For each folder, prioritizes .xyz files over .inp files for coordinates.
  Processes all folders without user interaction.

Extracts from each file:
  - atomic coordinates (from .inp or .xyz),
  - lattice vectors (from .inp &CELL or .xyz Lattice),
  - total energy (Hartree → eV) (from .log),
  - forces (Hartree/Bohr → eV/Å) (from .log),
  - stress tensor with automatic unit detection (Pa, MPa, GPa, atm, bar, kbar),
  - converts stress → GPa → virial (eV) using cell volume.

Skips frames missing any of {coordinates, energy, forces, stress}.
Outputs all valid frames to 'cp2k_selected.xyz' in extended XYZ format
compatible with ASE/NEP/MACE, with 'dirID' metadata for tracking.

Usage:
    python cp2k_log2xyz.py

Requirements:
  - .inp: &COORD and &CELL blocks.
  - .log: 'ENERGY| Total FORCE_EVAL', 'ATOMIC FORCES', and 'Analytical stress tensor [unit]'.
         Supported pressure units: Pa, MPa, GPa, atm, bar, kbar

Output fields:
  Lattice="Ax Ay Az Bx By Bz Cx Cy Cz"
  energy=<eV> virial="<9×eV>" pbc="T T T"
  Properties=species:S:1:pos:R:3:force:R:3
  dirID="<parent_dir>"

Conversion Summary:
  Energy:  Hartree × 27.2113838565563 = eV
  Forces:  Hartree/Bohr × 51.4220631857 = eV/Å
  Stress:  [Original unit] → GPa → Virial(eV)
           Virial = Stress(GPa) × Volume(Å³) / 160.2176634
"""

import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple
import numpy as np

def natural_sort_key(name: str):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', str(name))]

def find_log_for_inp(inp: Path) -> Optional[Path]:
    """Find corresponding .log file for .inp file.
    
    Tries multiple naming conventions:
    1. Same name with .log extension
    2. cp2k-{name}.log
    3. cp2k_{name}.log
    4. cp2k.log in same directory
    5. output.log in same directory
    6. {name}_output.log
    """
    candidates = [
        inp.with_suffix('.log'),
        inp.parent / f"cp2k-{inp.stem}.log",
        inp.parent / f"cp2k_{inp.stem}.log",
        inp.parent / "cp2k.log",
        inp.parent / "output.log",
        inp.parent / f"{inp.stem}_output.log",
    ]
    for log in candidates:
        if log.is_file():
            return log
    # If none found, return None (will prompt user to select)
    return None

def find_log_file_interactive(file_path: Path) -> Optional[Path]:
    """Interactively select a .log file for the given coordinate/inp file.
    
    First tries to find a matching log file automatically.
    If multiple or no matches found, prompts user to select.
    """
    # First, try automatic matching
    auto_log = find_log_for_inp(file_path.with_suffix('.inp'))
    if auto_log and auto_log.is_file():
        return auto_log
    
    # If not found, search for any .log files in the same directory
    log_candidates = list(file_path.parent.glob("*.log"))
    
    if not log_candidates:
        return None  # No log files found
    
    if len(log_candidates) == 1:
        return log_candidates[0]  # Only one log file, use it
    
    # Multiple log files found, let user choose
    print(f"\nFound {len(log_candidates)} .log files for: {file_path.name}")
    print("Please select which one to use:\n")
    for i, log in enumerate(sorted(log_candidates), 1):
        print(f"{i:3d}. {log.name}")
    
    while True:
        try:
            choice = input("\nYour choice (number): ").strip()
            idx = int(choice) - 1
            if 0 <= idx < len(log_candidates):
                selected_log = sorted(log_candidates)[idx]
                print(f"Selected: {selected_log.name}\n")
                return selected_log
            else:
                print(f"Please enter a number between 1 and {len(log_candidates)}")
        except ValueError:
            print("Invalid input. Please enter a number.")

def extract_lattice(content: str) -> List[float]:
    match = re.search(r'&CELL.*?A\s+([\d\.\-E\+\s]+)\s+B\s+([\d\.\-E\+\s]+)\s+C\s+([\d\.\-E\+\s]+).*?&END CELL', content, re.DOTALL | re.IGNORECASE)
    return [float(x) for group in match.groups() for x in group.split()] if match else [0.0] * 9

def compute_volume_from_lattice(lattice: List[float]) -> float:
    if len(lattice) != 9:
        return 1.0
    a = np.array(lattice[0:3])
    b = np.array(lattice[3:6])
    c = np.array(lattice[6:9])
    return abs(np.dot(a, np.cross(b, c)))

def extract_atoms(content: str) -> List[List]:
    match = re.search(r'&COORD\s*(.*?)&END COORD', content, re.DOTALL | re.IGNORECASE)
    if not match:
        return []
    atoms = []
    for line in match.group(1).strip().splitlines():
        line = line.strip()
        if line and not line.startswith('#'):
            parts = line.split()
            if len(parts) >= 4:
                atoms.append([parts[0]] + [float(x) for x in parts[1:4]])
    return atoms

def extract_energy(content: str) -> float:
    match = re.search(r'ENERGY\| Total FORCE_EVAL.*?:\s+([\d\.\-E\+]+)', content)
    return float(match.group(1)) * 27.2113838565563 if match else float('nan')

def extract_stress_in_gpa(content: str) -> Tuple[List[float], str]:
    """Extract stress tensor and identify its unit.
    
    Supports units: Pa, MPa, GPa, atm, bar, kbar
    Returns: (stress_values_in_gpa, unit_found)
    """
    # Try to find stress tensor with different units
    units_patterns = {
        'Pa': r'STRESS\|\s+Analytical stress tensor\s+\[Pa\]',
        'MPa': r'STRESS\|\s+Analytical stress tensor\s+\[MPa\]',
        'GPa': r'STRESS\|\s+Analytical stress tensor\s+\[GPa\]',
        'atm': r'STRESS\|\s+Analytical stress tensor\s+\[atm\]',
        'bar': r'STRESS\|\s+Analytical stress tensor\s+\[bar\]',
        'kbar': r'STRESS\|\s+Analytical stress tensor\s+\[kbar\]',
    }
    
    unit_found = None
    for unit, pattern in units_patterns.items():
        if re.search(pattern, content):
            unit_found = unit
            break
    
    if not unit_found:
        return [], 'unknown'
    
    # Extract stress values
    pattern = r'STRESS\|\s+[xyz]\s+([-\d\.E+\-]+)\s+([-\d\.E+\-]+)\s+([-\d\.E+\-]+)'
    matches = re.findall(pattern, content, re.IGNORECASE)
    if len(matches) < 3:
        return [], unit_found
    
    stress_values = [float(val) for row in matches[:3] for val in row]
    
    # Conversion factors to GPa
    conversions = {
        'Pa': 1e-9,
        'MPa': 1e-3,
        'GPa': 1.0,
        'atm': 0.0001013249848,  # 1 atm = 0.0001013249848 GPa
        'bar': 0.0001,  # 1 bar = 0.0001 GPa
        'kbar': 0.1,  # 1 kbar = 0.1 GPa
    }
    
    conversion_factor = conversions.get(unit_found, 1.0)
    stress_gpa = [s * conversion_factor for s in stress_values]
    
    return stress_gpa, unit_found

def stress_to_virial(stress_gpa: List[float], volume_ang3: float) -> List[float]:
    if len(stress_gpa) != 9 or volume_ang3 <= 0:
        return [0.0] * 9
    factor = volume_ang3 / 160.2176634
    return [s * factor for s in stress_gpa]

def extract_forces(content: str) -> List[List[float]]:
    match = re.search(r'ATOMIC FORCES in \[a\.u\.\]\n\n # Atom\s+Kind\s+Element\s+X\s+Y\s+Z\n(.*?)(?=\n SUM OF ATOMIC FORCES)', content, re.DOTALL)
    if not match:
        return []
    forces = []
    for line in match.group(1).strip().splitlines():
        parts = line.split()
        if len(parts) >= 6:
            fx, fy, fz = [float(parts[i]) * 51.4220631857 for i in range(3, 6)]
            forces.append([fx, fy, fz])
    return forces

def format_xyz_frame(lattice, atoms, energy, virial, forces, mat_id: str) -> str:
    lines = [f"{len(atoms)}\n"]
    lattice_str = ' '.join(f"{x:.10f}" for x in lattice)
    virial_str = ' '.join(f"{x:.10f}" for x in virial)
    info = [
        f'Lattice="{lattice_str}"',
        f'energy={energy:.10f}',
        f'virial="{virial_str}"',
        'pbc="T T T"',
        'Properties=species:S:1:pos:R:3:force:R:3',
        f'dirID="{mat_id}"'
    ]
    lines.append(' '.join(info) + '\n')
    for i, atom in enumerate(atoms):
        el, x, y, z = atom
        fx = fy = fz = 0.0
        if i < len(forces):
            fx, fy, fz = forces[i]
        lines.append(f"{el:2s} {x:14.8f} {y:14.8f} {z:14.8f} {fx:14.8f} {fy:14.8f} {fz:14.8f}\n")
    return ''.join(lines)

def select_files_interactive(inp_files: List[Path], source_type: str = 'inp') -> List[Path]:
    """Interactive file selection with appropriate prompts based on source type.
    
    Args:
        inp_files: List of files to select from
        source_type: 'inp' for input files or 'xyz' for coordinate files
    """
    file_type_name = 'input (.inp)' if source_type == 'inp' else 'coordinate (.xyz)'
    print(f"\nFound {len(inp_files)} {file_type_name} files:\n")
    for i, f in enumerate(inp_files, 1):
        print(f"{i:3d}. {f.relative_to(Path.cwd())}")
    print("\nEnter your choice:")
    print("  - Single number (e.g., 3)")
    print("  - Multiple numbers separated by commas (e.g., 1,4,7)")
    print("  - Range (e.g., 2-5)")
    print("  - 'all' to select all")
    print("  - 'quit' to exit")
    while True:
        choice = input("\nYour selection: ").strip()
        if choice.lower() in ('quit', 'q'):
            sys.exit(0)
        elif choice.lower() == 'all':
            return inp_files
        else:
            selected = []
            try:
                for part in choice.split(','):
                    part = part.strip()
                    if '-' in part:
                        start, end = map(int, part.split('-'))
                        selected.extend(range(start-1, end))
                    else:
                        selected.append(int(part)-1)
                return [inp_files[i] for i in selected if 0 <= i < len(inp_files)]
            except (ValueError, IndexError):
                print("Invalid input. Please try again.")

def select_coordinate_source():
    """Let user select coordinate source"""
    print("\nSelect processing mode:")
    print("  1. Manual selection: Choose files manually")
    print("  2. Auto batch mode: Auto-detect folders with .log files")
    while True:
        choice = input("\nYour choice (1-2): ").strip()
        if choice in ('1', '2'):
            return choice
        print("Invalid input. Please try again.")

def find_folders_with_logs(current: Path) -> List[Path]:
    """Find all folders containing .log files.
    
    Returns list of folder paths that contain .log files.
    """
    folders_with_logs = set()
    for log_file in current.rglob("*.log"):
        folders_with_logs.add(log_file.parent)
    return sorted(folders_with_logs, key=natural_sort_key)

def get_coord_files_for_folder(folder: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """Get coordinate files for a folder with .log file.
    
    Returns (xyz_file, inp_file) - prioritizes .xyz over .inp
    """
    xyz_files = list(folder.glob("*.xyz"))
    inp_files = list(folder.glob("*.inp"))
    
    # Prioritize .xyz files
    if xyz_files:
        # Return first xyz file found
        return sorted(xyz_files, key=natural_sort_key)[0], None
    elif inp_files:
        # Return first inp file found
        return None, sorted(inp_files, key=natural_sort_key)[0]
    else:
        return None, None

def process_batch_folders(current: Path) -> List[Tuple[Path, str]]:
    """Process all folders containing .log files in batch mode.
    
    Returns list of (file_path, source_type) tuples.
    """
    folders = find_folders_with_logs(current)
    
    if not folders:
        print("\nNo folders with .log files found.")
        return []
    
    print(f"\nFound {len(folders)} folder(s) with .log file(s):\n")
    
    selected_folders = []
    coord_files = []  # List of (file_path, source_type)
    
    for i, folder in enumerate(folders, 1):
        xyz_file, inp_file = get_coord_files_for_folder(folder)
        source_type = None
        source_file = None
        
        if xyz_file:
            source_type = 'xyz'
            source_file = xyz_file
            source_indicator = f"✓ {xyz_file.name} (XYZ)"
        elif inp_file:
            source_type = 'inp'
            source_file = inp_file
            source_indicator = f"✓ {inp_file.name} (INP)"
        else:
            source_indicator = "✗ No coordinate file found"
        
        rel_path = folder.relative_to(current)
        print(f"{i:3d}. {rel_path}")
        print(f"       {source_indicator}")
        
        if source_file:
            coord_files.append((source_file, source_type))
    
    return coord_files

def select_coordinate_source_manual():
    """Let user manually select coordinate source"""
    print("\nSelect coordinate source:")
    print("  1. Read coordinates from .inp file")
    print("  2. Read coordinates from .xyz file")
    while True:
        choice = input("\nYour choice (1-2): ").strip()
        if choice in ('1', '2'):
            return choice
        print("Invalid input. Please try again.")

def extract_atoms_from_xyz(xyz_content: str) -> List[List]:
    """Extract atomic information from XYZ file format"""
    lines = xyz_content.strip().split('\n')
    if len(lines) < 3:
        return []
    
    try:
        n_atoms = int(lines[0])
    except ValueError:
        return []
    
    atoms = []
    for i in range(2, min(2 + n_atoms, len(lines))):
        line = lines[i].strip()
        if line and not line.startswith('#'):
            parts = line.split()
            if len(parts) >= 4:
                try:
                    atoms.append([parts[0]] + [float(x) for x in parts[1:4]])
                except (ValueError, IndexError):
                    continue
    return atoms

def extract_lattice_from_xyz(xyz_content: str) -> List[float]:
    """Extract lattice vectors from XYZ file Lattice field"""
    match = re.search(r'Lattice="([^"]+)"', xyz_content)
    if match:
        try:
            return [float(x) for x in match.group(1).split()]
        except ValueError:
            pass
    return [0.0] * 9

def main():
    current = Path.cwd()
    
    # Select processing mode
    mode = select_coordinate_source()
    
    if mode == '2':
        # Batch auto-detection mode
        print("\n" + "=" * 70)
        print("AUTO BATCH MODE: Detecting folders with .log files...")
        print("=" * 70)
        
        coord_files = process_batch_folders(current)
        
        if not coord_files:
            print("\nNo coordinate files found in folders with .log files.")
            return
        
        print(f"\nAuto-selected {len(coord_files)} file(s) for processing")
        print(f"Priority: .xyz files preferred over .inp files\n")
        
        selected = [file_path for file_path, _ in coord_files]
        
    else:
        # Manual selection mode
        print("\n" + "=" * 70)
        print("MANUAL MODE: Please select coordinate source")
        print("=" * 70)
        
        coord_source = select_coordinate_source_manual()
        
        # Find files based on selection
        if coord_source == '1':
            # Only from .inp files
            inp_files = sorted(current.rglob("*.inp"), key=natural_sort_key)
            if not inp_files:
                print("No .inp files found in current directory or subdirectories.")
                return
            selected = select_files_interactive(inp_files, source_type='inp')
            coord_files = [(f, 'inp') for f in selected]
        else:  # '2'
            # Only from .xyz files
            xyz_files = sorted(current.rglob("*.xyz"), key=natural_sort_key)
            if not xyz_files:
                print("No .xyz files found in current directory or subdirectories.")
                return
            selected = select_files_interactive(xyz_files, source_type='xyz')
            coord_files = [(f, 'xyz') for f in selected]
    
    print(f"\nProcessing {len(selected)} selected file(s)...\n")

    success_count = 0
    error_counts = {
        'file: missing coordinates': 0,
        'log: missing file': 0,
        'log: missing energy': 0,
        'log: missing forces': 0,
        'log: missing stress': 0,
    }
    failed_details = []
    all_xyz_frames = []
    stress_units_found = {}  # Track stress units found and their conversion factors

    for file_path, source_type in coord_files:
        rel_path = file_path.relative_to(current)
        try:
            # Read coordinates based on source type
            if source_type == 'xyz':
                # Read from .xyz file
                xyz_text = file_path.read_text(encoding='utf-8', errors='ignore')
                atoms = extract_atoms_from_xyz(xyz_text)
                lattice = extract_lattice_from_xyz(xyz_text)
                
                # Find corresponding log file with interactive selection if needed
                log_file = find_log_file_interactive(file_path)
            else:  # 'inp'
                # Read from .inp file
                inp_text = file_path.read_text(encoding='utf-8', errors='ignore')
                atoms = extract_atoms(inp_text)
                lattice = extract_lattice(inp_text)
                
                # Find log file with interactive selection if needed
                log_file = find_log_file_interactive(file_path)
            
            if not atoms:
                error = 'file: missing coordinates'
                error_counts[error] += 1
                failed_details.append(f"{rel_path} → {error}")
                continue

            if not log_file:
                error = 'log: missing file'
                error_counts[error] += 1
                failed_details.append(f"{rel_path} → {error}")
                continue

            log_text = log_file.read_text(encoding='utf-8', errors='ignore')

            energy = extract_energy(log_text)
            forces = extract_forces(log_text)
            stress, stress_unit = extract_stress_in_gpa(log_text)

            # Validate presence
            missing = []
            if energy != energy:  # NaN check
                missing.append('energy')
            if not forces:
                missing.append('forces')
            if len(stress) != 9:
                missing.append('stress')

            if missing:
                for m in missing:
                    error_counts[f'log: missing {m}'] += 1
                failed_details.append(f"{rel_path} → log: missing {' / '.join(missing)}")
                continue

            # All valid → collect frame
            volume = compute_volume_from_lattice(lattice)
            virial = stress_to_virial(stress, volume)
            mat_id = file_path.parent.name
            frame = format_xyz_frame(lattice, atoms, energy, virial, forces, mat_id)
            all_xyz_frames.append(frame)
            # Track unit and conversion factor
            conversions = {
                'Pa': 1e-9,
                'MPa': 1e-3,
                'GPa': 1.0,
                'atm': 0.0001013249848,
                'bar': 0.0001,
                'kbar': 0.1,
            }
            if stress_unit not in stress_units_found:
                stress_units_found[stress_unit] = conversions.get(stress_unit, 1.0)
            success_count += 1

        except Exception as e:
            msg = f"{rel_path} → exception: {e}"
            failed_details.append(msg)
            print(msg, file=sys.stderr)

    # Write unified XYZ if any success
    output_file = current / "cp2k_selected.xyz"
    if all_xyz_frames:
        output_file.write_text(''.join(all_xyz_frames), encoding='utf-8')
        print(f"\nSuccessfully wrote {success_count} structures to: {output_file}")

    # Summary
    total = len(selected)
    failed = total - success_count
    print("\n" + "=" * 60)
    print(f"Summary: {success_count} succeeded, {failed} failed")
    for category, count in error_counts.items():
        if count > 0:
            print(f"  - {category}: {count}")
    if failed_details:
        print("\nFailed entries:")
        for d in failed_details:
            print(f"  {d}")
    
    # Print unit conversion coefficients
    print("\n" + "=" * 70)
    print("UNIT CONVERSION SUMMARY")
    print("=" * 70)
    
    print("\n1. ENERGY CONVERSION:")
    print("   Input unit:  Hartree (from CP2K)")
    print("   Output unit: eV (for NEP/MACE/ASE)")
    print("   Conversion factor: 1 Hartree = 27.2113838565563 eV")
    
    print("\n2. FORCE CONVERSION:")
    print("   Input unit:  Hartree/Bohr (from CP2K)")
    print("   Output unit: eV/Å (for NEP/MACE/ASE)")
    print("   Conversion factor: 1 Hartree/Bohr = 51.4220631857 eV/Å")
    
    print("\n3. STRESS CONVERSION:")
    print("   Detection: Automatically scans for pressure unit in .log file")
    print("   Supported input units: Pa, MPa, GPa, atm, bar, kbar")
    print("   Intermediate unit: GPa (for consistency)")
    print("   Final unit: Virial in eV (for NEP/MACE/ASE)")
    print("\n   Conversion factors to GPa:")
    print("     - Pa       → GPa: 1e-9")
    print("     - MPa      → GPa: 1e-3")
    print("     - bar      → GPa: 1e-4 (0.0001)")
    print("     - kbar     → GPa: 0.1")
    print("     - atm      → GPa: 0.0001013249848")
    print("     - GPa      → GPa: 1.0 (no conversion)")
    
    if stress_units_found:
        print(f"\n   Units detected in this run:")
        for unit, factor in sorted(stress_units_found.items()):
            print(f"     ✓ {unit:8s} → GPa: conversion factor = {factor}")
    
    print("\n4. VIRIAL CONVERSION:")
    print("   Formula: Virial(eV) = Stress(GPa) × Volume(Å³) / 160.2176634")
    print("   Where: 160.2176634 ≈ eV·Å³/GPa (unit conversion constant)")
    print("\n" + "=" * 70)

if __name__ == '__main__':
    main()
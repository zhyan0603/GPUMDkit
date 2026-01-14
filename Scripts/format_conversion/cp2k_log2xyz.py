#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Convert CP2K outputs (.log) together with structure files (.xyz or .inp) into an extended XYZ format
# containing energy, atomic forces, and virial tensor.
# Requirements: Each subdirectory must contain at least one .xyz/.inp and one .log file (cp2k.log preferred).
# Outputs: cp2k_exyz.xyz (extended XYZ trajectory), Logfile.txt (processing summary and diagnostics).

import re
from pathlib import Path
import numpy as np

HARTREE_TO_EV = 27.2113838565563
FORCE_AU_TO_EV_ANG = 51.4220631857
GPA_TO_EV_FACTOR = 160.2176634

def natural_sort_key(p): return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(p))]

def read_lines(p): return p.read_text(encoding='utf-8', errors='ignore').splitlines()

def extract_atoms_lattice_xyz(lines):
    if len(lines) < 3: return [], []
    try:
        n = int(lines[0])
        atoms = []
        for l in lines[2:2+n]:
            parts = l.split()
            if len(parts) >= 4:
                atoms.append([parts[0], float(parts[1]), float(parts[2]), float(parts[3])])
        lat = [0.0]*9
        if m := re.search(r'Lattice="([^"]+)"', '\n'.join(lines)):
            lat = [float(x) for x in m.group(1).split()]
        return atoms, lat
    except: return [], []

def extract_atoms_lattice_inp(text):
    atoms, lat = [], [0.0]*9
    if m := re.search(r'&COORD\s*(.*?)&END COORD', text, re.DOTALL | re.IGNORECASE):
        for line in m.group(1).strip().splitlines():
            if line.strip() and not line.startswith('#'):
                p = line.split()
                if len(p) >= 4:
                    atoms.append([p[0], float(p[1]), float(p[2]), float(p[3])])
    if m := re.search(r'&CELL.*?A\s+([\d\.\-E\+\s]+)\s+B\s+([\d\.\-E\+\s]+)\s+C\s+([\d\.\-E\+\s]+).*?&END CELL', text, re.DOTALL | re.IGNORECASE):
        lat = [float(x) for g in m.groups() for x in g.split()]
    return atoms, lat

def extract_mixed_atoms_lattice(folder):
    """Try to get atoms from .xyz and lattice from .inp; fall back to single-file mode."""
    xyzs = sorted(folder.glob("*.xyz"), key=natural_sort_key)
    inps = sorted(folder.glob("*.inp"), key=natural_sort_key)

    if xyzs and inps:
        try:
            atoms, _ = extract_atoms_lattice_xyz(read_lines(xyzs[0]))
            _, lat = extract_atoms_lattice_inp(inps[0].read_text(encoding='utf-8', errors='ignore'))
            if atoms and len(lat) == 9:
                return atoms, lat, f"{xyzs[0].name} (atoms) + {inps[0].name} (lattice)", "mixed"
        except Exception:
            pass

    if xyzs:
        try:
            atoms, lat = extract_atoms_lattice_xyz(read_lines(xyzs[0]))
            if atoms:
                return atoms, lat, xyzs[0].name, "xyz"
        except Exception:
            pass

    if inps:
        try:
            text = inps[0].read_text(encoding='utf-8', errors='ignore')
            atoms, lat = extract_atoms_lattice_inp(text)
            if atoms:
                return atoms, lat, inps[0].name, "inp"
        except Exception:
            pass

    return None, None, None, None

def extract_energy(log):
    """
    Extract total energy from CP2K log.
    Supports:
      - Old format: "ENERGY| Total FORCE_EVAL ... : <value>"
      - New format: "ENERGY| Total FORCE_EVAL (...) energy [hartree] <value>"
    Returns energy in eV.
    """
    # Try NEW format first (no colon, value after [hartree])
    if m := re.search(r'ENERGY\| Total FORCE_EVAL \(.*?\) energy \[hartree\]\s+([-\d\.E+\-]+)', log):
        return float(m.group(1)) * HARTREE_TO_EV

    # Try OLD format (with colon)
    if m := re.search(r'ENERGY\| Total FORCE_EVAL.*?:\s+([-\d\.E+\-]+)', log):
        return float(m.group(1)) * HARTREE_TO_EV

    return float('nan')

def extract_forces(log):
    """
    Extract atomic forces from CP2K log.
    Supports:
      - New format: "FORCES| Atomic forces [hartree/bohr]" (CP2K 2025+)
      - Old format: "ATOMIC FORCES in [a.u.]" with element table
    Returns forces in eV/angstrom as list of [fx, fy, fz].
    """
    forces = []

    # --- Try NEW format (CP2K 2025+) ---
    if 'FORCES| Atomic forces [hartree/bohr]' in log:
        lines = log.splitlines()
        in_block = False
        for line in lines:
            if line.strip().startswith('FORCES| Atomic forces [hartree/bohr]'):
                in_block = True
                continue
            if in_block:
                if not line.strip().startswith('FORCES|'):
                    break
                parts = line.split()
                if len(parts) >= 5 and parts[1].isdigit():
                    try:
                        fx = float(parts[2]) * FORCE_AU_TO_EV_ANG
                        fy = float(parts[3]) * FORCE_AU_TO_EV_ANG
                        fz = float(parts[4]) * FORCE_AU_TO_EV_ANG
                        forces.append([fx, fy, fz])
                    except (ValueError, IndexError):
                        continue
        if forces:
            return forces

    # --- Try OLD format ---
    if m := re.search(r'ATOMIC FORCES in \[a\.u\.\]\n\n # Atom\s+Kind\s+Element\s+X\s+Y\s+Z\n(.*?)(?=\n SUM OF ATOMIC FORCES)', log, re.DOTALL):
        for line in m.group(1).strip().splitlines():
            p = line.split()
            if len(p) >= 6:
                fx, fy, fz = [float(p[i]) * FORCE_AU_TO_EV_ANG for i in range(3,6)]
                forces.append([fx, fy, fz])
        if forces:
            return forces

    return []

def extract_stress_gpa(log):
    # Try GPa first
    if re.search(r'STRESS\|\s+Analytical stress tensor\s+\[GPa\]', log):
        matches = re.findall(r'STRESS\|\s+[xyz]\s+([-\d\.E+\-]+)\s+([-\d\.E+\-]+)\s+([-\d\.E+\-]+)', log, re.IGNORECASE)
        if len(matches) >= 3:
            return [float(v) for row in matches[:3] for v in row]
    
    # Try bar (CP2K 2025+)
    if re.search(r'STRESS\|\s+Analytical stress tensor\s+\[bar\]', log):
        matches = re.findall(r'STRESS\|\s+[xyz]\s+([-\d\.E+\-]+)\s+([-\d\.E+\-]+)\s+([-\d\.E+\-]+)', log, re.IGNORECASE)
        if len(matches) >= 3:
            stress_bar = [float(v) for row in matches[:3] for v in row]
            # Convert bar to GPa: 1 bar = 0.1 GPa
            return [s * 0.1 for s in stress_bar]
    
    return None

def compute_volume(lat):
    if len(lat) != 9: return 1.0
    a, b, c = map(np.array, [lat[:3], lat[3:6], lat[6:9]])
    return abs(np.dot(a, np.cross(b, c)))

def stress_to_virial(stress, vol):
    if len(stress) != 9 or vol <= 0: return [0.0]*9
    return [s * vol / GPA_TO_EV_FACTOR for s in stress]

def format_frame(lat, atoms, energy, virial, forces, dir_id):
    lines = [f"{len(atoms)}\n"]
    lat_str = ' '.join(f"{x:.10f}" for x in lat)
    vir_str = ' '.join(f"{x:.10f}" for x in virial)
    info = [
        f'Lattice="{lat_str}"',
        f'energy={energy:.10f}',
        f'virial="{vir_str}"',
        'pbc="T T T"',
        'Properties=species:S:1:pos:R:3:force:R:3',
        f'dirID="{dir_id}"'
    ]
    lines.append(' '.join(info) + '\n')
    for i, (el, x, y, z) in enumerate(atoms):
        fx, fy, fz = forces[i] if i < len(forces) else (0.0, 0.0, 0.0)
        lines.append(f"{el:2s} {x:14.8f} {y:14.8f} {z:14.8f} {fx:14.8f} {fy:14.8f} {fz:14.8f}\n")
    return ''.join(lines)

def main():
    cwd = Path.cwd()
    log_paths = sorted((p for p in cwd.rglob("*.log") if p.is_file() and p.stat().st_size > 0), key=natural_sort_key)
    folders = sorted({p.parent for p in log_paths}, key=natural_sort_key)

    if not folders:
        with open(cwd / "Logfile.txt", "w") as f:
            f.write("No .log files found.\n")
        print("No .log files found.")
        return

    all_frames = []
    success, failed, warnings = [], [], []

    for folder in folders:
        rel = str(folder.relative_to(cwd))

        atoms, lat, source_desc, mode = extract_mixed_atoms_lattice(folder)

        if atoms is None:
            failed.append(f"{rel} -> No valid structure files")
            continue

        xyzs = sorted(folder.glob("*.xyz"), key=natural_sort_key)
        inps = sorted(folder.glob("*.inp"), key=natural_sort_key)
        if len(xyzs) > 1:
            warnings.append(f"{rel} -> Multiple .xyz files; used one for atoms")
        if len(inps) > 1:
            warnings.append(f"{rel} -> Multiple .inp files; used one for lattice")

        # Choose log: prefer cp2k.log, else first .log
        logs = sorted(folder.glob("*.log"), key=natural_sort_key)
        log_file = folder / "cp2k.log" if (folder / "cp2k.log") in logs else logs[0]

        try:
            # Note: atoms and lat are already extracted above — no need to re-parse coord_file
            log_text = log_file.read_text(encoding='utf-8', errors='ignore')
            energy = extract_energy(log_text)
            forces = extract_forces(log_text)
            stress = extract_stress_gpa(log_text)

            if stress is None:
                failed.append(f"{rel} -> Stress not in [GPa] or missing")
                continue
            if energy != energy or not forces:
                missing = []
                if energy != energy: missing.append("energy")
                if not forces: missing.append("forces")
                failed.append(f"{rel} -> Missing: {' / '.join(missing)}")
                continue

            vol = compute_volume(lat)
            virial = stress_to_virial(stress, vol)
            frame = format_frame(lat, atoms, energy, virial, forces, folder.name)
            all_frames.append(frame)
            success.append(f"{rel} -> Used {source_desc} + {log_file.name}")

        except Exception as e:
            failed.append(f"{rel} -> Exception: {e}")

    # --- Write Logfile.txt ---
    with open(cwd / "Logfile.txt", "w") as f:
        w = lambda s="": f.write(s + "\n")
        w("-" * 60)
        w("UNIT CONVERSION CONSTANTS")
        w("-" * 60)
        w(f"Energy:   1 Hartree      = {HARTREE_TO_EV:.10f} eV")
        w(f"Forces:   1 Hartree/Bohr = {FORCE_AU_TO_EV_ANG:.10f} eV/angstrom")
        w(f"Virial(eV) = Stress(GPa) * Volume(Å³) / {GPA_TO_EV_FACTOR:.7f}")
        w("")
        w("-" * 60)
        w("SUMMARY STATISTICS")
        w("-" * 60)
        w(f"{'Item':<25} {'Count'}")
        w(f"{'-'*25} {'-----'}")
        w(f"{'Log folders processed':<25} {len(folders)}")
        w(f"{'Successfully converted':<25} {len(success)}")
        w(f"{'Failed conversions':<25} {len(failed)}")
        w(f"{'Warnings issued':<25} {len(warnings)}")
        w("")
        w("-" * 60)

        if failed:
            w("FAILED ENTRIES")
            w("-" * 60)
            for item in failed: w(item)
            w("")

        if warnings:
            w("WARNINGS")
            w("-" * 60)
            for item in warnings: w(item)
            w("")

        if success:
            w("SUCCESSFUL CONVERSIONS (used files)")
            w("-" * 60)
            for item in success: w(item)

    if all_frames:
        (cwd / "cp2k_exyz.xyz").write_text(''.join(all_frames))
        print("\nProcessing complete.")
        print("Log saved to: Logfile.txt")
        print("Output written to: cp2k_exyz.xyz")
    else:
        print("\nProcessing complete. No valid structures converted.")
        print("Log saved to: Logfile.txt")

if __name__ == '__main__':
    main()

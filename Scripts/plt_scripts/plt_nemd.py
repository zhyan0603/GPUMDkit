"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_nemd.py
Category:   Plot Scripts
Purpose:    Post-processing for NEMD (Non-Equilibrium Molecular Dynamics)
            thermal conductivity calculations, including temperature profile,
            thermostat energy, and optional SHC spectral analysis.
Usage:      gpumdkit.sh -plt nemd [real_length] [scale_eff_size] [cutoff_freq] [save] [quantum]
            python plt_nemd.py [real_length] [scale_eff_size] [cutoff_freq] [save] [quantum]
Arguments:
  real_length     Real length of heat transfer zone in nm (set to 'Auto' for auto)
  scale_eff_size  Scale factor for effective cross-sectional area (default: 1)
  cutoff_freq     Cutoff frequency for SHC in THz (default: 60)
  save            Save the plot as 'nemd.png' instead of displaying it
  quantum         Also apply the quantum correction to the SHC spectrum
Output:
  nemd.png    (if save is used, or if backend is non-interactive)
Author:     Xin Wu (xinwuchn97@gmial.com)
Last-modified: 2026-10-21
=============================================================================
"""

import ast
import operator as op
import os
import sys
from functools import lru_cache
from io import StringIO
from math import isfinite

import matplotlib.pyplot as plt
import numpy as np
from ase.io import read

FS = 12
TIME_UNIT_FS = 10.18051  # GPUMD 3.3.1 internal time unit.
trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
plt.rcParams.update({
    'font.size': FS, 'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'], 'axes.linewidth': 1.2,
})


def set_fig_properties(ax_list, tl=4, tw=1.2, tlm=4):
    """Apply matching tick properties to all axes."""
    for ax in ax_list:
        ax.tick_params(which='both', length=tl, width=tw, direction='in', right=True, top=True)
        ax.tick_params(which='minor', length=tlm)


def quantum_correction_factor(nu, T):
    """Return the mode heat-capacity ratio for frequency in THz and temperature in K."""
    h_over_kB = 6.62607015e-34 / 1.380649e-23  # K*s
    x = h_over_kB * np.asarray(nu, dtype=float) * 1e12 / T
    with np.errstate(divide='ignore', invalid='ignore'):
        factor = x**2 * np.exp(-x) / np.expm1(-x)**2
    return np.where(x == 0, 1.0, factor)


@lru_cache(maxsize=None)
def element_from_mass(mass):
    """Infer an element from its closest standard atomic mass."""
    from ase.data import atomic_masses, chemical_symbols

    return chemical_symbols[int(np.argmin(np.abs(atomic_masses[1:] - mass))) + 1]


def xyz_in_to_model(source, target):
    """Convert legacy GPUMD xyz.in to extended XYZ, preserving masses, velocities and groups."""
    with open(source, encoding="utf-8-sig") as inp:
        header = inp.readline().split()
        if len(header) != 6:
            raise ValueError("Line 1 must contain N MN rc triclinic has_velocity num_grouping_methods.")
        n, triclinic, has_velocity, ng = map(int, (header[0], *header[3:]))
        if n < 1 or triclinic not in (0, 1) or has_velocity not in (0, 1) or not 0 <= ng <= 10:
            raise ValueError("Invalid atom count, box flag, velocity flag or group count.")

        box = inp.readline().split()
        if len(box) != (12 if triclinic else 6) or any(v not in ("0", "1") for v in box[:3]):
            raise ValueError("Line 2 must contain three PBC flags followed by the cell.")
        if any(not isfinite(float(v)) for v in box[3:]):
            raise ValueError("Cell components must be finite.")
        lattice = box[3:] if triclinic else [box[3], "0", "0", "0", box[4], "0", "0", "0", box[5]]
        boundary = " ".join("T" if v == "1" else "F" for v in box[:3])
        properties = "species:S:1:pos:R:3:mass:R:1"
        if has_velocity:
            properties += ":vel:R:3"
        if ng:
            properties += f":group:I:{ng}"
        rows = [str(n), f'Lattice="{" ".join(lattice)}" Properties={properties} pbc="{boundary}"']

        type_elements = {}
        for index in range(n):
            fields = inp.readline().split()
            if len(fields) != 5 + 3 * has_velocity + ng:
                raise ValueError(f"Line {index + 3}: incorrect number of atom columns.")
            atom_type = int(fields[0])
            values = [float(v) for v in fields[1:5 + 3 * has_velocity]]
            if any(not isfinite(v) for v in values) or values[3] <= 0:
                raise ValueError(f"Line {index + 3}: invalid coordinates, mass or velocity.")
            labels = fields[5 + 3 * has_velocity:]
            if any(not 0 <= int(v) < n for v in labels):
                raise ValueError(f"Line {index + 3}: invalid group label.")
            element = element_from_mass(values[3])
            type_elements.setdefault(atom_type, set()).add((element, values[3]))
            row = [element, *fields[1:5]]
            if has_velocity:
                # Legacy xyz.in stores internal velocities; model.xyz uses A/fs.
                row += [f"{float(v) / TIME_UNIT_FS:.16g}" for v in fields[5:8]]
            rows.append(" ".join(row + labels))
        if inp.read().strip():
            raise ValueError("Extra data after the declared number of atoms.")

    with open(target, "w", encoding="utf-8") as out:
        out.write("\n".join(rows) + "\n")
    return n, type_elements


def prepare_model(paths):
    """Convert legacy xyz.in only when model.xyz is absent."""
    if os.path.exists(paths['model']) or not os.path.exists(paths['xyz_in']):
        return
    print("[INFO] model.xyz not found, converting xyz.in -> model.xyz ...")
    n, types = xyz_in_to_model(paths['xyz_in'], paths['model'])
    for atom_type, elements in sorted(types.items()):
        mapping = ", ".join(f"{element} (mass {mass:g})" for element, mass in sorted(elements))
        print(f"       type {atom_type} -> {mapping}")
    print(f"[INFO] Wrote {paths['model']}: {n} atoms.\n")


_ARITH_OPS = {
    ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul, ast.Div: op.truediv,
    ast.FloorDiv: op.floordiv, ast.Mod: op.mod, ast.Pow: op.pow,
    ast.UAdd: op.pos, ast.USub: op.neg,
}


def eval_number(expr):
    """Safely evaluate numeric literals and simple arithmetic expressions."""
    def _eval(node):
        if isinstance(node, ast.Expression):
            return _eval(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return node.value
        if isinstance(node, ast.BinOp) and type(node.op) in _ARITH_OPS:
            left, right = _eval(node.left), _eval(node.right)
            if isinstance(node.op, ast.Pow) and abs(right) > 20:
                # Guard against e.g. "10**99999999": CPython would try to compute
                # an exact, arbitrary-precision integer with millions of digits.
                raise ValueError(f"Exponent too large in {expr!r}")
            return _ARITH_OPS[type(node.op)](left, right)
        if isinstance(node, ast.UnaryOp) and type(node.op) in _ARITH_OPS:
            return _ARITH_OPS[type(node.op)](_eval(node.operand))
        raise ValueError(f"Unsupported expression: {expr!r}")
    try:
        tree = ast.parse(expr, mode='eval')
    except SyntaxError:
        raise ValueError(f"Unsupported expression: {expr!r}")
    try:
        return float(_eval(tree))
    except (ZeroDivisionError, OverflowError, TypeError) as e:
        raise ValueError(f"Could not evaluate {expr!r}: {e}")


def read_commands(path):
    """Read exact command tokens, excluding blank lines and comments."""
    with open(path, encoding="utf-8-sig") as source:
        return [tokens for line in source if (tokens := line.split("#", 1)[0].split())]


def last_command(commands, name):
    for command in reversed(commands):
        if command[0] == name:
            return command
    raise ValueError(f"Missing {name!r} command in run.in.")


def load_repeated(path, n_rows, names):
    """Load consecutive runs as (sample, repeat) arrays without a DataFrame."""
    if n_rows < 1:
        raise ValueError("The number of samples per run must be positive.")
    raw = np.loadtxt(path, ndmin=2)
    if raw.shape[1] != len(names):
        raise ValueError(f"{path}: expected {len(names)} columns, got {raw.shape[1]}.")
    if not len(raw) or len(raw) % n_rows:
        raise ValueError(f"{path}: incomplete data; expected a multiple of {n_rows} rows.")
    runs = raw.reshape(-1, n_rows, len(names)).transpose(1, 0, 2)
    return {name: runs[:, :, i] for i, name in enumerate(names)}


def add_statistics(values, sem=False):
    """Append the repeat mean and, optionally, its standard error."""
    columns = [values, values.mean(axis=1)]
    if sem:
        columns.append(values.std(axis=1) / np.sqrt(values.shape[1]))
    return np.column_stack(columns)


def load_shc(path, command):
    """Split each SHC block into correlation and spectral arrays."""
    n_corr = 2 * int(command[2]) - 1
    n_omega = int(command[4])
    raw = load_repeated(path, n_corr + n_omega, ('grid', 'in', 'out'))
    data = {'t': raw['grid'][:n_corr], 'omega': raw['grid'][n_corr:]}
    data['nu'] = data['omega'] / (2 * np.pi)
    for key, component in (('Ki', 'in'), ('Ko', 'out')):
        data[key] = raw[component][:n_corr]
        data['jw' + component[0]] = raw[component][n_corr:]
    return data, data['t'].shape[1]


def shc_geometry(model_path, command, scale_eff_size):
    """Apply the existing atom-fraction normalization for a selected group."""
    model = read(model_path)
    fraction = 1.0
    if 'group' in command:
        index = command.index('group')
        method, group_id = map(int, command[index + 1:index + 3])
        if group_id < 0:
            raise ValueError("SHC analysis requires a single group_id; group -1 outputs multiple groups.")
        labels = model.get_array('group')
        labels = labels if labels.ndim == 1 else labels[:, method]
        fraction = np.count_nonzero(labels == group_id) / len(labels)
        if not fraction:
            raise ValueError(f"SHC group {group_id} contains no atoms.")
    volume = model.get_volume() * fraction / scale_eff_size
    length = model.cell[int(command[3]), int(command[3])] * fraction
    return volume, length


def txt_block(name, unit, x_name, x_unit, x, arr, trailing_labels=("average",)):
    """Format a data block with the existing tab-separated export layout."""
    labels = [f"repeat_{i + 1}" for i in range(arr.shape[1] - len(trailing_labels))]
    header = [f"{x_name}({x_unit})", *labels, *trailing_labels]
    output = StringIO()
    np.savetxt(output, np.column_stack((x, arr)), delimiter="\t", fmt="%.6g")
    return [f"# ---- {name} ({unit}) ----", "# " + "\t".join(header),
            *output.getvalue().splitlines(), ""]


def print_usage():
    """Print usage instructions"""
    print("Usage: gpumdkit -plt nemd [real_length] [scale_eff_size] [cutoff_freq] [--save] [--save-data] [--quantum]")
    print("Params:")
    print("  real_length   : Real length of heat tranfer zone in nm (set to 'Auto', with auto-calculation)")
    print("  scale_eff_size: Optional, Scale factor for effective cross-sectional area (default: 1)")
    print("                   • For 3D bulk systems: use 1")
    print("                   • For low-dimensional systems with vacuum layer: S_box / S_eff")
    print("                     - S_box: box area perpendicular to heat transfer direction")
    print("                     - S_eff: real or effective area of the system")
    print("                   • May be a plain number or a simple arithmetic expression,")
    print("                     e.g. '100/2' (evaluates to 50)")
    print("  cutoff_freq   : Optional, Cutoff frequency for SHC calculation in THz (default: 60)")
    print("  --save        : Optional, save the plot as 'nemd.png' (default: show it interactively)")
    print("  --save-data   : Optional, also save the processed data as 'data_nemd.npz'/'data_shc.npz'")
    print("                  and their .txt equivalents (independent of --save; can be used on its own)")
    print("  --quantum     : Optional, also apply the quantum correction to the SHC spectral thermal")
    print("                  conductance, g_q(ω) = g(ω)·x²eˣ/(eˣ-1)², x = hν/(k_B·T), T from heat_lan in run.in")
    print("                  (requires shc.out)")
    print("  --save/--save-data/--quantum may appear anywhere on the command line, e.g.:")
    print("    gpumdkit -plt nemd --save-data              (defaults for everything else)")
    print("    gpumdkit -plt nemd Auto 1 60 --save --save-data --quantum")

class NEMD_Processor:
    def __init__(self, _directory, _real_length=None, _scale_eff_size=1, _cutoff_freq=60, _scale_vacuum=1,
                 _save=False, _save_data=False, _quantum=False):
        """Configure NEMD processing; lengths are in nm and frequencies in THz."""
        self.directory = _directory
        self.scale_vacuum = _scale_vacuum
        self.scale_eff_size = _scale_eff_size
        self.real_length = _real_length
        self.cutoff_freq = _cutoff_freq
        self.save = _save
        self.save_data = _save_data
        self.quantum = _quantum
        self.path = {key: os.path.join(_directory, filename) for key, filename in (
            ('run', 'run.in'), ('compute', 'compute.out'), ('model', 'model.xyz'),
            ('shc', 'shc.out'), ('xyz_in', 'xyz.in'),
        )}
        self.has_shc = os.path.exists(self.path['shc'])

    def prepare_model(self):
        """Convert a legacy xyz.in when model.xyz is absent."""
        prepare_model(self.path)

    def process_SHC(self, deltaT):
        """Compute classical and optional quantum-corrected spectral conductance."""
        commands = read_commands(self.path['run'])
        command = last_command(commands, 'compute_shc')
        data, n_repeat = load_shc(self.path['shc'], command)
        # NEMD historically scales this volume by scale_vacuum, not scale_eff_size.
        volume, length = shc_geometry(self.path['model'], command, self.scale_vacuum)
        convert = 1.602176634e7
        denom = volume * deltaT
        data['k_g_wi'] = data['jwi'] * convert / denom
        data['k_g_wo'] = data['jwo'] * convert / denom
        data['k_g_wt'] = data['k_g_wi'] + data['k_g_wo']
        # Preserve the original total-only cutoff and negative-value floor.
        data['k_g_wt'][data['nu'] > float(self.cutoff_freq)] = 0
        data['k_g_wt'][data['k_g_wt'] < 0] = 0.001
        data['Kt'] = data['Ki'] + data['Ko']
        data['L'] = length

        result_cols = [('in', 'k_g_wi'), ('out', 'k_g_wo'), ('tot', 'k_g_wt')]
        if self.quantum:
            heat = last_command([c for c in commands if c[:2] == ['ensemble', 'heat_lan']], 'ensemble')
            self.temperature = float(heat[2])
            factor = quantum_correction_factor(data['nu'], self.temperature)
            for key, col in list(result_cols):
                data[f'{col}_qc'] = data[col] * factor
                result_cols.append((f'{key}_qc', f'{col}_qc'))

        # Integration uses the original frequency spacing, including its first-bin convention.
        results = {}
        for key, col in result_cols:
            values = trap(data[col], dx=data['nu'][0, 0], axis=0)
            results[f'{key}_ave'] = values.mean()
            results[f'{key}_std'] = values.std() / np.sqrt(n_repeat)
        data = {key: value if key == 'L' else add_statistics(value, sem=True)
                for key, value in data.items()}
        data['Results'] = results

        if self.save_data:
            np.savez(os.path.join(self.directory, 'data_shc.npz'), **data)
            self._export_shc_txt(data)
        return data

    def process(self):
        """Process NEMD data, print results and plot the available panels."""
        self.prepare_model()
        direction = None
        found_nemd = False
        for command in read_commands(self.path['run']):
            if command[0] == 'time_step':
                time_step = float(command[1])
            elif command[:2] == ['ensemble', 'heat_lan']:
                n_temp_group = int(command[6]) + 1
            elif command[0] == 'compute' and 'temperature' in command[1:]:
                output_interval = int(command[2]) * int(command[3])
                found_nemd = True
            elif command[0] == 'compute_shc':
                direction = int(command[3])
            elif found_nemd and command[0] == 'run':
                nemd_run = int(command[1])

        n_rows = nemd_run // output_interval
        names = [f'temp_{i}' for i in range(n_temp_group)] + ['E_in', 'E_out']
        data = load_repeated(self.path['compute'], n_rows, names)
        n_repeat = data['E_in'].shape[1]
        time_upper = nemd_run * time_step * 1e-6  # ns
        for key, values in data.items():
            values = add_statistics(values)
            if key.startswith('temp_'):
                values = np.vstack((values, values[n_rows // 2:].mean(axis=0)))
            data[key] = values

        delta_T = (data['temp_1'][-1, :n_repeat]
                   - data[f'temp_{n_temp_group - 1}'][-1, :n_repeat]).reshape(1, -1)
        mid = n_rows // 2
        denom = (n_rows / 2) * (time_step * 0.001) * output_interval
        Q_in = (data['E_in'][mid, :n_repeat] - data['E_in'][-1, :n_repeat]) / denom
        Q_out = (data['E_out'][-1, :n_repeat] - data['E_out'][mid, :n_repeat]) / denom
        Q = ((Q_in + Q_out) / 2).reshape(1, -1)  # eV/ps

        model = read(self.path['model'])
        groups = model.get_array('group')
        group = groups if groups.ndim == 1 else groups[:, 0]
        heat_pos = model.positions[group == 1].mean(axis=0)
        cold_pos = model.positions[group == n_temp_group - 1].mean(axis=0)
        if direction is None:
            direction = int(np.argmax(np.abs(heat_pos - cold_pos)))
        coords_heat, coords_cold = heat_pos[direction], cold_pos[direction]
        xticks_length = np.linspace(coords_heat, coords_cold, n_temp_group - 1) * 0.1
        length = self.real_length if self.real_length is not None else abs(coords_heat - coords_cold) * 0.1
        data['L'] = length

        area = model.get_volume() / model.get_cell()[direction, direction] / self.scale_eff_size  # A^2
        G = 1.602176634e7 * Q / delta_T / area  # MW/(m^2 K)
        data['G'] = add_statistics(G, sem=True)
        data['k'] = add_statistics(G * length * 1e-3, sem=True)  # W/(m K)
        if self.save_data:
            np.savez(os.path.join(self.directory, 'data_nemd.npz'), **data)
            self._export_txt(data, n_rows, n_temp_group, xticks_length, time_upper, length)
        self._print_nemd_results(data, length)

        if self.has_shc:
            print("\n[INFO] SHC data detected, processing spectral heat current...")
            shc_data = self.process_SHC(deltaT=delta_T)
            self._print_shc_results(shc_data['Results'])
        else:
            print("\n[INFO] No SHC data found (shc.out not present), skipping SHC analysis.")
            if self.quantum:
                print("[WARNING] --quantum needs shc.out, quantum correction skipped.")
            shc_data = None
        self._plot_results(data, shc_data, time_upper, n_rows, n_repeat, xticks_length, length)

    def _print_nemd_results(self, Reformed_NEMD_data, Length):
        """Print NEMD thermal conductivity results"""
        print("\n" + "=" * 70)
        print("NEMD Thermal Conductivity Results")
        print("=" * 70)
        print(f"\nEffective length: {Length:.4f} nm")
        print(f"Scale_eff_size: {self.scale_eff_size}\n")
        print(f"Thermal conductance G = {Reformed_NEMD_data['G'][0, -2]:.6f} ± {Reformed_NEMD_data['G'][0, -1]:.6f} MW/(m²·K)")
        print(f"Thermal conductivity κ = {Reformed_NEMD_data['k'][0, -2]:.6f} ± {Reformed_NEMD_data['k'][0, -1]:.6f} W/(m·K)")
        print("=" * 70)

    def _print_shc_results(self, results):
        """Print SHC spectral thermal conductance results"""
        print("\n" + "=" * 70)
        print("SHC Spectral Thermal Conductance Results")
        print("=" * 70)
        print(f"\nCutoff frequency: {self.cutoff_freq} THz\n")
        print(f"G_in  (integrated) = {results['in_ave']:.6f} ± {results['in_std']:.6f} MW/(m²·K)")
        print(f"G_out (integrated) = {results['out_ave']:.6f} ± {results['out_std']:.6f} MW/(m²·K)")
        print(f"G_tot (integrated) = {results['tot_ave']:.6f} ± {results['tot_std']:.6f} MW/(m²·K)")
        if self.quantum:
            print(f"\nQuantum-corrected (T = {self.temperature:g} K):")
            print(f"G_in  (integrated) = {results['in_qc_ave']:.6f} ± {results['in_qc_std']:.6f} MW/(m²·K)")
            print(f"G_out (integrated) = {results['out_qc_ave']:.6f} ± {results['out_qc_std']:.6f} MW/(m²·K)")
            print(f"G_tot (integrated) = {results['tot_qc_ave']:.6f} ± {results['tot_qc_std']:.6f} MW/(m²·K)")
        print("=" * 70 + "\n")

    _txt_block = staticmethod(txt_block)

    def _export_txt(self, Reformed_NEMD_data, N_nemd_data, N_temp_group, xticks_length, Time_upper, Length):
        """Export NEMD data as a compact, Excel/Origin-friendly TAB-separated .txt file"""
        N_repeat = Reformed_NEMD_data['E_in'].shape[1] - 1
        lines = [
            "# ==== plt_nemd.py NEMD data export ====",
            f"# Directory: {self.directory}",
            f"# Effective length: {Length:.6g} nm",
            f"# Scale_eff_size: {self.scale_eff_size}",
            f"# N_repeat: {N_repeat}",
            "# Columns are TAB-separated; lines starting with '#' are comments/headers.",
            "",
        ]

        # Temperature profile at the converged state (what panel (a) plots): position vs T
        profile = np.array([Reformed_NEMD_data[f'temp_{i + 1}'][-1] for i in range(len(xticks_length))])
        lines += self._txt_block("temperature_profile", "K", "position", "nm", xticks_length, profile)

        # Raw per-timestep time series (the converged-average row above is omitted here to avoid duplication)
        time_ns = np.linspace(0, Time_upper, N_nemd_data)
        for i in range(N_temp_group):
            key = f'temp_{i}'
            lines += self._txt_block(key, "K", "t", "ns", time_ns, Reformed_NEMD_data[key][:N_nemd_data])
        lines += self._txt_block("E_in", "eV", "t", "ns", time_ns, Reformed_NEMD_data["E_in"])
        lines += self._txt_block("E_out", "eV", "t", "ns", time_ns, Reformed_NEMD_data["E_out"])

        lines.append("# ==== Summary results ====")
        lines.append("# quantity\taverage\tstd\tunit")
        lines.append(f"G\t{Reformed_NEMD_data['G'][0, -2]:.6g}\t{Reformed_NEMD_data['G'][0, -1]:.6g}\tMW/(m²·K)")
        lines.append(f"k\t{Reformed_NEMD_data['k'][0, -2]:.6g}\t{Reformed_NEMD_data['k'][0, -1]:.6g}\tW/(m·K)")

        with open(os.path.join(self.directory, 'data_nemd.txt'), 'w') as f:
            f.write("\n".join(lines) + "\n")

    def _export_shc_txt(self, Reformed_SHC_data):
        """Export NEMD's SHC data as a compact, Excel/Origin-friendly TAB-separated .txt file"""
        lines = [
            "# ==== plt_nemd.py SHC data export ====",
            f"# Directory: {self.directory}",
            f"# Cutoff frequency: {self.cutoff_freq} THz",
        ]
        if self.quantum:
            lines.append(f"# Quantum correction: T = {self.temperature:g} K ('*_qc' blocks)")
        lines += [
            "# Columns are TAB-separated; lines starting with '#' are comments/headers.",
            "",
        ]

        spectral_keys = ["k_g_wi", "k_g_wo", "k_g_wt"]
        result_keys = ["in", "out", "tot"]
        if self.quantum:
            spectral_keys += [f"{key}_qc" for key in spectral_keys]
            result_keys += [f"{key}_qc" for key in result_keys]

        lines += self._txt_block("Kt", "eV/ps", "t_corr", "ps", Reformed_SHC_data['t'][:, -2],
                                 Reformed_SHC_data['Kt'], trailing_labels=("average", "std"))
        for key in spectral_keys:
            lines += self._txt_block(key, "MW/(m²·K·THz)", "nu", "THz", Reformed_SHC_data['nu'][:, -2],
                                     Reformed_SHC_data[key], trailing_labels=("average", "std"))

        res = Reformed_SHC_data["Results"]
        lines.append("# ==== Summary results (frequency-integrated) ====")
        lines.append("# quantity\taverage\tstd\tunit")
        for key in result_keys:
            lines.append(f"G_{key}\t{res[key + '_ave']:.6g}\t{res[key + '_std']:.6g}\tMW/(m²·K)")

        with open(os.path.join(self.directory, 'data_shc.txt'), 'w') as f:
            f.write("\n".join(lines) + "\n")

    @staticmethod
    def _plot_temperature_profile(Reformed_NEMD_data, N_repeat, xticks_length, Length):
        """(a) Temperature profile along the heat transfer direction"""
        set_fig_properties([plt.gca()])
        profile = np.array([Reformed_NEMD_data[f'temp_{i + 1}'][-1]
                            for i in range(len(xticks_length))])
        plt.scatter(np.repeat(xticks_length, N_repeat), profile[:, :N_repeat].ravel(),
                    color="C0", alpha=0.3, s=30)
        plt.scatter(xticks_length, profile[:, -1], color="C1", s=60)
        plt.plot(xticks_length, profile[:, -1], '-', color="C0", alpha=0.4, linewidth=2)
        plt.text(0.95, 0.9, f'G={Reformed_NEMD_data["G"][0, -2]:.3f}±{Reformed_NEMD_data["G"][0, -1]:.2f} MW/(m$^2$·K)',
             ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(0.95, 0.8, f'κ={Reformed_NEMD_data["k"][0, -2]:.3f}±{Reformed_NEMD_data["k"][0, -1]:.3f} W/(m·K)',
             ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(0.08, 0.08, f'Effective length={Length:.2f} nm', ha='left', va='bottom',
             transform=plt.gca().transAxes)
        plt.xlabel("Length (nm)")
        plt.ylabel("Temperature (K)")
        plt.title("(a) Temperature profile")

    @staticmethod
    def _plot_thermostat_energy(Reformed_NEMD_data, Time_upper, N_nemd_data, N_repeat):
        """(b) Cumulative energy injected/extracted by the source/sink thermostats"""
        set_fig_properties([plt.gca()])
        Time = np.linspace(0, Time_upper, N_nemd_data)
        slope_in = abs(Reformed_NEMD_data["E_in"][-1, -1] - Reformed_NEMD_data["E_in"][0, -1]) / Time_upper
        slope_out = abs(Reformed_NEMD_data["E_out"][-1, -1] - Reformed_NEMD_data["E_out"][0, -1]) / Time_upper

        plt.plot(Time, Reformed_NEMD_data["E_in"][:, :N_repeat] / 1000, color="C0", alpha=0.3, lw=2)
        plt.plot(Time, Reformed_NEMD_data["E_out"][:, :N_repeat] / 1000, color="C1", alpha=0.3, lw=2)
        plt.plot(Time, Reformed_NEMD_data["E_in"][:, -1] / 1000, color="C0", lw=2, label="Source")
        plt.plot(Time, Reformed_NEMD_data["E_out"][:, -1] / 1000, color="C1", lw=2, label="Sink")

        plt.text(0.08, 0.10, f"slope$_{{in}}$={slope_in:.2f}", ha='left', va='bottom',
             transform=plt.gca().transAxes, fontsize=13, color="C0")
        plt.text(0.08, 0.88, f"slope$_{{out}}$={slope_out:.2f}", ha='left', va='top',
             transform=plt.gca().transAxes, fontsize=13, color="C1")

        plt.xlim(0, Time_upper)
        plt.xlabel("Time (ns)")
        plt.ylabel(r"Energy ($\times 10^3$ eV)")
        plt.legend(frameon=False, loc='center right')
        plt.title("(b) Thermostat energy")

    @staticmethod
    def _plot_kt_panel(Reformed_SHC_data):
        """(c) Force-virial correlation function K_tot(t)"""
        set_fig_properties([plt.gca()])
        plt.plot(Reformed_SHC_data['t'][:, -2], Reformed_SHC_data['Kt'][:, -2] / Reformed_SHC_data['L'], lw=2)
        plt.ylabel('K (eV/ps)')
        plt.xlabel('Correlation time (ps)')
        plt.title('(c) K$_{tot}$(t)')

    def _plot_shc_spectral_panel(self, Reformed_SHC_data):
        """(d) SHC spectral thermal conductance, split into in-/out-of-plane components"""
        set_fig_properties([plt.gca()])
        series = [("k_g_wi", "C1", "In-plane component", "-"),
                  ("k_g_wo", "C2", "Out-of-plane component", "-"),
                  ("k_g_wt", "C0", "Total", "-")]
        if self.quantum:
            series.append(("k_g_wt_qc", "C3", "Total (quantum-corrected)", "--"))
        for col, color, label, ls in series:
            plt.plot(Reformed_SHC_data['nu'][:, -2], Reformed_SHC_data[col][:, -2], linewidth=2, color=color, label=label,
                 linestyle=ls)
            plt.fill_between(Reformed_SHC_data['nu'][:, -2],
                             Reformed_SHC_data[col][:, -2] - Reformed_SHC_data[col][:, -1],
                             Reformed_SHC_data[col][:, -2] + Reformed_SHC_data[col][:, -1],
                             facecolor=color, alpha=0.3)
        plt.legend(frameon=False, fontsize=FS)
        plt.xlim(0, self.cutoff_freq)
        plt.ylabel(r'$g$($\omega$) (MW/(m$^2$·K·THz))')
        plt.xlabel(r'$\omega/2\pi$ (THz)')
        plt.title('(d) Spectral thermal conductance')

    def _plot_results(self, Reformed_NEMD_data, Reformed_SHC_data,
                      Time_upper, N_nemd_data, N_repeat,
                      xticks_length, Length):
        """Visualize NEMD and SHC results"""

        if self.has_shc:
            plt.figure(figsize=(10, 8))
            plt.subplot2grid((2, 5), (0, 0), colspan=3)
        else:
            plt.figure(figsize=(10, 4))
            plt.subplot(1, 2, 1)
        self._plot_temperature_profile(Reformed_NEMD_data, N_repeat, xticks_length, Length)

        if self.has_shc:
            plt.subplot2grid((2, 5), (0, 3), colspan=2)
        else:
            plt.subplot(1, 2, 2)
        self._plot_thermostat_energy(Reformed_NEMD_data, Time_upper, N_nemd_data, N_repeat)

        if self.has_shc:
            plt.subplot2grid((2, 4), (1, 0), colspan=1)
            self._plot_kt_panel(Reformed_SHC_data)

            plt.subplot2grid((2, 4), (1, 1), colspan=3)
            self._plot_shc_spectral_panel(Reformed_SHC_data)

            plt.subplots_adjust(wspace=1, hspace=0.3)
        else:
            plt.tight_layout()

        if self.save:
            plt.savefig('nemd.png', dpi=300, bbox_inches='tight')
        else:
            plt.show()


if __name__ == "__main__":

    argv = sys.argv[1:]

    if argv and argv[0] in ('-h', '--help', 'help'):
        print_usage()
        sys.exit(0)

    # --save/--save-data/--quantum are independent flags: they can appear anywhere,
    # in any combination, without disturbing the positional numeric arguments below.
    save = '--save' in argv or 'save' in argv
    save_data = '--save-data' in argv or 'save_data' in argv
    quantum = '--quantum' in argv or 'quantum' in argv
    positional = [a for a in argv if a not in ('--save', 'save', '--save-data', 'save_data', '--quantum', 'quantum')]

    try:
        if len(positional) > 0:
            if positional[0].lower() == "auto":
                real_length = None
            else:
                real_length = float(positional[0])
                if real_length <= 0:
                    raise ValueError
        else:
            real_length = None

        scale_eff_size = eval_number(positional[1]) if len(positional) > 1 else 1
        cutoff_freq = float(positional[2]) if len(positional) > 2 else 60
        if len(positional) > 3:
            raise ValueError
    except (ValueError, IndexError):
        print_usage()
        sys.exit(1)

    directory = os.getcwd()

    processor = NEMD_Processor(directory, real_length, scale_eff_size, cutoff_freq, _save=save, _save_data=save_data,
                               _quantum=quantum)
    processor.process()

"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_hnemd.py
Category:   Plot Scripts
Purpose:    Post-processing for HNEMD (Homogeneous Non-Equilibrium Molecular
            Dynamics) thermal conductivity calculations, including running
            average and optional SHC spectral analysis.
Usage:      gpumdkit.sh -plt hnemd [scale_eff_size] [cutoff_freq] [save] [quantum]
            python plt_hnemd.py [scale_eff_size] [cutoff_freq] [save] [quantum]
Arguments:
  scale_eff_size  Scale factor for effective cross-sectional area (default: 1)
  cutoff_freq     Cutoff frequency for SHC in THz (default: 60)
  save            Save the plot as 'hnemd.png' instead of displaying it
  quantum         Also apply the quantum correction to the SHC spectrum
Output:
  hnemd.png    (if save is used)
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
from scipy.integrate import cumulative_trapezoid

fs = 12
TIME_UNIT_FS = 10.18051  # Legacy GPUMD internal time unit.
trap = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
plt.rcParams.update({
    'font.size': fs, 'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'], 'axes.linewidth': 1.2,
})

def set_fig_properties(ax_list, tl=4, tw=1.2, tlm=4):
    """Set axis tick properties."""
    for ax in ax_list:
        ax.tick_params(which='both', length=tl, width=tw, direction='in', right=True, top=True)
        ax.tick_params(which='minor', length=tlm)

def quantum_correction_factor(nu, T):
    """Quantum/classical mode heat capacity; nu in THz and T in K."""
    h_over_kB = 6.62607015e-34 / 1.380649e-23  # K*s
    x = h_over_kB * np.asarray(nu, dtype=float) * 1e12 / T
    with np.errstate(divide='ignore', invalid='ignore'):
        factor = x**2 * np.exp(-x) / np.expm1(-x)**2
    return np.where(x == 0, 1.0, factor)

@lru_cache(maxsize=None)
def element_from_mass(mass):
    """Find the element with the closest standard atomic mass (amu)."""
    from ase.data import atomic_masses, chemical_symbols

    return chemical_symbols[int(np.argmin(np.abs(atomic_masses[1:] - mass))) + 1]

def xyz_in_to_model(source, target):
    """Convert legacy xyz.in, preserving masses, groups and velocities."""
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

_ARITH_OPS = {
    ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul, ast.Div: op.truediv,
    ast.FloorDiv: op.floordiv, ast.Mod: op.mod, ast.Pow: op.pow,
    ast.UAdd: op.pos, ast.USub: op.neg,
}

def eval_number(expr):
    """Evaluate numeric literals and basic arithmetic without executing code."""
    def _eval(node):
        if isinstance(node, ast.Expression):
            return _eval(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return node.value
        if isinstance(node, ast.BinOp) and type(node.op) in _ARITH_OPS:
            left, right = _eval(node.left), _eval(node.right)
            if isinstance(node.op, ast.Pow) and abs(right) > 20:
                # Avoid enormous integer powers.
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

HNEMD_KEYS = ('kx_in', 'kx_out', 'kx_tot', 'ky_in', 'ky_out', 'ky_tot', 'kz_tot')

def print_usage():
    """Print usage instructions"""
    print("Usage: gpumdkit.sh -plt hnemd [scale_eff_size] [cutoff_freq] [--save] [--save-data] [--quantum]")
    print("Params:")
    print("  scale_eff_size: Optional, Scale factor for effective cross-sectional area (default: 1)")
    print("                   • For 3D bulk systems: use 1")
    print("                   • For low-dimensional systems with vacuum layer: S_box / S_eff")
    print("                     - S_box: box area perpendicular to heat transfer direction")
    print("                     - S_eff: real or effective area of the system")
    print("                   • May be a plain number or a simple arithmetic expression,")
    print("                     e.g. '100/2' (evaluates to 50)")
    print("  cutoff_freq   : Optional, Cutoff frequency for SHC calculation in THz (default: 60)")
    print("  --save        : Optional, save the plot as 'hnemd.png' (default: show it interactively)")
    print("  --save-data   : Optional, also save the processed data as 'data_hnemd.npz'/'data_shc.npz'")
    print("                  and their .txt equivalents (independent of --save; can be used on its own)")
    print("  --quantum     : Optional, also apply the quantum correction to the SHC spectral thermal")
    print("                  conductivity (requires shc.out) ")
    print("  Legacy bare 'save', 'save_data', and 'quantum' are also accepted.")
    print("  --save/--save-data/--quantum may appear anywhere on the command line, e.g.:")
    print("    gpumdkit.sh -plt hnemd --save-data              (defaults for everything else)")
    print("    gpumdkit.sh -plt hnemd 1 60 --save --save-data --quantum")

class HNEMD_Processor:
    """Process HNEMD running averages and optional SHC spectra."""

    def __init__(self, _directory, _scale_eff_size=1, _cutoff_freq=60, _save=False, _save_data=False, _quantum=False):
        self.directory = _directory
        self.scale_eff_size = _scale_eff_size
        self.cutoff_freq = _cutoff_freq
        self.save = _save
        self.save_data = _save_data
        self.quantum = _quantum
        self.path = {
            key: os.path.join(self.directory, name)
            for key, name in (('run', 'run.in'), ('kappa', 'kappa.out'), ('shc', 'shc.out'),
                              ('model', 'model.xyz'), ('xyz_in', 'xyz.in'))
        }
        self.has_shc = os.path.exists(self.path['shc'])

    def prepare_model(self):
        """Convert legacy xyz.in when model.xyz is absent."""
        if os.path.exists(self.path['model']) or not os.path.exists(self.path['xyz_in']):
            return
        print("[INFO] model.xyz not found, converting xyz.in -> model.xyz ...")
        n, type_elements = xyz_in_to_model(self.path['xyz_in'], self.path['model'])
        for atom_type in sorted(type_elements):
            mapping = ", ".join(f"{element} (mass {mass:g})" for element, mass in sorted(type_elements[atom_type]))
            print(f"       type {atom_type} -> {mapping}")
        print(f"[INFO] Wrote {self.path['model']}: {n} atoms.\n")

    def process_SHC(self, Fe):
        """Compute classical and optional quantum-corrected SHC spectra."""
        commands = read_commands(self.path['run'])
        command = last_command(commands, 'compute_shc')
        thermostat = next(cmd for cmd in reversed(commands)
                          if len(cmd) > 2 and cmd[0] == 'ensemble' and cmd[1].startswith('nvt_'))
        temperature = int(thermostat[2])
        n_corr = 2 * int(command[2]) - 1
        n_omega = int(command[4])
        raw = load_repeated(self.path['shc'], n_corr + n_omega, ('grid', 'in', 'out'))
        data = {'t': raw['grid'][:n_corr], 'omega': raw['grid'][n_corr:],
                'nu': raw['grid'][n_corr:] / (2 * np.pi),
                'Ki': raw['in'][:n_corr], 'jwi': raw['in'][n_corr:],
                'Ko': raw['out'][:n_corr], 'jwo': raw['out'][n_corr:]}
        n_repeat = data['t'].shape[1]

        from ase.io import read

        model = read(self.path['model'])
        fraction = 1.0
        if 'group' in command:
            index = command.index('group')
            method, group_id = map(int, command[index + 1:index + 3])
            if group_id < 0:
                raise ValueError('SHC analysis requires a single group_id; group -1 outputs multiple groups.')
            labels = model.get_array('group')
            labels = labels if labels.ndim == 1 else labels[:, method]
            fraction = np.count_nonzero(labels == group_id) / len(labels)
            if not fraction:
                raise ValueError(f'SHC group {group_id} contains no atoms.')
        volume = model.get_volume() * fraction / self.scale_eff_size
        data['L'] = model.cell[int(command[3]), int(command[3])] * fraction
        convert = 1.602176634e3  # GPUMD spectral current -> W/(m·K·THz).
        denominator = Fe * temperature * volume
        data['k_g_wi'] = data['jwi'] * convert / denominator
        data['k_g_wo'] = data['jwo'] * convert / denominator
        data['k_g_wt'] = data['k_g_wi'] + data['k_g_wo']
        data['k_g_wt'][data['nu'] > float(self.cutoff_freq)] = 0
        data['k_g_wt'][data['k_g_wt'] < 0] = 0.001
        data['Kt'] = data['Ki'] + data['Ko']

        result_cols = [('in', 'k_g_wi'), ('out', 'k_g_wo'), ('tot', 'k_g_wt')]
        if self.quantum:
            self.temperature = temperature
            factor = quantum_correction_factor(data['nu'], temperature)
            for key, col in list(result_cols):
                data[f'{col}_qc'] = data[col] * factor
                result_cols.append((f'{key}_qc', f'{col}_qc'))

        data = {key: value if key == 'L' else add_statistics(value, sem=True)
                for key, value in data.items()}
        results = {}
        for key, col in result_cols:
            # Preserve the original constant-frequency-spacing integration.
            values = trap(data[col][:, :n_repeat], dx=data['nu'][0, 0], axis=0)
            results[f'{key}_ave'] = values.mean()
            results[f'{key}_std'] = values.std() / np.sqrt(n_repeat)
        data['Results'] = results
        if self.save_data:
            np.savez(os.path.join(self.directory, 'data_shc.npz'), **data)
            self._export_shc_txt(data)
        return data

    def process(self):
        """Load, average, export and plot HNEMD data."""
        self.prepare_model()
        commands = read_commands(self.path['run'])
        time_step = float(last_command(commands, 'time_step')[1])
        hnemd_index = max(i for i, cmd in enumerate(commands) if cmd[0] == 'compute_hnemd')
        command = commands[hnemd_index]
        output_interval = int(command[1])
        forces = [float(value) for value in command[2:5]]
        direction, force = next((axis, value) for axis, value in reversed(list(zip('xyz', forces)))
                                if value != 0)
        run_steps = int(last_command(commands[hnemd_index + 1:], 'run')[1])
        n_rows = run_steps // output_interval
        data = load_repeated(self.path['kappa'], n_rows,
                             ('kx_in', 'kx_out', 'ky_in', 'ky_out', 'kz_tot'))
        n_repeat = data['kx_in'].shape[1]
        time_upper = run_steps * time_step * 1e-6  # ns
        t = np.arange(1, n_rows + 1) * 1e-3  # Keep the established running-average time grid (ns).
        data = {key: cumulative_trapezoid(add_statistics(value) * self.scale_eff_size,
                                         t, axis=0, initial=0) / t[:, None]
                for key, value in data.items()}
        for prefix in ('kx', 'ky'):
            data[f'{prefix}_tot'] = data[f'{prefix}_in'] + data[f'{prefix}_out']

        results = {}
        for key in HNEMD_KEYS:
            values = data[key][-1, :-1]
            results[f'{key}_ave'] = values.mean()
            results[f'{key}_std'] = values.std() / np.sqrt(n_repeat)
        data['Results'] = results
        if self.save_data:
            np.savez(os.path.join(self.directory, 'data_hnemd.npz'), **data)
            self._export_txt(data, t)
        self._print_hnemd_results(results, direction)

        if self.has_shc:
            print("\n[INFO] SHC data detected, processing spectral heat current...")
            shc = self.process_SHC(Fe=force)
            self._print_shc_results(shc['Results'])
        else:
            print("\n[INFO] No SHC data found (shc.out not present), skipping SHC analysis.")
            if self.quantum:
                print("[WARNING] --quantum needs shc.out, quantum correction skipped.")
            shc = None
        self._plot_results(data, shc, results, direction, time_upper, n_repeat, t)

    def _print_hnemd_results(self, results, direction):
        """Print HNEMD thermal conductivity results"""
        print("\n" + "=" * 70)
        print("HNEMD Thermal Conductivity Results")
        print("=" * 70)
        print(f"\nDirection: {direction.lower()}")
        print(f"Scale_eff_size: {self.scale_eff_size}")

        if direction in ['x', 'y']:
            key_prefix = f'k{direction}'
            print(f"\nκ_in  = {results[key_prefix + '_in_ave']:.4f} ± {results[key_prefix + '_in_std']:.4f} W/(m·K)")
            print(f"κ_out = {results[key_prefix + '_out_ave']:.4f} ± {results[key_prefix + '_out_std']:.4f} W/(m·K)")
            print(f"κ_tot = {results[key_prefix + '_tot_ave']:.4f} ± {results[key_prefix + '_tot_std']:.4f} W/(m·K)")
        elif direction == 'z':
            print(f"\nκ = {results['kz_tot_ave']:.4f} ± {results['kz_tot_std']:.4f} W/(m·K)")

        print("=" * 70)

    def _print_shc_results(self, results):
        """Print SHC spectral thermal conductivity results"""
        print("\n" + "=" * 70)
        print("SHC Spectral Thermal Conductivity Results")
        print("=" * 70)
        print(f"\nCutoff frequency: {self.cutoff_freq} THz\n")
        print(f"κ_in  (integrated) = {results['in_ave']:.4f} ± {results['in_std']:.4f} W/(m·K)")
        print(f"κ_out (integrated) = {results['out_ave']:.4f} ± {results['out_std']:.4f} W/(m·K)")
        print(f"κ_tot (integrated) = {results['tot_ave']:.4f} ± {results['tot_std']:.4f} W/(m·K)")
        if self.quantum:
            print(f"\nQuantum-corrected (T = {self.temperature} K):")
            print(f"κ_in  (integrated) = {results['in_qc_ave']:.4f} ± {results['in_qc_std']:.4f} W/(m·K)")
            print(f"κ_out (integrated) = {results['out_qc_ave']:.4f} ± {results['out_qc_std']:.4f} W/(m·K)")
            print(f"κ_tot (integrated) = {results['tot_qc_ave']:.4f} ± {results['tot_qc_std']:.4f} W/(m·K)")
        print("=" * 70 + "\n")

    @staticmethod
    def _txt_block(name, unit, x_name, x_unit, x, arr, trailing_labels=("average",)):
        """Format a data block with the existing tab-separated export layout."""
        labels = [f"repeat_{i + 1}" for i in range(arr.shape[1] - len(trailing_labels))]
        header = [f"{x_name}({x_unit})", *labels, *trailing_labels]
        output = StringIO()
        np.savetxt(output, np.column_stack((x, arr)), delimiter="\t", fmt="%.6g")
        return [f"# ---- {name} ({unit}) ----", "# " + "\t".join(header),
                *output.getvalue().splitlines(), ""]

    def _export_txt(self, Reformed_HNEMD_data, t):
        """Export HNEMD data as tab-separated text."""
        res = Reformed_HNEMD_data['Results']
        lines = [
            "# ==== plt_hnemd.py HNEMD data export ====",
            f"# Directory: {self.directory}",
            f"# Scale_eff_size: {self.scale_eff_size}",
            f"# N_repeat: {Reformed_HNEMD_data['kz_tot'].shape[1] - 1}",
            "# Columns are TAB-separated; lines starting with '#' are comments/headers.",
            "",
        ]

        for key in HNEMD_KEYS:
            lines += self._txt_block(key, "W/(m·K)", "t", "ns", t, Reformed_HNEMD_data[key])

        lines.append("# ==== Summary results ====")
        lines.append("# quantity\taverage\tstd\tunit")
        for key in HNEMD_KEYS:
            lines.append(f"{key}\t{res[key + '_ave']:.6g}\t{res[key + '_std']:.6g}\tW/(m·K)")

        with open(os.path.join(self.directory, 'data_hnemd.txt'), 'w', encoding='utf-8') as f:
            f.write("\n".join(lines) + "\n")

    def _export_shc_txt(self, Reformed_SHC_data):
        """Export SHC data as tab-separated text."""
        lines = [
            "# ==== plt_hnemd.py SHC data export ====",
            f"# Directory: {self.directory}",
            f"# Cutoff frequency: {self.cutoff_freq} THz",
        ]
        if self.quantum:
            lines.append(f"# Quantum correction: T = {self.temperature} K ('*_qc' blocks)")
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
            lines += self._txt_block(key, "W/(m·K·THz)", "nu", "THz", Reformed_SHC_data['nu'][:, -2],
                                     Reformed_SHC_data[key], trailing_labels=("average", "std"))

        res = Reformed_SHC_data["Results"]
        lines.append("# ==== Summary results (frequency-integrated) ====")
        lines.append("# quantity\taverage\tstd\tunit")
        for key in result_keys:
            lines.append(f"kappa_{key}\t{res[key + '_ave']:.6g}\t{res[key + '_std']:.6g}\tW/(m·K)")

        with open(os.path.join(self.directory, 'data_shc.txt'), 'w', encoding='utf-8') as f:
            f.write("\n".join(lines) + "\n")

    @staticmethod
    def _plot_running_kappa_panel(t, N_repeat, Time_upper, curves, annotations, direction):
        """Plot repeated runs, their average and result labels."""
        ax = plt.gca()
        set_fig_properties([ax])
        for data, _, _ in curves:
            ax.plot(t, data[:, :N_repeat], color='k', alpha=0.3)
        for data, color, lw in curves:
            ax.plot(t, data[:, -1], color=color, lw=lw)
        for y, label, color in annotations:
            kwargs = {'color': color} if color is not None else {}
            ax.text(0.95, y, label, ha='right', va='top', transform=ax.transAxes, **kwargs)
        ax.set(xlim=(0, Time_upper), xlabel='time (ns)', ylabel=r'$\kappa$ (W/(m·K))',
               title=f'(a) Running average thermal conductivity: along {direction}')

    @staticmethod
    def _plot_kt_panel(data):
        """Plot the force-virial correlation function."""
        ax = plt.gca()
        set_fig_properties([ax])
        ax.plot(data['t'][:, -2], data['Kt'][:, -2] / data['L'], lw=2)
        ax.set(ylabel='K (eV/ps)', xlabel='Correlation time (ps)', title='(b) K$_{tot}$(t)')

    def _plot_shc_spectral_panel(self, data, series):
        """Plot spectral conductivities and standard-error bands."""
        ax = plt.gca()
        set_fig_properties([ax])
        nu = data['nu'][:, -2]
        for col, color, lw, label, *linestyle in series:
            average, error = data[col][:, -2:].T
            ax.plot(nu, average, linewidth=lw, color=color, label=label,
                    linestyle=linestyle[0] if linestyle else '-')
            ax.fill_between(nu, average - error, average + error, facecolor=color, alpha=0.3)
        if any(item[3] for item in series):
            ax.legend(frameon=False, fontsize=fs)
        ax.axhline(y=0, color='k', linestyle='--')
        ax.set(xlim=(0, self.cutoff_freq), ylabel=r'$\kappa$($\omega$) (W/(m·K·THz))',
               xlabel=r'$\nu$ (THz)', title='(c) Spectral thermal conductivity')

    def _plot_results(self, data, shc, results, direction, time_upper, n_repeat, t):
        """Plot HNEMD and available SHC results."""
        plt.figure(figsize=(10, 8) if self.has_shc else (8, 4))
        plt.subplot(2 if self.has_shc else 1, 1, 1)
        if direction in ('x', 'y'):
            prefix = f'k{direction}'
            components = (('in', 'C1'), ('out', 'C2'), ('tot', 'C0'))
            curves = [(data[f'{prefix}_{part}'], color, 3) for part, color in components]
            annotations = [
                (y, f"$\\kappa_{{\\mathrm{{{part}}}}}$ = {results[f'{prefix}_{part}_ave']:.3f} ± "
                 f"{results[f'{prefix}_{part}_std']:.2f} W/(m·K)", color)
                for y, (part, color) in zip((0.93, 0.83, 0.73), components)
            ]
            series = [('k_g_wi', 'C1', 2, 'In-plane component'),
                      ('k_g_wo', 'C2', 2, 'Out-of-plane component'),
                      ('k_g_wt', 'C0', 2, 'Total')]
            if self.quantum:
                series.append(('k_g_wt_qc', 'C3', 2, 'Total (quantum-corrected)', '--'))
        else:
            curves = [(data['kz_tot'], 'C0', 5)]
            annotations = [(0.9, fr"$\kappa_{{tot}}$ = {results['kz_tot_ave']:.3f} ± "
                            f"{results['kz_tot_std']:.2f} W/(m·K)", None)]
            series = [('k_g_wt', 'C0', 3, 'Classical' if self.quantum else None)]
            if self.quantum:
                series.append(('k_g_wt_qc', 'C3', 3, 'Quantum-corrected', '--'))
        self._plot_running_kappa_panel(t, n_repeat, time_upper, curves, annotations, direction)
        if self.has_shc:
            plt.subplot(2, 4, 5)
            self._plot_kt_panel(shc)
            plt.subplot2grid((2, 4), (1, 1), colspan=3)
            self._plot_shc_spectral_panel(shc, series)
        plt.tight_layout()
        if self.save:
            plt.savefig('hnemd.png', dpi=300, bbox_inches='tight')
        else:
            plt.show()


if __name__ == "__main__":

    argv = sys.argv[1:]

    if argv and argv[0] in ('-h', '--help', 'help'):
        print_usage()
        sys.exit(0)

    save = '--save' in argv or 'save' in argv
    save_data = '--save-data' in argv or 'save_data' in argv
    quantum = '--quantum' in argv or 'quantum' in argv
    positional = [a for a in argv if a not in ('--save', 'save', '--save-data', 'save_data', '--quantum', 'quantum')]

    try:
        scale_eff_size = eval_number(positional[0]) if len(positional) > 0 else 1
        cutoff_freq = float(positional[1]) if len(positional) > 1 else 60
        if len(positional) > 2:
            raise ValueError

    except (ValueError, IndexError):
        print_usage()
        sys.exit(1)

    directory = os.getcwd()

    processor = HNEMD_Processor(directory, scale_eff_size, cutoff_freq, _save=save, _save_data=save_data,
                                _quantum=quantum)
    processor.process()

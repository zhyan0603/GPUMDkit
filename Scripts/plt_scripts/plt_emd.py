"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_emd.py
Category:   Plot Scripts
Purpose:    Post-processing for EMD (Equilibrium Molecular Dynamics) thermal
            conductivity calculations, including HAC, running thermal
            conductivity, and results printing.
Usage:      gpumdkit.sh -plt emd <direction> [save]
            python plt_emd.py <x|y|z> [save]
Arguments:
  direction  Heat transfer direction: x, y, or z
  save       Save the plot as 'emd.png' instead of displaying it
Output:
  emd.png    (if save is used)
Author:     Xin Wu (xinwuchn97@gmial.com)
Last-modified: 2026-10-21
=============================================================================
"""

import os
import sys
from functools import lru_cache
from io import StringIO
from math import isfinite

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.size': 12, 'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'],
    'axes.linewidth': 1.2,
})
TIME_UNIT_FS = 10.18051  # Legacy GPUMD internal time unit.


def set_fig_properties(ax_list, tl=4, tw=1.2, tlm=4):
    """Set figure properties for axes"""
    for ax in ax_list:
        ax.tick_params(which='both', length=tl, width=tw, direction='in', right=True, top=True)
        ax.tick_params(which='minor', length=tlm)


@lru_cache(maxsize=None)
def element_from_mass(mass):
    """Return the element symbol whose standard atomic mass is closest to mass (amu)."""
    from ase.data import atomic_masses, chemical_symbols

    return chemical_symbols[int(np.argmin(np.abs(atomic_masses[1:] - mass))) + 1]


def xyz_in_to_model(source, target):
    """Convert legacy xyz.in, preserving masses, velocities and groups."""
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


def print_usage():
    print("""Usage: gpumdkit.sh -plt emd direction [--save] [--save-data]
Params:
  direction    : Heat transfer direction in lowercase (x, y, or z)
  --save       : Save the plot as 'emd.png' instead of displaying it
  --save-data  : Save 'data_emd.npz' and 'data_emd.txt' independently of --save
Flags may appear anywhere; 'save' and 'save_data' aliases are also accepted.
Examples:
  gpumdkit.sh -plt emd z --save-data
  gpumdkit.sh -plt emd z --save --save-data""")


class EMD_Processor:
    COLUMNS = ('time', 'jx_in', 'jx_out', 'jy_in', 'jy_out', 'jz_tot',
               'kx_in', 'kx_out', 'ky_in', 'ky_out', 'kz_tot')
    KAPPA_KEYS = ('kx_in', 'kx_out', 'kx_tot', 'ky_in', 'ky_out', 'ky_tot', 'kz_tot')

    def __init__(self, _directory, _direction='z', _save=False, _save_data=False):
        self.directory = _directory
        self.direction = _direction
        self.save = _save
        self.save_data = _save_data
        self.path = {key: os.path.join(_directory, filename) for key, filename in (
            ('run', 'run.in'), ('hac', 'hac.out'), ('model', 'model.xyz'), ('xyz_in', 'xyz.in'))}

    def prepare_model(self):
        if os.path.exists(self.path['model']) or not os.path.exists(self.path['xyz_in']):
            return
        print("[INFO] model.xyz not found, converting xyz.in -> model.xyz ...")
        n, types = xyz_in_to_model(self.path['xyz_in'], self.path['model'])
        for atom_type, elements in sorted(types.items()):
            mapping = ", ".join(f"{element} (mass {mass:g})" for element, mass in sorted(elements))
            print(f"       type {atom_type} -> {mapping}")
        print(f"[INFO] Wrote {self.path['model']}: {n} atoms.\n")

    def process(self):
        self.prepare_model()
        time_step = hac = None
        with open(self.path['run'], encoding='utf-8-sig') as source:
            for line in source:
                command = line.split('#', 1)[0].split()
                if not command:
                    continue
                if command[0] == 'time_step':
                    time_step = float(command[1])
                elif command[0] == 'compute_hac':
                    hac = list(map(int, command[1:4]))
        if time_step is None or hac is None:
            raise ValueError('run.in must define time_step and compute_hac.')
        sample_interval, correlation_steps, output_interval = hac
        n_rows = correlation_steps // output_interval
        raw = np.loadtxt(self.path['hac'], ndmin=2)
        if raw.shape[1] != len(self.COLUMNS):
            raise ValueError(f'hac.out must have {len(self.COLUMNS)} columns.')
        if n_rows < 1 or not len(raw) or len(raw) % n_rows:
            raise ValueError('Incomplete HAC data: check correlation steps and completed runs.')
        n_repeat = len(raw) // n_rows
        runs = raw.reshape(n_repeat, n_rows, len(self.COLUMNS)).transpose(1, 0, 2)
        data = {key: np.column_stack((values, values.mean(axis=1)))
                for key, values in zip(self.COLUMNS, runs.transpose(2, 0, 1))}
        for axis in 'xy':
            data[f'k{axis}_tot'] = data[f'k{axis}_in'] + data[f'k{axis}_out']
        results = {}
        for key in self.KAPPA_KEYS:
            # Preserve the existing estimate from the second half of the mean curve.
            values = data[key][n_rows // 2:, -1]
            results[key + '_ave'] = values.mean()
            results[key + '_std'] = values.std() / np.sqrt(n_repeat)
        data['Results'] = results
        if self.save_data:
            np.savez(os.path.join(self.directory, 'data_emd.npz'), **data)
            self._export_txt(data)
        self._print_results(results)
        time_upper = sample_interval * correlation_steps * time_step * 1e-6
        self._plot_results(data, time_upper, n_repeat)

    @staticmethod
    def _txt_block(name, unit, x_name, x_unit, x, arr, trailing_labels=('average',)):
        labels = [f'repeat_{i + 1}' for i in range(arr.shape[1] - len(trailing_labels))]
        header = [f'{x_name}({x_unit})', *labels, *trailing_labels]
        output = StringIO()
        np.savetxt(output, np.column_stack((x, arr)), delimiter='\t', fmt='%.6g')
        return [f'# ---- {name} ({unit}) ----', '# ' + '\t'.join(header),
                *output.getvalue().splitlines(), '']

    def _export_txt(self, data):
        lines = [
            '# ==== plt_emd.py EMD data export ====', f'# Directory: {self.directory}',
            f'# Direction: {self.direction}', f"# N_repeat: {data['jz_tot'].shape[1] - 1}",
            "# Columns are TAB-separated; lines starting with '#' are comments/headers.", '',
        ]
        for axis in 'xyz':
            components = ('in', 'out') if axis in 'xy' else ('tot',)
            keys = [f'j{axis}_{component}' for component in components]
            keys += [f'k{axis}_{component}' for component in components]
            if axis in 'xy':
                keys.append(f'k{axis}_tot')
            for key in keys:
                unit = 'a.u.' if key.startswith('j') else 'W/(m·K)'
                lines += self._txt_block(key, unit, 't_corr', 'ns', data['time'][:, -1] * 1e-3, data[key])
        lines += ['# ==== Summary results ====', '# quantity\taverage\tstd\tunit']
        for key in self.KAPPA_KEYS:
            lines.append(f"{key}\t{data['Results'][key + '_ave']:.6g}\t"
                         f"{data['Results'][key + '_std']:.6g}\tW/(m·K)")
        with open(os.path.join(self.directory, 'data_emd.txt'), 'w') as output:
            output.write('\n'.join(lines) + '\n')

    def _print_results(self, results):
        print('\n' + '=' * 60 + '\nEMD Thermal Conductivity Results\n' + '=' * 60)
        print(f'\nDirection: {self.direction}\n' if self.direction in 'xy' else 'Direction: z')
        components = ('in', 'out', 'tot') if self.direction in 'xy' else ('tot',)
        for component in components:
            key = f'k{self.direction}_{component}'
            label = f'κ_{component:<3}' if self.direction in 'xy' else 'κ'
            print(f"{label} = {results[key + '_ave']:.4f} ± {results[key + '_std']:.4f} W/(m·K)")
        print('=' * 60 + '\n')

    @staticmethod
    def _plot_hac_panel(time_data, series, N_repeat, Time_upper):
        set_fig_properties([plt.gca()])
        for data, color in series:
            normalized = data / data.max(axis=0)
            plt.loglog(time_data[:, :N_repeat], normalized[:, :N_repeat], color='k', alpha=0.3)
        for data, color in series:
            plt.loglog(time_data[:, -1], data[:, -1] / data[:, -1].max(), color=color)
        plt.xlim(1e-5, Time_upper)
        plt.xlabel('Correlation Time (ns)')
        plt.ylabel('Normalized HAC')
        plt.title('(a)')

    @staticmethod
    def _plot_kappa_panel(time_data, kappa_data, N_repeat, Time_upper, ave, std, ylabel_text, title_text):
        set_fig_properties([plt.gca()])
        plt.plot(time_data[:, :N_repeat], kappa_data[:, :N_repeat], color='k', alpha=0.3)
        plt.plot(time_data[:, -1], kappa_data[:, -1], color='C1', lw=3)
        plt.axhline(ave, color='C0', linestyle='--')
        plt.fill_between(time_data[:, -1], ave - std, ave + std, color='C0', alpha=0.2)
        plt.xlim(0, Time_upper)
        plt.xlabel('Correlation Time (ns)')
        plt.ylabel(ylabel_text)
        plt.title(title_text)

    def _plot_results(self, data, Time_upper, N_repeat):
        time_data = data['time'] * 1e-3
        in_plane = self.direction in ('x', 'y')
        components = ('in', 'out', 'tot') if in_plane else ('tot',)
        hac_series = [(data[f'j{self.direction}_{component}'], color)
                      for component, color in zip(components[:2] if in_plane else components, ('C0', 'C1'))]
        if not in_plane:
            hac_series[0] = (hac_series[0][0], 'C1')
        rows = 2 if in_plane else 1
        plt.figure(figsize=(10, 8 if in_plane else 4))
        plt.subplot(rows, 2, 1)
        self._plot_hac_panel(time_data, hac_series, N_repeat, Time_upper)
        for panel, component in enumerate(components, start=2):
            plt.subplot(rows, 2, panel)
            key = f'k{self.direction}_{component}'
            ave, std = data['Results'][key + '_ave'], data['Results'][key + '_std']
            symbol = rf'\kappa_{{\mathrm{{{component}}}}}' if in_plane else r'\kappa'
            title = f'${symbol}$ = {ave:.2f} ± {std:.2f} W/(m·K)'
            self._plot_kappa_panel(time_data, data[key], N_repeat, Time_upper, ave, std,
                                   f'${symbol}$ (W/(m·K))', title if in_plane else '(b) ' + title)
        plt.tight_layout()
        if self.save:
            plt.savefig('emd.png', dpi=300, bbox_inches='tight')
        else:
            plt.show()


if __name__ == '__main__':
    argv = sys.argv[1:]
    if not argv or argv[0] in ('-h', '--help', 'help'):
        print_usage()
        sys.exit(0)
    save = any(flag in argv for flag in ('--save', 'save'))
    save_data = any(flag in argv for flag in ('--save-data', 'save_data'))
    positional = [arg for arg in argv if arg not in ('--save', 'save', '--save-data', 'save_data')]
    if len(positional) != 1 or positional[0] not in ('x', 'y', 'z'):
        print_usage()
        sys.exit(1)
    EMD_Processor(os.getcwd(), positional[0], save, save_data).process()

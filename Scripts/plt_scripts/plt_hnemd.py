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
Usage:      gpumdkit.sh -plt hnemd [scale_eff_size] [cutoff_freq] [save]
            python plt_hnemd.py [scale_eff_size] [cutoff_freq] [save]
Arguments:
  scale_eff_size  Scale factor for effective cross-sectional area (default: 1)
  cutoff_freq     Cutoff frequency for SHC in THz (default: 60)
  save            Save the plot as 'hnemd.png' instead of displaying it
Output:
  hnemd.png    (if save is used)
Author:     Xin Wu (xinwuchn97@gmial.com)
Last-modified: 2026-05-16
=============================================================================
"""

import pandas as pd
from pylab import *
import numpy as np
import os
import sys
from scipy.integrate import cumulative_trapezoid
from ase.io import read

# Figure Properties
aw, fs = 1.2, 12
matplotlib.rc('font', size=fs)
matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', **{'sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans']})
matplotlib.rc('axes', linewidth=aw)

def set_fig_properties(ax_list, tl=4, tw=1.2, tlm=4):
    """Set figure properties for axes"""
    for ax in ax_list:
        ax.tick_params(which='both', length=tl, width=tw, direction='in', right=True, top=True)
        ax.tick_params(which='minor', length=tlm)

trap = np.trapezoid if hasattr(np, "trapezoid") else getattr(np, "trapz")

def print_usage():
    """Print usage instructions"""
    print("Usage: gpumdkit -plt hnemd [scale_eff_size] [cutoff_freq] [--save] [--save-data]")
    print("Params:")
    print("  scale_eff_size: Optional, Scale factor for effective cross-sectional area (default: 1)")
    print("                   • For 3D bulk systems: use 1")
    print("                   • For low-dimensional systems with vacuum layer: S_box / S_eff")
    print("                     - S_box: box area perpendicular to heat transfer direction")
    print("                     - S_eff: real or effective area of the system")
    print("  cutoff_freq   : Optional, Cutoff frequency for SHC calculation in THz (default: 60)")
    print("  --save        : Optional, save the plot as 'hnemd.png' (default: show it interactively)")
    print("  --save-data   : Optional, also save the processed data as 'data_hnemd.npz'/'data_shc.npz'")
    print("                  and their .txt equivalents (independent of --save; can be used on its own)")
    print("  --save/--save-data may appear anywhere on the command line, e.g.:")
    print("    gpumdkit -plt hnemd --save-data              (defaults for everything else)")
    print("    gpumdkit -plt hnemd 1 60 --save --save-data")

class HNEMD_Processor:
    def __init__(self, _directory, _scale_eff_size=1, _cutoff_freq=60, _save=False, _save_data=False):
        """
        Initialize HNEMD processor

        Parameters:
        -----------
        directory : str
        scale_eff_size : float
        cutoff_freq : float
        save : bool
            Save the plot as 'hnemd.png' instead of showing it interactively
        save_data : bool
            Also save the processed data as 'data_hnemd.npz'/'data_shc.npz' (and their .txt equivalents)
        """
        self.directory = _directory
        self.scale_eff_size = _scale_eff_size
        self.cutoff_freq = _cutoff_freq
        self.save = _save
        self.save_data = _save_data
        self.path = {
            'run': os.path.join(self.directory, 'run.in'),
            'kappa': os.path.join(self.directory, 'kappa.out'),
            'shc': os.path.join(self.directory, 'shc.out'),
            'model': os.path.join(self.directory, 'model.xyz')
        }
        self.has_shc = os.path.exists(self.path['shc'])

    def process_SHC(self, Fe):
        """
        Process spectral heat current (SHC) data

        Parameters:
        -----------
        Fe : float
            External driving force parameter

        Returns:
        --------
        dict : Reformed SHC data including spectral thermal conductivity
        """
        Reformed_SHC_data = {}

        col_shc_name = ['t_omega', 'Ki_jwi', 'Ko_jwo']
        raw_shc_data = np.loadtxt(self.path['shc'])
        raw_shc_data = pd.DataFrame(raw_shc_data, columns=col_shc_name)

        with open(self.path['run'], 'r') as file:
            for line in file:
                if line.lstrip().startswith('#'):
                    continue
                if 'compute_shc' in line:
                    Max_cor_step = int(line.split()[2])
                    N_omega = int(line.split()[4])
                    direction = int(line.split()[3])
                    if 'group' in line:
                        grouping_th = int(line.split()[7])
                        group_shc_th = int(line.split()[8])
                        part_ratio = -1
                    else:
                        part_ratio = 1
                if 'nvt_' in line:
                    Temp = int(line.split()[2])

        N_shc_data = 2 * Max_cor_step - 1 + N_omega
        N_repeat = len(raw_shc_data) // N_shc_data

        for col in raw_shc_data.columns:
            Reformed_SHC_data[col] = raw_shc_data[col].values.reshape(N_shc_data, N_repeat, order='F')

        Reformed_SHC_data = {
            "t": Reformed_SHC_data["t_omega"][:2 * Max_cor_step - 1, :],
            "omega": Reformed_SHC_data["t_omega"][2 * Max_cor_step - 1:, :],
            "nu": Reformed_SHC_data["t_omega"][2 * Max_cor_step - 1:, :] / (2 * np.pi),
            "Ki": Reformed_SHC_data["Ki_jwi"][:2 * Max_cor_step - 1, :],
            "jwi": Reformed_SHC_data["Ki_jwi"][2 * Max_cor_step - 1:, :],
            "Ko": Reformed_SHC_data["Ko_jwo"][:2 * Max_cor_step - 1, :],
            "jwo": Reformed_SHC_data["Ko_jwo"][2 * Max_cor_step - 1:, :]
        }

        # Calculate k(omega) from jw
        model = read(self.path['model'])
        if part_ratio == -1:
            group_arr = model.get_array('group')
            if group_arr.ndim == 1:
                group = group_arr
            else:
                group = group_arr[:, grouping_th]
            part_ratio = np.sum(group == group_shc_th) / group.size

        vol = model.get_volume() * part_ratio / self.scale_eff_size
        convert = 1.602176634e3  # ev*A/ps/THz * 1/A^3 *1/K * A ==> W/(m·K·THz)
        denom = Fe * Temp * vol

        Reformed_SHC_data["k_g_wi"] = Reformed_SHC_data["jwi"] * convert / denom
        Reformed_SHC_data["k_g_wo"] = Reformed_SHC_data["jwo"] * convert / denom
        Reformed_SHC_data["k_g_wt"] = Reformed_SHC_data["k_g_wi"] + Reformed_SHC_data["k_g_wo"]
        Reformed_SHC_data["k_g_wt"][Reformed_SHC_data["nu"] > float(self.cutoff_freq)] = 0
        Reformed_SHC_data["k_g_wt"][Reformed_SHC_data["k_g_wt"] < 0] = 0.001
        Reformed_SHC_data["Kt"] = Reformed_SHC_data["Ki"] + Reformed_SHC_data["Ko"]
        Reformed_SHC_data["L"] = model.get_cell()[direction, direction] * part_ratio

        # Average them and save in Reformed_SHC_data['Results']
        for key, value in Reformed_SHC_data.items():
            if key not in ["L"]:
                Reformed_SHC_data[key] = np.column_stack((value, value.mean(axis=1), value.std(axis=1) / sqrt(N_repeat)))

        if 'Results' not in Reformed_SHC_data:
            Reformed_SHC_data['Results'] = {}
        for key, col in zip(["in", "out", "tot"], ["k_g_wi", "k_g_wo", "k_g_wt"]):
            values = [trap(Reformed_SHC_data[col][:, i], dx=Reformed_SHC_data["nu"][0, 0]) for i in range(N_repeat)]
            Reformed_SHC_data['Results'][f"{key}_ave"] = np.mean(values)
            Reformed_SHC_data['Results'][f"{key}_std"] = np.std(values) / np.sqrt(N_repeat)

        if self.save_data:
            np.savez(os.path.join(self.directory, 'data_shc.npz'), **Reformed_SHC_data)
            self._export_shc_txt(Reformed_SHC_data)
        return Reformed_SHC_data

    def process(self):
        """Main processing function for HNEMD method"""
        # Initialization
        Reformed_HNEMD_data = {}
        col_HNEMD_name = ['kx_in', 'kx_out', 'ky_in', 'ky_out', 'kz_tot']
        raw_HNEMD_data = np.loadtxt(self.path['kappa'])
        raw_HNEMD_data = pd.DataFrame(raw_HNEMD_data, columns=col_HNEMD_name)

        # Get parameters from run.in
        with open(self.path['run'], 'r') as file:
            found_hnemd = False
            for line in file:
                if line.lstrip().startswith('#'):
                    continue
                if 'time_step' in line:
                    time_step = float(line.split()[1])
                if 'compute_hnemd' in line:
                    output_interval = int(line.split()[1])
                    parts = line.strip().split()
                    Fe_values = list(map(float, parts[2:5]))
                    directions = ['x', 'y', 'z']
                    for direction, value in zip(directions, Fe_values):
                        if value != 0:
                            HNEMD_direction = direction
                            Fe = value
                    found_hnemd = True
                if found_hnemd and 'run' in line:
                    HNEMD_run = int(line.split()[1])

        N_hnemd_data = HNEMD_run // output_interval
        N_repeat = len(raw_HNEMD_data) // N_hnemd_data
        if len(raw_HNEMD_data) % N_hnemd_data != 0:
            raise ValueError(f"The MD calculation seems to be not completed, please check it!")
        Time_upper = HNEMD_run * time_step * 1e-6  # ns

        # Classify and process data initially
        def running_ave(y, x):
            return cumulative_trapezoid(y, x, initial=0) / x

        for col in raw_HNEMD_data.columns:
            Reformed_HNEMD_data[col] = raw_HNEMD_data[col].values.reshape(N_hnemd_data, N_repeat, order='F')
            Reformed_HNEMD_data[col] = np.column_stack((Reformed_HNEMD_data[col], Reformed_HNEMD_data[col].mean(axis=1)))

        t = np.arange(1, Reformed_HNEMD_data['kx_in'].shape[0] + 1) * 1e-3  # ns
        for col in raw_HNEMD_data.columns:
            for i in range(N_repeat + 1):
                Reformed_HNEMD_data[col][:, i] = running_ave(Reformed_HNEMD_data[col][:, i] * self.scale_eff_size, t)

        for prefix in ["kx", "ky"]:
            Reformed_HNEMD_data[f"{prefix}_tot"] = Reformed_HNEMD_data[f"{prefix}_in"] + Reformed_HNEMD_data[f"{prefix}_out"]

        # Calculate the final results (average and standard error)
        def compute_mean_std(data_hnemd, keys_hnemd, N_hnemd):
            result = {}
            for key_ in keys_hnemd:
                values = data_hnemd[key_][-1, :-1]
                result[key_ + "_ave"] = values.mean()
                result[key_ + "_std"] = values.std() / np.sqrt(N_hnemd)
            return result

        keys = ["kx_in", "kx_out", "kx_tot", "ky_in", "ky_out", "ky_tot", "kz_tot"]
        Reformed_HNEMD_data['Results'] = compute_mean_std(Reformed_HNEMD_data, keys, N_repeat)
        res_h = Reformed_HNEMD_data['Results']

        if self.save_data:
            np.savez(os.path.join(self.directory, 'data_hnemd.npz'), **Reformed_HNEMD_data)
            self._export_txt(Reformed_HNEMD_data, t)

        # Print HNEMD results
        self._print_hnemd_results(res_h, HNEMD_direction)

        # Process SHC if available
        if self.has_shc:
            print("\n[INFO] SHC data detected, processing spectral heat current...")
            Reformed_SHC_data = self.process_SHC(Fe=Fe)
            res_s = Reformed_SHC_data['Results']
            self._print_shc_results(res_s)
        else:
            print("\n[INFO] No SHC data found (shc.out not present), skipping SHC analysis.")
            Reformed_SHC_data = None
            res_s = None

        # Visualization
        self._plot_results(Reformed_HNEMD_data, Reformed_SHC_data, res_h,
                           HNEMD_direction, Time_upper, N_repeat, t)

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
        print("=" * 70 + "\n")

    @staticmethod
    def _txt_block(name, unit, x_name, x_unit, x, arr, trailing_labels=("average",)):
        """Format one 2D array as a TAB-separated text block: x-column + repeat columns + trailing_labels"""
        n_repeat = arr.shape[1] - len(trailing_labels)
        header = [f"{x_name}({x_unit})"] + [f"repeat_{i + 1}" for i in range(n_repeat)] + list(trailing_labels)
        lines = [f"# ---- {name} ({unit}) ----", "# " + "\t".join(header)]
        for row in np.column_stack((x, arr)):
            lines.append("\t".join(f"{v:.6g}" for v in row))
        lines.append("")
        return lines

    def _export_txt(self, Reformed_HNEMD_data, t):
        """Export HNEMD data as a compact, Excel/Origin-friendly TAB-separated .txt file"""
        res = Reformed_HNEMD_data['Results']
        lines = [
            "# ==== plt_hnemd.py HNEMD data export ====",
            f"# Directory: {self.directory}",
            f"# Scale_eff_size: {self.scale_eff_size}",
            f"# N_repeat: {Reformed_HNEMD_data['kz_tot'].shape[1] - 1}",
            "# Columns are TAB-separated; lines starting with '#' are comments/headers.",
            "",
        ]

        for key in ["kx_in", "kx_out", "kx_tot", "ky_in", "ky_out", "ky_tot", "kz_tot"]:
            lines += self._txt_block(key, "W/(m·K)", "t", "ns", t, Reformed_HNEMD_data[key])

        lines.append("# ==== Summary results ====")
        lines.append("# quantity\taverage\tstd\tunit")
        for key in ["kx_in", "kx_out", "kx_tot", "ky_in", "ky_out", "ky_tot", "kz_tot"]:
            lines.append(f"{key}\t{res[key + '_ave']:.6g}\t{res[key + '_std']:.6g}\tW/(m·K)")

        with open(os.path.join(self.directory, 'data_hnemd.txt'), 'w') as f:
            f.write("\n".join(lines) + "\n")

    def _export_shc_txt(self, Reformed_SHC_data):
        """Export HNEMD's SHC data as a compact, Excel/Origin-friendly TAB-separated .txt file"""
        lines = [
            "# ==== plt_hnemd.py SHC data export ====",
            f"# Directory: {self.directory}",
            f"# Cutoff frequency: {self.cutoff_freq} THz",
            "# Columns are TAB-separated; lines starting with '#' are comments/headers.",
            "",
        ]

        lines += self._txt_block("Kt", "eV/ps", "t_corr", "ps", Reformed_SHC_data['t'][:, -2],
                                 Reformed_SHC_data['Kt'], trailing_labels=("average", "std"))
        for key in ("k_g_wi", "k_g_wo", "k_g_wt"):
            lines += self._txt_block(key, "W/(m·K·THz)", "nu", "THz", Reformed_SHC_data['nu'][:, -2],
                                     Reformed_SHC_data[key], trailing_labels=("average", "std"))

        res = Reformed_SHC_data["Results"]
        lines.append("# ==== Summary results (frequency-integrated) ====")
        lines.append("# quantity\taverage\tstd\tunit")
        for key in ["in", "out", "tot"]:
            lines.append(f"kappa_{key}\t{res[key + '_ave']:.6g}\t{res[key + '_std']:.6g}\tW/(m·K)")

        with open(os.path.join(self.directory, 'data_shc.txt'), 'w') as f:
            f.write("\n".join(lines) + "\n")

    @staticmethod
    def _plot_running_kappa_panel(t, N_repeat, Time_upper, curves, annotations, direction):
        """(a) Running-average thermal conductivity.

        curves: list of (data, color, lw) — each drawn as N_repeat faint gray traces
            plus one highlighted average trace.
        annotations: list of (y, text, color) result labels, color may be None for default.
        """
        set_fig_properties([gca()])
        for i in range(N_repeat):
            for data, _, _ in curves:
                plot(t, data[:, i], color='k', alpha=0.3)
        for data, color, lw in curves:
            plot(t, data[:, -1], color=color, lw=lw)

        for y, txt, color in annotations:
            kwargs = dict(ha='right', va='top', transform=plt.gca().transAxes)
            if color is not None:
                kwargs['color'] = color
            text(0.95, y, txt, **kwargs)

        xlim(0, Time_upper)
        xlabel('time (ns)')
        ylabel(r'$\kappa$ (W/(m·K))')
        title(f"(a) Running average thermal conductivity: along {direction}")

    @staticmethod
    def _plot_kt_panel(Reformed_SHC_data):
        """(b) Force-virial correlation function K_tot(t)"""
        set_fig_properties([gca()])
        plot(Reformed_SHC_data['t'][:, -2], Reformed_SHC_data['Kt'][:, -2] / Reformed_SHC_data['L'], lw=2)
        ylabel('K (eV/ps)')
        xlabel('Correlation time (ps)')
        title('(b) K$_{tot}$(t)')

    def _plot_shc_spectral_panel(self, Reformed_SHC_data, series):
        """(c) SHC spectral thermal conductivity.

        series: list of (col, color, lw, label) — label may be None to omit the legend entry.
        """
        set_fig_properties([gca()])
        for col, color, lw, label in series:
            plot(Reformed_SHC_data['nu'][:, -2], Reformed_SHC_data[col][:, -2], linewidth=lw, color=color, label=label)
            fill_between(Reformed_SHC_data['nu'][:, -2],
                         Reformed_SHC_data[col][:, -2] - Reformed_SHC_data[col][:, -1],
                         Reformed_SHC_data[col][:, -2] + Reformed_SHC_data[col][:, -1],
                         facecolor=color, alpha=0.3)
        if [label for *_, label in series if label]:
            legend(frameon=False, fontsize=fs)
        xlim(0, self.cutoff_freq)
        axhline(y=0, color='k', linestyle='--')
        ylabel(r'$\kappa$($\omega$) (W/(m·K·THz))')
        xlabel(r'$\nu$ (THz)')
        title('(c) Spectral thermal conductivity')

    def _plot_results(self, Reformed_HNEMD_data, Reformed_SHC_data, res_h,
                      HNEMD_direction, Time_upper, N_repeat, t):
        """Visualize HNEMD and SHC results"""

        if HNEMD_direction in ['x', 'y']:
            key = 'kx' if HNEMD_direction == 'x' else 'ky'

            if self.has_shc:
                figure(figsize=(10, 8))
                # (a) HNEMD running average
                subplot(2, 1, 1)
            else:
                figure(figsize=(8, 4))
                subplot(1, 1, 1)

            curves = [(Reformed_HNEMD_data[f"{key}_in"], 'C1', 3),
                      (Reformed_HNEMD_data[f"{key}_out"], 'C2', 3),
                      (Reformed_HNEMD_data[f"{key}_tot"], 'C0', 3)]
            annotations = [
                (0.93, f"$\\kappa_{{\\mathrm{{in}}}}$ = {res_h[f'{key}_in_ave']:.3f} ± {res_h[f'{key}_in_std']:.2f} W/(m·K)", 'C1'),
                (0.83, f"$\\kappa_{{\\mathrm{{out}}}}$ = {res_h[f'{key}_out_ave']:.3f} ± {res_h[f'{key}_out_std']:.2f} W/(m·K)", 'C2'),
                (0.73, f"$\\kappa_{{\\mathrm{{tot}}}}$ = {res_h[f'{key}_tot_ave']:.3f} ± {res_h[f'{key}_tot_std']:.2f} W/(m·K)", 'C0'),
            ]
            self._plot_running_kappa_panel(t, N_repeat, Time_upper, curves, annotations, HNEMD_direction)

            if self.has_shc:
                subplot(2, 4, 5)
                self._plot_kt_panel(Reformed_SHC_data)

                subplot2grid((2, 4), (1, 1), colspan=3)
                self._plot_shc_spectral_panel(Reformed_SHC_data, [
                    ('k_g_wi', 'C1', 2, 'In-plane component'),
                    ('k_g_wo', 'C2', 2, 'Out-of-plane component'),
                    ('k_g_wt', 'C0', 2, 'Total'),
                ])

        elif HNEMD_direction == "z":
            if self.has_shc:
                figure(figsize=(10, 8))
                subplot(2, 1, 1)
            else:
                figure(figsize=(8, 4))
                subplot(1, 1, 1)

            curves = [(Reformed_HNEMD_data["kz_tot"], 'C0', 5)]
            annotations = [(0.9, fr"$\kappa_{{tot}}$ = {res_h['kz_tot_ave']:.3f} ± {res_h['kz_tot_std']:.2f} W/(m·K)", None)]
            self._plot_running_kappa_panel(t, N_repeat, Time_upper, curves, annotations, HNEMD_direction)

            if self.has_shc:
                subplot(2, 4, 5)
                self._plot_kt_panel(Reformed_SHC_data)

                subplot2grid((2, 4), (1, 1), colspan=3)
                self._plot_shc_spectral_panel(Reformed_SHC_data, [('k_g_wt', 'C0', 3, None)])

        tight_layout()

        if self.save:
            savefig('hnemd.png', dpi=300, bbox_inches='tight')
        else:
            show()



if __name__ == "__main__":

    argv = sys.argv[1:]

    if argv and argv[0] in ('-h', '--help', 'help'):
        print_usage()
        sys.exit(0)

    # --save/--save-data are independent flags: they can appear anywhere,
    # in any combination, without disturbing the positional numeric arguments below.
    save = '--save' in argv or 'save' in argv
    save_data = '--save-data' in argv or 'save_data' in argv
    positional = [a for a in argv if a not in ('--save', 'save', '--save-data', 'save_data')]

    try:
        scale_eff_size = float(positional[0]) if len(positional) > 0 else 1
        cutoff_freq = float(positional[1]) if len(positional) > 1 else 60
        if len(positional) > 2:
            raise ValueError

    except (ValueError, IndexError):
        print_usage()
        sys.exit(1)

    directory = os.getcwd()

    processor = HNEMD_Processor(directory, scale_eff_size, cutoff_freq, _save=save, _save_data=save_data)
    processor.process()
    # python plt_hnemd.py [scale_eff_size] [cutoff_freq] [--save] [--save-data]
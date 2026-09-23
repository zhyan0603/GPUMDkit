"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     plt_emd.py
Category:   Plot Scripts
Purpose:    Post-processing for EMD (Equilibrium Molecular Dynamics) thermal
            conductivity calculations, including HAC, running thermal
            conductivity, results printing, and optional data export.
Usage:      gpumdkit.sh -plt emd <direction> [--save] [--save-data]
            python plt_emd.py <x|y|z> [--save] [--save-data]
Arguments:
  direction  Heat transfer direction: x, y, or z
  --save      Save the plot as 'emd.png' instead of displaying it
  --save-data Save processed data as 'data_emd.npz' and 'data_emd.txt'
Output:
  emd.png    (if --save is used)
  data_emd.* (if --save-data is used)
Author:     Xin Wu (xinwuchn97@gmial.com)
Last-modified: 2026-09-14
=============================================================================
"""

from pylab import *
import pandas as pd
import numpy as np
import os
import sys

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

def print_usage():
    """Print usage instructions"""
    print("Usage: gpumdkit.sh -plt emd direction [--save] [--save-data]")
    print("Params:")
    print("  direction    : Heat transfer direction in lowercase (x, y, or z)")
    print("  --save       : Optional, save the plot as 'emd.png' (default: show it interactively)")
    print("  --save-data  : Optional, also save the processed data as 'data_emd.npz' and 'data_emd.txt'")
    print("                 (independent of --save; can be used on its own)")
    print("  Legacy bare 'save' and 'save_data' are also accepted.")
    print("  --save/--save-data may appear anywhere on the command line, e.g.:")
    print("    gpumdkit.sh -plt emd z --save-data")
    print("    gpumdkit.sh -plt emd z --save --save-data")

class EMD_Processor:
    def __init__(self, _directory, _direction='z', _save=False, _save_data=False):
        """
        Initialize EMD processor

        Parameters:
        -----------
        directory : str
            Path to the directory containing GPUMD output files
        direction : str
            Heat transfer direction ('x', 'y', or 'z')
        save : bool
            Save the plot as 'emd.png' instead of showing it interactively
        save_data : bool
            Also save the processed data as 'data_emd.npz'/'data_emd.txt'
        """
        self.directory = _directory
        self.direction = _direction
        self.save = _save
        self.save_data = _save_data
        self.path = {
            'run': os.path.join(self.directory, 'run.in'),
            'hac': os.path.join(self.directory, 'hac.out')
        }

    def process(self):
        """Main processing function for EMD method"""
        # Initialization
        Reformed_EMD_data = {}
        col_name = ['time', 'jx_in', 'jx_out', 'jy_in', 'jy_out', 'jz_tot',
                    'kx_in', 'kx_out', 'ky_in', 'ky_out', 'kz_tot']
        raw_data = np.loadtxt(self.path['hac'])
        raw_data = pd.DataFrame(raw_data, columns=col_name)

        # Get how many repeats you did
        with open(self.path['run'], 'r') as file:
            for line in file:
                if line.lstrip().startswith('#'):
                    continue
                if 'time_step' in line:
                    time_step = float(line.split()[1])
                if 'compute_hac' in line:
                    N_hac_data = int(int(line.split()[2]) / int(line.split()[3]))
                    Max_cor_time = int(int(line.split()[1]) * int(line.split()[2]))

        N_repeat = len(raw_data) // N_hac_data
        if len(raw_data) % N_hac_data != 0:
            raise ValueError(f"The MD calculation seems to be not completed, please check it!")
        Time_upper = Max_cor_time * time_step * 1e-6  # ns

        # Classify and process data initially
        for col in raw_data.columns:
            Reformed_EMD_data[col] = raw_data[col].values.reshape(N_hac_data, N_repeat, order='F')
            Reformed_EMD_data[col] = np.column_stack((Reformed_EMD_data[col], Reformed_EMD_data[col].mean(axis=1)))

        Reformed_EMD_data["kx_tot"] = Reformed_EMD_data["kx_in"] + Reformed_EMD_data["kx_out"]
        Reformed_EMD_data["ky_tot"] = Reformed_EMD_data["ky_in"] + Reformed_EMD_data["ky_out"]

        # Calculate the final results (average and standard error)
        def compute_mean_std(data_emd, keys_emd, N_emd):
            result = {}
            for key in keys_emd:
                values = data_emd[key][len(data_emd[key]) // 2:, -1]
                result[key + "_ave"] = values.mean()
                result[key + "_std"] = values.std() / np.sqrt(N_emd)
            return result

        keys = ["kx_in", "kx_out", "kx_tot", "ky_in", "ky_out", "ky_tot", "kz_tot"]
        Reformed_EMD_data['Results'] = compute_mean_std(Reformed_EMD_data, keys, N_repeat)

        if self.save_data:
            np.savez(os.path.join(self.directory, 'data_emd.npz'), **Reformed_EMD_data)
            self._export_txt(Reformed_EMD_data)

        # Print results
        self._print_results(Reformed_EMD_data['Results'])

        # Visualization
        self._plot_results(Reformed_EMD_data, Time_upper, N_repeat)

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

    def _export_txt(self, Reformed_EMD_data):
        """Export EMD data as a compact, Excel/Origin-friendly TAB-separated .txt file"""
        time_ns = Reformed_EMD_data["time"][:, -1] * 1e-3  # ns (average column; the grid is identical across repeats)
        lines = [
            "# ==== plt_emd.py EMD data export ====",
            f"# Directory: {self.directory}",
            f"# Direction: {self.direction}",
            f"# N_repeat: {Reformed_EMD_data['jz_tot'].shape[1] - 1}",
            "# Columns are TAB-separated; lines starting with '#' are comments/headers.",
            "",
        ]

        hac_blocks = [
            ("jx_in", "a.u."), ("jx_out", "a.u."),
            ("kx_in", "W/(m·K)"), ("kx_out", "W/(m·K)"), ("kx_tot", "W/(m·K)"),
            ("jy_in", "a.u."), ("jy_out", "a.u."),
            ("ky_in", "W/(m·K)"), ("ky_out", "W/(m·K)"), ("ky_tot", "W/(m·K)"),
            ("jz_tot", "a.u."), ("kz_tot", "W/(m·K)"),
        ]
        for key, unit in hac_blocks:
            lines += self._txt_block(key, unit, "t_corr", "ns", time_ns, Reformed_EMD_data[key])

        res = Reformed_EMD_data["Results"]
        lines.append("# ==== Summary results ====")
        lines.append("# quantity\taverage\tstd\tunit")
        for key in ["kx_in", "kx_out", "kx_tot", "ky_in", "ky_out", "ky_tot", "kz_tot"]:
            lines.append(f"{key}\t{res[key + '_ave']:.6g}\t{res[key + '_std']:.6g}\tW/(m·K)")

        with open(os.path.join(self.directory, 'data_emd.txt'), 'w', encoding='utf-8') as f:
            f.write("\n".join(lines) + "\n")

    def _print_results(self, results):
        """Print thermal conductivity results"""
        print("\n" + "=" * 60)
        print("EMD Thermal Conductivity Results")
        print("=" * 60)

        if self.direction in ['x', 'y']:
            key_prefix = f'k{self.direction}'
            print(f"\nDirection: {self.direction.lower()}\n")
            print(f"κ_in  = {results[key_prefix + '_in_ave']:.4f} ± {results[key_prefix + '_in_std']:.4f} W/(m·K)")
            print(f"κ_out = {results[key_prefix + '_out_ave']:.4f} ± {results[key_prefix + '_out_std']:.4f} W/(m·K)")
            print(f"κ_tot = {results[key_prefix + '_tot_ave']:.4f} ± {results[key_prefix + '_tot_std']:.4f} W/(m·K)")
        elif self.direction == 'z':
            print(f"Direction: z")
            print(f"κ = {results['kz_tot_ave']:.4f} ± {results['kz_tot_std']:.4f} W/(m·K)")

        print("=" * 60 + "\n")

    @staticmethod
    def _plot_hac_panel(time_data, series, N_repeat, Time_upper):
        """Draw the normalized-HAC loglog panel for one or more (data, color) series"""
        set_fig_properties([gca()])
        for i in range(N_repeat):
            for data, _ in series:
                loglog(time_data[:, i], data[:, i] / data[:, i].max(), color='k', alpha=0.3)
        for data, color in series:
            loglog(time_data[:, -1], data[:, -1] / data[:, -1].max(), color=color)
        xlim([1e-5, Time_upper])
        xlabel('Correlation Time (ns)')
        ylabel('Normalized HAC')
        title('(a)')

    @staticmethod
    def _plot_kappa_panel(time_data, kappa_data, N_repeat, Time_upper, ave, std, ylabel_text, title_text):
        """Draw a single kappa-vs-correlation-time panel with its running average and error band"""
        set_fig_properties([gca()])
        for i in range(N_repeat):
            plot(time_data[:, i], kappa_data[:, i], color='k', alpha=0.3)
        plot(time_data[:, -1], kappa_data[:, -1], color='C1', lw=3)
        axhline(y=ave, color='C0', linestyle='--')
        fill_between(time_data[:, -1], ave - std, ave + std, color='C0', alpha=0.2)
        xlim([0, Time_upper])
        xlabel('Correlation Time (ns)')
        ylabel(ylabel_text)
        title(title_text)

    def _plot_results(self, Reformed_EMD_data, Time_upper, N_repeat):
        """Visualize EMD results"""
        time_data = Reformed_EMD_data["time"] * 1e-3  # ns
        res = Reformed_EMD_data['Results']

        if self.direction in ["x", "y"]:
            d = self.direction
            key_map = {k: f"{k[0]}{d}{k[1:]}" for k in ("j_in", "j_out", "k_in", "k_out", "k_tot")}
            symbols = {"k_in": "in", "k_out": "out", "k_tot": "tot"}

            figure(figsize=(10, 8))

            # (a) Plot the normalized HAC
            subplot(2, 2, 1)
            self._plot_hac_panel(
                time_data,
                [(Reformed_EMD_data[key_map["j_in"]], 'C0'), (Reformed_EMD_data[key_map["j_out"]], 'C1')],
                N_repeat, Time_upper)

            # (b)/(c)/(d) IN, OUT and TOT components
            for panel, key in zip((2, 3, 4), ("k_in", "k_out", "k_tot")):
                subplot(2, 2, panel)
                sym = symbols[key]
                ave, std = res[key_map[key] + "_ave"], res[key_map[key] + "_std"]
                self._plot_kappa_panel(
                    time_data, Reformed_EMD_data[key_map[key]], N_repeat, Time_upper, ave, std,
                    rf'$\kappa_{{\mathrm{{{sym}}}}}$ (W/(m·K))',
                    rf"$\kappa_{{\mathrm{{{sym}}}}}$ = {ave:.2f} ± {std:.2f} W/(m·K)")

        elif self.direction == "z":
            figure(figsize=(10, 4))

            # (a) only jz_tot for z direction
            subplot(1, 2, 1)
            self._plot_hac_panel(time_data, [(Reformed_EMD_data["jz_tot"], 'C1')], N_repeat, Time_upper)

            # (b) plot kz_tot
            subplot(1, 2, 2)
            ave, std = res["kz_tot_ave"], res["kz_tot_std"]
            self._plot_kappa_panel(
                time_data, Reformed_EMD_data["kz_tot"], N_repeat, Time_upper, ave, std,
                r'$\kappa$ (W/(m·K))', fr"(b) $\kappa$ = {ave:.2f} ± {std:.2f} W/(m·K)")

        tight_layout()

        if self.save:
            savefig('emd.png', dpi=300, bbox_inches='tight')
        else:
            show()


if __name__ == "__main__":

    argv = sys.argv[1:]

    if not argv or argv[0] in ('-h', '--help', 'help'):
        print_usage()
        sys.exit(0)

    # --save/--save-data are independent flags: they can appear anywhere,
    # in any combination, without disturbing the positional 'direction' argument.
    save = '--save' in argv or 'save' in argv
    save_data = '--save-data' in argv or 'save_data' in argv
    positional = [a for a in argv if a not in ('--save', 'save', '--save-data', 'save_data')]

    if len(positional) != 1 or positional[0] not in ("x", "y", "z"):
        print_usage()
        sys.exit(1)

    directory = os.getcwd()  # Change this to your directory
    direction = positional[0]  # Heat transfer direction: 'x', 'y', or 'z'

    processor = EMD_Processor(directory, direction, save, save_data)
    processor.process()
    # python plt_emd.py direction [--save] [--save-data]

"""
This script is part of GPUMDkit.
Repository: https://github.com/zhyan0603/GPUMDkit

Description:
    Plot equilibrium molecular dynamics results

Usage:
    python plt_emd.py [arguments]

Author: Zihan YAN
Contact: yanzihan@westlake.edu.cn
Last Modified: 2025-12-28
"""


from pylab import *
import pandas as pd
import numpy as np
import os

# Figure Properties
aw, fs = 1.2, 12
matplotlib.rc('font', size=fs)
matplotlib.rc('axes', linewidth=aw)

def set_fig_properties(ax_list, tl=4, tw=1.2, tlm=4):
    """Set figure properties for axes"""
    for ax in ax_list:
        ax.tick_params(which='both', length=tl, width=tw, direction='in', right=True, top=True)
        ax.tick_params(which='minor', length=tlm)

def print_usage():
    """Print usage instructions"""
    print("Usage: gpumdkit.sh -plt emd direction [save]")
    print("Params:")
    print("  direction: Heat transfer direction in lowercase (x, y, or z)")
    print("  save     : Optional, save the plot as 'emd.png'")

class EMD_Processor:
    def __init__(self, _directory, _direction='z'):
        """
        Initialize EMD processor

        Parameters:
        -----------
        directory : str
            Path to the directory containing GPUMD output files
        direction : str
            Heat transfer direction ('x', 'y', or 'z')
        """
        self.directory = _directory
        self.direction = _direction
        self.path = {
            'run': self.directory + '/run.in',
            'hac': self.directory + '/hac.out'
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
                if 'time_step' in line:
                    time_step = int(line.split()[1])
                if 'compute_hac' in line:
                    N_hac_data = int(int(line.split()[2]) / int(line.split()[3]))
                    Max_cor_time = int(int(line.split()[1]) * int(line.split()[2]))

        N_repeat = len(raw_data) // N_hac_data
        if len(raw_data) % N_hac_data != 0:
            raise ValueError(f"The MD calculation seems to be not completed, please check it!")
        Time_upper = Max_cor_time * 1e-6  # ns

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

        if len(sys.argv) > 3 and sys.argv[3] == 'save_data':
            np.savez(f'{self.directory}/data_emd.npz', **Reformed_EMD_data)

        # Print results
        self._print_results(Reformed_EMD_data['Results'])

        # Visualization
        self._plot_results(Reformed_EMD_data, Time_upper, N_repeat)

    def _print_results(self, results):
        """Print thermal conductivity results"""
        print("\n" + "=" * 60)
        print("EMD Thermal Conductivity Results")
        print("=" * 60)

        if self.direction in ['x', 'y']:
            key_prefix = f'k{self.direction}'
            print(f"\nDirection: {self.direction.lower()}\n")
            print(f"κ_in  = {results[key_prefix + '_in_ave']:.4f} ± {results[key_prefix + '_in_std']:.4f} W/mK")
            print(f"κ_out = {results[key_prefix + '_out_ave']:.4f} ± {results[key_prefix + '_out_std']:.4f} W/mK")
            print(f"κ_tot = {results[key_prefix + '_tot_ave']:.4f} ± {results[key_prefix + '_tot_std']:.4f} W/mK")
        elif self.direction == 'z':
            print(f"Direction: z")
            print(f"κ = {results['kz_tot_ave']:.4f} ± {results['kz_tot_std']:.4f} W/mK")

        print("=" * 60 + "\n")

    def _plot_results(self, Reformed_EMD_data, Time_upper, N_repeat):
        """Visualize EMD results"""
        time_data = Reformed_EMD_data["time"] * 1e-3  # ns
        res = Reformed_EMD_data['Results']

        if self.direction in ["x", "y"]:
            key_map = {
                "x": {"j_in": "jx_in", "j_out": "jx_out", "k_in": "kx_in", "k_out": "kx_out", "k_tot": "kx_tot"},
                "y": {"j_in": "jy_in", "j_out": "jy_out", "k_in": "ky_in", "k_out": "ky_out", "k_tot": "ky_tot"}
            }[self.direction]

            figure(figsize=(10, 8))

            # (a) Plot the normalized HAC
            subplot(2, 2, 1)
            set_fig_properties([gca()])
            for i in range(N_repeat):
                loglog(time_data[:, i], Reformed_EMD_data[key_map["j_in"]][:, i] / Reformed_EMD_data[key_map["j_in"]][:, i].max(),
                       color='k', alpha=0.3)
                loglog(time_data[:, i], Reformed_EMD_data[key_map["j_out"]][:, i] / Reformed_EMD_data[key_map["j_out"]][:, i].max(),
                       color='k', alpha=0.3)
            loglog(time_data[:, -1], Reformed_EMD_data[key_map["j_in"]][:, -1] / Reformed_EMD_data[key_map["j_in"]][:, -1].max(),
                   color='C0')
            loglog(time_data[:, -1], Reformed_EMD_data[key_map["j_out"]][:, -1] / Reformed_EMD_data[key_map["j_out"]][:, -1].max(),
                   color='C1')
            xlim([1e-5, Time_upper])
            # ylim([1e-2, 1])
            xlabel('Correlation Time (ns)')
            ylabel('Normalized HAC')
            title('(a)')

            # (b) IN component
            subplot(2, 2, 2)
            set_fig_properties([gca()])
            for i in range(N_repeat):
                plot(time_data[:, i], Reformed_EMD_data[key_map["k_in"]][:, i], color='k', alpha=0.3)
            plot(time_data[:, -1], Reformed_EMD_data[key_map["k_in"]][:, -1], color='C1', lw=3)
            axhline(y=res[key_map["k_in"] + "_ave"], color='C0', linestyle='--')
            fill_between(time_data[:, -1], res[key_map["k_in"] + "_ave"] - res[key_map["k_in"] + "_std"],
                         res[key_map["k_in"] + "_ave"] + res[key_map["k_in"] + "_std"], color='C0', alpha=0.2)
            xlim([0, Time_upper])
            xlabel('Correlation Time (ns)')
            ylabel(r'$\kappa^{in}$ (W/mK)')
            title(fr"(b) $\kappa_{{in}}$ = {res[key_map['k_in'] + '_ave']:.2f} ± {res[key_map['k_in'] + '_std']:.2f} W/mK")

            # (c) OUT component
            subplot(2, 2, 3)
            set_fig_properties([gca()])
            for i in range(N_repeat):
                plot(time_data[:, i], Reformed_EMD_data[key_map["k_out"]][:, i], color='k', alpha=0.3)
            plot(time_data[:, -1], Reformed_EMD_data[key_map["k_out"]][:, -1], color='C1', lw=3)
            axhline(y=res[key_map["k_out"] + "_ave"], color='C0', linestyle='--')
            fill_between(time_data[:, -1], res[key_map["k_out"] + "_ave"] - res[key_map["k_out"] + "_std"],
                         res[key_map["k_out"] + "_ave"] + res[key_map["k_out"] + "_std"], color='C0', alpha=0.2)
            xlim([0, Time_upper])
            xlabel('Correlation Time (ns)')
            ylabel(r'$\kappa^{out}$ (W/mK)')
            title(fr"(c) $\kappa_{{out}}$ = {res[key_map['k_out'] + '_ave']:.2f} ± {res[key_map['k_out'] + '_std']:.2f} W/mK")

            # (d) TOT
            subplot(2, 2, 4)
            set_fig_properties([gca()])
            for i in range(N_repeat):
                plot(time_data[:, i], Reformed_EMD_data[key_map["k_tot"]][:, i], color='k', alpha=0.3)
            plot(time_data[:, -1], Reformed_EMD_data[key_map["k_tot"]][:, -1], color='C1', lw=3)
            axhline(y=res[key_map["k_tot"] + "_ave"], color='C0', linestyle='--')
            fill_between(time_data[:, -1], res[key_map["k_tot"] + "_ave"] - res[key_map["k_tot"] + "_std"],
                         res[key_map["k_tot"] + "_ave"] + res[key_map["k_tot"] + "_std"], color='C0', alpha=0.2)
            xlim([0, Time_upper])
            xlabel('Correlation Time (ns)')
            ylabel(r'$\kappa^{tot}$ (W/mK)')
            title(fr"(d) $\kappa_{{tot}}$ = {res[key_map['k_tot'] + '_ave']:.2f} ± {res[key_map['k_tot'] + '_std']:.2f} W/mK")

        elif self.direction == "z":
            figure(figsize=(10, 4))

            # (a) only jz_tot for z direction
            subplot(1, 2, 1)
            set_fig_properties([gca()])
            for i in range(N_repeat):
                loglog(time_data[:, i], Reformed_EMD_data["jz_tot"][:, i] / Reformed_EMD_data["jz_tot"][:, i].max(),
                       color='k', alpha=0.3)
            loglog(time_data[:, -1], Reformed_EMD_data["jz_tot"][:, -1] / Reformed_EMD_data["jz_tot"][:, -1].max(),
                   color='C1')
            xlim([1e-5, Time_upper])
            # ylim([1e-2, 1])
            xlabel('Correlation Time (ns)')
            ylabel('Normalized HAC')
            title('(a)')

            # (b) plot kz_tot
            subplot(1, 2, 2)
            set_fig_properties([gca()])
            for i in range(N_repeat):
                plot(time_data[:, i], Reformed_EMD_data["kz_tot"][:, i], color='k', alpha=0.3)
            plot(time_data[:, -1], Reformed_EMD_data["kz_tot"][:, -1], color='C1', lw=3)
            axhline(y=res["kz_tot_ave"], color='C0', linestyle='--')
            fill_between(time_data[:, -1], res["kz_tot_ave"] - res["kz_tot_std"],
                         res["kz_tot_ave"] + res["kz_tot_std"], color='C0', alpha=0.2)
            xlim([0, Time_upper])
            xlabel('Correlation Time (ns)')
            ylabel(r'$\kappa$ (W/mK)')
            title(fr"(b) $\kappa$ = {res['kz_tot_ave']:.2f} ± {res['kz_tot_std']:.2f} W/mK")

        tight_layout()

        if len(sys.argv) > 2 and sys.argv[2] == 'save':
            savefig('emd.png', dpi=300, bbox_inches='tight')
        else:
            show()


if __name__ == "__main__":

    if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help', 'help']:
        print_usage()
        sys.exit(0)

    if len(sys.argv) < 2 or sys.argv[1] not in ("x", "y", "z"):
        print_usage()
        sys.exit(1)

    directory = os.getcwd()  # Change this to your directory
    direction = sys.argv[1]  # Heat transfer direction: 'x', 'y', or 'z'

    processor = EMD_Processor(directory, direction)
    processor.process()
    # python plt_emd.py [direction] [save]
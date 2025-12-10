"""
@Author   : Xin Wu
@Contact  : xinwuchn97@gmial.com
@Remark   : Post-processing script for HNEMD (Homogenous Non-Equilibrium Molecular Dynamics) thermal conductivity calculations
"""

import pandas as pd
from pylab import *
import numpy as np
import os
from scipy.integrate import cumulative_trapezoid
from ase.io import read

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
    print("Usage: gpumdkit -plt hnemd [scale_eff_size] [cutoff_freq] [save]")
    print("Params:")
    print("  scale_eff_size: Optional, Scale factor for effective cross-sectional area (default: 1)")
    print("                   • For 3D bulk systems: use 1")
    print("                   • For low-dimensional systems with vacuum layer: S_box / S_eff")
    print("                     - S_box: box area perpendicular to heat transfer direction")
    print("                     - S_eff: real or effective area of the system")
    print("  cutoff_freq   : Optional, Cutoff frequency for SHC calculation in THz (default: 60)")
    print("  save          : Optional, save the plot as 'hnemd.png'")
    print("  !!! Note !!!  : If no SHC data, set [scale_eff_size] and [cutoff_freq] to any number as placeholders when using 'save'.")

class HNEMD_Processor:
    def __init__(self, _directory, _scale_eff_size=1, _cutoff_freq=60):
        """
        Initialize HNEMD processor

        Parameters:
        -----------
        directory : str
        scale_eff_size : float
        cutoff_freq : float
        """
        self.directory = _directory
        self.scale_eff_size = _scale_eff_size
        self.cutoff_freq = _cutoff_freq
        self.path = {
            'run': self.directory + '/run.in',
            'kappa': self.directory + '/kappa.out'
        }
        self.has_shc = os.path.exists(self.directory + '/shc.out')

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
        self.path['shc'] = self.directory + "/shc.out"
        self.path['model'] = self.directory + "/model.xyz"

        col_shc_name = ['t_omega', 'Ki_jwi', 'Ko_jwo']
        raw_shc_data = np.loadtxt(self.path['shc'])
        raw_shc_data = pd.DataFrame(raw_shc_data, columns=col_shc_name)

        with open(self.path['run'], 'r') as file:
            for line in file:
                if 'compute_shc' in line:
                    Max_cor_step = int(line.split()[2])
                    N_omega = int(line.split()[4])
                    direction = int(line.split()[3])
                    grouping_th = int(line.split()[7])
                    group_shc_th = int(line.split()[8])
                if 'nvt_nhc' in line:
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
        if grouping_th == 0:
            group_index = 'group'
        elif grouping_th == 1:
            raise ValueError("Please check your grouping for SHC method, it seems not the Default setting !!!")
        group = model.get_array(group_index)
        part_ratio = np.sum(group == group_shc_th) / group.size

        vol = model.get_volume() * part_ratio / self.scale_eff_size
        convert = 1.602176634e3  # ev*A/ps/THz * 1/A^3 *1/K * A ==> W/m/K/THz
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
            values = [np.trapezoid(Reformed_SHC_data[col][:, i], dx=Reformed_SHC_data["nu"][0, 0]) for i in range(N_repeat)]
            Reformed_SHC_data['Results'][f"{key}_ave"] = np.mean(values)
            Reformed_SHC_data['Results'][f"{key}_std"] = np.std(values) / np.sqrt(N_repeat)

        if len(sys.argv) > 4 and sys.argv[4] == 'save_data':
            np.savez(f'{self.directory}/data_shc.npz', **Reformed_SHC_data)
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
                if 'time_step' in line:
                    time_step = int(line.split()[1])
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

        if len(sys.argv) > 4 and sys.argv[4] == 'save_data':
            np.savez(f'{self.directory}/data_hnemd.npz', **Reformed_HNEMD_data)

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
        self._plot_results(Reformed_HNEMD_data, Reformed_SHC_data, res_h, res_s,
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
            print(f"\nκ_in  = {results[key_prefix + '_in_ave']:.4f} ± {results[key_prefix + '_in_std']:.4f} W/mK")
            print(f"κ_out = {results[key_prefix + '_out_ave']:.4f} ± {results[key_prefix + '_out_std']:.4f} W/mK")
            print(f"κ_tot = {results[key_prefix + '_tot_ave']:.4f} ± {results[key_prefix + '_tot_std']:.4f} W/mK")
        elif direction == 'z':
            print(f"\nκ = {results['kz_tot_ave']:.4f} ± {results['kz_tot_std']:.4f} W/mK")

        print("=" * 70)

    def _print_shc_results(self, results):
        """Print SHC spectral thermal conductivity results"""
        print("\n" + "=" * 70)
        print("SHC Spectral Thermal Conductivity Results")
        print("=" * 70)
        print(f"Cutoff frequency: {self.cutoff_freq} THz\n")
        print(f"κ_in  (integrated) = {results['in_ave']:.4f} ± {results['in_std']:.4f} W/mK")
        print(f"κ_out (integrated) = {results['out_ave']:.4f} ± {results['out_std']:.4f} W/mK")
        print(f"κ_tot (integrated) = {results['tot_ave']:.4f} ± {results['tot_std']:.4f} W/mK")
        print("=" * 70 + "\n")

    def _plot_results(self, Reformed_HNEMD_data, Reformed_SHC_data, res_h, res_s,
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

            set_fig_properties([gca()])
            for i in range(N_repeat):
                plot(t, Reformed_HNEMD_data[f"{key}_in"][:, i], color='k', alpha=0.3)
                plot(t, Reformed_HNEMD_data[f"{key}_out"][:, i], color='k', alpha=0.3)
                plot(t, Reformed_HNEMD_data[f"{key}_tot"][:, i], color='k', alpha=0.3)
            plot(t, Reformed_HNEMD_data[f"{key}_in"][:, -1], color='C1', lw=3)
            plot(t, Reformed_HNEMD_data[f"{key}_out"][:, -1], color='C2', lw=3)
            plot(t, Reformed_HNEMD_data[f"{key}_tot"][:, -1], color='C0', lw=3)

            text(0.95, 0.93, f"$\kappa_{{in}}$ = {res_h[f'{key}_in_ave']:.3f} ± {res_h[f'{key}_in_std']:.2f} W/mK",
                 ha='right', va='top', transform=plt.gca().transAxes, color='C1')
            text(0.95, 0.83, f"$\kappa_{{out}}$ = {res_h[f'{key}_out_ave']:.3f} ± {res_h[f'{key}_out_std']:.2f} W/mK",
                 ha='right', va='top', transform=plt.gca().transAxes, color='C2')
            text(0.95, 0.73, f"$\kappa_{{tot}}$ = {res_h[f'{key}_tot_ave']:.3f} ± {res_h[f'{key}_tot_std']:.2f} W/mK",
                 ha='right', va='top', transform=plt.gca().transAxes, color='C0')
            xlim(0, Time_upper)
            xlabel('time (ns)')
            ylabel(r'$\kappa$ (W/mK)')
            title(f"(a) Running average thermal conductivity: along {HNEMD_direction}")

            if self.has_shc:
                # (b) Force-virial Correlation function
                subplot(2, 4, 5)
                set_fig_properties([gca()])
                plot(Reformed_SHC_data['t'][:, -2], Reformed_SHC_data['Kt'][:, -2] / Reformed_SHC_data['L'], lw=2)
                ylabel('K (eV/ps)')
                xlabel('Correlation time (ps)')
                title('(b) K$_{tot}$(t)')

                # (c) SHC
                subplot2grid((2, 4), (1, 1), colspan=3)
                set_fig_properties([gca()])
                plot(Reformed_SHC_data['nu'][:, -2], Reformed_SHC_data['k_g_wi'][:, -2], linewidth=2, color='C1')
                fill_between(Reformed_SHC_data['nu'][:, -2],
                             Reformed_SHC_data['k_g_wi'][:, -2] - Reformed_SHC_data['k_g_wi'][:, -1],
                             Reformed_SHC_data['k_g_wi'][:, -2] + Reformed_SHC_data['k_g_wi'][:, -1],
                             facecolor='C1', alpha=0.3)
                plot(Reformed_SHC_data['nu'][:, -2], Reformed_SHC_data['k_g_wo'][:, -2], linewidth=2, color='C2')
                fill_between(Reformed_SHC_data['nu'][:, -2],
                             Reformed_SHC_data['k_g_wo'][:, -2] - Reformed_SHC_data['k_g_wo'][:, -1],
                             Reformed_SHC_data['k_g_wo'][:, -2] + Reformed_SHC_data['k_g_wo'][:, -1],
                             facecolor='C2', alpha=0.3)
                plot(Reformed_SHC_data['nu'][:, -2], Reformed_SHC_data['k_g_wt'][:, -2], linewidth=2, color='C0')
                fill_between(Reformed_SHC_data['nu'][:, -2],
                             Reformed_SHC_data['k_g_wt'][:, -2] - Reformed_SHC_data['k_g_wt'][:, -1],
                             Reformed_SHC_data['k_g_wt'][:, -2] + Reformed_SHC_data['k_g_wt'][:, -1],
                             facecolor='C0', alpha=0.3)
                text(0.6, 0.9, f"$\\kappa_{{in}}$ = {res_s['in_ave']:.2f} ± {res_s['in_std']:.2f} W/mK",
                     ha='left', va='top', transform=gca().transAxes, color='C1')
                text(0.6, 0.8, f"$\\kappa_{{out}}$ = {res_s['out_ave']:.2f} ± {res_s['out_std']:.2f} W/mK",
                     ha='left', va='top', transform=gca().transAxes, color='C2')
                text(0.6, 0.7, f"$\\kappa_{{tot}}$ = {res_s['tot_ave']:.2f} ± {res_s['tot_std']:.2f} W/mK",
                     ha='left', va='top', transform=gca().transAxes, color='C0')
                xlim(0, self.cutoff_freq)
                axhline(y=0, color='k', linestyle='--')
                ylabel(r'$\kappa$($\omega$) (W/m/K/THz)')
                xlabel(r'$\nu$ (THz)')
                title('(c) Spectral thermal conductivity')

        elif HNEMD_direction == "z":
            if self.has_shc:
                figure(figsize=(10, 8))
                subplot(2, 1, 1)
            else:
                figure(figsize=(8, 4))
                subplot(1, 1, 1)

            set_fig_properties([gca()])
            for i in range(N_repeat):
                plot(t, Reformed_HNEMD_data["kz_tot"][:, i], color='k', alpha=0.3)
            plot(t, Reformed_HNEMD_data["kz_tot"][:, -1], color='C0', lw=5)

            text(0.95, 0.9, fr"$\kappa_{{tot}}$ = {res_h['kz_tot_ave']:.3f} ± {res_h['kz_tot_std']:.2f} W/mK",
                 ha='right', va='top', transform=plt.gca().transAxes)
            xlim(0, Time_upper)
            xlabel('time (ns)')
            ylabel(r'$\kappa$ (W/mK)')
            title(f"(a) Running average thermal conductivity: along {HNEMD_direction}")

            if self.has_shc:
                # (b) Force-virial Correlation function
                subplot(2, 4, 5)
                set_fig_properties([gca()])
                plot(Reformed_SHC_data['t'][:, -2], Reformed_SHC_data['Kt'][:, -2] / Reformed_SHC_data['L'], lw=2)
                ylabel('K (eV/ps)')
                xlabel('Correlation time (ps)')
                title('(b) K$_{tot}$(t)')

                # (c) SHC
                subplot2grid((2, 4), (1, 1), colspan=3)
                set_fig_properties([gca()])
                plot(Reformed_SHC_data['nu'][:, -2], Reformed_SHC_data['k_g_wt'][:, -2], linewidth=3, color='C0')
                fill_between(Reformed_SHC_data['nu'][:, -2],
                             Reformed_SHC_data['k_g_wt'][:, -2] - Reformed_SHC_data['k_g_wt'][:, -1],
                             Reformed_SHC_data['k_g_wt'][:, -2] + Reformed_SHC_data['k_g_wt'][:, -1],
                             facecolor='C0', alpha=0.3)
                text(0.6, 0.7, f"$\\kappa_{{tot}}$ = {res_s['tot_ave']:.2f} ± {res_s['tot_std']:.2f} W/mK",
                     ha='left', va='top', transform=gca().transAxes)
                xlim(0, self.cutoff_freq)
                axhline(y=0, color='k', linestyle='--')
                ylabel(r'$\kappa$($\omega$) (W/m/K/THz)')
                xlabel(r'$\nu$ (THz)')
                title('(c) Spectral thermal conductivity')

        tight_layout()
        if len(sys.argv) > 3 and sys.argv[3] == 'save':
            savefig(f'hnemd.png', dpi=300, bbox_inches='tight')
        show()



if __name__ == "__main__":

    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help', 'help']:
        print_usage()
        sys.exit(0)

    try:
        scale_eff_size = float(sys.argv[1]) if len(sys.argv) > 1 else 1
        cutoff_freq = float(sys.argv[2]) if len(sys.argv) > 2 else 60
        if len(sys.argv) > 3 and sys.argv[3] != 'save':
            raise ValueError
        if len(sys.argv) > 4 and sys.argv[4] != 'save_data':
            raise ValueError

    except (ValueError, IndexError):
        print_usage()
        sys.exit(1)

    directory = os.getcwd()

    processor = HNEMD_Processor(directory, scale_eff_size, cutoff_freq)
    processor.process()
    # python plt_hnemd.py [scale_eff_size] [cutoff_freq] [save]
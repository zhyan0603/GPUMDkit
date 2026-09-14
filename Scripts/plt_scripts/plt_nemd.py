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
Usage:      gpumdkit.sh -plt nemd [real_length] [scale_eff_size] [cutoff_freq] [save]
            python plt_nemd.py [real_length] [scale_eff_size] [cutoff_freq] [save]
Arguments:
  real_length     Real length of heat transfer zone in nm (set to 'Auto' for auto)
  scale_eff_size  Scale factor for effective cross-sectional area (default: 1)
  cutoff_freq     Cutoff frequency for SHC in THz (default: 60)
  save            Save the plot as 'nemd.png' instead of displaying it
Output:
  nemd.png    (if save is used, or if backend is non-interactive)
Author:     Xin Wu (xinwuchn97@gmial.com)
Last-modified: 2026-05-16
=============================================================================
"""

from pylab import *
import pandas as pd
import numpy as np
from ase.io import read
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

trap = np.trapezoid if hasattr(np, "trapezoid") else getattr(np, "trapz")

def print_usage():
    """Print usage instructions"""
    print("Usage: gpumdkit -plt nemd [real_length] [scale_eff_size] [cutoff_freq] [--save] [--save-data]")
    print("Params:")
    print("  real_length   : Real length of heat tranfer zone in nm (set to 'Auto', with auto-calculation)")
    print("  scale_eff_size: Optional, Scale factor for effective cross-sectional area (default: 1)")
    print("                   • For 3D bulk systems: use 1")
    print("                   • For low-dimensional systems with vacuum layer: S_box / S_eff")
    print("                     - S_box: box area perpendicular to heat transfer direction")
    print("                     - S_eff: real or effective area of the system")
    print("  cutoff_freq   : Optional, Cutoff frequency for SHC calculation in THz (default: 60)")
    print("  --save        : Optional, save the plot as 'nemd.png' (default: show it interactively)")
    print("  --save-data   : Optional, also save the processed data as 'data_nemd.npz'/'data_shc.npz'")
    print("                  and their .txt equivalents (independent of --save; can be used on its own)")
    print("  --save/--save-data may appear anywhere on the command line, e.g.:")
    print("    gpumdkit -plt nemd --save-data              (defaults for everything else)")
    print("    gpumdkit -plt nemd Auto 1 60 --save --save-data")

class NEMD_Processor:
    def __init__(self, _directory, _real_length=None, _scale_eff_size=1, _cutoff_freq=60, _scale_vacuum=1,
                 _save=False, _save_data=False):
        """
        Initialize NEMD processor

        Parameters:
        -----------
        directory : str
            Path to the directory containing GPUMD output files
        scale_vacuum : float
            Scale factor for vacuum along heat transfer direction (default: 1)
            For NEMD, box size along heat transfer may be longer than actual value
        scale_eff_size : float
            Scale factor for effective cross-sectional area (for low-dimensional systems)
            Scale = S_box / S_eff, where S_box is the box area perpendicular to heat transfer
        real_length : float or None
            Real/effective length in nm (if specified, overrides automatic calculation)
            Use this if N_bath != 1 or if automatic length is inaccurate
        cutoff_freq : float
            Cutoff frequency for SHC calculation (THz), material-dependent
        save : bool
            Save the plot as 'nemd.png' instead of showing it interactively
        save_data : bool
            Also save the processed data as 'data_nemd.npz'/'data_shc.npz' (and their .txt equivalents)
        """
        self.directory = _directory
        self.scale_vacuum = _scale_vacuum
        self.scale_eff_size = _scale_eff_size
        self.real_length = _real_length
        self.cutoff_freq = _cutoff_freq
        self.save = _save
        self.save_data = _save_data
        self.path = {
            'run': os.path.join(self.directory, 'run.in'),
            'compute': os.path.join(self.directory, 'compute.out'),
            'model': os.path.join(self.directory, 'model.xyz'),
            'shc': os.path.join(self.directory, 'shc.out')
        }
        self.has_shc = os.path.exists(self.path['shc'])

    def process_SHC(self, deltaT):
        """
        Process spectral heat current (SHC) data for NEMD

        Parameters:
        -----------
        deltaT : array
            Temperature difference across the system

        Returns:
        --------
        dict : Reformed SHC data including spectral thermal conductance
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
                # if 'nvt_' in line:
                #     Temp = int(line.split()[2])

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

        # Calculate g(omega) from jw for NEMD
        model = read(self.path['model'])
        if part_ratio == -1:
            group_arr = model.get_array('group')
            if group_arr.ndim == 1:
                group = group_arr
            else:
                group = group_arr[:, grouping_th]
            part_ratio = np.sum(group == group_shc_th) / group.size

        vol = model.get_volume() * part_ratio / self.scale_vacuum
        convert = 1.602176634e7  # ev*A/ps/THz * 1/A^3 *1/K * A ==> MW/(m·K·THz)
        denom = vol * deltaT

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
        """Main processing function for NEMD method"""
        raw_NEMD_data = np.loadtxt(self.path['compute'])
        Reformed_NEMD_data = {}

        # Get parameters from run.in
        direction = None  # heat transfer axis; set by compute_shc, else auto-detected below
        with open(self.path['run'], 'r') as file:
            found_nemd = False
            for line in file:
                if line.lstrip().startswith('#'):
                    continue
                if 'time_step' in line:
                    time_step = int(line.split()[1])
                if 'heat_lan' in line:
                    N_temp_group = int(int(line.split()[6]) + 1)
                if 'temperature' in line:
                    output_interval = int(int(line.split()[2]) * int(line.split()[3]))
                    found_nemd = True
                if 'compute_shc' in line:
                    direction = int(line.split()[3])
                if found_nemd and 'run' in line:
                    NEMD_run = int(line.split()[1])

        N_nemd_data = NEMD_run // output_interval
        N_repeat = len(raw_NEMD_data) // N_nemd_data
        col_name = [f"temp_{i}" for i in range(N_temp_group)] + ["E_in", "E_out"]
        raw_NEMD_data = pd.DataFrame(raw_NEMD_data, columns=col_name)

        if len(raw_NEMD_data) % N_nemd_data != 0:
            raise ValueError(f"The MD calculation seems to be not completed, please check it!")
        Time_upper = NEMD_run * time_step * 1e-6  # ns

        # Classify and process data initially
        for col in raw_NEMD_data.columns:
            Reformed_NEMD_data[col] = raw_NEMD_data[col].values.reshape(N_nemd_data, N_repeat, order='F')
            Reformed_NEMD_data[col] = np.column_stack((Reformed_NEMD_data[col], Reformed_NEMD_data[col].mean(axis=1)))
            if 'temp' in col:
                Reformed_NEMD_data[col] = np.vstack([Reformed_NEMD_data[col],
                                                     np.mean(Reformed_NEMD_data[col][N_nemd_data // 2:], axis=0)])

        delta_T = Reformed_NEMD_data['temp_1'][-1, :N_repeat] - Reformed_NEMD_data[f'temp_{N_temp_group - 1}'][-1, :N_repeat].reshape(1, -1)

        mid = int(N_nemd_data / 2)
        denom = (N_nemd_data / 2) * (time_step * 0.001) * output_interval
        Q_in = (Reformed_NEMD_data['E_in'][mid, :N_repeat] - Reformed_NEMD_data['E_in'][-1, :N_repeat]) / denom
        Q_out = (Reformed_NEMD_data['E_out'][-1, :N_repeat] - Reformed_NEMD_data['E_out'][mid, :N_repeat]) / denom
        Q = ((Q_in + Q_out) / 2).reshape(1, -1)  # eV/ps

        model = read(self.path['model'])
        group_arr = model.get_array('group')
        if group_arr.ndim == 1:
            group = group_arr
        else:
            group = group_arr[:, 0]

        if direction is None:
            # No compute_shc in run.in: infer the heat transfer axis as the
            # Cartesian direction where the source/sink groups are most separated
            heat_pos = model.positions[(group == 1)].mean(axis=0)
            cold_pos = model.positions[(group == N_temp_group - 1)].mean(axis=0)
            direction = int(np.argmax(np.abs(heat_pos - cold_pos)))

        coords_heat = model.positions[(group == 1), direction].mean()
        coords_cold = model.positions[(group == N_temp_group - 1), direction].mean()
        xticks_length = np.linspace(coords_heat, coords_cold, N_temp_group - 1) * 0.1

        if self.real_length is not None:
            Length = self.real_length
        else:
            Length = abs(coords_heat - coords_cold) * 0.1  # nm
        Reformed_NEMD_data['L'] = Length

        A = model.get_volume() / model.get_cell()[direction, direction] / self.scale_eff_size  # A^2
        convert = 1.602176634e7  # eV/ps * 1/A^2 * 1/K * A ==> MW/m^2/K
        G = convert * Q / delta_T / A
        k = G * Length * 1e-3  # W/(m·K)
        Reformed_NEMD_data['G'] = np.hstack((G, G.mean(axis=1, keepdims=True),
                                             G.std(axis=1, keepdims=True) / sqrt(N_repeat)))  # MW/m^2/K
        Reformed_NEMD_data['k'] = np.hstack((k, k.mean(axis=1, keepdims=True),
                                             k.std(axis=1, keepdims=True) / sqrt(N_repeat)))  # W/m/K

        if self.save_data:
            np.savez(os.path.join(self.directory, 'data_nemd.npz'), **Reformed_NEMD_data)
            self._export_txt(Reformed_NEMD_data, N_nemd_data, N_temp_group, xticks_length, Time_upper, Length)

        # Print NEMD results
        self._print_nemd_results(Reformed_NEMD_data, Length)

        # Process SHC if available
        if self.has_shc:
            print("\n[INFO] SHC data detected, processing spectral heat current...")
            Reformed_SHC_data = self.process_SHC(deltaT=delta_T)
            res_s = Reformed_SHC_data['Results']
            self._print_shc_results(res_s)
        else:
            print("\n[INFO] No SHC data found (shc.out not present), skipping SHC analysis.")
            Reformed_SHC_data = None
            res_s = None

        # Visualization
        self._plot_results(Reformed_NEMD_data, Reformed_SHC_data,
                           Time_upper, N_nemd_data, N_repeat,
                           xticks_length, Length)

    def _print_nemd_results(self, Reformed_NEMD_data, Length):
        """Print NEMD thermal conductivity results"""
        print("\n" + "=" * 70)
        print("NEMD Thermal Conductivity Results")
        print("=" * 70)
        print(f"\nEffective length: {Length:.4f} nm")
        # print(f"Scale_vacuum: {self.scale_vacuum}")
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
            if key in Reformed_NEMD_data:
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
            "# Columns are TAB-separated; lines starting with '#' are comments/headers.",
            "",
        ]

        lines += self._txt_block("Kt", "eV/ps", "t_corr", "ps", Reformed_SHC_data['t'][:, -2],
                                 Reformed_SHC_data['Kt'], trailing_labels=("average", "std"))
        for key, unit in (("k_g_wi", "MW/(m²·K·THz)"), ("k_g_wo", "MW/(m²·K·THz)"), ("k_g_wt", "MW/(m²·K·THz)")):
            lines += self._txt_block(key, unit, "nu", "THz", Reformed_SHC_data['nu'][:, -2],
                                     Reformed_SHC_data[key], trailing_labels=("average", "std"))

        res = Reformed_SHC_data["Results"]
        lines.append("# ==== Summary results (frequency-integrated) ====")
        lines.append("# quantity\taverage\tstd\tunit")
        for key in ["in", "out", "tot"]:
            lines.append(f"G_{key}\t{res[key + '_ave']:.6g}\t{res[key + '_std']:.6g}\tMW/(m²·K)")

        with open(os.path.join(self.directory, 'data_shc.txt'), 'w') as f:
            f.write("\n".join(lines) + "\n")

    @staticmethod
    def _plot_temperature_profile(Reformed_NEMD_data, N_repeat, xticks_length, Length):
        """(a) Temperature profile along the heat transfer direction"""
        set_fig_properties([gca()])
        x_values = []
        y_values = []
        for i in range(len(xticks_length)):
            key = f'temp_{i + 1}'
            if key in Reformed_NEMD_data:
                last_row = Reformed_NEMD_data[key][-1]
                for j in range(N_repeat):
                    plt.scatter(xticks_length[i], last_row[j], color="C0", alpha=0.3, s=30)
                plt.scatter(xticks_length[i], last_row[-1], color="C1", s=60)
                x_values.append(xticks_length[i])
                y_values.append(last_row[-1])
        plt.plot(x_values, y_values, '-', color="C0", alpha=0.4, linewidth=2)
        text(0.95, 0.9, f'G={Reformed_NEMD_data["G"][0, -2]:.3f}±{Reformed_NEMD_data["G"][0, -1]:.2f} MW/(m$^2$·K)',
             ha='right', va='top', transform=plt.gca().transAxes)
        text(0.95, 0.8, f'κ={Reformed_NEMD_data["k"][0, -2]:.3f}±{Reformed_NEMD_data["k"][0, -1]:.3f} W/(m·K)',
             ha='right', va='top', transform=plt.gca().transAxes)
        text(0.08, 0.08, f'Effective length={Length:.2f} nm', ha='left', va='bottom',
             transform=plt.gca().transAxes)
        xlabel("Length (nm)")
        ylabel("Temperature (K)")
        title("(a) Temperature profile")

    @staticmethod
    def _plot_thermostat_energy(Reformed_NEMD_data, Time_upper, N_nemd_data, N_repeat):
        """(b) Cumulative energy injected/extracted by the source/sink thermostats"""
        set_fig_properties([gca()])
        Time = np.linspace(0, Time_upper, N_nemd_data)
        slope_in = abs(Reformed_NEMD_data["E_in"][-1, -1] - Reformed_NEMD_data["E_in"][0, -1]) / Time_upper
        slope_out = abs(Reformed_NEMD_data["E_out"][-1, -1] - Reformed_NEMD_data["E_out"][0, -1]) / Time_upper

        for i in range(N_repeat):
            plt.plot(Time, Reformed_NEMD_data["E_in"][:, i] / 1000, color="C0", alpha=0.3, lw=2)
            plt.plot(Time, Reformed_NEMD_data["E_out"][:, i] / 1000, color="C1", alpha=0.3, lw=2)
        plt.plot(Time, Reformed_NEMD_data["E_in"][:, -1] / 1000, color="C0", lw=2, label="Source")
        plt.plot(Time, Reformed_NEMD_data["E_out"][:, -1] / 1000, color="C1", lw=2, label="Sink")

        text(0.08, 0.10, f"slope$_{{in}}$={slope_in:.2f}", ha='left', va='bottom',
             transform=plt.gca().transAxes, fontsize=13, color="C0")
        text(0.08, 0.88, f"slope$_{{out}}$={slope_out:.2f}", ha='left', va='top',
             transform=plt.gca().transAxes, fontsize=13, color="C1")

        xlim(0, Time_upper)
        xlabel("Time (ns)")
        ylabel(r"Energy ($\times 10^3$ eV)")
        legend(frameon=False, loc='center right')
        title("(b) Thermostat energy")

    @staticmethod
    def _plot_kt_panel(Reformed_SHC_data):
        """(c) Force-virial correlation function K_tot(t)"""
        set_fig_properties([gca()])
        plot(Reformed_SHC_data['t'][:, -2], Reformed_SHC_data['Kt'][:, -2] / Reformed_SHC_data['L'], lw=2)
        ylabel('K (eV/ps)')
        xlabel('Correlation time (ps)')
        title('(c) K$_{tot}$(t)')

    def _plot_shc_spectral_panel(self, Reformed_SHC_data):
        """(d) SHC spectral thermal conductance, split into in-/out-of-plane components"""
        set_fig_properties([gca()])
        for col, color, label in (("k_g_wi", "C1", "In-plane component"),
                                   ("k_g_wo", "C2", "Out-of-plane component"),
                                   ("k_g_wt", "C0", "Total")):
            plot(Reformed_SHC_data['nu'][:, -2], Reformed_SHC_data[col][:, -2], linewidth=2, color=color, label=label)
            plt.fill_between(Reformed_SHC_data['nu'][:, -2],
                             Reformed_SHC_data[col][:, -2] - Reformed_SHC_data[col][:, -1],
                             Reformed_SHC_data[col][:, -2] + Reformed_SHC_data[col][:, -1],
                             facecolor=color, alpha=0.3)
        legend(frameon=False, fontsize=fs)
        xlim(0, self.cutoff_freq)
        ylabel(r'$g$($\omega$) (MW/(m$^2$·K·THz))')
        xlabel(r'$\omega/2\pi$ (THz)')
        title('(d) Spectral thermal conductance')

    def _plot_results(self, Reformed_NEMD_data, Reformed_SHC_data,
                      Time_upper, N_nemd_data, N_repeat,
                      xticks_length, Length):
        """Visualize NEMD and SHC results"""

        if self.has_shc:
            figure(figsize=(10, 8))
            subplot2grid((2, 5), (0, 0), colspan=3)
        else:
            figure(figsize=(10, 4))
            subplot(1, 2, 1)
        self._plot_temperature_profile(Reformed_NEMD_data, N_repeat, xticks_length, Length)

        if self.has_shc:
            subplot2grid((2, 5), (0, 3), colspan=2)
        else:
            subplot(1, 2, 2)
        self._plot_thermostat_energy(Reformed_NEMD_data, Time_upper, N_nemd_data, N_repeat)

        if self.has_shc:
            subplot2grid((2, 4), (1, 0), colspan=1)
            self._plot_kt_panel(Reformed_SHC_data)

            subplot2grid((2, 4), (1, 1), colspan=3)
            self._plot_shc_spectral_panel(Reformed_SHC_data)

            plt.subplots_adjust(wspace=1, hspace=0.3)
        else:
            tight_layout()

        if self.save:
            savefig('nemd.png', dpi=300, bbox_inches='tight')
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
        if len(positional) > 0:
            if positional[0].lower() == "auto":
                real_length = None
            else:
                real_length = float(positional[0])
                if real_length <= 0:
                    raise ValueError
        else:
            real_length = None

        scale_eff_size = float(positional[1]) if len(positional) > 1 else 1
        cutoff_freq = float(positional[2]) if len(positional) > 2 else 60
        if len(positional) > 3:
            raise ValueError
    except (ValueError, IndexError):
        print_usage()
        sys.exit(1)

    directory = os.getcwd()

    processor = NEMD_Processor(directory, real_length, scale_eff_size, cutoff_freq, _save=save, _save_data=save_data)
    processor.process()
    # python plt_nemd.py [real_length] [scale_eff_size] [cutoff_freq] [--save] [--save-data]
"""
@Author   : Xin Wu
@Contact  : xinwuchn97@gmial.com
@Remark   : Post-processing script for NEMD (Non-Equilibrium Molecular Dynamics) thermal conductivity calculations
"""

from pylab import *
import pandas as pd
import numpy as np
from ase.io import read
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
    print("Usage: gpumdkit -plt nemd [real_length] [scale_eff_size] [cutoff_freq] [save]")
    print("Params:")
    print("  real_length   : Real length of heat tranfer zone in nm (set to 'Auto', with auto-calculation)")
    print("  scale_eff_size: Optional, Scale factor for effective cross-sectional area (default: 1)")
    print("                   • For 3D bulk systems: use 1")
    print("                   • For low-dimensional systems with vacuum layer: S_box / S_eff")
    print("                     - S_box: box area perpendicular to heat transfer direction")
    print("                     - S_eff: real or effective area of the system")
    print("  cutoff_freq   : Optional, Cutoff frequency for SHC calculation in THz (default: 60)")
    print("  save          : Optional, save the plot as 'nemd.png'")
    print("  !!! Note !!!  : If no SHC data, set [scale_eff_size] and [cutoff_freq] to any number as placeholders when using 'save'.")

class NEMD_Processor:
    def __init__(self, _directory, _real_length=None, _scale_eff_size=1, _cutoff_freq=60, _scale_vacuum=1):
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
        """
        self.directory = _directory
        self.scale_vacuum = _scale_vacuum
        self.scale_eff_size = _scale_eff_size
        self.real_length = _real_length
        self.cutoff_freq = _cutoff_freq
        self.path = {
            'run': self.directory + '/run.in',
            'compute': self.directory + '/compute.out',
            'model': self.directory + '/model.xyz'
        }
        self.has_shc = os.path.exists(self.directory + '/shc.out')

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
        self.path['shc'] = self.directory + "/shc.out"

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

        # Calculate g(omega) from jw for NEMD
        model = read(self.path['model'])
        if grouping_th == 0:
            group_index = 'group'
        elif grouping_th == 1:
            raise ValueError("Please check your grouping for SHC method, it seems not the Default setting !!!")
        group = model.get_array(group_index)
        part_ratio = np.sum(group == group_shc_th) / group.size

        vol = model.get_volume() * part_ratio / self.scale_vacuum
        convert = 1.602176634e7  # ev*A/ps/THz * 1/A^3 *1/K * A ==> MW/m/K/THz
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
            values = [np.trapz(Reformed_SHC_data[col][:, i], dx=Reformed_SHC_data["nu"][0, 0]) for i in range(N_repeat)]
            Reformed_SHC_data['Results'][f"{key}_ave"] = np.mean(values)
            Reformed_SHC_data['Results'][f"{key}_std"] = np.std(values) / np.sqrt(N_repeat)

        np.savez(f'{self.directory}/data_shc.npz', **Reformed_SHC_data)
        return Reformed_SHC_data

    def process(self):
        """Main processing function for NEMD method"""
        raw_NEMD_data = np.loadtxt(self.path['compute'])
        Reformed_NEMD_data = {}

        # Get parameters from run.in
        with open(self.path['run'], 'r') as file:
            found_nemd = False
            for line in file:
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
        group = model.get_array("group")
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
        k = G * Length * 1e-3  # W/mK
        Reformed_NEMD_data['G'] = np.hstack((G, G.mean(axis=1, keepdims=True),
                                             G.std(axis=1, keepdims=True) / sqrt(N_repeat)))  # MW/m^2/K
        Reformed_NEMD_data['k'] = np.hstack((k, k.mean(axis=1, keepdims=True),
                                             k.std(axis=1, keepdims=True) / sqrt(N_repeat)))  # W/m/K

        np.savez(f'{self.directory}/data_nemd.npz', **Reformed_NEMD_data)

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
        self._plot_results(Reformed_NEMD_data, Reformed_SHC_data, res_s,
                           Time_upper, N_nemd_data, N_repeat, N_temp_group,
                           xticks_length, Length)

    def _print_nemd_results(self, Reformed_NEMD_data, Length):
        """Print NEMD thermal conductivity results"""
        print("\n" + "=" * 70)
        print("NEMD Thermal Conductivity Results")
        print("=" * 70)
        print(f"\nEffective length: {Length:.4f} nm")
        # print(f"Scale_vacuum: {self.scale_vacuum}")
        print(f"Scale_eff_size: {self.scale_eff_size}\n")
        print(f"Thermal conductance G = {Reformed_NEMD_data['G'][0, -2]:.6f} ± {Reformed_NEMD_data['G'][0, -1]:.6f} MW/m²K")
        print(f"Thermal conductivity κ = {Reformed_NEMD_data['k'][0, -2]:.6f} ± {Reformed_NEMD_data['k'][0, -1]:.6f} W/mK")
        print("=" * 70)

    def _print_shc_results(self, results):
        """Print SHC spectral thermal conductance results"""
        print("\n" + "=" * 70)
        print("SHC Spectral Thermal Conductance Results")
        print("=" * 70)
        print(f"Cutoff frequency: {self.cutoff_freq} THz\n")
        print(f"G_in  (integrated) = {results['in_ave']:.6f} ± {results['in_std']:.6f} MW/m²K")
        print(f"G_out (integrated) = {results['out_ave']:.6f} ± {results['out_std']:.6f} MW/m²K")
        print(f"G_tot (integrated) = {results['tot_ave']:.6f} ± {results['tot_std']:.6f} MW/m²K")
        print("=" * 70 + "\n")

    def _plot_results(self, Reformed_NEMD_data, Reformed_SHC_data, res_s,
                      Time_upper, N_nemd_data, N_repeat, N_temp_group,
                      xticks_length, Length):
        """Visualize NEMD and SHC results"""

        if self.has_shc:
            figure(figsize=(10, 8))
            # (a) Temperature profile
            subplot2grid((2, 5), (0, 0), colspan=3)
        else:
            figure(figsize=(10, 4))
            subplot(1, 2, 1)

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
        text(0.95, 0.9, f'G={Reformed_NEMD_data["G"][0, -2]:.3f}±{Reformed_NEMD_data["G"][0, -1]:.2f} MW/m$^2$/K',
             ha='right', va='top', transform=plt.gca().transAxes)
        text(0.95, 0.8, f'κ={Reformed_NEMD_data["k"][0, -2]:.3f}±{Reformed_NEMD_data["k"][0, -1]:.3f} W/mK',
             ha='right', va='top', transform=plt.gca().transAxes)
        text(0.08, 0.08, f'Effective length={Length:.2f} nm', ha='left', va='bottom',
             transform=plt.gca().transAxes)
        xlabel("Length (nm)")
        ylabel("Temperature (K)")
        title("(a) Temperature profile")

        # (b) Thermostat energy
        if self.has_shc:
            subplot2grid((2, 5), (0, 3), colspan=2)
        else:
            subplot(1, 2, 2)

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

        if self.has_shc:
            # (c) Force-virial correlation function
            subplot2grid((2, 4), (1, 0), colspan=1)
            set_fig_properties([gca()])
            plot(Reformed_SHC_data['t'][:, -2], Reformed_SHC_data['Kt'][:, -2] / Reformed_SHC_data['L'], lw=2)
            ylabel('K (eV/ps)')
            xlabel('Correlation time (ps)')
            title('(c) K$_{tot}$(t)')

            # (d) SHC spectral conductance
            subplot2grid((2, 4), (1, 1), colspan=3)
            set_fig_properties([gca()])
            plot(Reformed_SHC_data['nu'][:, -2], Reformed_SHC_data['k_g_wi'][:, -2], linewidth=2, color='C1')
            plt.fill_between(Reformed_SHC_data['nu'][:, -2],
                             Reformed_SHC_data['k_g_wi'][:, -2] - Reformed_SHC_data['k_g_wi'][:, -1],
                             Reformed_SHC_data['k_g_wi'][:, -2] + Reformed_SHC_data['k_g_wi'][:, -1],
                             facecolor='C1', alpha=0.3)
            plot(Reformed_SHC_data['nu'][:, -2], Reformed_SHC_data['k_g_wo'][:, -2], linewidth=2, color='C2')
            plt.fill_between(Reformed_SHC_data['nu'][:, -2],
                             Reformed_SHC_data['k_g_wo'][:, -2] - Reformed_SHC_data['k_g_wo'][:, -1],
                             Reformed_SHC_data['k_g_wo'][:, -2] + Reformed_SHC_data['k_g_wo'][:, -1],
                             facecolor='C2', alpha=0.3)
            plot(Reformed_SHC_data['nu'][:, -2], Reformed_SHC_data['k_g_wt'][:, -2], linewidth=2, color='C0')
            plt.fill_between(Reformed_SHC_data['nu'][:, -2],
                             Reformed_SHC_data['k_g_wt'][:, -2] - Reformed_SHC_data['k_g_wt'][:, -1],
                             Reformed_SHC_data['k_g_wt'][:, -2] + Reformed_SHC_data['k_g_wt'][:, -1],
                             facecolor='C0', alpha=0.3)
            text(0.45, 0.9, f"$G_{{in}}$ = {res_s['in_ave']:.3f} ± {res_s['in_std']:.2f} MW/m$^2$K",
                 ha='left', va='top', transform=plt.gca().transAxes, color="C1")
            text(0.45, 0.8, f"$G_{{out}}$ = {res_s['out_ave']:.3f} ± {res_s['out_std']:.2f} MW/m$^2$K",
                 ha='left', va='top', transform=plt.gca().transAxes, color="C2")
            text(0.45, 0.7, f"$G_{{tot}}$ = {res_s['tot_ave']:.3f} ± {res_s['tot_std']:.2f} MW/m$^2$K",
                 ha='left', va='top', transform=plt.gca().transAxes, color="C0")
            xlim(0, self.cutoff_freq)
            ylabel(r'$g$($\omega$) (MW/m$^2$/K/THz)')
            xlabel(r'$\nu$ (THz)')
            title('(d) Spectral thermal conductance')

            plt.subplots_adjust(wspace=1, hspace=0.3)
        else:
            tight_layout()

        if len(sys.argv) > 4 and sys.argv[4] == 'save':
            savefig('nemd.png', dpi=300, bbox_inches='tight')
            show()
        else:
            show()



if __name__ == "__main__":

    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help', 'help']:
        print_usage()
        sys.exit(0)

    try:
        if len(sys.argv) > 1:
            if sys.argv[1].lower() == "auto":
                real_length = None
            else:
                real_length = float(sys.argv[1])
                if real_length <= 0:
                    raise ValueError
        else:
            real_length = None

        scale_eff_size = float(sys.argv[2]) if len(sys.argv) > 2 else 1
        cutoff_freq = float(sys.argv[3]) if len(sys.argv) > 3 else 60
        if len(sys.argv) > 4 and sys.argv[4] != 'save':
            raise ValueError
        if len(sys.argv) > 5 and sys.argv[5] != 'save_data':
            raise ValueError
    except (ValueError, IndexError):
        print_usage()
        sys.exit(1)

    directory = os.getcwd()

    processor = NEMD_Processor(directory, real_length, scale_eff_size, cutoff_freq)
    processor.process()
    # python plt_nemd.py [real_length] [scale_eff_size] [cutoff_freq] [save]
"""
@Author   : Ziyang Wang
@Email    : m15566605404@163.com
@Remark   : Post-processing script for VAC (Velocity Autocorrelation) and PDOS (Phonon Density of States)
            Calculates normalized VAC, PDOS, and Heat Capacity (Cv) from GPUMD outputs.
"""

import pandas as pd
from pylab import *
import numpy as np
import os
import sys
from scipy.integrate import cumulative_trapezoid
from ase.io import read

# Figure Properties
aw, lw, fs = 1.2, 2, 12
matplotlib.rc('font', size=fs)
matplotlib.rc('axes', linewidth=aw)

def set_fig_properties(ax_list, tl=4, tw=1.2, tlm=4):
    """Set figure properties for axes to match previous style"""
    for ax in ax_list:
        ax.tick_params(which='both', length=tl, width=tw, direction='in', right=True, top=True)
        ax.tick_params(which='minor', length=tlm)

trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

def print_usage():
    print("Usage: python plt_vac_dos.py [save]")
    print("  save : Optional. If provided, saves figures (.png) and data (.xlsx).")

class VAC_DOS_Processor:
    def __init__(self, directory):
        self.directory = directory
        self.files = {
            'mvac': os.path.join(directory, 'mvac.out'),
            'dos': os.path.join(directory, 'dos.out'),
            'run': os.path.join(directory, 'run.in'),
            'model': os.path.join(directory, 'model.xyz')
        }
        
        # Check files
        if not os.path.exists(self.files['mvac']) or not os.path.exists(self.files['dos']):
            print("[Error] 'mvac.out' or 'dos.out' not found in current directory.")
            sys.exit(1)

    def get_simulation_params(self):
        """Parse N (correlation steps) from run.in and num_atoms from model.xyz"""
        N = None
        num_atoms = None

        # 1. Get N from run.in
        if os.path.exists(self.files['run']):
            with open(self.files['run'], 'r') as f:
                for line in f:
                    # Look for 'compute_dos' or 'compute_vac'
                    # Format: compute_dos interval N sample_interval ...
                    if 'compute_dos' in line and not line.strip().startswith('#'):
                        parts = line.split()
                        try:
                            N = int(parts[2])
                        except IndexError:
                            pass
        
        if N is None:
            print("[Warning] Could not parse N from run.in. Using default N=500 (from original script).")
            N = 500

        # 2. Get num_atoms from model.xyz
        if os.path.exists(self.files['model']):
            try:
                atoms = read(self.files['model'])
                num_atoms = len(atoms)
                print(f"[Info] Auto-detected num_atoms = {num_atoms} from model.xyz")
            except Exception as e:
                print(f"[Warning] Failed to read model.xyz: {e}")
        
        if num_atoms is None:
            # Fallback to the value in your original script if file missing
            print("[Warning] model.xyz not found. Using default num_atoms = 1152.")
            num_atoms = 1152

        return N, num_atoms

    def process(self, save_mode=False):
        # 1. Initialization
        N, num_atoms = self.get_simulation_params()
        
        # Load Data
        raw_mvac = np.loadtxt(self.files['mvac'])
        raw_dos = np.loadtxt(self.files['dos'])

        # Determine number of runs (M)
        # Total rows = N * M
        if len(raw_mvac) % N != 0:
            print(f"[Warning] mvac.out length ({len(raw_mvac)}) is not divisible by N ({N}). Checking logical consistency...")
        
        M = len(raw_mvac) // N
        print(f"[Info] Data loaded. Steps (N)={N}, Independent runs (M)={M}")

        # 2. Reshape and Average (Ensemble Average)
        # Reshape to (N, M) then mean over axis 1 (columns)
        # Note: GPUMD outputs are typically just concatenated.
        
        # Extract Time and Frequency (Assume same for all runs, take first N)
        t = raw_mvac[:N, 0]             # Correlation time (ps)
        nu = raw_dos[:N, 0] / (2*np.pi) # omega -> nu (THz)

        # Process VAC (Cols: t, vac_x, vac_y, vac_z)
        vac_x = raw_mvac[:, 1].reshape(N, M, order='F').mean(axis=1)
        vac_y = raw_mvac[:, 2].reshape(N, M, order='F').mean(axis=1)
        vac_z = raw_mvac[:, 3].reshape(N, M, order='F').mean(axis=1)

        # Process DOS (Cols: omega, dos_x, dos_y, dos_z)
        dos_x = raw_dos[:, 1].reshape(N, M, order='F').mean(axis=1)
        dos_y = raw_dos[:, 2].reshape(N, M, order='F').mean(axis=1)
        dos_z = raw_dos[:, 3].reshape(N, M, order='F').mean(axis=1)

        # 3. Normalization Check
        # PDOS should integrate to num_atoms
        norm_x = trap(dos_x, x=nu) / num_atoms
        norm_y = trap(dos_y, x=nu) / num_atoms
        norm_z = trap(dos_z, x=nu) / num_atoms
        
        print("\n" + "="*50)
        print("PDOS Normalization Check (Should be ~1.0)")
        print("="*50)
        print(f"X: {norm_x:.4f}")
        print(f"Y: {norm_y:.4f}")
        print(f"Z: {norm_z:.4f}")
        print("="*50 + "\n")

        # 4. Heat Capacity Calculation
        print("[Info] Calculating Heat Capacity...")
        dT = 100
        NT = 50
        temperatures = np.arange(1, NT + 1) * dT # 100, 200, ... 5000
        
        # Constants
        k_B = 1.380649e-23 # J/K
        h = 6.62607015e-34 # J s
        THz_to_Hz = 1.0e12

        # Prepare arrays for broadcasting
        # T: (NT, 1), nu: (1, N)
        T_grid = temperatures[:, np.newaxis]
        nu_grid = nu[np.newaxis, :] * THz_to_Hz # Convert to Hz
        
        # x = h*nu / k_B*T
        # Avoid division by zero if T=0 (but we start at 100)
        x = (h * nu_grid) / (k_B * T_grid)
        
        # Modal Heat Capacity (Einstein term), in units of k_B
        # Formula: x^2 * exp(x) / (exp(x) - 1)^2
        # Handle x close to 0 to avoid singularity (though nu=0 usually dos=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            modal_cv = (x**2 * np.exp(x)) / ((np.exp(x) - 1)**2)
        modal_cv[np.isnan(modal_cv)] = 1.0 # Limit x->0 is 1 (classical limit)
        
        # Integrate: int(DOS(nu) * Cv(nu) dnu) / num_atoms
        # Result is array of shape (NT,)
        C_x = trap(dos_x * modal_cv, x=nu, axis=1) / num_atoms
        C_y = trap(dos_y * modal_cv, x=nu, axis=1) / num_atoms
        C_z = trap(dos_z * modal_cv, x=nu, axis=1) / num_atoms

        # 5. Plotting
        self._plot_vac_pdos(t, nu, vac_x, vac_y, vac_z, dos_x, dos_y, dos_z, save_mode)
        self._plot_heat_capacity(temperatures, C_x, C_y, C_z, save_mode)

        # 6. Save Data to Excel
        if save_mode:
            self._save_excel(nu, dos_x, dos_y, dos_z, temperatures, C_x, C_y, C_z)

    def _plot_vac_pdos(self, t, nu, vx, vy, vz, dx, dy, dz, save_mode):
        """Plot Figure 1: VAC and PDOS components and averages"""
        fig = figure(figsize=(10, 8))
        
        # (a) VAC components normalized
        ax1 = subplot(2, 2, 1)
        plot(t, vx/vx[0], 'C0-', lw=lw, label='x')
        plot(t, vy/vy[0], 'C1--', lw=lw, label='y')
        plot(t, vz/vz[0], 'C2-.', lw=lw, label='z')
        xlabel('Correlation time (ps)')
        ylabel('VAC (Normalized)')
        #xlim(0, 3)
        legend(loc='upper right', frameon=False)
        title('(a)')
        set_fig_properties([ax1])

        # (b) PDOS components
        ax2 = subplot(2, 2, 2)
        plot(nu, dx, 'C0-', lw=lw, label='x')
        plot(nu, dy, 'C1--', lw=lw, label='y')
        plot(nu, dz, 'C2-.', lw=lw, label='z')
        xlabel(r'$\nu$ (THz)')
        ylabel('PDOS (1/THz)')
        #xlim(0, 25)
        #ylim(0, 300)
        legend(loc='upper right', frameon=False)
        title('(b)')
        set_fig_properties([ax2])

        # (c) Average VAC
        ax3 = subplot(2, 2, 3)
        vac_tot = (vx + vy + vz)
        plot(t, vac_tot/vac_tot[0], 'k-', lw=lw)
        xlabel('Correlation time (ps)')
        ylabel('VAC (Normalized)')
        #xlim(0, 3)
        title('(c)')
        set_fig_properties([ax3])

        # (d) Average PDOS
        ax4 = subplot(2, 2, 4)
        dos_ave = (dx + dy + dz) / 3
        plot(nu, dos_ave, 'k-', lw=lw)
        xlabel(r'$\nu$ (THz)')
        ylabel('PDOS (1/THz)')
        #xlim(0, 25)
        #ylim(0, 300)
        title('(d)')
        set_fig_properties([ax4])

        tight_layout()
        
        if save_mode:
            savefig('vac_pdos.png', dpi=300, bbox_inches='tight')
            print("[Info] Plot saved: vac_pdos.png")
        else:
            show()

    def _plot_heat_capacity(self, T, cx, cy, cz, save_mode):
        """Plot Figure 2: Heat Capacity"""
        fig = figure(figsize=(5, 4))
        ax = gca()
        
        plot(T, cx, 'd-', color='C0', lw=lw, label='x', markersize=4)
        plot(T, cy, 's-', color='C1', lw=lw, label='y', markersize=4)
        plot(T, cz, 'o-', color='C2', lw=lw, label='z', markersize=4)
        
        xlabel(r'Temperature (K)')
        ylabel(r'Heat capacity ($k_{\rm B}$/atom)')
        xlim(0, 5100)
        ylim(0, 1.1)
        legend(frameon=False, loc='lower right')
        set_fig_properties([ax])
        
        tight_layout()
        
        if save_mode:
            savefig('heat_capacity.png', dpi=300, bbox_inches='tight')
            print("[Info] Plot saved: heat_capacity.png")
        else:
            show()

    def _save_excel(self, nu, dx, dy, dz, T, cx, cy, cz):
        """Save computed data to Excel"""
        filename = 'dos_data.xlsx'
        print(f"[Info] Saving data to {filename} ...")
        
        ave_dos = (dx + dy + dz) / 3
        
        with pd.ExcelWriter(filename) as writer:
            # Sheet 1: DOS Data
            df_dos = pd.DataFrame({
                'nu (THz)': nu,
                'dos_x (1/THz)': dx,
                'dos_y (1/THz)': dy,
                'dos_z (1/THz)': dz,
                'ave_dos (1/THz)': ave_dos
            })
            df_dos.to_excel(writer, sheet_name='PDOS_Data', index=False)
            
            # Sheet 2: Heat Capacity Data
            df_cv = pd.DataFrame({
                'Temperature (K)': T,
                'Cx (kB/atom)': cx,
                'Cy (kB/atom)': cy,
                'Cz (kB/atom)': cz,
                'Ctot_ave (kB/atom)': (cx+cy+cz)/3
            })
            df_cv.to_excel(writer, sheet_name='Heat_Capacity', index=False)
        
        print("[Info] Data save complete.")

if __name__ == "__main__":
    save_mode = False
    if len(sys.argv) > 1:
        if sys.argv[1] == 'save':
            save_mode = True
        else:
            print_usage()
            sys.exit(0)

    directory = os.getcwd()
    processor = VAC_DOS_Processor(directory)
    processor.process(save_mode)
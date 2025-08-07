import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext, LogLocator
from matplotlib.gridspec import GridSpec
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# Configuration
FONT_SIZE = 10
BINS = 200
COLORMAP = 'viridis' 
# preferred colormaps: 'rainbow', 'plasma', 'inferno', 'magma', 'cividis', 'viridis'

# Load data
energy_data = np.loadtxt('energy_train.out')
force_data = np.loadtxt('force_train.out')
stress_data = np.loadtxt('stress_train.out')

# === Metric Function ===
def compute_metrics(y_true, y_pred):
    rmse = np.sqrt(mean_squared_error(y_true, y_pred)) * 1000
    mae = mean_absolute_error(y_true, y_pred) * 1000
    r2 = r2_score(y_true, y_pred)
    return rmse, mae, r2

# === Plot Function ===
def parity_density_plot(ax, x, y, xlabel, ylabel, unit_str, cmap, bins):
    xmin, xmax = np.min([x, y]), np.max([x, y])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(xmin, xmax)

    hb = ax.hexbin(x, y, gridsize=bins, cmap=cmap,
                   norm=LogNorm(), mincnt=1, linewidths=0)

    ax.plot([xmin, xmax], [xmin, xmax], color='gray', linestyle='-', linewidth=1)
    ax.set_xlabel(xlabel, fontsize=FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=FONT_SIZE)
    ax.tick_params(labelsize=FONT_SIZE)

    # Metrics
    rmse, mae, r2 = compute_metrics(y, x)
    ax.text(0.05, 0.95,
            f"RMSE = {rmse:.1f} meV {unit_str}\n"
            f"MAE = {mae:.1f} meV {unit_str}\n"
            r"$R^2$" + f" = {r2:.5f}",
            transform=ax.transAxes,
            fontsize=FONT_SIZE,
            va='top', ha='left')

    return hb

# === Layout Setup ===
fig = plt.figure(figsize=(12, 4.2), dpi=120)
gs = GridSpec(2, 3, height_ratios=[4, 0.2], hspace=0.35)
cmap = plt.get_cmap(COLORMAP)

# === Plot Creation Function ===
def create_plot(ax, cb_ax, x, y, xlabel, ylabel, unit_str, cmap, bins):
    hb = parity_density_plot(ax, x, y, xlabel, ylabel, unit_str, cmap, bins)
    cb = fig.colorbar(hb, cax=cb_ax, orientation='horizontal')
    cb.set_label("Data density", fontsize=FONT_SIZE)
    cb.ax.tick_params(labelsize=FONT_SIZE)
    cb.locator = LogLocator(base=10.0)  
    cb.formatter = LogFormatterMathtext(base=10, labelOnlyBase=True)  
    cb.update_ticks()

# === Energy ===
ax1 = fig.add_subplot(gs[0, 0])
cb1_ax = fig.add_subplot(gs[1, 0])
create_plot(ax1, cb1_ax, energy_data[:, 1], energy_data[:, 0],
            "DFT energy (eV atom$^{-1}$)",
            "NEP energy (eV atom$^{-1}$)",
            "atom$^{-1}$", cmap, BINS)

# === Force ===
ax2 = fig.add_subplot(gs[0, 1])
cb2_ax = fig.add_subplot(gs[1, 1])
create_plot(ax2, cb2_ax, force_data[:, 3:6].reshape(-1), force_data[:, 0:3].reshape(-1),
            "DFT force (eV Å$^{-1}$)",
            "NEP force (eV Å$^{-1}$)",
            "Å$^{-1}$", cmap, BINS)

# === Stress ===
ax3 = fig.add_subplot(gs[0, 2])
cb3_ax = fig.add_subplot(gs[1, 2])
create_plot(ax3, cb3_ax, stress_data[:, 6:12].reshape(-1), stress_data[:, 0:6].reshape(-1),
            "DFT stress (GPa)",
            "NEP stress (GPa)",
            "GPa", cmap, BINS)

plt.subplots_adjust(top=0.968, bottom=0.122, left=0.073, right=0.983, hspace=0.2, wspace=0.286)

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('parity_plot_density.png', dpi=300)
else:
    # Check if the current backend is non-interactive
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'parity_plot_density.png'.")
        plt.savefig('parity_plot_density.png', dpi=300)
    else:
        plt.show()        
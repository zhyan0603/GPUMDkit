#!/usr/bin/env python3
"""
Plot GPUMD viscosity.out: stress autocorrelation and viscosity components.
https://gpumd.org/dev/gpumd/output_files/viscosity_out.html#viscosity-out

Columns: time, S_diag(3), S_off(6), eta_diag(3), eta_off(6).
Derived: η_L = (1/3)(η_xx+η_yy+η_zz), η_S = (1/3)(η_xy+η_xz+η_yz), η_B = η_L - (4/3)η_S.

Usage:
    python plot_viscosity.py [viscosity.out [viscosity.png]]
    python plot_viscosity.py --tex   # use external LaTeX if available
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import Axes

# Column layout (0-based): time, then S and eta components
TIME = 0
S_DIAG = slice(1, 4)    # S_xx, S_yy, S_zz
S_OFF = slice(4, 10)    # S_xy, S_xz, S_yz, S_yx, S_zx, S_zy
ETA_DIAG = slice(10, 13)
ETA_OFF = slice(13, 19)

# LaTeX labels
LABELS_S_DIAG = [r"$S_{xx}$", r"$S_{yy}$", r"$S_{zz}$"]
LABELS_S_OFF = [r"$S_{xy}$", r"$S_{xz}$", r"$S_{yz}$", r"$S_{yx}$", r"$S_{zx}$", r"$S_{zy}$"]
LABELS_ETA_DIAG = [r"$\eta_{xx}$", r"$\eta_{yy}$", r"$\eta_{zz}$"]
LABELS_ETA_OFF = [r"$\eta_{xy}$", r"$\eta_{xz}$", r"$\eta_{yz}$", r"$\eta_{yx}$", r"$\eta_{zx}$", r"$\eta_{zy}$"]

COLORS = ["#1f77b4", "#2ca02c", "#d62728", "#ff7f0e", "#9467bd", "#8c564b"]


def load_viscosity(path: Path) -> np.ndarray:
    """Load viscosity.out; return array of shape (n_steps, 20)."""
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def derived_viscosities(data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (η_L, η_S, η_B) from eta diagonal and off-diagonal components."""
    eta_xx, eta_yy, eta_zz = data[:, 10], data[:, 11], data[:, 12]
    eta_xy, eta_xz, eta_yz = data[:, 13], data[:, 14], data[:, 15]
    eta_L = (eta_xx + eta_yy + eta_zz) / 3
    eta_S = (eta_xy + eta_xz + eta_yz) / 3
    eta_B = eta_L - (4 / 3) * eta_S
    return eta_L, eta_S, eta_B


def plot_panel(ax: Axes, t: np.ndarray, curves: list[np.ndarray], labels: list[str], ylabel: str, title: str) -> None:
    """Plot multiple curves on ax with consistent styling."""
    for y, label, color in zip(curves, labels, COLORS):
        ax.plot(t, y, lw=1.5, color=color, label=label)
    ax.set_xlabel(r"Correlation time (ps)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(style="scientific", axis="y", scilimits=(-2, 2))


def main() -> None:
    args = sys.argv[1:]
    use_tex = "--tex" in args
    if use_tex:
        args = [a for a in args if a != "--tex"]

    # Input: default to viscosity.out in the current working directory,
    # or use the provided path.
    if len(args) < 1:
        input_path = Path.cwd() / "viscosity.out"
    else:
        input_path = Path(args[0])

    # Output: always save alongside the input file with .png extension,
    # so viscosity.png shares the same path as viscosity.out.
    output_path = input_path.with_suffix(".png")

    if not input_path.exists():
        print(f"Missing: {input_path}", file=sys.stderr)
        sys.exit(1)

    plt.rcParams["text.usetex"] = use_tex
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.size"] = 10

    data = load_viscosity(input_path)
    t = data[:, TIME]

    fig, axes = plt.subplots(3, 2, figsize=(10, 10), sharex=True)
    axes = axes.flatten()

    # 1) S diagonal (3)
    plot_panel(
        axes[0],
        t,
        [data[:, i] for i in range(1, 4)],
        LABELS_S_DIAG,
        r"$\langle S_{\alpha\alpha}(0)\,S_{\alpha\alpha}(t)\rangle$ (eV${}^2$)",
        r"$S$ diagonal (3)",
    )

    # 2) S non-diagonal (6)
    plot_panel(
        axes[1],
        t,
        [data[:, i] for i in range(4, 10)],
        LABELS_S_OFF,
        r"$\langle S_{\alpha\beta}(0)\,S_{\alpha\beta}(t)\rangle$ (eV${}^2$)",
        r"$S$ non-diagonal (6)",
    )

    # 3) η diagonal (3)
    plot_panel(
        axes[2],
        t,
        [data[:, i] for i in range(10, 13)],
        LABELS_ETA_DIAG,
        r"$\eta_{\alpha\alpha}$ (Pa$\cdot$s)",
        r"$\eta$ diagonal (3)",
    )

    # 4) η non-diagonal (6)
    plot_panel(
        axes[3],
        t,
        [data[:, i] for i in range(13, 19)],
        LABELS_ETA_OFF,
        r"$\eta_{\alpha\beta}$ (Pa$\cdot$s)",
        r"$\eta$ non-diagonal (6)",
    )

    # 5) η_L, η_S, η_B
    eta_L, eta_S, eta_B = derived_viscosities(data)
    plot_panel(
        axes[4],
        t,
        [eta_L, eta_S, eta_B],
        [r"$\eta_{\mathrm{L}}$", r"$\eta_{\mathrm{S}}$", r"$\eta_{\mathrm{B}}$"],
        r"$\eta$ (Pa$\cdot$s)",
        r"$\eta_{\mathrm{L}}$, $\eta_{\mathrm{S}}$, $\eta_{\mathrm{B}}$",
    )

    axes[5].set_visible(False)
    fig.suptitle(r"GPUMD viscosity.out", fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()

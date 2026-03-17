#!/usr/bin/env python3
"""
Plot GPUMD cohesive.out: isotropic scaling factor vs cohesive energy.

File format:
- column 1: scaling factor (dimensionless)
- column 2: potential energy (eV) after minimization

Usage:
  conda activate pymatgen_env
  python plot_cohesive.py                # uses ./cohesive.out -> cohesive.png
  python plot_cohesive.py path/to/cohesive.out
  python plot_cohesive.py in.out out.png
"""

import os
import sys
from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def configure_backend() -> str:
    """
    Choose a sensible matplotlib backend.

    - If running without a DISPLAY (e.g. on a cluster), force Agg.
    - Otherwise respect the current/default backend.
    """
    backend = matplotlib.get_backend()
    if not os.environ.get("DISPLAY") and "agg" not in str(backend).lower():
        matplotlib.use("Agg", force=True)
        backend = matplotlib.get_backend()
    # Optional: uncomment if you want to see which backend is used
    # print(f\"Using matplotlib backend: {backend}\")
    return str(backend)

_BACKEND = configure_backend()


def load_cohesive(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError(f"{path} must have at least 2 columns")
    scale = data[:, 0]
    energy = data[:, 1]
    return scale, energy


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    args = sys.argv[1:]

    if len(args) == 0:
        in_path = script_dir / "cohesive.out"
        out_path = script_dir / "cohesive.png"
    elif len(args) == 1:
        in_path = Path(args[0])
        out_path = in_path.with_suffix(".png")
    else:
        in_path = Path(args[0])
        out_path = Path(args[1])

    if not in_path.exists():
        print(f"Missing input file: {in_path}", file=sys.stderr)
        sys.exit(1)

    # LaTeX-style math via mathtext (no external TeX required)
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.size"] = 11

    scale, energy = load_cohesive(in_path)

    # Fit a parabola E(a) ~ a*s^2 + b*s + c near the minimum
    # and also mark the true discrete minimum from the data.
    x = scale
    y = energy
    x0_theory = 1.0
    mask = np.isfinite(x) & np.isfinite(y)
    x_fit = x[mask]
    y_fit = y[mask]
    # Prefer points near 1.0 if available.
    if x_fit.size >= 5:
        # Keep the central 50% of points by |x-1|, heuristically focusing around the minimum.
        order = np.argsort(np.abs(x_fit - x0_theory))
        k = max(5, x_fit.size // 2)
        idx = order[:k]
        x_fit = x_fit[idx]
        y_fit = y_fit[idx]

    # True minimum from data
    x_min_true = None
    y_min_true = None
    if x_fit.size > 0:
        j = np.argmin(y_fit)
        x_min_true = x_fit[j]
        y_min_true = y_fit[j]

    x_min_fit = None
    y_min_fit = None
    if x_fit.size >= 3:
        coeffs = np.polyfit(x_fit, y_fit, deg=2)
        a, b, c = coeffs
        if a != 0.0:
            x_min_fit = -b / (2.0 * a)
            y_min_fit = np.polyval(coeffs, x_min_fit)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(scale, energy, marker="o", ms=4, lw=1.5, color="#1f77b4", label="data")

    # Plot fitted parabola for visual guidance.
    if x_fit.size >= 3 and x_min_fit is not None:
        x_dense = np.linspace(x.min(), x.max(), 200)
        y_dense = np.polyval(coeffs, x_dense)
        ax.plot(x_dense, y_dense, color="#ff7f0e", lw=1.2, ls="--", label="parabolic fit")

        # Vertical line at theoretical minimum (1.0) and fitted minimum.
        ax.axvline(x0_theory, color="gray", lw=1.0, ls=":", label=r"$a=1.0$")
        ax.axvline(x_min_fit, color="#2ca02c", lw=1.0, ls="--", label=r"$a_{\min}$")

        # Mark fitted minimum.
        ax.scatter([x_min_fit], [y_min_fit], color="#d62728", zorder=3)
        ax.annotate(
            f"a_min={x_min_fit:.4f}",
            xy=(x_min_fit, y_min_fit),
            xytext=(5, -15),
            textcoords="offset points",
            fontsize=9,
        )

    # Mark true discrete minimum (may differ slightly from fit)
    if x_min_true is not None and y_min_true is not None:
        ax.scatter([x_min_true], [y_min_true], color="#9467bd", zorder=3)
        ax.annotate(
            f"a_min,data={x_min_true:.4f}",
            xy=(x_min_true, y_min_true),
            xytext=(5, 10),
            textcoords="offset points",
            fontsize=9,
        )

    ax.set_xlabel(r"Isotropic scaling factor")
    ax.set_ylabel(r"Potential energy $E$ (eV)")
    ax.set_title(r"GPUMD cohesive energy vs scaling")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()

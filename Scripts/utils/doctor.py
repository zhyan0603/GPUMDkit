"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     doctor.py
Category:   Utilities
Purpose:    Check the GPUMDkit Python environment and report available common
            and feature-specific packages without importing them.
Usage:      gpumdkit.sh -doctor
            python3 doctor.py
Arguments:
  None
Output:
  Terminal environment report. No files are created or modified.
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-07-11
=============================================================================
"""

import importlib.metadata
import importlib.util
import os
import platform
import sys


COMMON_PACKAGES = (
    ("numpy", "numpy", "numpy"),
    ("scipy", "scipy", "scipy"),
    ("matplotlib", "matplotlib", "matplotlib"),
    ("ase", "ase", "ase"),
    ("pandas", "pandas", "pandas"),
    ("pymatgen", "pymatgen", "pymatgen"),
)

OPTIONAL_PACKAGES = (
    ("dpdata", "dpdata", "dpdata", "DeepMD data conversion and perturbation"),
    ("calorine", "calorine", "calorine", "NEP calculations and model utilities"),
    ("NepTrain", "NepTrain", "NepTrain", "NepTrain structure sampling"),
    ("ovito", "ovito", "ovito", "OVITO-based conversion and analysis"),
    ("seaborn", "seaborn", "seaborn", "Selected plots"),
    ("scikit-learn", "sklearn", "scikit-learn", "PCA and parity-density plots"),
    ("umap-learn", "umap", "umap-learn", "UMAP descriptor plots"),
    ("tqdm", "tqdm", "tqdm", "Progress display in selected calculations"),
    ("ferrodispcalc", "ferrodispcalc", "ferrodispcalc", "Polar-material analysis"),
)


def package_version(distribution):
    """Return an installed distribution version or an empty string."""
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return ""


def package_status(module, distribution):
    """Return whether a module is discoverable and its distribution version."""
    try:
        available = importlib.util.find_spec(module) is not None
    except (ImportError, ModuleNotFoundError, ValueError):
        available = False
    return available, package_version(distribution) if available else ""


def print_package(name, module, distribution, note=""):
    """Print one aligned package status row."""
    available, version = package_status(module, distribution)
    status = "[OK]" if available else "[MISS]"
    detail = version or note
    print(f"   {status:<6} {name:<16} {detail}")


def main():
    """Print the GPUMDkit environment report."""
    gpumdkit_path = os.environ.get("GPUMDkit_path", "")
    path_ok = bool(gpumdkit_path and os.path.isdir(gpumdkit_path))

    print(" GPUMDkit environment check")
    print("")
    print(" Environment")
    print(f"   {'[OK]' if path_ok else '[MISS]':<6} {'GPUMDkit_path':<16} {gpumdkit_path or 'not set'}")
    print(f"   {'[OK]':<6} {'Python':<16} {platform.python_version()}")
    bash_version = os.environ.get("GPUMDKIT_BASH_VERSION", "unknown")
    print(f"   {'[OK]':<6} {'Bash':<16} {bash_version}")
    print("")
    print(" Common packages")
    for name, module, distribution in COMMON_PACKAGES:
        print_package(name, module, distribution)
    print("")
    print(" Optional packages")
    for name, module, distribution, note in OPTIONAL_PACKAGES:
        print_package(name, module, distribution, note)
    print("")
    print(" Packages are not required by every GPUMDkit function.")
    print(" Install a missing package when the corresponding function is needed.")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] in ("-h", "--help"):
            print(" Usage: gpumdkit.sh -doctor")
            print("    or: python3 doctor.py")
            print("")
            print(" Checks the GPUMDkit path, Python and Bash versions, and Python packages.")
            sys.exit(0)
        print(" Error: -doctor does not accept arguments.")
        sys.exit(1)
    main()

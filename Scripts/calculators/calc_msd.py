"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     calc_msd.py
Category:   Calculator Scripts
Purpose:    Compute directional MSD (x, y, z) from an extxyz trajectory via
            the Wiener-Khinchin FFT algorithm for a target element.
            For wrapped trajectories in a fixed cell, periodic jumps are removed
            with the minimum-image convention along periodic directions.
Usage:      gpumdkit.sh -calc msd <extxyz_file> <element_symbol> <dt_fs> [max_corr_steps]
            python calc_msd.py <extxyz_file> <element_symbol> <dt_fs> [max_corr_steps]
Arguments:
  extxyz_file       Path to the input extxyz trajectory file
  element_symbol    Chemical symbol (e.g., Li, O, Na)
  dt_fs             Time step between consecutive frames (fs)
  max_corr_steps    Max correlation lag steps (optional)
Output:
  msd.out  (4-column: Time/ps, MSD_x, MSD_y, MSD_z)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-11-10
=============================================================================
"""

import sys
import os
import time
import numpy as np


# ==============================================================================
# Argument parsing
# ==============================================================================

def parse_args():
    """
    Parse and validate positional command-line arguments.

    Returns
    -------
    extxyz_file    : str
    element_symbol : str
    dt_fs          : float
    max_corr_steps : int or None
    """
    args = sys.argv[1:]

    if len(args) == 0 or args[0] in ("-h", "--help"):
        print(" Usage: gpumdkit.sh -calc msd <extxyz_file> <element_symbol> <dt_fs> [max_corr_steps]")
        print("    or: python calc_msd.py <extxyz_file> <element_symbol> <dt_fs> [max_corr_steps]")
        print("")
        print(" Arguments:")
        print("   extxyz_file       Path to the input extxyz trajectory file")
        print("   element_symbol    Chemical symbol (e.g., Li, O, Na)")
        print("   dt_fs             Time step between consecutive frames (fs)")
        print("   max_corr_steps    Max correlation lag steps (optional)")
        print("   Default: use unwrapped_position if present, otherwise unwrap pos by minimum image.")
        print("")
        print(" Example: gpumdkit.sh -calc msd dump.xyz Li 1.0 500")
        print("")
        sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

    if len(args) not in (3, 4):
        print(" Error: Expected 3 or 4 positional arguments, got %d." % len(args))
        print(" Usage: gpumdkit.sh -calc msd <extxyz_file> <element_symbol> <dt_fs> [max_corr_steps]")
        print("    or: python calc_msd.py <extxyz_file> <element_symbol> <dt_fs> [max_corr_steps]")
        sys.exit(1)

    extxyz_file    = args[0]
    element_symbol = args[1]
    dt_str         = args[2]
    max_corr_str   = args[3] if len(args) == 4 else None

    # Validate file
    if not os.path.isfile(extxyz_file):
        print(" Error: File not found: '%s'" % extxyz_file)
        sys.exit(1)

    if not extxyz_file.lower().endswith((".xyz", ".extxyz")):
        print(" Warning: '%s' does not have a recognised .xyz / .extxyz extension. "
              " Proceeding anyway." % extxyz_file)

    # Validate element symbol
    if not element_symbol.isalpha() or not element_symbol[0].isupper():
        print(" Warning: '%s' may not be a valid element symbol. "
              " Expected format: e.g. 'Li', 'O', 'Na'." % element_symbol)

    # Validate dt
    try:
        dt_fs = float(dt_str)
    except ValueError:
        print(" Error: dt_fs must be a number, got: '%s'" % dt_str)
        sys.exit(1)
    if not np.isfinite(dt_fs) or dt_fs <= 0:
        print(" Error: dt_fs must be positive, got: %s" % dt_fs)
        sys.exit(1)

    # Validate max_corr_steps (optional)
    max_corr_steps = None
    if max_corr_str is not None:
        try:
            max_corr_steps = int(max_corr_str)
        except ValueError:
            print(" Error: max_corr_steps must be a positive integer, got: '%s'" % max_corr_str)
            sys.exit(1)
        if max_corr_steps <= 0:
            print(" Error: max_corr_steps must be a positive integer, got: %d" % max_corr_steps)
            sys.exit(1)

    return extxyz_file, element_symbol, dt_fs, max_corr_steps


# ==============================================================================
# Trajectory I/O
# ==============================================================================

def read_extxyz_file(filename, element_symbol):
    """
    Read coordinates of the specified element from an extxyz file.

    Parameters
    ----------
    filename       : str
    element_symbol : str

    Returns
    -------
    all_positions : ndarray, shape (n_frames, n_atoms, 3)
    cells         : ndarray, shape (n_frames, 3, 3)
    pbc           : ndarray, shape (n_frames, 3)
    position_key  : str
    """
    try:
        from ase.io import read
    except ImportError:
        print(" Error: ASE is required but not installed.")
        print("        Install it with: pip install ase")
        sys.exit(1)

    try:
        atoms_list = read(filename, index=":")
    except Exception as exc:
        print(" Error: Failed to read '%s': %s" % (filename, exc))
        sys.exit(1)

    if not atoms_list:
        print(" Error: No frames found in '%s'." % filename)
        sys.exit(1)

    position_key = "positions"
    if all("unwrapped_position" in atoms.arrays for atoms in atoms_list):
        position_key = "unwrapped_position"

    all_positions = []
    cells = []
    pbc = []
    n_target = None

    for frame_idx, atoms in enumerate(atoms_list):
        cells.append(atoms.cell.array)
        pbc.append(atoms.pbc)

        symbols = atoms.get_chemical_symbols()
        element_indices = [i for i, s in enumerate(symbols) if s == element_symbol]

        if not element_indices:
            print(" Error: Element '%s' not found in frame %d. "
                  " Check the element symbol and trajectory file."
                  % (element_symbol, frame_idx))
            sys.exit(1)

        if n_target is None:
            n_target = len(element_indices)
        elif len(element_indices) != n_target:
            print(" Error: frame %d contains %d %s atoms, but frame 0 contains %d."
                  % (frame_idx, len(element_indices), element_symbol, n_target))
            print(" The number of target atoms must be the same in every frame.")
            sys.exit(1)

        all_positions.append(atoms.arrays[position_key][element_indices])

    all_positions = np.array(all_positions)   # shape: (n_frames, n_atoms, 3)
    cells = np.array(cells)                   # shape: (n_frames, 3, 3)
    return all_positions, cells, np.asarray(pbc), position_key


def unwrap_coordinates(positions, cell, pbc):
    """Accumulate the shortest periodic image of each consecutive displacement."""
    from ase.geometry import find_mic

    delta = np.diff(positions, axis=0)
    vectors = delta.reshape(-1, 3)
    for start in range(0, len(vectors), 8192):
        block = vectors[start:start + 8192]
        block[:] = find_mic(block, cell, pbc=pbc)[0]

    unwrapped = np.empty_like(positions)
    unwrapped[0] = positions[0]
    unwrapped[1:] = positions[0] + np.cumsum(delta, axis=0)
    return unwrapped


# ==============================================================================
# MSD calculation
# ==============================================================================

def msd_fft(positions):
    """
    Compute directional MSD per axis via the Wiener-Khinchin FFT algorithm.

    Uses the identity:
        MSD(dt) = 2*<r^2> - 2*C(dt)
    where C is the position autocorrelation function evaluated via FFT,
    giving O(N log N) complexity.

    Parameters
    ----------
    positions : ndarray, shape (n_frames, n_atoms, 3)

    Returns
    -------
    msd_xyz : ndarray, shape (n_frames - 2, 3)
        Columns correspond to MSD_x, MSD_y, MSD_z, each averaged over all
        atoms. Lag index starts from 1 (first row = dt * 1 frame).
    """
    axis_time = 0
    axis_atom = 1
    n_time    = positions.shape[axis_time]
    sq_pos    = np.square(positions)

    # Autocorrelation via FFT (Wiener-Khinchin theorem)
    fft_pos = np.fft.fft(positions, n=2 * n_time, axis=axis_time)
    S2 = np.fft.ifft(np.abs(fft_pos) ** 2, axis=axis_time).real
    S2 = S2[:n_time]                              # keep first half only

    # Direct-sum (S1) term
    sq_aug   = np.concatenate(
        [sq_pos, np.zeros((1,) + sq_pos.shape[1:])], axis=axis_time
    )                                             # shape: (n_time+1, n_atoms, 3)
    total    = 2.0 * np.sum(sq_aug, axis=axis_time, keepdims=True)
    forward  = np.concatenate(
        [np.zeros((1,) + sq_pos.shape[1:]), sq_pos], axis=axis_time
    )
    backward = np.flip(sq_aug, axis=axis_time)
    S1       = (total - np.cumsum(forward + backward, axis=axis_time))[:n_time]

    # Per-atom MSD: shape (n_time, n_atoms, 3)
    counts   = (n_time - np.arange(n_time)).reshape(-1, 1, 1)
    atom_msd = (S1 - 2.0 * S2) / counts

    # Drop lag=0 and the last lag (statistically unreliable)
    lag_range = np.arange(1, n_time - 1)
    atom_msd  = atom_msd[lag_range]              # shape: (n_time-2, n_atoms, 3)

    # Average over atoms, keep x/y/z separate -> shape: (n_time-2, 3)
    msd_xyz = atom_msd.mean(axis=axis_atom)
    return msd_xyz


# ==============================================================================
# Main workflow
# ==============================================================================

def calculate_msd_from_extxyz(extxyz_file, element_symbol, dt_fs, max_corr_steps):
    """
    Full MSD calculation pipeline:
      1. Read trajectory
      2. Unwrap coordinates (minimum image between consecutive frames)
      3. Compute directional MSD via FFT
      4. Truncate to max_corr_steps if specified
      5. Write msd.out (time, MSD_x, MSD_y, MSD_z)
    """
    print(" ")
    print(" "+ "=" * 60)
    print("  Trajectory file  : %s" % extxyz_file)
    print("  Target element   : %s" % element_symbol)
    print("  Time step (dt)   : %.4f fs" % dt_fs)
    if max_corr_steps is not None:
        print("  Max corr. steps  : %d  (%.4f ps)"
              % (max_corr_steps, max_corr_steps * dt_fs / 1000.0))
    else:
        print("  Max corr. steps  : all frames")
    print(" "+ "-" * 60)

    # ------------------------------------------------------------------
    # Read trajectory
    # ------------------------------------------------------------------
    print(" ")
    print(" Reading file: %s" % extxyz_file)
    all_positions, cells, pbc, position_key = read_extxyz_file(extxyz_file, element_symbol)
    n_frames, n_atoms, _ = all_positions.shape
    print("   Loaded %d frames, %d %s atoms per frame"
          % (n_frames, n_atoms, element_symbol))

    if n_frames < 3:
        print(" Error: At least 3 frames are required for MSD, "
              " but only %d found." % n_frames)
        sys.exit(1)

    if position_key == "unwrapped_position":
        print(" Coordinate mode: using the unwrapped_position column.")
    elif not pbc.any():
        print(" Coordinate mode: Cartesian coordinates as stored.")
    else:
        if np.any(pbc != pbc[0]) or not np.allclose(cells, cells[0], rtol=0, atol=1e-6):
            print(" Error: wrapped MSD requires a fixed cell and fixed periodic directions.")
            print(" For variable-cell trajectories, provide an unwrapped_position column.")
            sys.exit(1)
        periodic_cell = cells[0][pbc[0]]
        if np.linalg.matrix_rank(periodic_cell) != len(periodic_cell):
            print(" Error: invalid cell vectors for the periodic directions.")
            sys.exit(1)
        # Always unwrap: for coordinates that are already continuous, the
        # minimum-image step leaves them unchanged as long as consecutive
        # displacements are shorter than half the cell. Deciding from the
        # coordinate range instead would miss partially wrapped trajectories.
        print(" Coordinate mode: minimum-image unwrapping in periodic directions.")
        all_positions = unwrap_coordinates(all_positions, cells[0], pbc[0])

    print(" Calculating MSD (Wiener-Khinchin FFT)...")
    msd_xyz = msd_fft(all_positions)             # shape: (n_frames-2, 3)

    # Build lag-time axis: lag 1, 2, ..., n_frames-2
    n_lags  = msd_xyz.shape[0]
    time_ps = np.arange(1, n_lags + 1) * dt_fs / 1000.0   # fs -> ps

    # Truncate to max_corr_steps if requested
    if max_corr_steps is not None:
        if max_corr_steps > n_lags:
            print(" Warning: max_corr_steps (%d) exceeds available lag steps (%d). "
                  " Using all available steps." % (max_corr_steps, n_lags))
        else:
            time_ps = time_ps[:max_corr_steps]
            msd_xyz = msd_xyz[:max_corr_steps]

    print("   MSD computed for %d lag times" % len(time_ps))
    print("   Lag-time range  : %.6f ps  ->  %.6f ps"
          % (time_ps[0], time_ps[-1]))

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    print(" ")
    output_file = "msd.out"
    header = ("msd.out by GPUMDkit\n"
               "Element: %s  |  dt: %.4f fs  |  Frames: %d  |  Atoms: %d\n"
               "%12s  %14s  %14s  %14s"
               % (element_symbol, dt_fs, n_frames, n_atoms,
                  "Time(ps)", "MSD_x(Ang^2)", "MSD_y(Ang^2)", "MSD_z(Ang^2)"))

    np.savetxt(
        output_file,
        np.column_stack((time_ps, msd_xyz)),
        fmt="%14.6f",
        header=header,
    )

    print(" MSD calculation completed. Results saved to '%s'." % output_file)
    print(" You can use [gpumdkit.sh -plt msd] to visualize the results.")

# ==============================================================================
# Entry point
# ==============================================================================

if __name__ == "__main__":
    t0 = time.perf_counter()

    calculate_msd_from_extxyz(*parse_args())

    print(" Total elapsed time: %.2f s" % (time.perf_counter() - t0))
    print(" "+ "=" * 60)

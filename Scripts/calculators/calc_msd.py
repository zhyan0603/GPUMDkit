"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     calc_msd.py
Category:   Calculator Scripts
Purpose:    Compute directional MSD (x, y, z) from an extxyz trajectory via
            the Wiener-Khinchin FFT algorithm for a target element.
Usage:      python calc_msd.py <extxyz_file> <element_symbol> <dt_fs> [max_corr_steps]
Arguments:
  extxyz_file       Path to the input extxyz trajectory file
  element_symbol    Chemical symbol (e.g., Li, O, Na)
  dt_fs             Time step between consecutive frames (fs)
  max_corr_steps    Max correlation lag steps (optional)
Output:
  msd.out  (4-column: Time/ps, MSD_x, MSD_y, MSD_z)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-05-16
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
        print(__doc__)
        sys.exit(0)

    if len(args) not in (3, 4):
        print(" Error: Expected 3 or 4 positional arguments, got %d." % len(args))
        print("        Run with -h for usage information.")
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
    if dt_fs <= 0:
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
    cell_info     : ndarray, shape (3, 3), or None
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

    all_positions = []
    cell_info = None

    for frame_idx, atoms in enumerate(atoms_list):
        if frame_idx == 0:
            cell_info = np.array(atoms.get_cell())

        symbols = atoms.get_chemical_symbols()
        element_indices = [i for i, s in enumerate(symbols) if s == element_symbol]

        if not element_indices:
            print(" Error: Element '%s' not found in frame %d. "
                  " Check the element symbol and trajectory file."
                  % (element_symbol, frame_idx))
            sys.exit(1)

        all_positions.append(atoms.positions[element_indices])

    all_positions = np.array(all_positions)   # shape: (n_frames, n_atoms, 3)
    return all_positions, cell_info


# ==============================================================================
# Coordinate wrapping utilities
# ==============================================================================

def check_wrapped_coordinates(positions, cell):
    """
    Heuristically check whether Cartesian coordinates are PBC-wrapped.

    Returns True if any coordinate is negative, exceeds 1.5x the cell
    length, or if the per-dimension coordinate span exceeds the cell length.
    """
    if cell is None or positions.size == 0:
        return False

    cell_lengths = np.array([np.linalg.norm(cell[i]) for i in range(3)])

    for dim in range(3):
        col = positions[:, :, dim]
        if np.any(col < 0) or np.any(col > cell_lengths[dim] * 1.5):
            return True
        if col.max() - col.min() > cell_lengths[dim] * 1.5:
            return True

    return False


def unwrap_coordinates(positions, cell):
    """
    Unwrap PBC-wrapped coordinates using the minimum-image convention.

    For each atom and each consecutive frame pair, if the displacement
    along any axis exceeds half the cell length, all subsequent frames
    are shifted by one full cell vector to remove the jump.

    Parameters
    ----------
    positions : ndarray, shape (n_frames, n_atoms, 3)
    cell      : ndarray, shape (3, 3)

    Returns
    -------
    ndarray, shape (n_frames, n_atoms, 3)
    """
    if positions.shape[0] < 2 or cell is None:
        return positions

    cell_lengths = np.array([np.linalg.norm(cell[i]) for i in range(3)])
    unwrapped = np.copy(positions)

    for atom_idx in range(positions.shape[1]):
        for frame in range(1, positions.shape[0]):
            displacement = positions[frame, atom_idx] - positions[frame - 1, atom_idx]
            for dim in range(3):
                half = cell_lengths[dim] / 2.0
                if displacement[dim] > half:
                    unwrapped[frame:, atom_idx, dim] -= cell_lengths[dim]
                elif displacement[dim] < -half:
                    unwrapped[frame:, atom_idx, dim] += cell_lengths[dim]

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
      2. Check / unwrap coordinates
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
    all_positions, cell = read_extxyz_file(extxyz_file, element_symbol)
    n_frames, n_atoms, _ = all_positions.shape
    print("   Loaded %d frames, %d %s atoms per frame"
          % (n_frames, n_atoms, element_symbol))

    if cell is not None:
        cell_lengths = np.array([np.linalg.norm(cell[i]) for i in range(3)])
        print("   Cell lengths (a, b, c) : %.4f  %.4f  %.4f  Ang"
              % (cell_lengths[0], cell_lengths[1], cell_lengths[2]))
    else:
        print(" Warning: No cell information found. PBC unwrapping will be skipped.")

    # ------------------------------------------------------------------
    # Coordinate wrapping check
    # ------------------------------------------------------------------
    print(" ")
    # print(" "+ "-" * 60)
    if cell is not None:
        is_wrapped = check_wrapped_coordinates(all_positions, cell)
        if is_wrapped:
            print(" Coordinate status : Wrapped (PBC detected)")
            print(" Unwrapping coordinates...")
            all_positions = unwrap_coordinates(all_positions, cell)
            print(" Coordinate unwrapping completed.")
        else:
            print(" Coordinate status : Unwrapped (no action needed)")
    else:
        print(" Coordinate status : Skipped (no cell information)")

    # ------------------------------------------------------------------
    # MSD calculation
    # ------------------------------------------------------------------
    print(" ")
    # print(" "+ "-" * 60)
    if n_frames < 3:
        print(" Error: At least 3 frames are required for MSD, "
              " but only %d found." % n_frames)
        sys.exit(1)

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
    # print(" "+ "-" * 60)
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

    extxyz_file, element_symbol, dt_fs, max_corr_steps = parse_args()
    calculate_msd_from_extxyz(extxyz_file, element_symbol, dt_fs, max_corr_steps)

    print(" Total elapsed time: %.2f s" % (time.perf_counter() - t0))
    print(" "+ "=" * 60)
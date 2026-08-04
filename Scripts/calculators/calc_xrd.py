#!/usr/bin/env python3
"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     calc_xrd.py
Category:   Calculator
Purpose:    Calculate and average LAMMPS-compatible XRD intensities from
            extended XYZ trajectories.
Usage:      gpumdkit.sh -> 4) Calculators -> 413) X-ray diffraction (XRD)
            python3 Scripts/calculators/calc_xrd.py
Arguments:  None. All values are collected through the interactive prompts.
Interactive inputs:
  The input trajectory, output file, wavelength, two-theta range, bin count,
  atom selection, and optional CPU worker count are entered interactively by
  the script.
Output:
  <output>   Histogrammed XRD data with intensity and normalized intensity.
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-07-31
=============================================================================
"""

from __future__ import annotations

from collections import deque
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
import math
import os
from pathlib import Path
import re
import sys
import time
from typing import Iterator

import numpy as np

PI = math.pi
MAX_PACKED_PHASE_BYTES = 128 * 1024**2

# Atomic scattering factors copied from LAMMPS src/DIFFRACTION/compute_xrd_consts.h.
# The table is embedded so this script does not depend on a LAMMPS installation.
ASF_LABELS = (
    "H", "He1-", "He", "Li", "Li1+", "Be", "Be2+", "B",
    "C", "Cval", "N", "O", "O1-", "F", "F1-", "Ne",
    "Na", "Na1+", "Mg", "Mg2+", "Al", "Al3+", "Si", "Sival",
    "Si4+", "P", "S", "Cl", "Cl1-", "Ar", "K", "Ca",
    "Ca2+", "Sc", "Sc3+", "Ti", "Ti2+", "Ti3+", "Ti4+", "V",
    "V2+", "V3+", "V5+", "Cr", "Cr2+", "Cr3+", "Mn", "Mn2+",
    "Mn3+", "Mn4+", "Fe", "Fe2+", "Fe3+", "Co", "Co2+", "Co",
    "Ni", "Ni2+", "Ni3+", "Cu", "Cu1+", "Cu2+", "Zn", "Zn2+",
    "Ga", "Ga3+", "Ge", "Ge4+", "As", "Se", "Br", "Br1-",
    "Kr", "Rb", "Rb1+", "Sr", "Sr2+", "Y", "Y3+", "Zr",
    "Zr4+", "Nb", "Nb3+", "Nb5+", "Mo", "Mo3+", "Mo5+", "Mo6+",
    "Tc", "Ru", "Ru3+", "Ru4+", "Rh", "Rh3+", "Rh4+", "Pd",
    "Pd2+", "Pd4+", "Ag", "Ag1+", "Ag2+", "Cd", "Cd2+", "In",
    "In3+", "Sn", "Sn2+", "Sn4+", "Sb", "Sb3+", "Sb5+", "Te",
    "I", "I1-", "Xe", "Cs", "Cs1+", "Ba", "Ba2+", "La",
    "La3+", "Ce", "Ce3+", "Ce4+", "Pr", "Pr3+", "Pr4+", "Nd",
    "Nd3+", "Pm", "Pm3+", "Sm", "Sm3+", "Eu", "Eu2+", "Eu3+",
    "Gd", "Gd3+", "Tb", "Tb3+", "Dy", "Dy3+", "Ho", "Ho3+",
    "Er", "Er3+", "Tm", "Tm3+", "Yb", "Yb2+", "Yb3+", "Lu",
    "Lu3+", "Hf", "Hf4+", "Ta", "Ta5+", "W", "W6+", "Re",
    "Os", "Os4+", "Ir", "Ir3+", "Ir4+", "Pt", "Pt2+", "Pt4+",
    "Au", "Au1+", "Au3+", "Hg", "Hg1+", "Hg2+", "Tl", "Tl1+",
    "Tl3+", "Pb", "Pb2+", "Pb4+", "Bi", "Bi3+", "Bi5+", "Po",
    "At", "Rn", "Fr", "Ra", "Ra2+", "Ac", "Ac3+", "Th",
    "Th4+", "Pa", "U", "U3+", "U4+", "U6+", "Np", "Np3+",
    "Np4+", "Np6+", "Pu", "Pu3+", "Pu4+", "Pu6+", "Am", "Cm",
    "Bk", "Cf",
)
ASF_COEFFS = np.array((
    (0.489918, 20.6593, 0.262003, 7.74039, 0.196767, 49.5519, 0.049879, 2.20159, 0.001305),
    (0.897661, 53.1368, 0.565616, 15.187, 0.415815, 186.576, 0.116973, 3.56709, 0.002389),
    (0.8734, 9.1037, 0.6309, 3.3568, 0.3112, 22.9276, 0.178, 0.9821, 0.0064),
    (1.1282, 3.9546, 0.7508, 1.0524, 0.6175, 85.3905, 0.4653, 68.261, 0.0377),
    (0.6968, 4.6237, 0.7888, 1.9557, 0.3414, 0.6316, 0.7029, 0.542, 0.0167),
    (1.5919, 43.6427, 1.1278, 1.8623, 0.5391, 103.483, 0.7029, 0.542, 0.0385),
    (6.2603, 0.0027, 0.8849, 0.8313, 0.7993, 2.2758, 0.1647, 5.1146, -6.1092),
    (2.0545, 23.2185, 1.3326, 1.021, 1.0979, 60.3498, 0.7068, 0.1403, -0.1932),
    (2.31, 20.8439, 1.02, 10.2075, 1.5886, 0.5687, 0.865, 51.6512, 0.2156),
    (2.26069, 22.6907, 1.56165, 0.656665, 1.05075, 9.75618, 0.839259, 55.5949, 0.286977),
    (12.2126, 0.0057, 3.1322, 9.8933, 2.0125, 28.9975, 1.1663, 0.5826, -11.529),
    (3.0485, 13.2771, 2.2868, 5.7011, 1.5463, 0.3239, 0.867, 32.9089, 0.2508),
    (4.1916, 12.8573, 1.63969, 4.17236, 1.52673, -47.0179, 20.307, -0.01404, 21.9412),
    (3.5392, 10.2825, 2.6412, 4.2944, 1.517, 0.2615, 1.0243, 26.1476, 0.2776),
    (3.6322, 5.27756, 3.51057, 14.7353, 1.26064, 0.442258, 0.940706, 47.3437, 0.653396),
    (3.9553, 8.4042, 3.1125, 3.4262, 1.4546, 0.2306, 1.1251, 21.7184, 0.3515),
    (4.7626, 3.285, 3.1736, 8.8422, 1.2674, 0.3136, 1.1128, 129.424, 0.676),
    (3.2565, 2.6671, 3.9362, 6.1153, 1.3998, 0.2001, 1.0032, 14.039, 0.404),
    (5.4204, 2.8275, 2.1735, 79.2611, 1.2269, 0.3808, 2.3073, 7.1937, 0.8584),
    (3.4988, 2.1676, 3.8378, 4.7542, 1.3284, 0.185, 0.8497, 10.1411, 0.4853),
    (6.4202, 3.0387, 1.9002, 0.7426, 1.5936, 31.5472, 1.9646, 85.0886, 1.1151),
    (4.17448, 1.93816, 3.3876, 4.14553, 1.20296, 0.228753, 0.528137, 8.28524, 0.706786),
    (6.2915, 2.4386, 3.0353, 32.3337, 1.9891, 0.6785, 1.541, 81.6937, 1.1407),
    (5.66269, 2.6652, 3.07164, 38.6634, 2.62446, 0.916946, 1.3932, 93.5458, 1.24707),
    (4.43918, 1.64167, 3.20345, 3.43757, 1.19453, 0.2149, 0.41653, 6.65365, 0.746297),
    (6.4345, 1.9067, 4.1791, 27.157, 1.78, 0.526, 1.4908, 68.1645, 1.1149),
    (6.9053, 1.4679, 5.2034, 22.2151, 1.4379, 0.2536, 1.5863, 56.172, 0.8669),
    (11.4604, 0.0104, 7.1964, 1.1662, 6.2556, 18.5194, 1.6455, 47.7784, -9.5574),
    (18.2915, 0.0066, 7.2084, 1.1717, 6.5337, 19.5424, 2.3386, 60.4486, -16.378),
    (7.4845, 0.9072, 6.7723, 14.8407, 0.6539, 43.8983, 1.6442, 33.3929, 1.4445),
    (8.2186, 12.7949, 7.4398, 0.7748, 1.0519, 213.187, 0.8659, 41.6841, 1.4228),
    (8.6266, 10.4421, 7.3873, 0.6599, 1.5899, 85.7484, 1.0211, 178.437, 1.3751),
    (15.6348, -0.0074, 7.9518, 0.6089, 8.4372, 10.3116, 0.8537, 25.9905, -14.875),
    (9.189, 9.0213, 7.3679, 0.5729, 1.6409, 136.108, 1.468, 51.3531, 1.3329),
    (13.4008, 0.29854, 8.0273, 7.9629, 1.65943, -0.28604, 1.57936, 16.0662, -6.6667),
    (9.7595, 7.8508, 7.3558, 0.5, 1.6991, 35.6338, 1.9021, 116.105, 1.2807),
    (9.11423, 7.5243, 7.62174, 0.457585, 2.2793, 19.5361, 0.087899, 61.6558, 0.897155),
    (17.7344, 0.22061, 8.73816, 7.04716, 5.25691, -0.15762, 1.92134, 15.9768, -14.652),
    (19.5114, 0.178847, 8.23473, 6.67018, 2.01341, -0.29263, 1.5208, 12.9464, -13.28),
    (10.2971, 6.8657, 7.3511, 0.4385, 2.0703, 26.8938, 2.0571, 102.478, 1.2199),
    (10.106, 6.8818, 7.3541, 0.4409, 2.2884, 20.3004, 0.0223, 115.122, 1.2298),
    (9.43141, 6.39535, 7.7419, 0.383349, 2.15343, 15.1908, 0.016865, 63.969, 0.656565),
    (15.6887, 0.679003, 8.14208, 5.40135, 2.03081, 9.97278, -9.576, 0.940464, 1.7143),
    (10.6406, 6.1038, 7.3537, 0.392, 3.324, 20.2626, 1.4922, 98.7399, 1.1832),
    (9.54034, 5.66078, 7.7509, 0.344261, 3.58274, 13.3075, 0.509107, 32.4224, 0.616898),
    (9.6809, 5.59463, 7.81136, 0.334393, 2.87603, 12.8288, 0.113575, 32.8761, 0.518275),
    (11.2819, 5.3409, 7.3573, 0.3432, 3.0193, 17.8674, 2.2441, 83.7543, 1.0896),
    (10.8061, 5.2796, 7.362, 0.3435, 3.5268, 14.343, 0.2184, 41.3235, 1.0874),
    (9.84521, 4.91797, 7.87194, 0.294393, 3.56531, 10.8171, 0.323613, 24.1281, 0.393974),
    (9.96253, 4.8485, 7.97057, 0.283303, 2.76067, 10.4852, 0.054447, 27.573, 0.251877),
    (11.7695, 4.7611, 7.3573, 0.3072, 3.5222, 15.3535, 2.3045, 76.8805, 1.0369),
    (11.0424, 4.6538, 7.374, 0.3053, 4.1346, 12.0546, 0.4399, 31.2809, 1.0097),
    (11.1764, 4.6147, 7.3863, 0.3005, 3.3948, 11.6729, 0.0724, 38.5566, 0.9707),
    (12.2841, 4.2791, 7.3409, 0.2784, 4.0034, 13.5359, 2.3488, 71.1692, 1.0118),
    (11.2296, 4.1231, 7.3883, 0.2726, 4.7393, 10.2443, 0.7108, 25.6466, 0.9324),
    (10.338, 3.90969, 7.88173, 0.238668, 4.76795, 8.35583, 0.725591, 18.3491, 0.286667),
    (12.8376, 3.8785, 7.292, 0.2565, 4.4438, 12.1763, 2.38, 66.3421, 1.0341),
    (11.4166, 3.6766, 7.4005, 0.2449, 5.3442, 8.873, 0.9773, 22.1626, 0.8614),
    (10.7806, 3.5477, 7.75868, 0.22314, 5.22746, 7.64468, 0.847114, 16.9673, 0.386044),
    (13.338, 3.5828, 7.1676, 0.247, 5.6158, 11.3966, 1.6735, 64.8126, 1.191),
    (11.9475, 3.3669, 7.3573, 0.2274, 6.2455, 8.6625, 1.5578, 25.8487, 0.89),
    (11.8168, 3.37484, 7.11181, 0.244078, 5.78135, 7.9876, 1.14523, 19.897, 1.14431),
    (14.0743, 3.2655, 7.0318, 0.2333, 5.1652, 10.3163, 2.41, 58.7097, 1.3041),
    (11.9719, 2.9946, 7.3862, 0.2031, 6.4668, 7.0826, 1.394, 18.0995, 0.7807),
    (15.2354, 3.0669, 6.7006, 0.2412, 4.3591, 10.7805, 2.9623, 61.4135, 1.7189),
    (12.692, 2.81262, 6.69883, 0.22789, 6.06692, 6.36441, 1.0066, 14.4122, 1.53545),
    (16.0816, 2.8509, 6.3747, 0.2516, 3.7068, 11.4468, 3.683, 54.7625, 2.1313),
    (12.9172, 2.53718, 6.70003, 0.205855, 6.06791, 5.47913, 0.859041, 11.603, 1.45572),
    (16.6723, 2.6345, 6.0701, 0.2647, 3.4313, 12.9479, 4.2779, 47.7972, 2.531),
    (17.0006, 2.4098, 5.8196, 0.2726, 3.9731, 15.2372, 4.3543, 43.8163, 2.8409),
    (17.1789, 2.1723, 5.2358, 16.5796, 5.6377, 0.2609, 3.9851, 41.4328, 2.9557),
    (17.1718, 2.2059, 6.3338, 19.3345, 5.5754, 0.2871, 3.7272, 58.1535, 3.1776),
    (17.3555, 1.9384, 6.7286, 16.5623, 5.5493, 0.2261, 3.5375, 39.3972, 2.825),
    (17.1784, 1.7888, 9.6435, 17.3151, 5.1399, 0.2748, 1.5292, 164.934, 3.4873),
    (17.5816, 1.7139, 7.6598, 14.7957, 5.8981, 0.1603, 2.7817, 31.2087, 2.0782),
    (17.5663, 1.5564, 9.8184, 14.0988, 5.422, 0.1664, 2.6694, 132.376, 2.5064),
    (18.0874, 1.4907, 8.1373, 12.6963, 2.5654, -24.5651, 34.193, -0.0138, 41.4025),
    (17.776, 1.4029, 10.2946, 12.8006, 5.72629, 0.125599, 3.26588, 104.354, 1.91213),
    (17.9268, 1.35417, 9.1531, 11.2145, 1.76795, -22.6599, 33.108, -0.01319, 40.2602),
    (17.8765, 1.27618, 10.948, 11.916, 5.41732, 0.117622, 3.65721, 87.6627, 2.06929),
    (18.1668, 1.2148, 10.0562, 10.1483, 1.01118, 21.6054, -2.6479, -0.10276, 9.41454),
    (17.6142, 1.18865, 12.0144, 11.766, 4.04183, 0.204785, 3.53346, 69.7957, 3.75591),
    (19.8812, 0.019175, 18.0653, 1.13305, 11.0177, 10.1621, 1.94715, 28.3389, -12.912),
    (17.9163, 1.12446, 13.3417, 0.028781, 10.799, 9.28206, 0.337905, 25.7228, -6.3934),
    (3.7025, 0.2772, 17.2356, 1.0958, 12.8876, 11.004, 3.7429, 61.6584, 4.3875),
    (21.1664, 0.014734, 18.2017, 1.03031, 11.7423, 9.53659, 2.30951, 26.6307, -14.421),
    (21.0149, 0.014345, 18.0992, 1.02238, 11.4632, 8.78809, 0.740625, 23.3452, -14.316),
    (17.8871, 1.03649, 11.175, 8.48061, 6.57891, 0.058881, 0, 0, 0.344941),
    (19.1301, 0.864132, 11.0948, 8.14487, 4.64901, 21.5707, 2.71263, 86.8472, 5.40428),
    (19.2674, 0.80852, 12.9182, 8.43467, 4.86337, 24.7997, 1.56756, 94.2928, 5.37814),
    (18.5638, 0.847329, 13.2885, 8.37164, 9.32602, 0.017662, 3.00964, 22.887, -3.1892),
    (18.5003, 0.844582, 13.1787, 8.12534, 4.71304, 0.36495, 2.18535, 20.8504, 1.42357),
    (19.2957, 0.751536, 14.3501, 8.21758, 4.73425, 25.8749, 1.28918, 98.6062, 5.328),
    (18.8785, 0.764252, 14.1259, 7.84438, 3.32515, 21.2487, -6.1989, -0.01036, 11.8678),
    (18.8545, 0.760825, 13.9806, 7.62436, 2.53464, 19.3317, -5.6526, -0.0102, 11.2835),
    (19.3319, 0.698655, 15.5017, 7.98929, 5.29537, 25.2052, 0.605844, 76.8986, 5.26593),
    (19.1701, 0.696219, 15.2096, 7.55573, 4.32234, 22.5057, 0, 0, 5.2916),
    (19.2493, 0.683839, 14.79, 7.14833, 2.89289, 17.9144, -7.9492, 0.005127, 13.0174),
    (19.2808, 0.6446, 16.6885, 7.4726, 4.8045, 24.6605, 1.0463, 99.8156, 5.179),
    (19.1812, 0.646179, 15.9719, 7.19123, 5.27475, 21.7326, 0.357534, 66.1147, 5.21572),
    (19.1643, 0.645643, 16.2456, 7.18544, 4.3709, 21.4072, 0, 0, 5.21404),
    (19.2214, 0.5946, 17.6444, 6.9089, 4.461, 24.7008, 1.6029, 87.4825, 5.0694),
    (19.1514, 0.597922, 17.2535, 6.80639, 4.47128, 20.2521, 0, 0, 5.11937),
    (19.1624, 0.5476, 18.5596, 6.3776, 4.2948, 25.8499, 2.0396, 92.8029, 4.9391),
    (19.1045, 0.551522, 18.1108, 6.3247, 3.78897, 17.3595, 0, 0, 4.99635),
    (19.1889, 5.8303, 19.1005, 0.5031, 4.4585, 26.8909, 2.4663, 83.9571, 4.7821),
    (19.1094, 0.5036, 19.0548, 5.8378, 4.5648, 23.3752, 0.487, 62.2061, 4.7861),
    (18.9333, 5.764, 19.7131, 0.4655, 3.4182, 14.0049, 0.0193, -0.7583, 3.9182),
    (19.6418, 5.3034, 19.0455, 0.4607, 5.0371, 27.9074, 2.6827, 75.2825, 4.5909),
    (18.9755, 0.467196, 18.933, 5.22126, 5.10789, 19.5902, 0.288753, 55.5113, 4.69626),
    (19.8685, 5.44853, 19.0302, 0.467973, 2.41253, 14.1259, 0, 0, 4.69263),
    (19.9644, 4.81742, 19.0138, 0.420885, 6.14487, 28.5284, 2.5239, 70.8403, 4.352),
    (20.1472, 4.347, 18.9949, 0.3814, 7.5138, 27.766, 2.2735, 66.8776, 4.0712),
    (20.2332, 4.3579, 18.997, 0.3815, 7.8069, 29.5259, 2.8868, 84.9304, 4.0714),
    (20.2933, 3.9282, 19.0298, 0.344, 8.9767, 26.4659, 1.99, 64.2658, 3.7118),
    (20.3892, 3.569, 19.1062, 0.3107, 10.662, 24.3879, 1.4953, 213.904, 3.3352),
    (20.3524, 3.552, 19.1278, 0.3086, 10.2821, 23.7128, 0.9615, 59.4565, 3.2791),
    (20.3361, 3.216, 19.297, 0.2756, 10.888, 20.2073, 2.6959, 167.202, 2.7731),
    (20.1807, 3.21367, 19.1136, 0.28331, 10.9054, 20.0558, 0.77634, 51.746, 3.02902),
    (20.578, 2.94817, 19.599, 0.244475, 11.3727, 18.7726, 3.28719, 133.124, 2.14678),
    (20.2489, 2.9207, 19.3763, 0.250698, 11.6323, 17.8211, 0.336048, 54.9453, 2.4086),
    (21.1671, 2.81219, 19.7695, 0.226836, 11.8513, 17.6083, 3.33049, 127.113, 1.86264),
    (20.8036, 2.77691, 19.559, 0.23154, 11.9369, 16.5408, 0.612376, 43.1692, 2.09013),
    (20.3235, 2.65941, 19.8186, 0.21885, 12.1233, 15.7992, 0.144583, 62.2355, 1.5918),
    (22.044, 2.77393, 19.6697, 0.222087, 12.3856, 16.7669, 2.82428, 143.644, 2.0583),
    (21.3727, 2.6452, 19.7491, 0.214299, 12.1329, 15.323, 0.97518, 36.4065, 1.77132),
    (20.9413, 2.54467, 20.0539, 0.202481, 12.4668, 14.8137, 0.296689, 45.4643, 1.24285),
    (22.6845, 2.66248, 19.6847, 0.210628, 12.774, 15.885, 2.85137, 137.903, 1.98486),
    (21.961, 2.52722, 19.9339, 0.199237, 12.12, 14.1783, 1.51031, 30.8717, 1.47588),
    (23.3405, 2.5627, 19.6095, 0.202088, 13.1235, 15.1009, 2.87516, 132.721, 2.02876),
    (22.5527, 2.4174, 20.1108, 0.185769, 12.0671, 13.1275, 2.07492, 27.4491, 1.19499),
    (24.0042, 2.47274, 19.4258, 0.196451, 13.4396, 14.3996, 2.89604, 128.007, 2.20963),
    (23.1504, 2.31641, 20.2599, 0.174081, 11.9202, 12.1571, 2.71488, 24.8242, 0.954586),
    (24.6274, 2.3879, 19.0886, 0.1942, 13.7603, 13.7546, 2.9227, 123.174, 2.5745),
    (24.0063, 2.27783, 19.9504, 0.17353, 11.8034, 11.6096, 3.87243, 26.5156, 1.36389),
    (23.7497, 2.22258, 20.3745, 0.16394, 11.8509, 11.311, 3.26503, 22.9966, 0.759344),
    (25.0709, 2.25341, 19.0798, 0.181951, 13.8518, 12.9331, 3.54545, 101.398, 2.4196),
    (24.3466, 2.13553, 20.4208, 0.155525, 11.8708, 10.5782, 3.7149, 21.7029, 0.645089),
    (25.8976, 2.24256, 18.2185, 0.196143, 14.3167, 12.6648, 2.95354, 115.362, 3.58324),
    (24.9559, 2.05601, 20.3271, 0.149525, 12.2471, 10.0499, 3.773, 21.2773, 0.691967),
    (26.507, 2.1802, 17.6383, 0.202172, 14.5596, 12.1899, 2.96577, 111.874, 4.29728),
    (25.5395, 1.9804, 20.2861, 0.143384, 11.9812, 9.34972, 4.50073, 19.581, 0.68969),
    (26.9049, 2.07051, 17.294, 0.19794, 14.5583, 11.4407, 3.63837, 92.6566, 4.56796),
    (26.1296, 1.91072, 20.0994, 0.139358, 11.9788, 8.80018, 4.93676, 18.5908, 0.852795),
    (27.6563, 2.07356, 16.4285, 0.223545, 14.9779, 11.3604, 2.98233, 105.703, 5.92046),
    (26.722, 1.84659, 19.7748, 0.13729, 12.1506, 8.36225, 5.17379, 17.8974, 1.17613),
    (28.1819, 2.02859, 15.8851, 0.238849, 15.1542, 10.9975, 2.98706, 102.961, 1.63929),
    (27.3083, 1.78711, 19.332, 0.136974, 12.3339, 7.96778, 5.38348, 17.2922, 1.63929),
    (28.6641, 1.9889, 15.4345, 0.257119, 15.3087, 10.6647, 2.98963, 100.417, 7.56672),
    (28.1209, 1.78503, 17.6817, 0.15997, 13.3335, 8.18304, 5.14657, 20.39, 3.70983),
    (27.8917, 1.73272, 18.7614, 0.13879, 12.6072, 7.64412, 5.47647, 16.8153, 2.26001),
    (28.9476, 1.90182, 15.2208, 9.98519, 15.1, 0.261033, 3.71601, 84.3298, 7.97628),
    (28.4628, 1.68216, 18.121, 0.142292, 12.8429, 7.33727, 5.59415, 16.3535, 2.97573),
    (29.144, 1.83262, 15.1726, 9.5999, 14.7586, 0.275116, 4.30013, 72.029, 8.58154),
    (28.8131, 1.59136, 18.4601, 0.128903, 12.7285, 6.76232, 5.59927, 14.0366, 2.39699),
    (29.2024, 1.77333, 15.2293, 9.37046, 14.5135, 0.295977, 4.76492, 63.3644, 9.24354),
    (29.1587, 1.50711, 18.8407, 0.116741, 12.8268, 6.31524, 5.38695, 12.4244, 1.78555),
    (29.0818, 1.72029, 15.43, 9.2259, 14.4327, 0.321703, 5.11982, 57.056, 9.8875),
    (29.4936, 1.42755, 19.3763, 0.104621, 13.0544, 5.93667, 5.06412, 11.1972, 1.01074),
    (28.7621, 1.67191, 15.7189, 9.09227, 14.5564, 0.3505, 5.44174, 52.0861, 10.472),
    (28.1894, 1.62903, 16.155, 8.97948, 14.9305, 0.382661, 5.67589, 48.1647, 11.0005),
    (30.419, 1.37113, 15.2637, 6.84706, 14.7458, 0.165191, 5.06795, 18.003, 6.49804),
    (27.3049, 1.59279, 16.7296, 8.86553, 15.6115, 0.417916, 5.83377, 45.0011, 11.4722),
    (30.4156, 1.34323, 15.862, 7.10909, 13.6145, 0.204633, 5.82008, 20.3254, 8.27903),
    (30.7058, 1.30923, 15.5512, 6.71983, 14.2326, 0.167252, 5.53672, 17.4911, 6.96824),
    (27.0059, 1.51293, 17.7639, 8.81174, 15.7131, 0.424593, 5.7837, 38.6103, 11.6883),
    (29.8429, 1.32927, 16.7224, 7.38979, 13.2153, 0.263297, 6.35234, 22.9426, 9.85329),
    (30.9612, 1.24813, 15.9829, 6.60834, 13.7348, 0.16864, 5.92034, 16.9392, 7.39534),
    (16.8819, 0.4611, 18.5913, 8.6216, 25.5582, 1.4826, 5.86, 36.3956, 12.0658),
    (28.0109, 1.35321, 17.8204, 7.7395, 14.3359, 0.356752, 6.58077, 26.4043, 11.2299),
    (30.6886, 1.2199, 16.9029, 6.82872, 12.7801, 0.212867, 6.52354, 18.659, 9.0968),
    (20.6809, 0.545, 19.0417, 8.4484, 21.6575, 1.5729, 5.9676, 38.3246, 12.6089),
    (25.0853, 1.39507, 18.4973, 7.65105, 16.8883, 0.443378, 6.48216, 28.2262, 12.0205),
    (29.5641, 1.21152, 18.06, 7.05639, 12.8374, 0.284738, 6.89912, 20.7482, 10.6268),
    (27.5446, 0.65515, 19.1584, 8.70751, 15.538, 1.96347, 5.52593, 45.8149, 13.1746),
    (21.3985, 1.4711, 20.4723, 0.517394, 18.7478, 7.43463, 6.82847, 28.8482, 12.5258),
    (30.8695, 1.1008, 18.3481, 6.53852, 11.9328, 0.219074, 7.00574, 17.2114, 9.8027),
    (31.0617, 0.6902, 13.0637, 2.3576, 18.442, 8.618, 5.9696, 47.2579, 13.4118),
    (21.7886, 1.3366, 19.5682, 0.488383, 19.1406, 6.7727, 7.01107, 23.8132, 12.4734),
    (32.1244, 1.00566, 18.8003, 6.10926, 12.0175, 0.147041, 6.96886, 14.714, 8.08428),
    (33.3689, 0.704, 12.951, 2.9238, 16.5877, 8.7937, 6.4692, 48.0093, 13.5782),
    (21.8053, 1.2356, 19.5026, 6.24149, 19.1053, 0.469999, 7.10295, 20.3185, 12.4711),
    (33.5364, 0.91654, 25.0946, 0.39042, 19.2497, 5.71414, 6.91555, 12.8285, 6.7994),
    (34.6726, 0.700999, 15.4733, 3.55078, 13.1138, 9.55642, 7.02588, 47.0045, 13.677),
    (35.3163, 0.68587, 19.0211, 3.97458, 9.49887, 11.3824, 7.42518, 45.4715, 13.7108),
    (35.5631, 0.6631, 21.2816, 4.0691, 8.0037, 14.0422, 7.4433, 44.2473, 13.6905),
    (35.9299, 0.646453, 23.0547, 4.17619, 12.1439, 23.1052, 2.11253, 150.645, 13.7247),
    (35.763, 0.616341, 22.9064, 3.87135, 12.4739, 19.9887, 3.21097, 142.325, 13.6211),
    (35.215, 0.604909, 21.67, 3.5767, 7.91342, 12.601, 7.65078, 29.8436, 13.5431),
    (35.6597, 0.589092, 23.1032, 3.65155, 12.5977, 18.599, 4.08655, 117.02, 13.5266),
    (35.1736, 0.579689, 22.1112, 3.41437, 8.19216, 12.9187, 7.05545, 25.9443, 13.4637),
    (35.5645, 0.563359, 23.4219, 3.46204, 12.7473, 17.8309, 4.80703, 99.1722, 13.4314),
    (35.1007, 0.555054, 22.4418, 3.24498, 9.78554, 13.4661, 5.29444, 23.9533, 13.376),
    (35.8847, 0.547751, 23.2948, 3.41519, 14.1891, 16.9235, 4.17287, 105.251, 13.4287),
    (36.0228, 0.5293, 23.4128, 3.3253, 14.9491, 16.0927, 4.188, 100.613, 13.3966),
    (35.5747, 0.52048, 22.5259, 3.12293, 12.2165, 12.7148, 5.37073, 26.3394, 13.3092),
    (35.3715, 0.516598, 22.5326, 3.05053, 12.0291, 12.5723, 4.7984, 23.4582, 13.2671),
    (34.8509, 0.507079, 22.7584, 2.8903, 14.0099, 13.1767, 1.21457, 25.2017, 13.1665),
    (36.1874, 0.511929, 23.5964, 3.25396, 15.6402, 15.3622, 4.1855, 97.4908, 13.3573),
    (35.7074, 0.502322, 22.613, 3.03807, 12.9898, 12.1449, 5.43227, 25.4928, 13.2544),
    (35.5103, 0.498626, 22.5787, 2.96627, 12.7766, 11.9484, 4.92159, 22.7502, 13.2116),
    (35.0136, 0.48981, 22.7286, 2.81099, 14.3884, 12.33, 1.75669, 22.6581, 13.113),
    (36.5254, 0.499384, 23.8083, 3.26371, 16.7707, 14.9455, 3.47947, 105.98, 13.3812),
    (35.84, 0.484938, 22.7169, 2.96118, 13.5807, 11.5331, 5.66016, 24.3992, 13.1991),
    (35.6493, 0.481422, 22.646, 2.8902, 13.3595, 11.316, 5.18831, 21.8301, 13.1555),
    (35.1736, 0.473204, 22.7181, 2.73848, 14.7635, 11.553, 2.28678, 20.9303, 13.0582),
    (36.6706, 0.483629, 24.0992, 3.20647, 17.3415, 14.3136, 3.49331, 102.273, 13.3592),
    (36.6488, 0.465154, 24.4096, 3.08997, 17.399, 13.4346, 4.21665, 88.4834, 13.2887),
    (36.7881, 0.451018, 24.7736, 3.04619, 17.8919, 12.8946, 4.23284, 86.003, 13.2754),
    (36.9185, 0.437533, 25.1995, 3.00775, 18.3317, 12.4044, 4.24391, 83.7881, 13.2674),
), dtype=np.float64)
SCATTERING_FACTORS = {
    "h": ASF_COEFFS[0],
    "he1-": ASF_COEFFS[1],
    "he": ASF_COEFFS[2],
    "li": ASF_COEFFS[3],
    "li1+": ASF_COEFFS[4],
    "be": ASF_COEFFS[5],
    "be2+": ASF_COEFFS[6],
    "b": ASF_COEFFS[7],
    "c": ASF_COEFFS[8],
    "cval": ASF_COEFFS[9],
    "n": ASF_COEFFS[10],
    "o": ASF_COEFFS[11],
    "o1-": ASF_COEFFS[12],
    "f": ASF_COEFFS[13],
    "f1-": ASF_COEFFS[14],
    "ne": ASF_COEFFS[15],
    "na": ASF_COEFFS[16],
    "na1+": ASF_COEFFS[17],
    "mg": ASF_COEFFS[18],
    "mg2+": ASF_COEFFS[19],
    "al": ASF_COEFFS[20],
    "al3+": ASF_COEFFS[21],
    "si": ASF_COEFFS[22],
    "sival": ASF_COEFFS[23],
    "si4+": ASF_COEFFS[24],
    "p": ASF_COEFFS[25],
    "s": ASF_COEFFS[26],
    "cl": ASF_COEFFS[27],
    "cl1-": ASF_COEFFS[28],
    "ar": ASF_COEFFS[29],
    "k": ASF_COEFFS[30],
    "ca": ASF_COEFFS[31],
    "ca2+": ASF_COEFFS[32],
    "sc": ASF_COEFFS[33],
    "sc3+": ASF_COEFFS[34],
    "ti": ASF_COEFFS[35],
    "ti2+": ASF_COEFFS[36],
    "ti3+": ASF_COEFFS[37],
    "ti4+": ASF_COEFFS[38],
    "v": ASF_COEFFS[39],
    "v2+": ASF_COEFFS[40],
    "v3+": ASF_COEFFS[41],
    "v5+": ASF_COEFFS[42],
    "cr": ASF_COEFFS[43],
    "cr2+": ASF_COEFFS[44],
    "cr3+": ASF_COEFFS[45],
    "mn": ASF_COEFFS[46],
    "mn2+": ASF_COEFFS[47],
    "mn3+": ASF_COEFFS[48],
    "mn4+": ASF_COEFFS[49],
    "fe": ASF_COEFFS[50],
    "fe2+": ASF_COEFFS[51],
    "fe3+": ASF_COEFFS[52],
    "co": ASF_COEFFS[53],
    "co2+": ASF_COEFFS[54],
    "co": ASF_COEFFS[55],
    "ni": ASF_COEFFS[56],
    "ni2+": ASF_COEFFS[57],
    "ni3+": ASF_COEFFS[58],
    "cu": ASF_COEFFS[59],
    "cu1+": ASF_COEFFS[60],
    "cu2+": ASF_COEFFS[61],
    "zn": ASF_COEFFS[62],
    "zn2+": ASF_COEFFS[63],
    "ga": ASF_COEFFS[64],
    "ga3+": ASF_COEFFS[65],
    "ge": ASF_COEFFS[66],
    "ge4+": ASF_COEFFS[67],
    "as": ASF_COEFFS[68],
    "se": ASF_COEFFS[69],
    "br": ASF_COEFFS[70],
    "br1-": ASF_COEFFS[71],
    "kr": ASF_COEFFS[72],
    "rb": ASF_COEFFS[73],
    "rb1+": ASF_COEFFS[74],
    "sr": ASF_COEFFS[75],
    "sr2+": ASF_COEFFS[76],
    "y": ASF_COEFFS[77],
    "y3+": ASF_COEFFS[78],
    "zr": ASF_COEFFS[79],
    "zr4+": ASF_COEFFS[80],
    "nb": ASF_COEFFS[81],
    "nb3+": ASF_COEFFS[82],
    "nb5+": ASF_COEFFS[83],
    "mo": ASF_COEFFS[84],
    "mo3+": ASF_COEFFS[85],
    "mo5+": ASF_COEFFS[86],
    "mo6+": ASF_COEFFS[87],
    "tc": ASF_COEFFS[88],
    "ru": ASF_COEFFS[89],
    "ru3+": ASF_COEFFS[90],
    "ru4+": ASF_COEFFS[91],
    "rh": ASF_COEFFS[92],
    "rh3+": ASF_COEFFS[93],
    "rh4+": ASF_COEFFS[94],
    "pd": ASF_COEFFS[95],
    "pd2+": ASF_COEFFS[96],
    "pd4+": ASF_COEFFS[97],
    "ag": ASF_COEFFS[98],
    "ag1+": ASF_COEFFS[99],
    "ag2+": ASF_COEFFS[100],
    "cd": ASF_COEFFS[101],
    "cd2+": ASF_COEFFS[102],
    "in": ASF_COEFFS[103],
    "in3+": ASF_COEFFS[104],
    "sn": ASF_COEFFS[105],
    "sn2+": ASF_COEFFS[106],
    "sn4+": ASF_COEFFS[107],
    "sb": ASF_COEFFS[108],
    "sb3+": ASF_COEFFS[109],
    "sb5+": ASF_COEFFS[110],
    "te": ASF_COEFFS[111],
    "i": ASF_COEFFS[112],
    "i1-": ASF_COEFFS[113],
    "xe": ASF_COEFFS[114],
    "cs": ASF_COEFFS[115],
    "cs1+": ASF_COEFFS[116],
    "ba": ASF_COEFFS[117],
    "ba2+": ASF_COEFFS[118],
    "la": ASF_COEFFS[119],
    "la3+": ASF_COEFFS[120],
    "ce": ASF_COEFFS[121],
    "ce3+": ASF_COEFFS[122],
    "ce4+": ASF_COEFFS[123],
    "pr": ASF_COEFFS[124],
    "pr3+": ASF_COEFFS[125],
    "pr4+": ASF_COEFFS[126],
    "nd": ASF_COEFFS[127],
    "nd3+": ASF_COEFFS[128],
    "pm": ASF_COEFFS[129],
    "pm3+": ASF_COEFFS[130],
    "sm": ASF_COEFFS[131],
    "sm3+": ASF_COEFFS[132],
    "eu": ASF_COEFFS[133],
    "eu2+": ASF_COEFFS[134],
    "eu3+": ASF_COEFFS[135],
    "gd": ASF_COEFFS[136],
    "gd3+": ASF_COEFFS[137],
    "tb": ASF_COEFFS[138],
    "tb3+": ASF_COEFFS[139],
    "dy": ASF_COEFFS[140],
    "dy3+": ASF_COEFFS[141],
    "ho": ASF_COEFFS[142],
    "ho3+": ASF_COEFFS[143],
    "er": ASF_COEFFS[144],
    "er3+": ASF_COEFFS[145],
    "tm": ASF_COEFFS[146],
    "tm3+": ASF_COEFFS[147],
    "yb": ASF_COEFFS[148],
    "yb2+": ASF_COEFFS[149],
    "yb3+": ASF_COEFFS[150],
    "lu": ASF_COEFFS[151],
    "lu3+": ASF_COEFFS[152],
    "hf": ASF_COEFFS[153],
    "hf4+": ASF_COEFFS[154],
    "ta": ASF_COEFFS[155],
    "ta5+": ASF_COEFFS[156],
    "w": ASF_COEFFS[157],
    "w6+": ASF_COEFFS[158],
    "re": ASF_COEFFS[159],
    "os": ASF_COEFFS[160],
    "os4+": ASF_COEFFS[161],
    "ir": ASF_COEFFS[162],
    "ir3+": ASF_COEFFS[163],
    "ir4+": ASF_COEFFS[164],
    "pt": ASF_COEFFS[165],
    "pt2+": ASF_COEFFS[166],
    "pt4+": ASF_COEFFS[167],
    "au": ASF_COEFFS[168],
    "au1+": ASF_COEFFS[169],
    "au3+": ASF_COEFFS[170],
    "hg": ASF_COEFFS[171],
    "hg1+": ASF_COEFFS[172],
    "hg2+": ASF_COEFFS[173],
    "tl": ASF_COEFFS[174],
    "tl1+": ASF_COEFFS[175],
    "tl3+": ASF_COEFFS[176],
    "pb": ASF_COEFFS[177],
    "pb2+": ASF_COEFFS[178],
    "pb4+": ASF_COEFFS[179],
    "bi": ASF_COEFFS[180],
    "bi3+": ASF_COEFFS[181],
    "bi5+": ASF_COEFFS[182],
    "po": ASF_COEFFS[183],
    "at": ASF_COEFFS[184],
    "rn": ASF_COEFFS[185],
    "fr": ASF_COEFFS[186],
    "ra": ASF_COEFFS[187],
    "ra2+": ASF_COEFFS[188],
    "ac": ASF_COEFFS[189],
    "ac3+": ASF_COEFFS[190],
    "th": ASF_COEFFS[191],
    "th4+": ASF_COEFFS[192],
    "pa": ASF_COEFFS[193],
    "u": ASF_COEFFS[194],
    "u3+": ASF_COEFFS[195],
    "u4+": ASF_COEFFS[196],
    "u6+": ASF_COEFFS[197],
    "np": ASF_COEFFS[198],
    "np3+": ASF_COEFFS[199],
    "np4+": ASF_COEFFS[200],
    "np6+": ASF_COEFFS[201],
    "pu": ASF_COEFFS[202],
    "pu3+": ASF_COEFFS[203],
    "pu4+": ASF_COEFFS[204],
    "pu6+": ASF_COEFFS[205],
    "am": ASF_COEFFS[206],
    "cm": ASF_COEFFS[207],
    "bk": ASF_COEFFS[208],
    "cf": ASF_COEFFS[209],
}

def resolve_labels(labels: list[str], table: dict[str, np.ndarray]) -> np.ndarray:
    """Return ASF coefficients using LAMMPS case-insensitive label lookup."""
    missing = [label for label in labels if label.lower() not in table]
    if missing:
        available = ", ".join(sorted(table))
        missing_labels = ", ".join(missing)
        raise ValueError(
            "Unknown LAMMPS XRD scattering-factor label(s): "
            f"{missing_labels}. Available labels include: {available}"
        )
    return np.asarray([table[label.lower()] for label in labels], dtype=np.float64)


# Extended-XYZ reader
@dataclass(frozen=True)
class Frame:
    """One extended-XYZ structure."""

    natoms: int
    species: np.ndarray
    positions: np.ndarray
    cell: np.ndarray | None
    pbc: tuple[bool, bool, bool]
    comment: str


def _parse_comment(comment: str) -> tuple[np.ndarray | None, tuple[bool, bool, bool]]:
    """Extract the cell matrix and periodic flags from an extxyz comment."""
    lattice_match = re.search(r'Lattice="([^"]+)"', comment)
    cell = None
    if lattice_match:
        values = np.fromstring(lattice_match.group(1), sep=" ", dtype=np.float64)
        if values.size != 9:
            raise ValueError("Lattice must contain exactly 9 numbers")
        cell = values.reshape(3, 3)

    pbc_match = re.search(r'pbc="([^"]+)"', comment)
    if pbc_match:
        pbc_values = pbc_match.group(1).split()
        if len(pbc_values) != 3:
            raise ValueError("pbc must contain three flags")
        pbc = tuple(value.lower() in {"t", "true", "1"} for value in pbc_values)
    else:
        pbc = (True, True, True)
    return cell, pbc


def _read_frame(stream) -> Frame | None:
    """Read one extxyz frame from a binary stream."""
    first = stream.readline()
    while first and not first.strip():
        first = stream.readline()
    if not first:
        return None
    try:
        natoms = int(first.decode("utf-8").strip())
    except ValueError as exc:
        raise ValueError(f"Invalid XYZ atom-count line: {first!r}") from exc
    if natoms <= 0:
        raise ValueError(f"XYZ atom count must be positive, got {natoms}")

    comment_bytes = stream.readline()
    if not comment_bytes:
        raise ValueError("Unexpected end of file while reading XYZ comment")
    comment = comment_bytes.decode("utf-8").rstrip("\r\n")
    cell, pbc = _parse_comment(comment)

    species = []
    positions = np.empty((natoms, 3), dtype=np.float64)
    for index in range(natoms):
        line = stream.readline()
        if not line:
            raise ValueError(f"Unexpected end of file in atom line {index + 1}/{natoms}")
        fields = line.decode("utf-8").split()
        if len(fields) < 4:
            raise ValueError(f"Invalid XYZ atom line: {line!r}")
        species.append(fields[0])
        try:
            positions[index] = [float(fields[1]), float(fields[2]), float(fields[3])]
        except ValueError as exc:
            raise ValueError(f"Invalid XYZ coordinates: {line!r}") from exc

    return Frame(natoms, np.asarray(species, dtype="U16"), positions, cell, pbc, comment)


def iter_frames(path: str | Path) -> Iterator[Frame]:
    """Yield frames without materializing the trajectory."""
    with Path(path).open("rb") as stream:
        while True:
            frame = _read_frame(stream)
            if frame is None:
                return
            yield frame


# XRD numerical core
@dataclass(frozen=True)
class XRDConfig:
    """Validated settings for one XRD calculation."""

    wavelength: float
    two_theta_min: float
    two_theta_max: float
    bins: int = 250
    selection: tuple[str, ...] | None = None
    chunk_size: int = 1024


@dataclass(frozen=True)
class Mesh:
    """Reciprocal points and their two-theta coordinates for one frame."""

    q: np.ndarray
    two_theta: np.ndarray
    spacing: tuple[float, float, float]


@dataclass(frozen=True)
class FrameResult:
    """Histogrammed intensity produced for one frame."""

    intensity: np.ndarray
    intensity_over_total: np.ndarray
    count: int
    mesh_points: int


_MESH_CACHE: dict[tuple[float, ...], Mesh] = {}


def _validate_config(config: XRDConfig) -> None:
    """Validate the small set of user-facing XRD settings."""
    numeric_values = (
        config.wavelength,
        config.two_theta_min,
        config.two_theta_max,
    )
    if not all(math.isfinite(value) for value in numeric_values):
        raise ValueError("all numeric XRD parameters must be finite")
    if config.wavelength <= 0:
        raise ValueError("wavelength must be greater than zero")
    if not 0 < config.two_theta_min < config.two_theta_max < 180:
        raise ValueError("two-theta range must satisfy 0 < min < max < 180 degrees")
    if config.bins <= 0:
        raise ValueError("bins must be positive")
    if config.chunk_size <= 0:
        raise ValueError("chunk-size must be positive")


def _selection_mask(species: np.ndarray, selection: tuple[str, ...] | None) -> np.ndarray:
    """Build a case-insensitive mask for all atoms or selected elements."""
    if selection is None:
        return np.ones(species.size, dtype=bool)
    selected = {value.lower() for value in selection}
    return np.asarray(
        [value.lower() in selected for value in species],
        dtype=bool,
    )


def _inverse_lengths(frame: Frame) -> np.ndarray:
    """Return reciprocal cell lengths under the LAMMPS-compatible convention."""
    if frame.cell is None:
        raise ValueError("Input frame has no Lattice data")
    cell = np.asarray(frame.cell, dtype=np.float64)
    if cell.shape != (3, 3):
        raise ValueError("cell must be a 3x3 matrix")
    diagonal = np.diag(cell)
    if np.any(diagonal <= 0):
        raise ValueError("cell lengths must be positive")
    if not np.allclose(cell, np.diag(diagonal), rtol=1e-10, atol=1e-10):
        raise ValueError(
            "triclinic/non-orthogonal cells are not supported, matching LAMMPS compute xrd"
        )
    periodic = np.asarray(frame.pbc, dtype=bool)
    if not np.any(periodic):
        raise ValueError("at least one periodic direction is required")
    inverse = np.zeros(3, dtype=np.float64)
    inverse[periodic] = 1.0 / diagonal[periodic]
    average_inverse = float(np.mean(inverse[periodic]))
    inverse[~periodic] = average_inverse
    return inverse


def build_mesh(frame: Frame, config: XRDConfig) -> Mesh:
    """Build the reciprocal mesh using the same convention as LAMMPS."""
    _validate_config(config)
    spacing = _inverse_lengths(frame)
    cache_key = (
        *tuple(float(value) for value in spacing),
        float(config.wavelength),
        float(config.two_theta_min),
        float(config.two_theta_max),
    )
    cached = _MESH_CACHE.get(cache_key)
    if cached is not None:
        return cached

    theta_max = math.radians(config.two_theta_max / 2.0)
    kmax = 2.0 * math.sin(theta_max) / config.wavelength
    maximum = np.ceil(kmax / spacing).astype(int)
    axes = [np.arange(-limit, limit + 1, dtype=np.int32) for limit in maximum]
    i, j, k = np.meshgrid(*axes, indexing="ij")
    q = np.column_stack(
        (
            i.ravel().astype(np.float64) * spacing[0],
            j.ravel().astype(np.float64) * spacing[1],
            k.ravel().astype(np.float64) * spacing[2],
        )
    )
    q_norm = np.linalg.norm(q, axis=1)
    valid = config.wavelength * q_norm <= 2.0 + 1e-13
    theta = np.empty_like(q_norm)
    theta[valid] = np.arcsin(np.clip(0.5 * config.wavelength * q_norm[valid], 0.0, 1.0))
    valid &= theta >= math.radians(config.two_theta_min / 2.0)
    valid &= theta <= theta_max
    q = q[valid]
    two_theta = np.degrees(2.0 * theta[valid])
    if q.shape[0] == 0:
        raise ValueError(
            "no reciprocal points fall within the requested 2theta range; "
            "adjust the wavelength, range, or reciprocal spacing"
        )
    mesh = Mesh(q=q, two_theta=two_theta, spacing=tuple(float(x) for x in spacing))
    _MESH_CACHE[cache_key] = mesh
    return mesh


def _bin_indices(values: np.ndarray, config: XRDConfig) -> tuple[np.ndarray, np.ndarray]:
    """Map reciprocal points into bins spanning the requested 2theta range."""
    width = (config.two_theta_max - config.two_theta_min) / config.bins
    inside = (values >= config.two_theta_min) & (values <= config.two_theta_max)
    indices = np.floor(
        (values[inside] - config.two_theta_min) / width
    ).astype(np.int64)
    indices = np.minimum(indices, config.bins - 1)
    return indices, inside


def _atomic_form_factors(sin_theta_over_lambda: np.ndarray, coefficients: np.ndarray) -> np.ndarray:
    """Evaluate f(q) for every reciprocal point and element."""
    a = coefficients[:, 0:8:2]
    b = coefficients[:, 1:8:2]
    c = coefficients[:, 8]
    s2 = sin_theta_over_lambda * sin_theta_over_lambda
    return np.sum(
        a[None, :, :] * np.exp(-s2[:, None, None] * b[None, :, :]), axis=2
    ) + c[None, :]


def compute_frame(
    frame: Frame, config: XRDConfig, table: dict[str, np.ndarray]
) -> FrameResult:
    """Calculate and histogram one frame."""
    mesh = build_mesh(frame, config)
    mask = _selection_mask(frame.species, config.selection)
    positions = np.asarray(frame.positions[mask], dtype=np.float64)
    selected_species = frame.species[mask]
    if positions.shape[0] == 0:
        selected_text = ", ".join(config.selection or ())
        raise ValueError(
            "element selection contains no atoms"
            + (f": {selected_text}" if selected_text else "")
        )

    unique_species, species_indices = np.unique(
        selected_species, return_inverse=True
    )
    coefficients = resolve_labels(unique_species.tolist(), table)
    counts = np.bincount(species_indices, minlength=unique_species.size)
    order = np.argsort(species_indices, kind="stable")
    grouped_positions = positions[order]
    group_starts = np.concatenate(
        (np.array([0], dtype=np.int64), np.cumsum(counts[:-1], dtype=np.int64))
    )
    packed_phase_bytes = (
        config.chunk_size
        * positions.shape[0]
        * np.dtype(np.complex128).itemsize
    )
    use_packed_structure_factor = packed_phase_bytes <= MAX_PACKED_PHASE_BYTES
    intensity = np.empty(mesh.q.shape[0], dtype=np.float64)
    norm = np.linalg.norm(mesh.q, axis=1)

    for start in range(0, mesh.q.shape[0], config.chunk_size):
        stop = min(start + config.chunk_size, mesh.q.shape[0])
        q_chunk = mesh.q[start:stop]
        factors = _atomic_form_factors(0.5 * norm[start:stop], coefficients)
        if use_packed_structure_factor:
            phases = 2j * PI * (q_chunk @ grouped_positions.T)
            structure_factors = np.add.reduceat(
                np.exp(phases), group_starts, axis=1
            )
            amplitude = np.sum(factors * structure_factors, axis=1)
        else:
            real = np.zeros(stop - start, dtype=np.float64)
            imaginary = np.zeros(stop - start, dtype=np.float64)
            for species_index, start_atom in enumerate(group_starts):
                stop_atom = start_atom + counts[species_index]
                species_positions = grouped_positions[start_atom:stop_atom]
                phases = 2.0 * PI * (q_chunk @ species_positions.T)
                real += factors[:, species_index] * np.sum(
                    np.cos(phases), axis=1
                )
                imaginary += factors[:, species_index] * np.sum(
                    np.sin(phases), axis=1
                )
            amplitude = real + 1j * imaginary
        theta = 0.5 * np.radians(mesh.two_theta[start:stop])
        lp = (1.0 + np.cos(2.0 * theta) ** 2) / (
            np.cos(theta) * np.sin(theta) ** 2
        )
        intensity[start:stop] = (
            lp * (amplitude.real * amplitude.real + amplitude.imag * amplitude.imag)
            / positions.shape[0]
        )

    indices, valid = _bin_indices(mesh.two_theta, config)
    selected_intensity = intensity[valid]
    if selected_intensity.size == 0:
        raise ValueError(
            "no reciprocal points fall within the output 2theta range; adjust the range"
        )
    binned = np.bincount(
        indices,
        weights=selected_intensity,
        minlength=config.bins,
    ).astype(np.float64)
    count = int(np.count_nonzero(valid))
    total_intensity = float(np.sum(selected_intensity))
    if total_intensity:
        intensity_over_total = binned / total_intensity
    else:
        intensity_over_total = np.zeros_like(binned)
    return FrameResult(
        intensity=binned,
        intensity_over_total=intensity_over_total,
        count=count,
        mesh_points=mesh.q.shape[0],
    )


def bin_centers(config: XRDConfig) -> np.ndarray:
    width = (config.two_theta_max - config.two_theta_min) / config.bins
    return config.two_theta_min + (np.arange(config.bins) + 0.5) * width


# Interactive inputs, threaded processing, progress, and output


class Progress:
    """Lightweight progress reporter for streaming frame processing."""

    def __init__(self, total: int | None, enabled: bool) -> None:
        """Initialize a lightweight stderr progress reporter."""
        self.total = total
        self.enabled = enabled
        self.start = time.monotonic()
        self.done = 0

    def update(self, done: int) -> None:
        """Print the current progress when progress reporting is enabled."""
        if not self.enabled:
            return
        self.done = done
        width = 30
        elapsed = time.monotonic() - self.start
        if self.total is None:
            line = f"\r XRD processed {done} frame(s); elapsed {elapsed:7.1f}s"
        else:
            fraction = done / self.total if self.total else 1.0
            filled = int(width * fraction)
            rate = done / elapsed if elapsed > 0 else 0.0
            remaining = (self.total - done) / rate if rate > 0 else 0.0
            line = (
                f"\r XRD [{('=' * filled).ljust(width, '.')}] "
                f"{done}/{self.total} {fraction * 100:6.2f}% "
                f"elapsed {elapsed:7.1f}s"
            )
            if done < self.total:
                line += f" eta {remaining:7.1f}s"
        print(line, end="", file=sys.stderr, flush=True)

    def close(self) -> None:
        """Finish the progress line."""
        if self.enabled:
            if self.total is not None and self.done < self.total:
                self.update(self.total)
            print(file=sys.stderr)


def _write_output(
    path: str | Path,
    config: XRDConfig,
    nframes: int,
    intensity: np.ndarray,
    intensity_over_total: np.ndarray,
) -> None:
    """Write averaged histogrammed intensity data to a text file."""
    coordinates = bin_centers(config)
    if coordinates.size != intensity.size:
        raise ValueError("internal error: output coordinate and intensity sizes differ")
    with Path(path).open("w", encoding="utf-8") as output:
        output.write("# Histogrammed data for standalone LAMMPS-compatible XRD\n")
        output.write(f"# Frames {nframes}\n")
        elements = "all" if config.selection is None else ",".join(config.selection)
        output.write(f"# Elements {elements}\n")
        output.write(
            f"# wavelength {config.wavelength:g} 2theta {config.two_theta_min:g} "
            f"{config.two_theta_max:g} LP 1\n"
        )
        output.write("# Bin Coord Count Count/Total\n")
        for index, (coord, value, normalized) in enumerate(
            zip(coordinates, intensity, intensity_over_total), start=1
        ):
            output.write(f"{index} {coord:.10g} {value:.12g} {normalized:.12g}\n")


def _accumulate(results) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Average frame results and count processed reciprocal points."""
    sum_intensity = None
    sum_normalized = None
    total_points = 0
    nframes = 0
    for result in results:
        if sum_intensity is None:
            sum_intensity = np.zeros_like(result.intensity)
            sum_normalized = np.zeros_like(result.intensity_over_total)
        if result.intensity.shape != sum_intensity.shape:
            raise ValueError("different frame results have different bin counts")
        sum_intensity += result.intensity
        sum_normalized += result.intensity_over_total
        total_points += result.mesh_points
        nframes += 1
    if nframes == 0 or sum_intensity is None or sum_normalized is None:
        raise ValueError("no frames were found in the input")
    return sum_intensity / nframes, sum_normalized / nframes, nframes, total_points


def calculate(
    input_path: str | Path,
    output_path: str | Path,
    config: XRDConfig,
    workers: int,
) -> tuple[int, int]:
    """Stream, calculate, and average all frames in an extended XYZ file."""
    _validate_config(config)
    input_path = str(Path(input_path).expanduser().resolve())
    if workers == 0:
        workers = os.cpu_count() or 1
    if workers < 1:
        raise ValueError("workers must be positive, or zero for automatic selection")

    print(" Reading the first trajectory frame...", flush=True)
    frames = iter_frames(input_path)
    try:
        first_frame = next(frames)
    except StopIteration as exc:
        raise ValueError("input trajectory contains no frames") from exc

    print(" Preparing reciprocal mesh...", flush=True)
    mesh = build_mesh(first_frame, config)
    print(
        f" Reciprocal mesh ready: {mesh.q.shape[0]:,} points; "
        f"using {workers} CPU worker(s).",
        flush=True,
    )

    progress = Progress(None, True)
    results = []
    try:
        if workers == 1:
            results.append(compute_frame(first_frame, config, SCATTERING_FACTORS))
            progress.update(1)
            for index, frame in enumerate(frames, start=2):
                results.append(compute_frame(frame, config, SCATTERING_FACTORS))
                progress.update(index)
        else:
            pending = deque()
            with ThreadPoolExecutor(max_workers=workers) as pool:
                pending.append(
                    pool.submit(compute_frame, first_frame, config, SCATTERING_FACTORS)
                )
                for frame in frames:
                    pending.append(
                        pool.submit(compute_frame, frame, config, SCATTERING_FACTORS)
                    )
                    if len(pending) >= workers:
                        results.append(pending.popleft().result())
                        progress.update(len(results))
                while pending:
                    results.append(pending.popleft().result())
                    progress.update(len(results))
    finally:
        progress.close()

    average, normalized, nframes, total_points = _accumulate(results)
    _write_output(output_path, config, nframes, average, normalized)
    return nframes, total_points


def _read_input(prompt: str, default: str | None = None) -> str:
    """Read one interactive value and apply an optional displayed default."""
    suffix = f" [{default}]" if default is not None else ""
    try:
        value = input(f" {prompt}{suffix}: ").strip()
    except EOFError as exc:
        raise RuntimeError("input closed") from exc
    except KeyboardInterrupt as exc:
        raise RuntimeError("operation canceled") from exc
    return value if value else (default or "")


def _prompt_text(prompt: str, default: str | None = None) -> str:
    """Read a non-empty text value."""
    while True:
        value = _read_input(prompt, default)
        if value:
            return value
        print(" Error: input cannot be empty.")


def _prompt_float(
    prompt: str,
    default: float | None = None,
    minimum: float | None = None,
) -> float:
    """Read a floating-point value with an optional lower bound."""
    default_text = None if default is None else f"{default:g}"
    while True:
        value = _read_input(prompt, default_text)
        try:
            number = float(value)
        except ValueError:
            print(" Error: please enter a number.")
            continue
        if not math.isfinite(number):
            print(" Error: value must be a finite number.")
            continue
        if minimum is not None and number <= minimum:
            print(f" Error: value must be greater than {minimum:g}.")
            continue
        return number


def _prompt_int(
    prompt: str,
    default: int | None = None,
    minimum: int | None = None,
) -> int:
    """Read an integer value with an optional lower bound."""
    default_text = None if default is None else str(default)
    while True:
        value = _read_input(prompt, default_text)
        try:
            number = int(value)
        except ValueError:
            print(" Error: please enter an integer.")
            continue
        if minimum is not None and number < minimum:
            print(f" Error: value must be at least {minimum}.")
            continue
        return number


def _prompt_float_tuple(
    prompt: str, count: int, default: tuple[float, ...] | None = None
) -> tuple[float, ...]:
    """Read a whitespace-separated tuple of floating-point values."""
    default_text = (
        None if default is None else " ".join(f"{value:g}" for value in default)
    )
    while True:
        value = _read_input(prompt, default_text)
        try:
            numbers = tuple(float(item) for item in value.replace(",", " ").split())
        except ValueError:
            print(" Error: please enter numbers separated by spaces.")
            continue
        if len(numbers) != count:
            print(f" Error: please enter exactly {count} numbers.")
            continue
        if not all(math.isfinite(number) for number in numbers):
            print(" Error: values must be finite numbers.")
            continue
        return numbers


def _prompt_two_theta() -> tuple[float, float]:
    """Read and validate the calculated two-theta range."""
    while True:
        values = _prompt_float_tuple("2theta range (degrees; min max)", 2)
        if 0 < values[0] < values[1] < 180:
            return values
        print(" Error: 2theta must satisfy 0 < min < max < 180 degrees.")


def _parse_selection(value: str) -> tuple[str, ...] | None:
    """Parse all atoms or a comma/space-separated element selection."""
    if value.strip().lower() == "all":
        return None
    selection = tuple(
        item for item in re.split(r"[,\s]+", value.strip()) if item
    )
    if not selection:
        raise ValueError("enter all or at least one element symbol")
    return selection


def _prompt_selection() -> tuple[str, ...] | None:
    """Read the optional element selection."""
    while True:
        value = _read_input("Elements to include (all or comma-separated)", "all")
        try:
            return _parse_selection(value)
        except ValueError as exc:
            print(f" Error: {exc}.")


def _prompt_input_path() -> str:
    """Read the path of an existing input trajectory."""
    while True:
        value = _prompt_text("Input extended XYZ trajectory")
        path = Path(value).expanduser()
        if path.is_file():
            return str(path)
        print(f" Error: input file '{value}' does not exist.")


def interactive_inputs() -> tuple[str, str, XRDConfig, int]:
    """Collect the small set of XRD settings that users normally need."""
    print(" XRD calculation from an extended XYZ trajectory")
    print(" Defaults: all atoms, input Lattice/pbc, standard scattering factors, LP on.")

    input_path = _prompt_input_path()
    output_path = str(Path(_prompt_text("Output XRD file")).expanduser())
    wavelength = _prompt_float("X-ray wavelength (Angstrom)", minimum=0.0)
    two_theta_min, two_theta_max = _prompt_two_theta()
    bins = _prompt_int("Number of bins in this 2theta range", minimum=1)
    selection = _prompt_selection()
    workers = _prompt_int("CPU workers (0 means automatic)", default=0, minimum=0)

    config = XRDConfig(
        wavelength=wavelength,
        two_theta_min=two_theta_min,
        two_theta_max=two_theta_max,
        bins=bins,
        selection=selection,
    )
    return input_path, output_path, config, workers


def _print_help() -> None:
    """Print the interactive-only usage information."""
    print(" Usage: gpumdkit.sh -> 4) Calculators -> 413) Calc XRD")
    print("    or: python3 Scripts/calculators/calc_xrd.py")
    print("")
    print(" This calculator has no command-line calculation interface.")
    print(" It asks for the trajectory, output, wavelength, 2theta range, bins,")
    print(" element selection, and CPU workers interactively.")


def main() -> int:
    """Run the interactive XRD workflow and report user-facing errors."""
    if len(sys.argv) > 1:
        if sys.argv[1] in {"-h", "--help"}:
            _print_help()
            return 0
        print(" Error: XRD calculation is interactive only; remove command-line arguments.")
        return 1
    try:
        input_path, output_path, config, workers = interactive_inputs()
        nframes, total_points = calculate(input_path, output_path, config, workers)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f" Error: {exc}.")
        return 1
    print(
        f" Wrote {output_path} from {nframes} frame(s); "
        f"processed {total_points} reciprocal points."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

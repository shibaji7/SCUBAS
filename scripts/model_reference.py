# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""The formulas of the Models pages, written from the papers and independent of ``scubas``.

Nothing here imports SCUBAS. Each function names the paper and the equation it follows, so that ``tests/test_docs_model.py`` can compare the
code of the package with the published theory, and so that the figure script can overlay the theory on what the code computes.

Papers (all read for the Models pages):

* BP2003: Boteler and Pirjola (2003), Pure and Applied Geophysics 160, 1695-1716, "Magnetic and electric fields produced in the sea during
  geomagnetic disturbances".
* C2022: Chakraborty et al. (2022), Frontiers in Physics 10, 1022475, "Modeling geomagnetic induction in submarine cables".
* B2023: Boteler et al. (2023), International Journal of Geosciences 14, 767-791, "Transmission line modelling of geomagnetic induction in the
  Ocean/Earth conductivity structure".

Units are SI unless a name says otherwise: metres, ohm, siemens, volt, hertz, V/m. Line constants are per metre (the papers quote them per
kilometre, which is a factor of 1000).
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy.constants import mu_0

# --------------------------------------------------------------------------------------------------------------------------------------
# The seafloor field: a plane wave in a seawater layer over a layered Earth (BP2003)
# --------------------------------------------------------------------------------------------------------------------------------------


def seawater_k_and_z(rho_s: float, f: np.ndarray | float) -> tuple:
    """Propagation constant k and characteristic impedance Z of the seawater, BP2003 Eqs 3 and 4 (time factor exp(i w t))."""
    omega = 2.0 * np.pi * np.asarray(f, dtype=float)
    sigma = 1.0 / rho_s
    return np.sqrt(1j * omega * mu_0 * sigma), np.sqrt(1j * omega * mu_0 / sigma)


def surface_impedance(
    rho: Sequence[float], d: Sequence[float], f: np.ndarray | float
) -> np.ndarray:
    """Impedance E/H at the top of a stack of layers (``rho`` ohm m, ``d`` m, top first; the last layer is a half-space).

    The standard recursion for a layered Earth (Wait 1981; Weaver 1994, as used by BP2003 and C2022): starting from the half-space
    impedance, each layer j above it gives ``Z_j (Z_below + Z_j tanh(k_j d_j)) / (Z_j + Z_below tanh(k_j d_j))``.
    """
    omega = 2.0 * np.pi * np.asarray(f, dtype=float)
    k = [np.sqrt(1j * omega * mu_0 / r) for r in rho]
    z_in = 1j * omega * mu_0 / k[-1]
    for j in range(len(rho) - 2, -1, -1):
        z_j = 1j * omega * mu_0 / k[j]
        t = np.tanh(k[j] * d[j])
        z_in = z_j * (z_in + z_j * t) / (z_j + z_in * t)
    return z_in


def ed_over_h0(
    rho_s: float, d: float, z_d: np.ndarray, f: np.ndarray | float
) -> np.ndarray:
    """E at the seafloor over H at the sea surface, BP2003 Eq 10 (C2022 Eq 2 is the same expression divided by mu_0)."""
    k, z = seawater_k_and_z(rho_s, f)
    return (
        z * 2.0 / ((1.0 + z / z_d) * np.exp(k * d) - (1.0 - z / z_d) * np.exp(-k * d))
    )


def ed_over_h0_cosh(
    rho_s: float, d: float, z_d: np.ndarray, f: np.ndarray | float
) -> np.ndarray:
    """The same ratio in its hyperbolic form, BP2003 Eq 11."""
    k, z = seawater_k_and_z(rho_s, f)
    return z_d / (np.cosh(k * d) + (z_d / z) * np.sinh(k * d))


def bd_over_b0(
    rho_s: float, d: float, z_d: np.ndarray, f: np.ndarray | float
) -> np.ndarray:
    """B at the seafloor over B at the sea surface, BP2003 Eq 13."""
    k, z = seawater_k_and_z(rho_s, f)
    return 1.0 / (np.cosh(k * d) + (z_d / z) * np.sinh(k * d))


def ed_over_e0(
    rho_s: float, d: float, z_d: np.ndarray, f: np.ndarray | float
) -> np.ndarray:
    """E at the seafloor over E at the sea surface, BP2003 Eq 15."""
    k, z = seawater_k_and_z(rho_s, f)
    return 1.0 / (np.cosh(k * d) + (z / z_d) * np.sinh(k * d))


def transfer_function(
    rho: Sequence[float], d: Sequence[float], f: np.ndarray | float
) -> np.ndarray:
    """T_x = E_d / B_0 in mV/km/nT for seawater (``rho[0]``, ``d[0]``) over the layers below, C2022 Eq 2.

    The factor 1e-3 turns (V/m)/T into (mV/km)/nT.
    """
    z_d = surface_impedance(rho[1:], d[1:], f)
    return ed_over_h0(rho[0], d[0], z_d, f) / mu_0 * 1e-3


def field_ratios_with_depth(
    rho_s: float, d: float, z_d: complex, f: float, z: np.ndarray
) -> tuple:
    """E(z)/E_0 and B(z)/B_0 inside the ocean, BP2003 Eqs 5 to 8 (alpha is the reflected-to-incident ratio, Eq 8)."""
    k, zs = seawater_k_and_z(rho_s, f)
    alpha = ((z_d - zs) / (z_d + zs)) * np.exp(-2.0 * k * d)
    e_ratio = (np.exp(-k * z) + alpha * np.exp(k * z)) / (1.0 + alpha)
    b_ratio = (np.exp(-k * z) - alpha * np.exp(k * z)) / (1.0 - alpha)
    return e_ratio, b_ratio


def skin_depth(rho_s: float, f: float) -> float:
    """Skin depth in seawater, BP2003 Eq 16."""
    return float(np.sqrt(2.0 * rho_s / (2.0 * np.pi * f * mu_0)))


# --------------------------------------------------------------------------------------------------------------------------------------
# The transmission-line (DSTL) model (B2023)
# --------------------------------------------------------------------------------------------------------------------------------------


def line_constants(
    sigma_top: Sequence[float],
    d_top: Sequence[float],
    rho_bottom: Sequence[float],
    d_bottom: Sequence[float],
    width: float = 1.0,
) -> dict:
    """Z, Y, gamma, Z_0 per metre from the thin double layer, B2023 Eqs 1 to 5.

    The conducting layers (seawater and sediments) carry the current along the line in parallel: C = w sum(sigma_i d_i), Z = 1/C. The
    resistive layers (crust and mantle lithosphere) leak it downward in series: R = sum(rho_i d_i) / w, Y = 1/R.
    """
    c = width * float(np.sum(np.asarray(sigma_top) * np.asarray(d_top)))
    r = float(np.sum(np.asarray(rho_bottom) * np.asarray(d_bottom))) / width
    z, y = 1.0 / c, 1.0 / r
    return {
        "C": c,
        "R": r,
        "Z": z,
        "Y": y,
        "gamma": np.sqrt(z * y),
        "Z0": np.sqrt(z / y),
    }


def end_potential(e_field: float, length: float, gamma: float) -> float:
    """Potential at the start of an isolated section (open ends), C2022 Eq 11 in the form ``-(E/gamma) tanh(gamma L / 2)``.

    Written as ``-(E/gamma) sinh(gamma L) / (cosh(gamma L) + 1)`` in C2022 Eq 11; the two are equal because ``sinh x / (cosh x + 1) =
    tanh(x / 2)``. The potential at the other end has the opposite sign (Eq 12).
    """
    return (
        -(e_field / gamma) * np.sinh(gamma * length) / (np.cosh(gamma * length) + 1.0)
    )


def section_potential(
    x: np.ndarray, length: float, gamma: float, v_i: float, v_k: float
) -> np.ndarray:
    """Potential along a section from its end potentials, B2023 Eq 19 (C2022 Eq 10)."""
    x = np.asarray(x, dtype=float)
    den = np.exp(gamma * length) - np.exp(-gamma * length)
    return (v_k * np.exp(gamma * length) - v_i) / den * np.exp(
        -gamma * (length - x)
    ) + (v_i * np.exp(gamma * length) - v_k) / den * np.exp(-gamma * x)


def equivalent_pi(z: float, y: float, length: float, e_field: float) -> dict:
    """The equivalent-pi circuit of a section: B2023 Eqs 29 to 34 (``Ye`` is Y_E, ``Yp2`` is Y'/2, ``Ie`` is I_E)."""
    gamma, z0 = np.sqrt(z * y), np.sqrt(z / y)
    sh, ch = np.sinh(gamma * length), np.cosh(gamma * length)
    y_e = 1.0 / (z0 * sh)
    return {
        "Zp": z0 * sh,
        "Yp2": (ch - 1.0) * y_e,
        "Ep": e_field / gamma * sh,
        "Ye": y_e,
        "Ie": e_field / z,
        "gamma": gamma,
        "Z0": z0,
    }


def thevenin_of_long_line(e_field: float, z: float, y: float) -> dict:
    """The Thevenin and Norton equivalents of an electrically long line seen from one end, B2023 Eqs 25 to 28."""
    gamma, z0 = np.sqrt(z * y), np.sqrt(z / y)
    return {"VTh": -e_field / gamma, "ZTh": z0, "JN": e_field / z, "YN": 1.0 / z0}


def boundary_potential(v_th1: float, z_th1: float, v_th2: float, z_th2: float) -> tuple:
    """Current across and potential at the boundary of two electrically long regions, B2023 Eqs 41 and 42."""
    current = (v_th1 + v_th2) / (z_th1 + z_th2)
    return current, v_th1 - current * z_th1


def network_solve(
    sections: Sequence[dict], left: dict | None = None, right: dict | None = None
) -> tuple:
    """Nodal potentials of a chain of equivalent-pi sections, B2023 Eqs 43 to 55 (C2022 Eqs 6 to 9).

    ``sections`` holds the ``equivalent_pi`` dictionaries from the west end to the east end. ``left`` and ``right`` are the Norton circuits of
    the land at the two ends (``{"J": ..., "Y": ...}``), or None for an open end. The currents entering the nodes are
    ``J_1 = j_W - j_A``, ``J_2 = j_A - j_B``, ..., ``J_n = j_last - j_E`` (Eqs 46 to 52), and the admittance matrix has the sums of the pi
    admittances on its diagonal and ``-y_series`` off it (Eqs 43 to 45 and 53).

    Returns the node potentials, the admittance matrix and the current vector.
    """
    n = len(sections) + 1
    y_mat = np.zeros((n, n))
    j_vec = np.zeros(n)
    for i, s in enumerate(sections):
        y_mat[i, i] += s["Ye"] + s["Yp2"]
        y_mat[i + 1, i + 1] += s["Ye"] + s["Yp2"]
        y_mat[i, i + 1] -= s["Ye"]
        y_mat[i + 1, i] -= s["Ye"]
        j_vec[i] -= s["Ie"]
        j_vec[i + 1] += s["Ie"]
    if left is not None:
        y_mat[0, 0] += left["Y"]
        j_vec[0] += left["J"]
    if right is not None:
        y_mat[-1, -1] += right["Y"]
        j_vec[-1] -= right["J"]
    return np.linalg.solve(y_mat, j_vec), y_mat, j_vec


# --------------------------------------------------------------------------------------------------------------------------------------
# From the magnetic record to the electric field, and the voltage on the cable (C2022 Section 2)
# --------------------------------------------------------------------------------------------------------------------------------------


def synthetic_field(
    t: np.ndarray,
    amplitude_nt: Sequence[float],
    period_min: Sequence[float],
    phase_deg: Sequence[float],
) -> np.ndarray:
    """The synthetic magnetic record of C2022 Eq 24: a sum of sines, ``t`` in seconds."""
    return sum(
        a * np.sin(2 * np.pi * t / (p * 60.0) + np.radians(ph))
        for a, p, ph in zip(amplitude_nt, period_min, phase_deg, strict=True)
    )


def analytic_efield(
    t: np.ndarray,
    amplitude_nt: Sequence[float],
    period_min: Sequence[float],
    phase_deg: Sequence[float],
    tf_values: Sequence[complex],
) -> np.ndarray:
    """The analytic seafloor field of C2022 Eq 25: each sine scaled and phase-shifted by the transfer function at its own frequency."""
    return sum(
        a * abs(tf) * np.sin(2 * np.pi * t / (p * 60.0) + np.radians(ph) + np.angle(tf))
        for a, p, ph, tf in zip(
            amplitude_nt, period_min, phase_deg, tf_values, strict=True
        )
    )


def cable_voltage(emf: float, u_a: float, u_b: float) -> float:
    """V_C = emf + U_A - U_B, C2022 Eq 27 and B2023 Eq 60."""
    return emf + u_a - u_b


def km_lengths(lat_deg: float, dlat_deg: float, dlon_deg: float) -> tuple:
    """North and east lengths in km of a section between two points, C2022 Eqs 4 and 5 (WGS84 series; ``lat_deg`` is the mean latitude)."""
    lat = np.radians(lat_deg)
    north = (111.133 - 0.56 * np.cos(2 * lat)) * dlat_deg
    east = (111.5065 - 0.1872 * np.cos(2 * lat)) * np.cos(lat) * dlon_deg
    return north, east

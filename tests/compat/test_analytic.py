# SPDX-License-Identifier: MIT
# Copyright (c) 2023 Shibaji Chakraborty
r"""Closed-form expectations from ``docs/tutorial/theory.md`` and ``elsc.md``, checked against the legacy code as it is.

These agree with the code today to the tolerances below; they are pins, so a refactor that moves them has changed the physics.
"""

import numpy as np
import pandas as pd
import pytest
from loguru import logger
from scipy import constants as C

from scubas.cables import Cable, TransmissionLine
from scubas.datasets import PROFILES, Site
from scubas.utils import GreatCircle

TIMES = pd.date_range("2020-01-01", periods=3, freq="60s", name="Time")


@pytest.fixture(autouse=True)
def _quiet():
    logger.disable("scubas")
    yield
    logger.enable("scubas")


def single_section(length_km, profile, e_mv_per_km=100.0):
    efield = pd.DataFrame({"X": [e_mv_per_km] * 3}, index=TIMES)
    tl = TransmissionLine(
        "s", directed_length={"length_north": length_km}, elec_params={"site": profile}
    )
    tl.compute_eqv_pi_circuit(Efield=efield, components=["X"])
    return tl, Cable([tl], components=["X"])


def test_electrically_long_cable_end_potential_is_minus_E_over_gamma():
    """V_i = -V_k = -E / gamma once gamma L >> 1 (Boteler et al.)."""
    tl, cable = single_section(20000.0, PROFILES.DO)
    e = 100.0e-6  # V/m
    assert tl.gma * 20000.0e3 > 15
    assert cable.V["X"][0, 0] == pytest.approx(-e / tl.gma, rel=1e-6)
    assert cable.V["X"][1, 0] == pytest.approx(e / tl.gma, rel=1e-6)


def test_electrically_short_cable_end_potential_is_minus_L_E_over_2():
    """V_i = -V_k = -L E / 2 once gamma L << 1."""
    _, cable = single_section(1.0, PROFILES.DO)
    assert cable.V["X"][0, 0] == pytest.approx(-0.5 * 1.0e3 * 100.0e-6, rel=1e-5)
    assert cable.V["X"][1, 0] == pytest.approx(0.5 * 1.0e3 * 100.0e-6, rel=1e-5)


def test_the_equivalent_pi_branch_admittance_of_a_long_line():
    """Y_E = 2 / (Z0 e^{gamma L}) for a long line, and the source current is I_E = E / Z."""
    tl, _ = single_section(20000.0, PROFILES.DO)
    length_m = 20000.0e3
    assert tl.Ye["X"] == pytest.approx(
        1.0 / (tl.Z0 * np.sinh(tl.gma * length_m)), rel=1e-12
    )
    assert tl.Ye["X"] == pytest.approx(
        2.0 / (tl.Z0 * np.exp(tl.gma * length_m)), rel=1e-6
    )
    np.testing.assert_allclose(tl.Ie["X"], 100.0e-6 / tl.Z, rtol=1e-12)


def test_line_parameters_follow_from_the_layers():
    """gamma = sqrt(Z Y), Z0 = sqrt(Z / Y) with Z = 1/C and Y = 1/R, C and R from the layer thicknesses and resistivities."""
    tl = TransmissionLine(
        "p",
        directed_length={"length": 10.0},
        elec_params={"site": PROFILES.CS, "width": 2.0},
    )
    s = PROFILES.CS
    c_expected = 2.0 * (
        s.get_thicknesses(1) / s.get_resistivities(1)
        + s.get_thicknesses(0) / s.get_resistivities(0)
    )
    r_expected = (
        s.get_thicknesses(2) * s.get_resistivities(2)
        + s.get_thicknesses(3) * s.get_resistivities(3)
    ) / 2.0
    assert (
        pytest.approx(c_expected, rel=1e-12) == tl.C
        and pytest.approx(r_expected, rel=1e-12) == tl.R
    )
    assert tl.gma == pytest.approx(np.sqrt(1.0 / (tl.C * tl.R)), rel=1e-12)
    assert pytest.approx(np.sqrt(tl.R / tl.C), rel=1e-12) == tl.Z0


def test_half_space_impedance_is_sqrt_i_omega_mu_rho():
    """A uniform Earth: Z = sqrt(i w mu rho), here scaled to mV/km per nT as the code returns it."""
    site = Site.init(
        [0.1, 0.1, 0.1], [10e3, 20e3, np.inf], ["a", "b", "c"], "uniform", "uniform"
    )
    f = np.logspace(-4, -1, 7)
    w = 2 * np.pi * f
    expected = np.sqrt(1j * w * C.mu_0 * 10.0) * (1.0e-3 / C.mu_0)
    np.testing.assert_allclose(site.calcZ(f)[1], expected, rtol=1e-12)
    np.testing.assert_allclose(site.calcZ(f)[2], -expected, rtol=1e-12)


def test_skin_depth_is_the_impedance_over_i_omega_mu():
    site = Site.init(
        [0.1, 0.1, 0.1], [10e3, 20e3, np.inf], ["a", "b", "c"], "uniform", "uniform"
    )
    f = np.logspace(-4, -1, 7)
    z = site.calcZ(f)[1]
    w = 2 * np.pi * f
    np.testing.assert_allclose(
        site.calcP(f), z / (1j * w * C.mu_0) / (1.0e-3 / C.mu_0), rtol=1e-12
    )


def test_great_circle_and_haversine_agree():
    for a, b in [
        ((0.0, 0.0), (0.0, 90.0)),
        ((39.6, -74.33), (50.79, -4.55)),
        ((-33.9, 151.2), (51.5, -0.1)),
    ]:
        from types import SimpleNamespace

        gc = GreatCircle(
            SimpleNamespace(lat=a[0], lon=a[1]), SimpleNamespace(lat=b[0], lon=b[1])
        )
        assert gc.great_circle() == pytest.approx(gc.haversine(), rel=1e-9)


def test_subdividing_a_uniform_cable_does_not_change_the_end_potentials():
    """A property of the nodal solution: two equal halves of a uniform cable give the same end potentials as the whole."""
    whole = single_section(600.0, PROFILES.CS)[1]
    efield = pd.DataFrame({"X": [100.0] * 3}, index=TIMES)
    halves = []
    for i in range(2):
        tl = TransmissionLine(
            f"h{i}",
            directed_length={"length_north": 300.0},
            elec_params={"site": PROFILES.CS},
        )
        tl.compute_eqv_pi_circuit(Efield=efield, components=["X"])
        halves.append(tl)
    split = Cable(halves, components=["X"])
    np.testing.assert_allclose(split.V["X"][[0, 2]], whole.V["X"], rtol=1e-9)


def test_the_response_is_linear_in_the_electric_field():
    a = single_section(600.0, PROFILES.CS, 100.0)[1].V["X"]
    b = single_section(600.0, PROFILES.CS, 300.0)[1].V["X"]
    np.testing.assert_allclose(b, 3.0 * a, rtol=1e-12)

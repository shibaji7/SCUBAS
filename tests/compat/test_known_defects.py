# SPDX-License-Identifier: MIT
# Copyright (c) 2023 Shibaji Chakraborty
"""Known defects of SCUBAS 0.1.9, pinned as they are.

Each test documents behaviour that looks wrong but that users and published results may depend on. They are NOT a request to
keep the behaviour forever: a fix is a separate, approved change that updates the one test here and says so in the changelog
(an output-changing fix goes behind a new function or option, never silently).
"""

import numpy as np
import pandas as pd
import pytest
from loguru import logger

from scubas.cables import Cable, TransmissionLine
from scubas.datasets import PROFILES
from scubas.sources.current_sources import Current, Direction


@pytest.fixture(autouse=True)
def _quiet():
    logger.disable("scubas")
    yield
    logger.enable("scubas")


def small_cable():
    index = pd.date_range("2020-01-01", periods=4, freq="60s", name="Time")
    lines = []
    for i in range(2):
        tl = TransmissionLine(f"S{i}", directed_length={"length_north": 300.0})
        tl.compute_eqv_pi_circuit(
            Efield=pd.DataFrame({"X": [137.7, 51.3, -20.9, 8.8]}, index=index),
            components=["X"],
        )
        lines.append(tl)
    return Cable(lines, components=["X"])


def test_direction_enumerate_key_cannot_be_called():
    """It compares the class with a member, so no branch matches and ``key`` is never bound."""
    with pytest.raises(UnboundLocalError):
        Direction.enumerate_key("B")


def test_compute_b_field_rescales_the_callers_array_in_place():
    """``x`` is multiplied by the unit factor (km to m) in the array the caller passed."""
    cur = Current(site=PROFILES.Quebec)
    x = np.linspace(-100.0, 100.0, len(cur.omega))
    original = x.copy()
    cur.compute_B_field(x=x, length_unit="km")
    np.testing.assert_allclose(x, original * 1.0e3)


def test_end_potentials_are_rounded_to_two_decimals_and_feed_the_total():
    """``_pot_end*`` round to 0.01 V, and ``Vt(v) = U0 - U1 + V(v)`` uses the rounded values."""
    cable = small_cable()
    u0, u1 = cable._pot_end_cable_()
    np.testing.assert_array_equal(u0, np.round(cable.V["X"][0, :], 2))
    np.testing.assert_array_equal(u1, np.round(cable.V["X"][-1, :], 2))
    np.testing.assert_allclose(
        cable.tot_params["Vt(v)"], u0 - u1 + cable.tot_params["V(v)"], rtol=1e-12
    )
    assert not np.allclose(
        u0, cable.V["X"][0, :], rtol=0, atol=1e-6
    )  # the rounding is visible


def test_unit_argument_of_the_end_potentials_is_v_or_anything_else_means_millivolts():
    cable = small_cable()
    u_v, _ = cable._pot_end_cable_(unit="V")
    u_other, _ = cable._pot_end_cable_(unit="whatever")
    np.testing.assert_allclose(u_other, np.round(cable.V["X"][0, :] * 1000.0, 2))
    assert not np.allclose(u_v, u_other)


def test_layer_thicknesses_are_in_metres_although_the_docstring_says_km():
    """Seawater of the continental shelf is 100 (m); the ocean model uses them as metres."""
    assert PROFILES.CS.get_thicknesses(0) == 100
    assert "km" in PROFILES.CS.get_thicknesses.__doc__


def test_a_section_without_an_active_termination_flag_before_the_efield_is_computed():
    tl = TransmissionLine("S", directed_length={"length": 1.0})
    with pytest.raises(RuntimeError, match="E-field"):
        tl.add_active_termination()


def test_potential_profile_before_the_nodal_solution_fails_with_an_attribute_error():
    """The guard written for this (``RuntimeError: ... run nodal analysis first``) never fires: ``end_pot`` already exists, empty."""
    tl = TransmissionLine("S", directed_length={"length": 1.0})
    with pytest.raises(AttributeError, match="Vi"):
        tl._pot_alongCS_()

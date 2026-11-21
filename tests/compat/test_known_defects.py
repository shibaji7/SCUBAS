# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
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


def test_ifft_of_fft_returns_twice_the_signal():
    """``ifft`` multiplies by ``2 n`` (it was ``n`` before the commit of 2023-02-20), so the round trip doubles the signal.

    This is why ``OceanModel.to_Efields`` returns twice the analytic field of Chakraborty et al. (2022), Eq 25 (see
    ``tests/test_docs_model.py``).
    """
    from scubas.utils import fft, ifft

    t = np.arange(1000)
    x = np.sin(2 * np.pi * t / 50.0)
    spectrum, _ = fft(x, 60.0)
    assert np.abs(ifft(spectrum)).max() / np.abs(x).max() == pytest.approx(
        2.0, abs=0.01
    )


def test_the_taper_ramp_depends_on_the_cadence_of_the_record():
    """The cosine ramp uses elapsed seconds over a sample count: smooth for one-second data, oscillating for one-minute data."""
    from scubas.models import Preprocess

    n = 4320
    ramp_1s = Preprocess(np.arange(n) * 1.0, np.zeros(n)).get_tapering_function(0.1)[
        :216
    ]
    ramp_60s = Preprocess(np.arange(n) * 60.0, np.zeros(n)).get_tapering_function(0.1)[
        :216
    ]
    assert np.all(np.diff(ramp_1s) >= 0)
    assert not np.all(np.diff(ramp_60s) >= 0)  # it goes up and down

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""The Models pages against the code and the papers.

Every equation of ``docs/model`` is written out independently in ``scripts/model_reference.py`` from the papers (BP2003, C2022, B2023), and
the tests compare it with the package. The numbers printed in the papers (their tables, and values read off their figures where the paper
gives no table) are compared with the package too. ``docs/model/figures/numbers.json``, which the pages quote, must be what the scripts
compute now.
"""

import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.constants import mu_0

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import model_numbers  # noqa: E402
import model_reference as ref  # noqa: E402
from model_cases import (  # noqa: E402
    atlantic_site,
    efield_frame,
    network_case,
    potentials,
    single_section,
)
from scubas.cables import CableSection  # noqa: E402
from scubas.datasets import PROFILES, Layer, Site  # noqa: E402
from scubas.models import OceanModel  # noqa: E402
from scubas.sources.current_sources import Current  # noqa: E402

SYNTH_AMPLITUDE = [200, 90, 30, 17, 8, 3.5]
SYNTH_PERIOD = [180, 80, 36, 15, 8, 3]
SYNTH_PHASE = [10, 20, 30, 40, 50, 60]
ATLANTIC_RHO = [0.3, 3.0, 3000.0, 1000.0, 100.0, 10.0, 1.0]
ATLANTIC_D = [4000.0, 2000.0, 10000.0, 140000.0, 254000.0, 250000.0, 340000.0]  # m


def tf(site: Site, freqs, kind: str = "E2B") -> np.ndarray:
    return np.asarray(OceanModel(site, flim=[1e-6, 1.0]).calcTF(freqs=freqs)[kind])


def half_space(rho_sea: float, depth_m: float, rho_earth: float = 10.0) -> Site:
    return Site(
        [Layer("Sea", depth_m, 1.0 / rho_sea), Layer("Earth", 1e6, 1.0 / rho_earth)],
        "BP2003",
        "BP2003",
    )


# --------------------------------------------------------------------------------------------------------------------------------------
# The seafloor field (BP2003, C2022 Section 2.2)
# --------------------------------------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("depth", [10.0, 100.0, 1000.0, 4000.0, 5000.0])
def test_the_transfer_function_is_eq_10_of_bp2003_for_the_atlantic_model(depth):
    freqs = np.logspace(-6, -1, 25)
    site = atlantic_site(depth / 1e3)
    rho = [r for _, _, r in [(0, 0, x) for x in ATLANTIC_RHO]]
    d = [depth, *ATLANTIC_D[1:]]
    np.testing.assert_allclose(
        tf(site, freqs), ref.transfer_function(rho, d, freqs), rtol=1e-9
    )


def test_the_two_forms_of_the_ratio_agree():
    """BP2003 Eq 10 (exponentials) and Eq 11 (hyperbolic functions) are the same function."""
    f = np.logspace(-5, -1, 9)
    z_d = ref.surface_impedance(ATLANTIC_RHO[1:], ATLANTIC_D[1:], f)
    np.testing.assert_allclose(
        ref.ed_over_h0(0.3, 4000.0, z_d, f),
        ref.ed_over_h0_cosh(0.3, 4000.0, z_d, f),
        rtol=1e-9,
    )


def test_the_b_transfer_function_is_eq_13_of_bp2003():
    f = np.logspace(-4, -2, 7)
    site = half_space(0.25, 5000.0)
    z_d = ref.surface_impedance([10.0], [1e6], f)
    np.testing.assert_allclose(
        tf(site, f, "B2B"), ref.bd_over_b0(0.25, 5000.0, z_d, f), rtol=1e-9
    )


def test_the_skin_depth_of_seawater_is_the_4_36_km_of_bp2003():
    assert ref.skin_depth(0.25, 1.0 / 300.0) == pytest.approx(4.36e3, rel=2e-3)


@pytest.mark.parametrize("m", range(6))
def test_table_4_of_c2022_follows_from_the_code(m):
    """Amplitude of the analytic field and phase (phi + theta) of the six frequencies of the synthetic test."""
    paper_amplitude = [10.35, 4.86, 1.65, 0.95, 0.45, 0.19]
    paper_phase = [20.66, 25.0, 31.22, 36.88, 41.93, 36.37]
    value = tf(atlantic_site(4.0), [1.0 / (SYNTH_PERIOD[m] * 60.0)])[0]
    assert SYNTH_AMPLITUDE[m] * abs(value) == pytest.approx(
        paper_amplitude[m], abs=0.006
    )
    assert SYNTH_PHASE[m] + np.degrees(np.angle(value)) == pytest.approx(
        paper_phase[m], abs=0.01
    )


@pytest.mark.parametrize(
    "depth, expected_b, expected_e",
    [(100.0, 0.87, 0.9999), (1000.0, 0.35, 0.96), (5000.0, 0.085, 0.575)],
)
def test_figure_1_of_bp2003_ocean_depths(depth, expected_b, expected_e):
    """Values read off BP2003 Figure 1 (5 min period, seawater 0.25 ohm m, Earth 10 ohm m): B and E at the seafloor over their surface values."""
    f = 1.0 / 300.0
    z_d = ref.surface_impedance([10.0], [1e6], f)
    e_ratio, b_ratio = ref.field_ratios_with_depth(
        0.25, depth, z_d, f, np.array([depth])
    )
    assert abs(b_ratio[0]) == pytest.approx(expected_b, abs=0.02)
    assert abs(e_ratio[0]) == pytest.approx(expected_e, abs=0.02)
    code_b = abs(tf(half_space(0.25, depth), [f], "B2B")[0])
    assert code_b == pytest.approx(abs(b_ratio[0]), rel=1e-9)


@pytest.mark.parametrize(
    "depth, f, expected",
    [
        (100.0, 1e-4, 0.065),
        (100.0, 1e-2, 0.52),
        (5000.0, 1e-4, 0.027),
        (5000.0, 1e-3, 0.034),
        (5000.0, 1e-2, 0.028),
    ],
)
def test_figure_5_of_bp2003_seafloor_e_over_surface_b(depth, f, expected):
    """Values read off BP2003 Figure 5, within the width of a line on its log axis."""
    assert abs(tf(half_space(0.25, depth), [f])[0]) == pytest.approx(expected, rel=0.12)


def test_figure_4_of_bp2003_the_deep_ocean_magnetic_field_at_2_min():
    assert abs(tf(half_space(0.25, 5000.0), [1.0 / 120.0], "B2B")[0]) == pytest.approx(
        0.045, rel=0.12
    )  # the curve ends at the 2 min period


# --------------------------------------------------------------------------------------------------------------------------------------
# The line constants (B2023 Eqs 1 to 5)
# --------------------------------------------------------------------------------------------------------------------------------------

PAPER_TABLE_2 = {  # B2023 Table 2: Z (ohm/km), Y (S/km), gamma (1/km), Z0 (ohm), adjustment distance (km)
    "CS": (7.5e-1, 5.0e-6, 1.936e-3, 3.873e2, 516.4),
    "SO": (2.5e-1, 1.0e-5, 1.58e-3, 1.58e2, 632.5),
    "DO": (7.143e-2, 1.0e-5, 8.45e-4, 8.45e1, 1183.2),
}


@pytest.mark.parametrize("name", sorted(PAPER_TABLE_2))
def test_table_2_of_b2023(name):
    line = single_section(getattr(PROFILES, name), 100.0, 1.0)
    z, y, gamma, z0, adjustment = PAPER_TABLE_2[name]
    assert float(line.Z) * 1e3 == pytest.approx(z, rel=1e-3)
    assert float(line.Y) * 1e3 == pytest.approx(y, rel=1e-3)
    assert float(line.gma) * 1e3 == pytest.approx(gamma, rel=3e-3)
    assert float(line.Z0) == pytest.approx(z0, rel=3e-3)
    assert 1.0 / (float(line.gma) * 1e3) == pytest.approx(adjustment, rel=1e-3)


def test_table_2_of_b2023_land_row():
    line = TransmissionLine_for_land()
    assert float(line.Z) * 1e3 == pytest.approx(3.0, rel=1e-3)
    assert float(line.Y) * 1e3 == pytest.approx(5.0e-6, rel=1e-3)
    assert float(line.gma) * 1e3 == pytest.approx(3.873e-3, rel=1e-3)
    assert float(line.Z0) == pytest.approx(7.746e2, rel=1e-3)
    assert 1.0 / (float(line.gma) * 1e3) == pytest.approx(258.2, rel=1e-3)


def TransmissionLine_for_land():
    return model_numbers.land_line()


def test_table_3_of_c2022_atlantic_model():
    line = single_section(atlantic_site(4.0), 100.0, 1.0)
    assert float(line.Z) * 1e3 == pytest.approx(7.14e-2, rel=1e-3)
    assert float(line.Y) * 1e3 == pytest.approx(5.88e-6, rel=1e-3)
    assert float(line.Z0) == pytest.approx(110.0, rel=3e-3)
    assert float(line.gma) * 1e3 == pytest.approx(6.48e-4, rel=1e-3)
    assert 1.0 / (float(line.gma) * 1e3) == pytest.approx(1543.21, rel=5e-4)


@pytest.mark.parametrize("name", ["CS", "SO", "DO"])
def test_the_line_constants_are_the_thin_double_layer_of_eqs_1_to_4(name):
    site = getattr(PROFILES, name)
    sigma = [site.get_conductivities(i) for i in (0, 1)]
    d = [site.get_thicknesses(i) for i in (0, 1)]
    rho = [site.get_resistivities(i) for i in (2, 3)]
    dd = [site.get_thicknesses(i) for i in (2, 3)]
    expected = ref.line_constants(sigma, d, rho, dd)
    line = single_section(site, 100.0, 1.0)
    for key, attribute in [("Z", "Z"), ("Y", "Y"), ("gamma", "gma"), ("Z0", "Z0")]:
        assert float(getattr(line, attribute)) == pytest.approx(
            expected[key], rel=1e-12
        )


# --------------------------------------------------------------------------------------------------------------------------------------
# One section: the end potentials, the equivalent-pi circuit, the potential along the section (C2022 Section 2.4, B2023 Section 3)
# --------------------------------------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("length", [10.0, 100.0, 600.0, 1500.0, 4000.0, 10000.0])
def test_the_end_potential_of_a_section_is_eq_11_of_c2022(length):
    """A section with open ends; E = 1 V/km. The package rounds the end potentials to 0.01 V (a known defect), the path does not."""
    site = atlantic_site(4.0)
    line = single_section(site, length, 1000.0)
    gamma = float(line.gma) * 1e3
    values, _ = potentials(
        network_case([(site, length, 1000.0)], land_e_mv_per_km=None)
    )
    assert values[0] == pytest.approx(ref.end_potential(1.0, length, gamma), rel=1e-9)
    assert values[-1] == pytest.approx(-ref.end_potential(1.0, length, gamma), rel=1e-9)


@pytest.mark.parametrize("x", [1e-3, 0.1, 0.5, 1.0, 2.0, 4.0, 10.0])
def test_sinh_over_cosh_plus_one_is_tanh_of_half(x):
    assert np.sinh(x) / (np.cosh(x) + 1.0) == pytest.approx(np.tanh(x / 2.0), rel=1e-12)


def test_the_two_limits_of_c2022_eqs_17_and_22():
    site = atlantic_site(4.0)
    gamma = float(single_section(site, 10.0, 1.0).gma) * 1e3
    adjustment = 1.0 / gamma
    long_l, short_l = 10.0 * adjustment, 0.01 * adjustment
    assert ref.end_potential(1.0, long_l, gamma) == pytest.approx(
        -1.0 / gamma, rel=2e-4
    )  # -E/gamma
    assert ref.end_potential(1.0, short_l, gamma) == pytest.approx(
        -short_l / 2.0, rel=1e-4
    )  # -LE/2
    # at the limit of the "electrically long" criterion, L = 4/gamma, the potential is tanh(2) = 96.4 % of E/gamma: the criterion is a few percent
    assert ref.end_potential(1.0, 4.0 * adjustment, gamma) / (
        -1.0 / gamma
    ) == pytest.approx(np.tanh(2.0), rel=1e-12)


def test_the_potential_along_a_section_is_eq_19_of_b2023():
    site = atlantic_site(4.0)
    length = 1500.0
    line = single_section(site, length, 1000.0)
    gamma = float(line.gma) * 1e3
    values, distance = potentials(
        network_case([(site, length, 1000.0)], land_e_mv_per_km=None)
    )
    u_a = ref.end_potential(1.0, length, gamma)
    np.testing.assert_allclose(
        values,
        ref.section_potential(distance, length, gamma, u_a, -u_a),
        rtol=1e-7,
        atol=1e-7,
    )


@pytest.mark.parametrize(
    "name, length", [("CS", 100.0), ("SO", 1800.0), ("DO", 7800.0)]
)
def test_the_equivalent_pi_circuit_is_eqs_29_to_34_of_b2023(name, length):
    site = getattr(PROFILES, name)
    line = single_section(site, length, 300.0)
    z, y = float(line.Z), float(line.Y)
    expected = ref.equivalent_pi(z, y, length * 1e3, 3e-4)
    assert float(line.Ye["X"]) == pytest.approx(expected["Ye"], rel=1e-9)
    assert float(line.Yp2["X"]) == pytest.approx(expected["Yp2"], rel=1e-9)
    assert float(line.Ie["X"][0]) == pytest.approx(expected["Ie"], rel=1e-9)
    # the voltage-source form: the same circuit with Z' = Z0 sinh and E' = E sinh / gamma, so Y_E = 1/Z' and I_E Z' = E'
    assert expected["Ye"] == pytest.approx(1.0 / expected["Zp"], rel=1e-12)
    assert expected["Ie"] / expected["Ye"] == pytest.approx(expected["Ep"], rel=1e-12)


def test_the_short_section_limit_of_the_pi_circuit_eqs_35_to_37():
    z, y, length, e = 0.75e-3, 5e-9, 1e3, 3e-4  # 1 km of the shelf line: gamma L = 2e-3
    pi = ref.equivalent_pi(z, y, length, e)
    assert pi["Zp"] == pytest.approx(z * length, rel=1e-5)
    assert pi["Yp2"] / pi["Ye"] < 1e-5
    assert pi["Ep"] == pytest.approx(e * length, rel=1e-5)


def test_the_cable_voltage_is_c2022_eq_27():
    """Vt(v) = U0 - U1 + V(v) in the code is V_C = emf + U_A - U_B, with the emf the section sum of E L."""
    cable = network_case([(PROFILES.CS, 600.0, 300.0)])
    row = cable.tot_params.iloc[0]
    emf = 0.3 * 600.0  # V: 0.3 V/km over 600 km
    assert row["V(v)"] == pytest.approx(emf, rel=1e-9)
    assert row["Vt(v)"] == pytest.approx(
        ref.cable_voltage(emf, row["U0"], row["U1"]), abs=1e-9
    )
    assert row["Vt(v)"] == pytest.approx(emf - 2 * 48.2, abs=0.05)


# --------------------------------------------------------------------------------------------------------------------------------------
# Land at the ends, and the network (B2023 Sections 3.3, 4, 5 and 6)
# --------------------------------------------------------------------------------------------------------------------------------------


def test_the_land_norton_circuit_is_eqs_27_and_28_of_b2023():
    land = model_numbers.land_line()
    z, y = float(land.Z), float(land.Y)
    expected = ref.thevenin_of_long_line(3e-4, z, y)
    termination = (
        network_case([(PROFILES.CS, 100.0, 300.0)])
        .cable_sections[0]
        .term_params["left"]
    )
    assert float(termination.Jn["X"][0]) == pytest.approx(expected["JN"], rel=1e-12)
    assert float(termination.Yn) == pytest.approx(expected["YN"], rel=1e-12)


def test_table_3_and_eqs_56_to_59_of_b2023():
    numbers = model_numbers.collect()
    assert float(numbers["dstl_vth_do"]) == pytest.approx(118.32, abs=0.01)
    assert float(numbers["dstl_vth_cs"]) == pytest.approx(154.92, abs=0.01)
    assert float(numbers["dstl_vth_ld"]) == pytest.approx(77.46, abs=0.01)
    assert float(numbers["dstl_zth_do"]) == pytest.approx(84.5, abs=0.05)
    assert float(numbers["dstl_zth_cs"]) == pytest.approx(387.3, abs=0.05)
    assert float(numbers["dstl_zth_ld"]) == pytest.approx(774.6, abs=0.05)
    assert float(numbers["dstl_i_cs_ld"]) == pytest.approx(0.2, abs=0.0005)
    assert float(numbers["dstl_uc_limit"]) == pytest.approx(77.46, abs=0.005)
    assert float(numbers["dstl_i_do_cs"]) == pytest.approx(0.579, abs=0.0005)
    assert float(numbers["dstl_ue_limit"]) == pytest.approx(69.38, abs=0.005)


def test_the_boundary_current_and_potential_are_eqs_41_and_42():
    current, potential = ref.boundary_potential(154.92, 387.3, 77.46, 774.6)
    assert current == pytest.approx(0.2, abs=1e-3)
    assert potential == pytest.approx(77.46, abs=0.01)


def test_a_wide_shelf_gives_the_limiting_coast_potential_of_b2023_eq_57():
    values, _ = potentials(network_case([(PROFILES.CS, 20000.0, 300.0)]))
    assert abs(values[0]) == pytest.approx(77.46, abs=0.01)
    assert abs(values[-1]) == pytest.approx(77.46, abs=0.01)


@pytest.mark.parametrize(
    "sections, expected",
    [
        ([(PROFILES.CS, 600.0, 300.0)], 48.0),  # B2023 Figure 10a
        (
            [
                (PROFILES.CS, 100.0, 300.0),
                (PROFILES.SO, 1800.0, 100.0),
                (PROFILES.CS, 100.0, 300.0),
            ],
            48.0,
        ),  # Figure 10b
        (
            [
                (PROFILES.CS, 100.0, 300.0),
                (PROFILES.DO, 7800.0, 100.0),
                (PROFILES.CS, 100.0, 300.0),
            ],
            105.0,
        ),  # Figure 10c
    ],
)
def test_figure_10_of_b2023_end_potentials(sections, expected):
    """Read off the figure, to the accuracy of its axes."""
    values, _ = potentials(network_case(sections))
    assert abs(values[0]) == pytest.approx(expected, abs=1.5)
    assert values[0] == pytest.approx(-values[-1], rel=1e-9)


def test_a_zero_width_shelf_gives_the_99_v_of_the_b2023_discussion():
    values, _ = potentials(
        network_case(
            [
                (PROFILES.CS, 1.0, 300.0),
                (PROFILES.DO, 10000.0, 100.0),
                (PROFILES.CS, 1.0, 300.0),
            ]
        )
    )
    assert abs(values[0]) == pytest.approx(99.0, abs=1.0)


def test_the_coast_potential_of_a_wide_shelf_peaks_near_200_km_figure_12_of_b2023():
    widths = [50.0, 150.0, 250.0, 500.0, 1000.0, 2000.0]
    uc = [
        abs(
            potentials(
                network_case(
                    [
                        (PROFILES.CS, w, 300.0),
                        (PROFILES.DO, 10000.0, 100.0),
                        (PROFILES.CS, w, 300.0),
                    ]
                )
            )[0][0]
        )
        for w in widths
    ]
    assert uc.index(max(uc)) in (1, 2)  # the maximum is between 100 and 300 km
    assert uc[-1] < uc[0]


@pytest.mark.parametrize("which", ["shallow", "ocean", "deep"])
def test_the_network_of_the_code_is_the_network_of_eqs_43_to_55(which):
    cs, so, do = PROFILES.CS, PROFILES.SO, PROFILES.DO
    sections = {
        "shallow": [(cs, 600.0, 300.0)],
        "ocean": [(cs, 100.0, 300.0), (so, 1800.0, 100.0), (cs, 100.0, 300.0)],
        "deep": [(cs, 100.0, 300.0), (do, 7800.0, 100.0), (cs, 100.0, 300.0)],
    }[which]
    cable = network_case(sections)
    pis = []
    for site, length, e in sections:
        line = single_section(site, length, e)
        pis.append(
            ref.equivalent_pi(float(line.Z), float(line.Y), length * 1e3, e * 1e-6)
        )
    land = model_numbers.land_line()
    norton = ref.thevenin_of_long_line(3e-4, float(land.Z), float(land.Y))
    solved, _, _ = ref.network_solve(
        pis,
        left={"J": norton["JN"], "Y": norton["YN"]},
        right={"J": norton["JN"], "Y": norton["YN"]},
    )
    np.testing.assert_allclose(cable.V["X"][:, 0], solved, rtol=1e-9)


# --------------------------------------------------------------------------------------------------------------------------------------
# From the record to the field (C2022 Section 2.2 and 2.4.2)
# --------------------------------------------------------------------------------------------------------------------------------------


def synthetic_run():
    n = 72 * 60
    t = np.arange(n) * 60.0
    b = ref.synthetic_field(t, SYNTH_AMPLITUDE, SYNTH_PERIOD, SYNTH_PHASE)
    site = atlantic_site(4.0)
    values = [tf(site, [1.0 / (p * 60.0)])[0] for p in SYNTH_PERIOD]
    analytic = ref.analytic_efield(
        t, SYNTH_AMPLITUDE, SYNTH_PERIOD, SYNTH_PHASE, values
    )
    frame = pd.DataFrame(
        {"X": 0.0 * b, "Y": b},
        index=pd.date_range("2020-01-01", periods=n, freq="60s", name="Time"),
    )
    model = OceanModel(site, flim=[1e-6, 1.0])
    model.to_Efields(Bfield=frame, components=["X", "Y"], p=0.1)
    keep = slice(n // 10, -n // 10)
    return analytic[keep], np.asarray(model.Efield["X"])[keep]


def test_the_synthetic_test_reproduces_the_shape_of_the_analytic_field():
    analytic, numeric = synthetic_run()
    assert np.corrcoef(analytic, numeric)[0, 1] > 0.9999


def test_the_field_the_code_returns_is_twice_the_analytic_field():
    """Known defect, recorded: ``ifft(fft(x))`` is 2 x since the commit of 2023-02-20, so ``to_Efields`` is twice C2022 Eq 25 (docs/dev/testing.md)."""
    analytic, numeric = synthetic_run()
    assert np.polyfit(analytic, numeric, 1)[0] == pytest.approx(2.0, abs=0.01)


def test_the_north_and_east_lengths_are_eqs_4_and_5_of_c2022():
    section = CableSection(
        "A",
        directed_length={
            "edge_locations": {
                "initial": {"lat": 40.0, "lon": -70.0},
                "final": {"lat": 43.0, "lon": -60.0},
            }
        },
    )
    north, east = ref.km_lengths(41.5, 3.0, 10.0)
    assert section.length_north == pytest.approx(north, rel=1e-12)
    assert section.length_east == pytest.approx(east, rel=1e-12)


def test_the_lengths_do_not_change_sign_with_the_direction_of_the_section():
    """The package takes absolute differences, so a section laid south-west has the same lengths as one laid north-east."""
    north_east = CableSection(
        "A",
        directed_length={
            "edge_locations": {
                "initial": {"lat": 40.0, "lon": -70.0},
                "final": {"lat": 43.0, "lon": -60.0},
            }
        },
    )
    south_west = CableSection(
        "B",
        directed_length={
            "edge_locations": {
                "initial": {"lat": 43.0, "lon": -60.0},
                "final": {"lat": 40.0, "lon": -70.0},
            }
        },
    )
    assert south_west.length_north == pytest.approx(north_east.length_north, rel=1e-12)
    assert south_west.length_east == pytest.approx(north_east.length_east, rel=1e-12)


# --------------------------------------------------------------------------------------------------------------------------------------
# The source of the electrojet page: the line current of ``scubas.sources.Current``
# --------------------------------------------------------------------------------------------------------------------------------------

MU0 = mu_0  # scipy's value; BP2003 write 4 pi 1e-7, which differs by 1e-10


def line_current(p_m: float, site=None, omega: float = 2 * np.pi / 300.0) -> Current:
    return Current(
        h=110.0,
        a=0.0,
        I=1.0,
        p=np.array([p_m / 1e3]),
        omega=np.array([omega]),
        site=site,
        current_unit="MA",
    )


def test_the_primary_field_of_the_line_current_is_amperes_law():
    """B_x = mu_0 I h / (2 pi (h^2 + x^2)) and |B_z| = mu_0 I |x| / (2 pi (h^2 + x^2)) for a line current I at height h."""
    x_km = np.linspace(-1000.0, 1000.0, 41)
    current = line_current(0.0)
    current.compute_B_field(x=x_km.copy())
    bx, bz = current.B_fields
    h, x = 110e3, x_km * 1e3
    assert (
        bx["direction"].name == "GEO_NORTH_SOUTH"
        and bz["direction"].name == "GEO_VERT_DOWN"
    )
    np.testing.assert_allclose(
        bx["external"], MU0 * 1e6 * h / (2 * np.pi * (h**2 + x**2)), rtol=1e-12
    )
    np.testing.assert_allclose(
        abs(bz["external"]),
        MU0 * 1e6 * abs(x) / (2 * np.pi * (h**2 + x**2)),
        rtol=1e-12,
    )


def test_a_perfectly_conducting_earth_doubles_the_horizontal_field_and_cancels_the_vertical_one():
    """Zero skin depth: the image is a mirror image at depth h, the total horizontal field is twice the primary one, B_z and E vanish."""
    x_km = np.linspace(-500.0, 500.0, 21)
    current = line_current(0.0)
    current.compute_B_field(x=x_km.copy())
    bx, bz = current.B_fields
    np.testing.assert_allclose(bx["total"], 2.0 * bx["external"], rtol=1e-12)
    np.testing.assert_allclose(bz["total"], 0.0, atol=1e-18)
    np.testing.assert_allclose(current.E_fields[0]["total"], 0.0, atol=1e-18)


def test_the_induced_field_is_that_of_an_image_at_complex_depth_h_plus_2p():
    """The image of the current lies at z = h + 2p (Pirjola and Viljanen 1998, Eq 1 and the text after Eq 2); p = Z / (i omega mu_0) comes from Site.calcP."""
    site = PROFILES.Quebec
    omega = 2 * np.pi / 300.0
    current = line_current(0.0, site=site, omega=omega)
    z_surface = site.calcZ(np.array([omega / (2 * np.pi)]), layer=0)[1][0] * (
        mu_0 / 1e-3
    )  # back to ohm
    p = z_surface / (1j * omega * MU0)
    np.testing.assert_allclose(current.p, p, rtol=1e-9)
    x_km = np.array([0.0, 150.0, 600.0])
    current.compute_B_field(x=x_km.copy())
    h, x = 110e3, x_km * 1e3
    image = MU0 * 1e6 / (2 * np.pi) * (h + 2 * p) / ((h + 2 * p) ** 2 + x**2)
    np.testing.assert_allclose(current.B_fields[0]["internal"], image, rtol=1e-12)


# --------------------------------------------------------------------------------------------------------------------------------------
# The numbers the pages quote
# --------------------------------------------------------------------------------------------------------------------------------------


def test_numbers_json_is_what_the_scripts_compute_now():
    """Run ``python scripts/make_model_figures.py`` after a change to the code or the scripts, and commit the new file."""
    committed = json.loads(
        (ROOT / "docs/model/figures/numbers.json").read_text(encoding="utf-8")
    )
    assert committed == model_numbers.collect()


def test_every_figure_the_pages_use_exists():
    text = "".join(p.read_text() for p in (ROOT / "docs/model").glob("*.md"))
    import re

    for name in set(re.findall(r"figures/([\w\-]+\.(?:png|svg))", text)):
        assert (ROOT / "docs/model/figures" / name).exists() or (
            ROOT / "docs/figures" / name
        ).exists(), name


def test_every_number_a_page_substitutes_is_in_numbers_json():
    import re

    keys = set(
        json.loads(
            (ROOT / "docs/model/figures/numbers.json").read_text(encoding="utf-8")
        )
    )
    for page in (ROOT / "docs/model").glob("*.md"):
        used = set(re.findall(r"\{\{\s*([\w]+)\s*\}\}", page.read_text()))
        assert used <= keys, (page.name, sorted(used - keys))


def test_the_efield_frame_helper_is_in_millivolts_per_kilometre():
    assert float(efield_frame(300.0)["X"].iloc[0]) == 300.0


def test_the_python_blocks_of_the_pages_run():
    """The code shown on the Models pages is executed."""
    import re

    for page in (ROOT / "docs/model").glob("*.md"):
        for block in re.findall(r"```python\n(.*?)```", page.read_text(), flags=re.S):
            exec(compile(block, str(page), "exec"), {})

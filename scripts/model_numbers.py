# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""The numbers that the Models pages quote, computed from the SCUBAS code and written to ``docs/model/figures/numbers.json``.

``docs/conf.py`` loads that file as MyST substitutions, so a number in the text of a page is written ``{{ dstl_do_gamma }}`` and cannot differ
from what the package computes. ``tests/test_docs_model.py`` checks the file against the formulas of ``scripts/model_reference.py`` and
against the values printed in the papers, and fails when the file is stale (``python scripts/make_model_figures.py`` rewrites it).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

from model_cases import atlantic_site, network_case, potentials, single_section
from model_reference import analytic_efield, synthetic_field
from scubas.cables import TransmissionLine
from scubas.datasets import PROFILES
from scubas.models import OceanModel

SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(x: float, digits: int = 4) -> str:
    """``7.143 × 10⁻²``: the form of the tables of the papers (a plain number between 0.01 and 10000)."""
    x = float(x)
    if x == 0:
        return "0"
    if 1e-2 <= abs(x) < 1e4:
        return f"{x:.{digits}g}"
    mantissa, exponent = f"{x:.{digits - 1}e}".split("e")
    return f"{float(mantissa):g} × 10{str(int(exponent)).translate(SUPERSCRIPT)}"


def line_numbers(prefix: str, line: TransmissionLine) -> dict:
    """Z (ohm/km), Y (S/km), gamma (1/km), Z0 (ohm) and the adjustment distance (km) of a section."""
    z, y, gamma, z0 = (
        float(line.Z) * 1e3,
        float(line.Y) * 1e3,
        float(line.gma) * 1e3,
        float(line.Z0),
    )
    return {
        f"{prefix}_Z": sci(z),
        f"{prefix}_Y": sci(y),
        f"{prefix}_gamma": sci(gamma),
        f"{prefix}_Z0": sci(z0),
        f"{prefix}_adj": sci(1.0 / gamma, 5),
    }


def land_line() -> TransmissionLine:
    """The land profile as the package treats it (its own branch in ``calc_trasmission_line_parameters``)."""
    return TransmissionLine(
        "land", directed_length={"length_north": 1.0}, elec_params={"site": PROFILES.LD}
    )


def collect() -> dict:
    """Every number the pages quote."""
    n: dict[str, str] = {}
    # Line constants: B2023 Table 2 and C2022 Table 3.
    for key, site in [
        ("cs", PROFILES.CS),
        ("so", PROFILES.SO),
        ("do", PROFILES.DO),
        ("atl", atlantic_site(4.0)),
    ]:
        n.update(line_numbers(f"dstl_{key}", single_section(site, 100.0, 1.0)))
    n.update(line_numbers("dstl_ld", land_line()))

    # Thevenin circuits and the two limiting boundary potentials: B2023 Table 3 and Eqs 56 to 59.
    land = land_line()
    vth, zth = {}, {}
    for key, e_v_per_km, line in [
        ("do", 0.1, single_section(PROFILES.DO, 100.0, 1.0)),
        ("cs", 0.3, single_section(PROFILES.CS, 100.0, 1.0)),
        ("ld", 0.3, land),
    ]:
        vth[key], zth[key] = e_v_per_km / (float(line.gma) * 1e3), float(line.Z0)
        n[f"dstl_vth_{key}"], n[f"dstl_zth_{key}"] = sci(vth[key], 5), sci(zth[key], 4)
    i_cs_ld = (vth["cs"] + vth["ld"]) / (zth["cs"] + zth["ld"])
    n["dstl_i_cs_ld"] = sci(i_cs_ld, 3)
    n["dstl_uc_limit"] = f"{vth['cs'] - i_cs_ld * zth['cs']:.2f}"
    i_do_cs = (vth["do"] + vth["cs"]) / (zth["do"] + zth["cs"])
    n["dstl_i_do_cs"] = sci(i_do_cs, 3)
    n["dstl_ue_limit"] = f"{vth['do'] - i_do_cs * zth['do']:.2f}"

    # The worked network: the deep ocean of B2023 Figure 9c (two shelves of 100 km and 7800 km of ocean, land at both ends).
    sections = [
        (PROFILES.CS, 100.0, 300.0),
        (PROFILES.DO, 7800.0, 100.0),
        (PROFILES.CS, 100.0, 300.0),
    ]
    cable = network_case(sections)
    land = cable.cable_sections[0].term_params["left"]
    for i, section in enumerate(cable.cable_sections):
        tag = "abc"[i]
        n[f"net_{tag}_L"] = f"{float(section.length):g}"
        n[f"net_{tag}_gL"] = sci(float(section.gma) * float(section.length) * 1e3, 4)
        n[f"net_{tag}_Ye"] = sci(float(section.Ye["X"]), 4)
        n[f"net_{tag}_Yp2"] = sci(float(section.Yp2["X"]), 4)
        n[f"net_{tag}_Ie"] = sci(float(section.Ie["X"][0]), 4)
    n["net_land_Yn"] = sci(float(land.Yn), 4)
    n["net_land_Jn"] = sci(float(land.Jn["X"][0]), 4)
    for i in range(4):
        node = cable.nodes[i]["X"]
        n[f"net_J{i + 1}"] = sci(float(node.Ji[0]), 4)
        for j in range(4):
            if abs(node.Yii[j]) > 0:
                n[f"net_Y{i + 1}{j + 1}"] = sci(float(node.Yii[j]), 4)
    volts = np.linalg.solve(
        np.array([cable.nodes[i]["X"].Yii for i in range(4)]),
        np.array([cable.nodes[i]["X"].Ji[0] for i in range(4)]),
    )
    for i, v in enumerate(volts):
        n[f"net_U{i + 1}"] = f"{v:.2f}"
    values, _ = potentials(cable)
    n["net_deep_end"] = f"{abs(values[0]):.1f}"
    shallow, _ = potentials(network_case([(PROFILES.CS, 600.0, 300.0)]))
    n["net_shallow_end"] = f"{abs(shallow[0]):.1f}"
    so_values, _ = potentials(
        network_case(
            [
                (PROFILES.CS, 100.0, 300.0),
                (PROFILES.SO, 1800.0, 100.0),
                (PROFILES.CS, 100.0, 300.0),
            ]
        )
    )
    n["net_so_end"] = f"{abs(so_values[0]):.1f}"

    # The synthetic test (C2022 Section 2.4.2, Table 4).
    site = atlantic_site(4.0)
    model = OceanModel(site, flim=[1e-6, 1.0])
    amp, per, phi = (
        [200, 90, 30, 17, 8, 3.5],
        [180, 80, 36, 15, 8, 3],
        [10, 20, 30, 40, 50, 60],
    )
    for m, (a, p, ph) in enumerate(zip(amp, per, phi, strict=True), start=1):
        tf = complex(np.asarray(model.calcTF(freqs=[1.0 / (p * 60.0)])["E2B"])[0])
        n[f"synth_amp_{m}"] = f"{a * abs(tf):.2f}"
        n[f"synth_phase_{m}"] = f"{ph + np.degrees(np.angle(tf)):.2f}"
    nsamp = 72 * 60
    t = np.arange(nsamp) * 60.0
    tfs = [
        complex(np.asarray(model.calcTF(freqs=[1.0 / (p * 60.0)])["E2B"])[0])
        for p in per
    ]
    e_an = analytic_efield(t, amp, per, phi, tfs)
    frame = pd.DataFrame(
        {"X": 0.0 * t, "Y": synthetic_field(t, amp, per, phi)},
        index=pd.date_range("2020-01-01", periods=nsamp, freq="60s", name="Time"),
    )
    model.to_Efields(Bfield=frame, components=["X", "Y"], p=0.1)
    e_num = np.asarray(model.Efield["X"])
    keep = slice(nsamp // 10, -nsamp // 10)
    n["synth_slope"] = f"{np.polyfit(e_an[keep], e_num[keep], 1)[0]:.3f}"
    n["synth_corr"] = f"{np.corrcoef(e_an[keep], e_num[keep])[0, 1]:.5f}"
    return n


def write_numbers(path: Path) -> None:
    """Write ``numbers.json`` (sorted keys, so that a re-run gives the same bytes)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(collect(), indent=1, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print("wrote", path.name)

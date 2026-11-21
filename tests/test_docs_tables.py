# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Numbers printed in the documentation must agree with the code and with themselves."""

import pathlib
import re

import numpy as np
import pytest

from scubas.cables import TransmissionLine
from scubas.datasets import PROFILES

ROOT = pathlib.Path(__file__).resolve().parents[1]
SUPERSCRIPT = str.maketrans("⁻⁰¹²³⁴⁵⁶⁷⁸⁹", "-0123456789")
PROFILE_OF = {
    "Continental shelf": "CS",
    "Shallow ocean": "SO",
    "Deep ocean": "DO",
    "Land": "LD",
}


def number(cell: str) -> float:
    """Read "7.5 × 10⁻¹" (or a plain number) from a table cell."""
    match = re.fullmatch(r"([\d.]+)(?: × 10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+))?", cell.strip())
    assert match, cell
    exponent = int(match.group(2).translate(SUPERSCRIPT)) if match.group(2) else 0
    return float(match.group(1)) * 10.0**exponent


def conduct_rows() -> dict[str, list[float]]:
    """The rows of the transmission-line parameter table of the conduct tutorial."""
    text = (ROOT / "docs/tutorial/conduct.md").read_text()
    rows = {}
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if line.startswith("| ") and cells[0] in PROFILE_OF:
            rows[cells[0]] = [number(c) for c in cells[1:]]
    return rows


def test_the_table_has_the_four_sections():
    assert sorted(conduct_rows()) == sorted(PROFILE_OF)


@pytest.mark.parametrize("section", sorted(PROFILE_OF))
def test_each_row_follows_from_its_own_z_and_y(section):
    _, z, y, gamma, z0, adjustment = conduct_rows()[section]
    assert gamma == pytest.approx(np.sqrt(z * y), rel=2e-3)
    assert z0 == pytest.approx(np.sqrt(z / y), rel=2e-3)
    assert adjustment == pytest.approx(1 / gamma, rel=2e-3)


@pytest.mark.parametrize("section", sorted(PROFILE_OF))
def test_z_and_y_are_those_of_the_profile_in_the_code(section):
    _, z, y, *_ = conduct_rows()[section]
    line = TransmissionLine(
        section,
        directed_length={"length_north": 100.0},
        elec_params={"site": getattr(PROFILES, PROFILE_OF[section])},
    )
    assert z == pytest.approx(float(line.Z) * 1e3, rel=2e-3)  # ohm/m to ohm/km
    assert y == pytest.approx(float(line.Y) * 1e3, rel=2e-3)  # S/m to S/km


def test_the_readthedocs_badge_uses_the_lowercase_project_name():
    """shields.io looks the project up case-sensitively: "SCUBAS" says "no finished builds"."""
    for path in [ROOT / "README.md", *(ROOT / "docs").rglob("*.md")]:
        assert "readthedocs/SCUBAS" not in path.read_text(), path

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""The test cases of the papers, set up with the SCUBAS classes, for the Models pages, their figures and their tests.

* ``atlantic_site``: the Ocean-Earth model of C2022 Table 2 ("Atlantic"), with the seawater depth as a parameter.
* ``network_case``: a cable of several sections with the land at the two ends as Norton circuits, as in B2023 Section 5 (the shallow sea, the
  shallow ocean and the deep ocean of its Figure 9) and Section 6 (the widths of the sea and of the continental shelf).

The terminations are set from a given land field and not from a magnetic record: ``TransmissionLine.add_active_termination`` computes the
land field from a B-field file, which is the normal use, but the papers specify the land by a uniform field of 0.3 V/km, so the same Norton
circuit (J = E/Z, Y = 1/Z0 of the land profile) is attached here by hand through the ``term_params`` attribute that the package itself fills.
"""

from __future__ import annotations

from collections.abc import Sequence
from types import SimpleNamespace

import numpy as np
import pandas as pd

from scubas.cables import Cable, TransmissionLine
from scubas.datasets import PROFILES, Layer, Site

# C2022 Table 2: the Atlantic model (thickness in km, resistivity in ohm m), top first.
ATLANTIC = [
    ("Seawater", 4.0, 0.3),
    ("Sediments", 2.0, 3.0),
    ("Crust", 10.0, 3000.0),
    ("Mantle Lithosphere", 140.0, 1000.0),
    ("Upper Mantle", 254.0, 100.0),
    ("Transition Zone", 250.0, 10.0),
    ("Lower Mantle", 340.0, 1.0),
]


def site_from_rows(name: str, rows: Sequence[tuple]) -> Site:
    """A ``Site`` from ``(name, thickness_km, resistivity_ohm_m)`` rows. The package keeps thicknesses in metres."""
    return Site([Layer(n, t * 1e3, 1.0 / r) for n, t, r in rows], name, name)


def atlantic_site(ocean_depth_km: float = 4.0) -> Site:
    """The Atlantic model of C2022 Table 2 with the given seawater depth (the 4 km of the table by default)."""
    rows = [(ATLANTIC[0][0], ocean_depth_km, ATLANTIC[0][2]), *ATLANTIC[1:]]
    return site_from_rows("Atlantic", rows)


def efield_frame(value_mv_per_km: float) -> pd.DataFrame:
    """A one-sample E-field table in mV/km for the X (north) component."""
    return pd.DataFrame({"X": [value_mv_per_km]}, index=pd.RangeIndex(1, name="Time"))


def single_section(
    site: Site, length_km: float, e_mv_per_km: float
) -> TransmissionLine:
    """One transmission-line section with a uniform northward field and no terminations."""
    line = TransmissionLine(
        "S", directed_length={"length_north": length_km}, elec_params={"site": site}
    )
    line.compute_eqv_pi_circuit(Efield=efield_frame(e_mv_per_km), components=["X"])
    return line


def land_norton(e_mv_per_km: float) -> SimpleNamespace:
    """The Norton circuit of the land profile ``PROFILES.LD`` for a uniform field: ``Jn = E / Z`` and ``Yn = 1 / Z0``, as the package sets it."""
    land = TransmissionLine(
        "land", directed_length={"length_north": 1.0}, elec_params={"site": PROFILES.LD}
    )
    return SimpleNamespace(
        Yn=1.0 / land.Z0,
        Jn={"X": np.array([e_mv_per_km * 1e-6 / land.Z])},
        Efield=efield_frame(e_mv_per_km),
    )


def network_case(
    sections: Sequence[tuple], land_e_mv_per_km: float | None = 300.0
) -> Cable:
    """A ``Cable`` of sections ``(site, length_km, e_mv_per_km)`` from west to east, land (Norton) at both ends when a land field is given."""
    lines = []
    for i, (site, length_km, e_field) in enumerate(sections):
        line = TransmissionLine(
            f"s{i}",
            directed_length={"length_north": length_km},
            elec_params={"site": site},
        )
        line.compute_eqv_pi_circuit(Efield=efield_frame(e_field), components=["X"])
        lines.append(line)
    if land_e_mv_per_km is not None:
        for line in (lines[0], lines[-1]):
            line.has_active_term = True
            line.term_params = {"left": None, "right": None}
        lines[0].term_params["left"] = land_norton(land_e_mv_per_km)
        lines[-1].term_params["right"] = land_norton(land_e_mv_per_km)
    return Cable(lines, components=["X"])


def potentials(cable: Cable) -> tuple:
    """Potential in volts along the whole cable and the distance in km, from the package (``Cable._pot_along_cable_``)."""
    values, distance = cable._pot_along_cable_(timestamp=0)
    return np.asarray(values, dtype=float), np.asarray(distance, dtype=float)

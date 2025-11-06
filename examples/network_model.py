"""Multisection network modelling demonstrations for the conductance tutorial.

This script recreates the three canonical scenarios described in
``docs/tutorial/conduct.md`` – continental shelf, shallow ocean, and deep ocean –
using the modern SCUBAS API.  For each configuration it:

* configures Matplotlib styling via ``scubas.plotlib.update_rc_params``;
* assembles one or more ``TransmissionLine`` sections with the appropriate
  conductivity profiles; and
* solves for the induced cable potentials, exporting publication-ready figures
  used directly in the documentation.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

os.environ["MPLBACKEND"] = "Agg"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import matplotlib.pyplot as plt
import pandas as pd

from scubas.cables import Cable, TransmissionLine
from scubas.datasets import PROFILES
from scubas.plotlib import cable_potential, update_rc_params


def create_induced_field(millivolts_per_km: float) -> pd.DataFrame:
    """Return a single-sample electric-field dataframe (X component only)."""
    return pd.DataFrame(
        {"X": [millivolts_per_km]},
        index=pd.RangeIndex(1, name="Time"),
    )


def build_transmission_line(
    *,
    sec_id: str,
    length_km: float,
    profile: Any,
    efield: pd.DataFrame,
    terminations: Optional[Mapping[str, Mapping[str, Any]]] = None,
) -> TransmissionLine:
    """Instantiate and parameterise a transmission-line section."""
    transmission_line = TransmissionLine(
        sec_id=sec_id,
        directed_length={"length_north": length_km},
        elec_params={"site": profile, "width": 1.0, "flim": [1e-6, 1.0]},
        active_termination=terminations,
    )
    transmission_line.compute_eqv_pi_circuit(Efield=efield, components=["X"])
    return transmission_line


def simulate_case(
    *,
    title: str,
    slug: str,
    sections: Sequence[Mapping[str, Any]],
    figures_dir: Path,
    ylim: Sequence[float],
) -> None:
    """Solve the cable potentials for a multi-section configuration."""
    transmission_lines = [
        build_transmission_line(
            sec_id=section["sec_id"],
            length_km=section["length_km"],
            profile=section["profile"],
            efield=section["efield"],
            terminations=section.get("terminations"),
        )
        for section in sections
    ]
    cable = Cable(transmission_lines, components=["X"])
    potentials, distances = cable._pot_along_cable_(timestamp=0)

    artifacts = cable_potential(potentials, distances, ylim=ylim)
    total_length = sum(section["length_km"] for section in sections)
    artifacts.axes.set_title(title)
    artifacts.axes.text(
        0.05,
        0.85,
        rf"$L_{{\mathrm{{total}}}}$ = {total_length:.0f} km",
        ha="left",
        va="center",
        transform=artifacts.axes.transAxes,
    )
    artifacts.figure.savefig(
        figures_dir / f"{slug}.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(artifacts.figure)


def main() -> None:
    """Generate figures for the three canonical conductivity case studies."""
    figures_dir = (
        Path(__file__).resolve().parent.parent / "docs" / "tutorial" / "figures"
    )
    figures_dir.mkdir(parents=True, exist_ok=True)

    update_rc_params(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"],
        },
        science=True,
        latex=False,
    )

    # Shared electric-field definitions (values in mV/km).
    efield_cs = create_induced_field(300.0)  # 0.3 V/km
    efield_ocean = create_induced_field(100.0)  # 0.1 V/km

    simulate_case(
        title="Case 1 – Continental Shelf",
        slug="conduct_case1_continental_shelf",
        sections=[
            {
                "sec_id": "CS",
                "length_km": 600.0,
                "profile": PROFILES.CS,
                "efield": efield_cs,
            },
        ],
        figures_dir=figures_dir,
        ylim=(-200.0, 200.0),
    )

    simulate_case(
        title="Case 2 – Shallow Ocean",
        slug="conduct_case2_shallow_ocean",
        sections=[
            {
                "sec_id": "CS-west",
                "length_km": 100.0,
                "profile": PROFILES.CS,
                "efield": efield_cs,
                "terminations": {"left": {"site": PROFILES.LD, "width": 1.0}},
            },
            {
                "sec_id": "SO",
                "length_km": 1800.0,
                "profile": PROFILES.SO,
                "efield": efield_ocean,
            },
            {
                "sec_id": "CS-east",
                "length_km": 100.0,
                "profile": PROFILES.CS,
                "efield": efield_cs,
                "terminations": {"right": {"site": PROFILES.LD, "width": 1.0}},
            },
        ],
        figures_dir=figures_dir,
        ylim=(-60.0, 60.0),
    )

    simulate_case(
        title="Case 3 – Deep Ocean",
        slug="conduct_case3_deep_ocean",
        sections=[
            {
                "sec_id": "CS-west",
                "length_km": 100.0,
                "profile": PROFILES.CS,
                "efield": efield_cs,
                "terminations": {"left": {"site": PROFILES.LD, "width": 1.0}},
            },
            {
                "sec_id": "DO",
                "length_km": 7800.0,
                "profile": PROFILES.DO,
                "efield": efield_ocean,
            },
            {
                "sec_id": "CS-east",
                "length_km": 100.0,
                "profile": PROFILES.CS,
                "efield": efield_cs,
                "terminations": {"right": {"site": PROFILES.LD, "width": 1.0}},
            },
        ],
        figures_dir=figures_dir,
        ylim=(-120.0, 120.0),
    )


if __name__ == "__main__":
    main()

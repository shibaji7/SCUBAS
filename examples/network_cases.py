"""Electrically short cable demonstration matching the tutorial narrative.

The accompanying ``docs/tutorial/elsc.md`` file introduces two limiting
transmission-line cases.  This script reproduces the electrically short
configuration with the refactored SCUBAS API:

* configure Matplotlib styling through ``scubas.plotlib.update_rc_params``;
* build an ``OceanModel`` for the deep-ocean profile and inspect its
  transfer function; and
* evaluate the voltage induced along a short continental shelf cable.
"""

import os
from pathlib import Path

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
from scubas.models import OceanModel
from scubas.plotlib import (
    cable_potential,
    plot_transfer_function,
    update_rc_params,
)


def save_transfer_function(figures_dir: Path) -> None:
    """Render and persist the deep-ocean transfer function plot."""
    ocean_model = OceanModel(PROFILES.DO_3)
    transfer_function = ocean_model.get_TFs()
    artifacts = plot_transfer_function(transfer_function)
    artifacts.figure.suptitle("Deep Ocean Transfer Function")
    artifacts.figure.savefig(
        figures_dir / "electrically_cable_transfer_function.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(artifacts.figure)


def simulate_case(
    *,
    length_km: float,
    sec_id: str,
    case_slug: str,
    figures_dir: Path,
    y_limits: tuple[float, float] = (-200.0, 200.0),
) -> None:
    """Compute potentials for a cable section and persist the resulting figure."""
    induced_e_field = pd.DataFrame(
        {"X": [300.0]},
        index=pd.RangeIndex(1, name="Time"),
    )

    transmission_line = TransmissionLine(
        sec_id=sec_id,
        directed_length={"length_north": length_km},
        elec_params={"site": PROFILES.CS, "width": 1.0, "flim": [1e-6, 1.0]},
    )
    transmission_line.compute_eqv_pi_circuit(
        Efield=induced_e_field,
        components=["X"],
    )

    cable = Cable([transmission_line], components=["X"])
    potentials, distances = cable._pot_along_cable_(timestamp=0)

    plot = cable_potential(potentials, distances, ylim=y_limits)
    plot.axes.text(
        0.05,
        0.85,
        rf"$L_{{cs}}$={length_km:.0f} km",
        ha="left",
        va="center",
        transform=plot.axes.transAxes,
    )
    plot.figure.savefig(
        figures_dir / f"{case_slug}_cable_potential.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(plot.figure)


def main() -> None:
    """Generate plots for electrically long and short cable cases."""

    # Match the documentation styling so the generated figures are consistent.
    figures_dir = Path(__file__).resolve().parent.parent / "docs" / "tutorial" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    update_rc_params(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"],
        },
        science=True,
        latex=False,
    )

    save_transfer_function(figures_dir)

    ####################################################################
    # Electrically long cable:
    #   * Induced electric field of 0.3 V/km along the X component
    #   * Continental shelf 600 km in length, 100 m depth (profile CS)
    ####################################################################
    simulate_case(
        length_km=600.0,
        sec_id="CS-long",
        case_slug="electrically_long",
        figures_dir=figures_dir,
    )

    ####################################################################
    # Electrically short cable:
    #   * Induced electric field of 0.3 V/km along the X component
    #   * Continental shelf 4000 km in length, 100 m depth (profile CS)
    ####################################################################
    simulate_case(
        length_km=4000.0,
        sec_id="CS-short",
        case_slug="electrically_short",
        figures_dir=figures_dir,
    )


if __name__ == "__main__":
    main()

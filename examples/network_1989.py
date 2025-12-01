"""Reproduce the March 1989 event study using the refreshed SCUBAS API.

The tutorial walkthrough in ``docs/tutorial/1989.md`` highlights the TAT-8
trans-Atlantic cable response during the March 1989 geomagnetic storm.  This
script mirrors that workflow end-to-end with modern helpers:

* ingest geomagnetic station data from ``examples/datasets/1989``;
* build transmission-line sections that span the continental shelves, deep
  ocean segments, and the Mid-Atlantic Ridge; and
* export publication-ready plots of the induced electric fields and potentials
  directly into the documentation figures directory.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, List, Mapping, Sequence

os.environ["MPLBACKEND"] = "Agg"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import dates as mdates
from matplotlib.dates import DateFormatter

from scubas.cables import Cable, TransmissionLine
from scubas.datasets import PROFILES
from scubas.plotlib import update_rc_params

DATASET_DIR = Path(__file__).resolve().parent / "datasets" / "1989"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "docs" / "tutorial" / "figures"

EVENT_START = pd.Timestamp("1989-03-13 00:00:00")
EVENT_END = pd.Timestamp("1989-03-14 12:00:00")


SECTION_DEFINITIONS: Sequence[Mapping[str, Any]] = [
    {
        "sec_id": "CS-W",
        "profile": PROFILES.CS_W,
        "station": "FRD",
        "dataset": DATASET_DIR / "FRD.csv",
        "terminations": {"left": {"site": PROFILES.LD, "width": 1.0}},
        "edge_locations": {
            "initial": {"lat": 39.60, "lon": -74.33},
            "final": {"lat": 38.79, "lon": -72.62},
        },
    },
    {
        "sec_id": "DO-1",
        "profile": PROFILES.DO_1,
        "station": "FRD",
        "dataset": DATASET_DIR / "FRD.csv",
        "edge_locations": {
            "initial": {"lat": 38.79, "lon": -72.62},
            "final": {"lat": 37.11, "lon": -68.94},
        },
    },
    {
        "sec_id": "DO-2",
        "profile": PROFILES.DO_2,
        "station": "STJ",
        "dataset": DATASET_DIR / "STJ.csv",
        "edge_locations": {
            "initial": {"lat": 37.11, "lon": -68.94},
            "final": {"lat": 39.80, "lon": -48.20},
        },
    },
    {
        "sec_id": "DO-3",
        "profile": PROFILES.DO_3,
        "station": "STJ",
        "dataset": DATASET_DIR / "STJ.csv",
        "edge_locations": {
            "initial": {"lat": 39.80, "lon": -48.20},
            "final": {"lat": 40.81, "lon": -45.19},
        },
    },
    {
        "sec_id": "DO-4",
        "profile": PROFILES.DO_4,
        "station": "STJ",
        "dataset": DATASET_DIR / "STJ.csv",
        "edge_locations": {
            "initial": {"lat": 40.81, "lon": -45.19},
            "final": {"lat": 43.15, "lon": -39.16},
        },
    },
    {
        "sec_id": "DO-5",
        "profile": PROFILES.DO_5,
        "station": "STJ",
        "dataset": DATASET_DIR / "STJ.csv",
        "edge_locations": {
            "initial": {"lat": 43.15, "lon": -39.16},
            "final": {"lat": 44.83, "lon": -34.48},
        },
    },
    {
        "sec_id": "MAR",
        "profile": PROFILES.MAR,
        "station": "STJ",
        "dataset": DATASET_DIR / "STJ.csv",
        "edge_locations": {
            "initial": {"lat": 44.83, "lon": -34.48},
            "final": {"lat": 46.51, "lon": -22.43},
        },
    },
    {
        "sec_id": "DO-6",
        "profile": PROFILES.DO_6,
        "station": "HAD",
        "dataset": DATASET_DIR / "HAD.csv",
        "edge_locations": {
            "initial": {"lat": 46.51, "lon": -22.43},
            "final": {"lat": 47.85, "lon": -9.05},
        },
    },
    {
        "sec_id": "CS-E",
        "profile": PROFILES.CS_E,
        "station": "HAD",
        "dataset": DATASET_DIR / "HAD.csv",
        "terminations": {"right": {"site": PROFILES.LD, "width": 1.0}},
        "edge_locations": {
            "initial": {"lat": 47.85, "lon": -9.05},
            "final": {"lat": 50.79, "lon": -4.55},
        },
    },
]


def build_section(config: Mapping[str, Any]) -> TransmissionLine:
    """Instantiate, populate, and return a transmission-line section."""
    transmission_line = TransmissionLine(
        sec_id=config["sec_id"],
        directed_length={"edge_locations": config["edge_locations"]},
        elec_params={"site": config["profile"], "width": 1.0, "flim": [1e-6, 1.0]},
        active_termination=config.get("terminations"),
    )
    transmission_line.compile_oml(
        bfield_data_files=[config["dataset"]],
        csv_file_date_name="Date",
    )
    transmission_line.compute_eqv_pi_circuit(components=["X", "Y"])
    return transmission_line


def subset_parameters(frame: pd.DataFrame) -> pd.DataFrame:
    """Focus the parameter frame on the peak of the March 1989 storm."""
    if EVENT_START in frame.index and EVENT_END in frame.index:
        return frame.loc[EVENT_START:EVENT_END]
    return frame


def plot_e_fields(
    params: pd.DataFrame,
    section_labels: Sequence[str],
    station_labels: Sequence[str],
    figures_dir: Path,
) -> None:
    """Overlay E-field estimates for each section and station."""
    plotted_params = subset_parameters(params)
    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        dpi=200,
        figsize=(8, 6),
        sharex=True,
    )
    colors = plt.cm.tab10(np.linspace(0, 1, len(section_labels)))
    offsets = np.linspace(len(section_labels) - 1, 0, len(section_labels)) * 150.0

    for idx, (sec, station, offset, color) in enumerate(
        zip(section_labels, station_labels, offsets, colors)
    ):
        ex_column = f"E.X.{idx:02d}"
        ey_column = f"E.Y.{idx:02d}"
        if ex_column not in plotted_params or ey_column not in plotted_params:
            continue
        ex_series = plotted_params[ex_column] - plotted_params[ex_column].mean()
        ey_series = plotted_params[ey_column] - plotted_params[ey_column].mean()

        axes[0].plot(
            plotted_params.index,
            ex_series + offset,
            color=color,
            linewidth=0.9,
            ls="-",
        )
        axes[1].plot(
            plotted_params.index,
            ey_series + offset,
            color=color,
            linewidth=0.9,
            ls="-",
            label=f"{sec} [{station}]",
        )

    for ax, label in zip(axes, [r"$E_x$ offset (mV/km)", r"$E_y$ offset (mV/km)"]):
        ax.set_ylabel(label)
        ax.set_ylim(
            ax.get_ylim()[0] - 100,
            ax.get_ylim()[1] + 100,
        )
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
        ax.xaxis.set_minor_locator(mdates.HourLocator())
        ax.xaxis.set_major_formatter(DateFormatter("%H"))

    axes[0].set_title("De-trended electric field estimates along TAT-8")
    axes[1].legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=False,
    )
    axes[1].set_xlabel("Time (UT)")
    fig.tight_layout()
    fig.savefig(
        figures_dir / "event_1989_e_fields.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_potentials(params: pd.DataFrame, figures_dir: Path) -> None:
    """Visualise cumulative potentials and end-point voltages."""
    plotted_params = subset_parameters(params)
    fig, ax = plt.subplots(dpi=200, figsize=(8, 4))

    ax.plot(
        plotted_params.index,
        -1 * plotted_params["V(v)"],
        label="Section potential sum",
        color="tab:green",
        ls="-",
        linewidth=1.2,
    )
    if "Vt(v)" in plotted_params.columns:
        ax.plot(
            plotted_params.index,
            -1 * plotted_params["Vt(v)"],
            label=r"$V_{TAT-8}$",
            color="tab:purple",
            ls="-",
            linewidth=1.0,
        )
    ax.plot(
        plotted_params.index,
        -1 * plotted_params["U0"],
        label=r"$U_W$ (west termination)",
        color="tab:red",
        ls="-",
        linewidth=1.0,
    )
    ax.plot(
        plotted_params.index,
        -1 * plotted_params["U1"],
        label=r"$U_E$ (east termination)",
        color="tab:blue",
        ls="-",
        linewidth=1.0,
    )

    ax.set_ylabel("Potential (V)")
    ax.set_xlabel("Time (UT)")
    ax.set_title("TAT-8 induced potentials during March 1989 storm")
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
    ax.xaxis.set_minor_locator(mdates.HourLocator())
    ax.xaxis.set_major_formatter(DateFormatter("%H"))
    ax.grid(alpha=0.3, linewidth=0.5)
    ax.set_ylim(-700, 500)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(
        figures_dir / "event_1989_potentials.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def main() -> None:
    """Run the March 1989 reproduction and save all artefacts."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    update_rc_params(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"],
        },
        science=True,
        latex=False,
    )

    sections: List[TransmissionLine] = [
        build_section(section_config) for section_config in SECTION_DEFINITIONS
    ]
    cable = Cable(sections, components=["X", "Y"])
    params = cable.tot_params

    section_labels = [config["sec_id"] for config in SECTION_DEFINITIONS]
    station_labels = [config["station"] for config in SECTION_DEFINITIONS]
    plot_e_fields(params, section_labels, station_labels, FIGURES_DIR)
    plot_potentials(params, FIGURES_DIR)


if __name__ == "__main__":
    main()

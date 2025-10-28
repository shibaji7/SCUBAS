"""
Plotting helpers for transfer functions and cable potentials leveraged across SCUBAS.
"""

from __future__ import annotations
__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np

from scubas.utils import frexp102str


@dataclass(frozen=True)
class PlotArtifacts:
    """
    Container for plot objects returned to callers.
    """

    figure: plt.Figure
    axes: Any


def update_rc_params(
    params: Optional[Mapping[str, Any]] = None,
    science: bool = False,
) -> None:
    """
    Update matplotlib rcParams, optionally enabling the ``science`` style.
    """
    if science:
        plt.style.use(["science", "ieee"])
        plt.rcParams.update(
            {
                "figure.figsize": np.array([8, 6]),
                "text.usetex": True,
                "font.family": "sans-serif",
                "font.sans-serif": [
                    "Tahoma",
                    "DejaVu Sans",
                    "Lucida Grande",
                    "Verdana",
                ],
                "font.size": 10,
            }
        )
    if params:
        plt.rcParams.update(params)


def plot_transfer_function(
    Tx: Any,
    xlim: Sequence[float] = (1e-6, 1e-2),
    ylim: Sequence[float] = (1e-3, 1e0),
    science: bool = False,
) -> PlotArtifacts:
    """
    Plot amplitude and phase of an ``E2B`` transfer function.
    """
    update_rc_params(science=science)
    fig = plt.figure(dpi=180, figsize=(3, 2.5))
    ax_amp = fig.add_subplot(111)
    ax_amp.loglog(Tx.freq, np.abs(Tx.E2B), "r", lw=0.6, ls="-")
    ax_amp.set_xlabel(r"Frequency [Hz]")
    ax_amp.set_ylabel(r"Amplitude [mV/km/nT]", color="r")
    ax_amp.set_xlim(xlim)
    ax_amp.set_ylim(ylim)

    ax_phase = ax_amp.twinx()
    ax_phase.semilogx(Tx.freq, np.angle(Tx.E2B, deg=True), "b", lw=0.6, ls="-")
    ax_phase.set_ylabel(r"Phase [degree, $^\circ$]", color="b")
    ax_phase.set_ylim(-90, 90)
    ax_phase.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax_phase.set_yticklabels([-90, -60, -30, 0, 30, 60, 90])
    ax_phase.set_xlim(xlim)
    return PlotArtifacts(figure=fig, axes=ax_phase)


def potential_along_section(
    V: Sequence[float],
    x: Sequence[float],
    sec: Optional[int] = None,
    Vi: Optional[float] = None,
    Vk: Optional[float] = None,
    Z: Optional[float] = None,
    Y: Optional[float] = None,
    gma: Optional[float] = None,
    Z0: Optional[float] = None,
    science: bool = False,
) -> PlotArtifacts:
    """
    Visualise along-section potential with optional annotations.
    """
    update_rc_params(
        {"xtick.labelsize": 12, "ytick.labelsize": 12, "font.size": 12},
        science,
    )
    fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        dpi=150,
        figsize=(6, 3),
        sharex="all",
        sharey="all",
    )
    ax.set_ylabel("Voltage, V")
    ax.set_xlabel("Cable Length, km")
    ax.plot(x, V, "k", lw=0.8, ls="-")

    Z_scaled = Z * 1e3 if Z is not None else None
    Y_scaled = Y * 1e3 if Y is not None else None
    gma_scaled = gma * 1e3 if gma is not None else None

    details = []
    if sec is not None:
        details.append(f"Along: Bin{sec:02d}")
    if Vi is not None and Vk is not None:
        details.append(rf"$V_i,V_k\sim {Vi:.1f} V, {Vk:.1f} V$")
    if Z_scaled is not None and Y_scaled is not None:
        details.append(
            rf"$Z,Y\sim$ {frexp102str(Z_scaled)} $\Omega/km$, "
            rf"{frexp102str(Y_scaled)} $\mho/km$"
        )
    if gma_scaled is not None and Z0 is not None:
        details.append(
            rf"$\gamma,Z_0\sim$ {frexp102str(gma_scaled)} /km, {frexp102str(Z0)} $\Omega$"
        )
    details.append(f"L={np.max(x):.0f} km")

    ax.text(
        0.05,
        0.95,
        "\n".join(details),
        ha="left",
        va="top",
        transform=ax.transAxes,
        fontsize="small",
    )
    ax.set_xlim(x[0], x[-1])
    return PlotArtifacts(figure=fig, axes=ax)


def cable_potential(
    V: Sequence[float],
    x: Sequence[float],
    science: bool = False,
    ylim: Sequence[float] = (-50, 50),
) -> PlotArtifacts:
    """
    Plot potential along the full cable.
    """
    update_rc_params(
        {"xtick.labelsize": 12, "ytick.labelsize": 12, "font.size": 12},
        science,
    )
    fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        dpi=150,
        figsize=(6, 3),
        sharex="all",
        sharey="all",
    )
    ax.set_ylabel("Earth potential, V")
    ax.set_xlabel("Distance, km")
    ax.plot(x, V, "k", lw=0.8, ls="-")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(ylim)
    return PlotArtifacts(figure=fig, axes=ax)


#!/usr/bin/env python

"""
    plotlib.py: This is a data holder module with various datatypes
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import matplotlib.pyplot as plt
import numpy as np

from scubas.utils import frexp102str


def update_rc_params(params=dict(), science=False):
    """ """
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
    if len(params.keys()) > 0:
        plt.rcParams.update(params)
    return


def plot_transfer_function(Tx, xlim=[1e-6, 1e-2], ylim=[1e-3, 1e0], science=False):
    """ """
    update_rc_params(science=science)
    fig = plt.figure(dpi=180, figsize=(3, 2.5))
    ax = fig.add_subplot(111)
    ax.loglog(Tx.freq, np.abs(Tx.E2B), "r", lw=0.6, ls="-")
    ax.set_xlabel(r"Frequency [Hz]")
    ax.set_ylabel(r"Amplitude [mV/km/nT]", color="r")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax = ax.twinx()
    ax.semilogx(Tx.freq, np.angle(Tx.E2B, deg=True), "b", lw=0.6, ls="-")
    ax.set_ylabel(r"Phase [degree, $^\circ$]", color="b")
    ax.set_ylim(-90, 90)
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_yticklabels([-90, -60, -30, 0, 30, 60, 90])
    _ = ax.set_xlim(xlim)
    return {"figure": fig, "axes": ax}


def potential_along_section(
    V,
    x,
    sec=None,
    Vi=None,
    Vk=None,
    Z=None,
    Y=None,
    gma=None,
    Z0=None,
    science=False,
):
    update_rc_params(
        {"xtick.labelsize": 12, "ytick.labelsize": 12, "font.size": 12}, science
    )
    fig, axes = plt.subplots(
        nrows=1, ncols=1, dpi=150, figsize=(6, 3), sharex="all", sharey="all"
    )
    ax = axes
    ax.set_ylabel("Voltage, V")
    ax.set_xlabel("Cable Length, km")
    ax.plot(x, V, "k", lw=0.8, ls="-")
    if Z is not None:
        Z *= 1e3
    if Y is not None:
        Y *= 1e3
    if gma is not None:
        gma *= 1e3
    txt = ""
    if sec is not None:
        txt += "Along: Bin%02d\n" % (sec)
    if (Vi is not None) and (Vk is not None):
        txt += r"$V_i,V_k\sim %.1f V, %.1f V$" % (Vi, Vk) + "\n"
    if (Z is not None) and (Y is not None):
        txt += (
            r"$Z,Y\sim$ %s $\Omega/km$, %s $\mho/km$" % (frexp102str(Z), frexp102str(Y))
            + "\n"
        )
    if (gma is not None) and (Z0 is not None):
        txt += (
            r"$\gamma,Z_0\sim$ %s /km, %s $\Omega$"
            % (frexp102str(gma), frexp102str(Z0))
            + "\n"
        )
    txt += "L=%d km" % np.max(x)
    ax.text(
        0.05, 0.95, txt, ha="left", va="top", transform=ax.transAxes, fontsize="small"
    )
    ax.set_xlim(x[0], x[-1])
    return {"figure": fig, "axes": ax}


def cable_potential(V, x, science=False):
    update_rc_params(
        {"xtick.labelsize": 12, "ytick.labelsize": 12, "font.size": 12}, science
    )
    fig, axes = plt.subplots(
        nrows=1, ncols=1, dpi=150, figsize=(6, 3), sharex="all", sharey="all"
    )
    ax = axes
    ax.set_ylabel("Voltage, V")
    ax.set_xlabel("Cable Length, km")
    ax.plot(x, V, "k", lw=0.8, ls="-")
    ax.set_xlim(x[0], x[-1])
    txt = ""
    ax.text(
        0.05, 0.95, txt, ha="left", va="top", transform=ax.transAxes, fontsize="small"
    )
    return {"figure": fig, "axes": ax}

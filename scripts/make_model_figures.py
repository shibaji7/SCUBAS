#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Regenerate the figures of the Models pages in ``docs/model/figures`` from the SCUBAS code.

Every plot curve is computed by the package (``OceanModel``, ``TransmissionLine``, ``Cable``); where a paper gives a closed form, the formula
of ``scripts/model_reference.py`` is overlaid as markers, so that agreement is visible and not asserted. The schematics are drawn with
Matplotlib. Run ``python scripts/make_model_figures.py [OUTDIR]`` (default ``docs/model/figures``; ``make model-figures``). The same numbers
that the pages quote are written to ``numbers.json`` next to the figures (see ``scripts/model_numbers.py``).

Style: the SCUBAS figure palette (deep blue, red, golden yellow, grey on white), axes labelled with units, 150 dpi PNG.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

from model_cases import (  # noqa: E402
    atlantic_site,
    network_case,
    potentials,
    single_section,
)
from model_reference import (  # noqa: E402
    analytic_efield,
    bd_over_b0,
    end_potential,
    field_ratios_with_depth,
    section_potential,
    surface_impedance,
    synthetic_field,
)
from scubas.datasets import PROFILES, Layer, Site  # noqa: E402
from scubas.models import OceanModel  # noqa: E402
from scubas.sources.current_sources import Current  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "docs" / "model" / "figures"

PALETTE = {
    "deep_blue": "#003366",
    "golden_yellow": "#FFCC00",
    "white": "#FFFFFF",
    "light_grey": "#F0F0F0",
    "accent_red": "#990000",
    "accent_grey": "#555555",
}
BLUE, RED, YELLOW, GREY = (
    PALETTE[k] for k in ("deep_blue", "accent_red", "golden_yellow", "accent_grey")
)

# C2022 Table 4: the synthetic magnetic record (amplitude nT, period min, phase deg).
SYNTH_A = [200, 90, 30, 17, 8, 3.5]
SYNTH_T = [180, 80, 36, 15, 8, 3]
SYNTH_PHI = [10, 20, 30, 40, 50, 60]


def apply_style() -> None:
    """The SCUBAS figure style (palette, light grid, no legend frame)."""
    mpl.rcParams.update(
        {
            "figure.facecolor": PALETTE["white"],
            "axes.facecolor": PALETTE["white"],
            "axes.edgecolor": GREY,
            "axes.grid": True,
            "grid.color": PALETTE["light_grey"],
            "grid.linewidth": 0.8,
            "axes.prop_cycle": mpl.cycler(color=[BLUE, RED, YELLOW, GREY]),
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.frameon": False,
            "savefig.dpi": 150,
            "savefig.facecolor": PALETTE["white"],
            "lines.linewidth": 2.0,
        }
    )


def save(fig: plt.Figure, out: Path, name: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / name, bbox_inches="tight")
    plt.close(fig)


def tf_of(site: Site, freqs: np.ndarray, kind: str = "E2B") -> np.ndarray:
    """A transfer function of ``OceanModel`` on the given frequencies."""
    return np.asarray(OceanModel(site, flim=[1e-6, 1.0]).calcTF(freqs=freqs)[kind])


def half_space_site(rho_sea: float, depth_m: float, rho_earth: float) -> Site:
    """Seawater of the given depth over a uniform half-space, the model of BP2003 Section 'Seafloor/Surface Transfer Functions'."""
    return Site(
        [
            Layer("Seawater", depth_m, 1.0 / rho_sea),
            Layer("Earth", 1.0e6, 1.0 / rho_earth),
        ],
        "BP2003 half-space",
        "BP2003",
    )


# --------------------------------------------------------------------------------------------------------------------------------------
# Schematics
# --------------------------------------------------------------------------------------------------------------------------------------


def _box(ax, xy, w, h, text, fc, ec=BLUE, size=7.5, color=BLUE):
    ax.add_patch(
        FancyBboxPatch(
            xy, w, h, boxstyle="round,pad=0.02,rounding_size=0.06", fc=fc, ec=ec, lw=1.4
        )
    )
    ax.text(
        xy[0] + w / 2,
        xy[1] + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=size,
        color=color,
    )


def _arrow(ax, a, b, color=GREY, lw=1.4):
    ax.add_patch(
        FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=12, color=color, lw=lw)
    )


def fig_chain(out: Path) -> None:
    """The chain from the magnetic record to the cable voltage, with the class that does each step."""
    fig, ax = plt.subplots(figsize=(9.0, 4.9))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5.4)
    ax.axis("off")
    ax.grid(False)
    soft, warm, ec2, c2 = "#E8EEF5", "#FFF6CC", "#8A6D00", "#5A4600"
    xs = [0.05, 2.1, 4.15, 6.2, 8.25]
    w = 1.7
    top = [
        "Magnetic record\nB(t), observatory\n(IAGA or CSV)",
        "Detrend, taper,\nFFT\n(Preprocess, fft)",
        "× transfer function\nT_x(f) per section\n(OceanModel)",
        "inverse FFT\nE(t), seafloor\n(ifft)",
        "emf of a section\nE_N L_N + E_E L_E\n(CableSection)",
    ]
    bottom = [
        "Layered Earth\nper section\n(Site, LITHO1.0)",
        "Z, Y, γ, Z0\n(thin double layer)\n(DSTL)",
        "equivalent-π\ncircuit\n(Ye, Yp2, Ie)",
        "nodal network\n[J] = [Y][U]\n(Cable)",
        "end potentials\nU_A, U_B\n(U0, U1)",
    ]
    for x, t in zip(xs, top, strict=True):
        _box(ax, (x, 3.8), w, 1.1, t, soft, size=8.5)
    for x, t in zip(xs, bottom, strict=True):
        _box(ax, (x, 1.75), w, 1.1, t, warm, ec=ec2, size=8.5, color=c2)
    for i in range(4):
        _arrow(ax, (xs[i] + w, 4.35), (xs[i + 1], 4.35))
        _arrow(ax, (xs[i] + w, 2.3), (xs[i + 1], 2.3))
    _arrow(ax, (xs[3] + w / 2, 3.8), (xs[2] + w / 2, 2.85), color=BLUE)
    ax.text(5.55, 3.32, "E/Z: current sources", fontsize=8.5, color=BLUE, ha="center")
    _box(
        ax,
        (2.6, 0.1),
        4.8,
        0.95,
        "Cable voltage   V_C = Σ emf + U_A − U_B\n(tot_params['Vt(v)'])",
        "#DCE8F4",
        ec=BLUE,
        size=10,
    )
    _arrow(ax, (xs[4] + w / 2, 3.8), (7.4, 0.85), color=BLUE)
    _arrow(ax, (xs[4] + w / 2 - 0.3, 1.75), (7.4, 0.5), color=ec2)
    ax.text(
        0.05,
        0.6,
        "blue: the field along the cable\nyellow: the Earth under the cable",
        fontsize=8.5,
        color=GREY,
        va="center",
    )
    save(fig, out, "chain.png")


def fig_voltage_assembly(out: Path) -> None:
    """How the cable voltage is made: the emf in the cable and the Earth potentials at the two ends."""
    fig, ax = plt.subplots(figsize=(8.6, 3.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(-1.7, 2.2)
    ax.axis("off")
    ax.grid(False)
    ax.fill_between([0, 10], -1.7, -0.05, color="#DCE8F4")
    ax.text(
        5, -1.35, "sea water and the Earth below", ha="center", color=BLUE, fontsize=9
    )
    ax.fill_between([0, 0.9], -0.05, 0.25, color="#CBB88A")
    ax.fill_between([9.1, 10], -0.05, 0.25, color="#CBB88A")
    ax.text(0.45, 0.45, "coast A", ha="center", fontsize=9, color=GREY)
    ax.text(9.55, 0.45, "coast B", ha="center", fontsize=9, color=GREY)
    ax.plot([0.45, 1.2, 8.8, 9.55], [0.0, -0.05, -0.05, 0.0], color=BLUE, lw=3)
    ax.text(
        5,
        -0.35,
        "cable (copper conductor of the power feed)",
        ha="center",
        fontsize=9,
        color=BLUE,
    )
    _arrow(ax, (2.0, 0.55), (8.0, 0.55), color=RED, lw=2.2)
    ax.text(
        5,
        0.8,
        "induced electric field E along the cable:  emf = Σ E·L",
        ha="center",
        color=RED,
        fontsize=10,
    )
    for x, lab, col in [(0.45, "U_A", GREY), (9.55, "U_B", GREY)]:
        ax.plot([x, x], [0.0, 1.55], color=col, lw=1.2, ls=":")
        ax.text(x, 1.7, f"Earth potential\n{lab}", ha="center", fontsize=9, color=col)
    ax.text(
        5,
        1.8,
        "V_C = emf + U_A − U_B",
        ha="center",
        fontsize=12,
        color=BLUE,
        weight="bold",
    )
    save(fig, out, "voltage_assembly.png")


# --------------------------------------------------------------------------------------------------------------------------------------
# Seafloor field (SCUBAS page)
# --------------------------------------------------------------------------------------------------------------------------------------


def fig_attenuation_with_depth(out: Path) -> None:
    """E and B against depth in the sea for 5 min (BP2003 Figure 1), the seafloor values from the code marked."""
    f, rho_s, rho_e = 1.0 / 300.0, 0.25, 10.0
    fig, axes = plt.subplots(3, 1, figsize=(6.2, 7.2), sharex=True)
    for ax, d in zip(axes, [100.0, 1000.0, 5000.0], strict=True):
        z = np.linspace(0, d, 200)
        z_d = surface_impedance([rho_e], [1e6], f)
        e_ratio, b_ratio = field_ratios_with_depth(rho_s, d, z_d, f, z)
        ax.plot(abs(e_ratio), -z, color=BLUE, label="electric field E(z)/E₀")
        ax.plot(abs(b_ratio), -z, color=RED, ls="--", label="magnetic field B(z)/B₀")
        code_b = abs(
            OceanModel(half_space_site(rho_s, d, rho_e), flim=[1e-6, 1.0]).calcTF(
                ["B2B"], freqs=[f]
            )["B2B"][0]
        )
        ax.plot(
            [code_b],
            [-d],
            "o",
            color=YELLOW,
            mec=GREY,
            ms=7,
            label="B at the seafloor from OceanModel" if d == 100.0 else None,
        )
        ax.set_ylabel("Depth (m)")
        ax.set_xlim(0, 1.2)
        ax.text(
            0.03,
            0.12,
            f"ocean depth {d / 1000:g} km",
            transform=ax.transAxes,
            color=GREY,
        )
    axes[0].legend(loc="lower left", bbox_to_anchor=(0.0, 0.25))
    axes[-1].set_xlabel("Amplitude ratio to the sea surface")
    fig.suptitle(
        "Fields in the sea for a 5 min period (seawater 0.25 Ω m, Earth 10 Ω m)",
        y=0.995,
        fontsize=10,
    )
    fig.tight_layout()
    save(fig, out, "attenuation_with_depth.png")


def fig_seafloor_ratios(out: Path) -> None:
    """B_d/B_0 and E_d/B_0 against frequency for 100 m and 5 km (BP2003 Figures 4 and 5)."""
    f = np.logspace(-4, -2, 120)
    rho_s, rho_e = 0.25, 10.0
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9))
    for d, col in [(100.0, BLUE), (5000.0, RED)]:
        site = half_space_site(rho_s, d, rho_e)
        b = abs(tf_of(site, f, "B2B"))
        e = abs(tf_of(site, f, "E2B"))
        axes[0].loglog(f, b, color=col, label=f"{d / 1000:g} km")
        axes[1].loglog(f, e, color=col, label=f"{d / 1000:g} km")
        axes[0].loglog(
            f[::12],
            abs(
                bd_over_b0(
                    rho_s, d, surface_impedance([rho_e], [1e6], f[::12]), f[::12]
                )
            ),
            "o",
            color=col,
            ms=4,
            mfc="none",
        )
    z_d = surface_impedance([rho_e], [1e6], f)
    axes[1].loglog(
        f,
        abs(z_d) * 1e-3 / (4e-7 * np.pi),
        color=GREY,
        ls="--",
        label="zero depth (the seafloor response)",
    )
    axes[0].set(
        xlabel="Frequency (Hz)",
        ylabel="B at the seafloor / B at the sea surface",
        ylim=(0.01, 1.2),
    )
    axes[1].set(
        xlabel="Frequency (Hz)",
        ylabel="E at the seafloor / B at the sea surface (mV/km/nT)",
        ylim=(0.01, 4.0),
    )
    axes[0].set_title("magnetic field (markers: BP2003 Eq 13)")
    axes[1].set_title("electric field")
    axes[0].legend(title="ocean depth", loc="lower left")
    axes[1].legend(title="ocean depth", loc="upper left")
    fig.tight_layout()
    save(fig, out, "seafloor_ratios.png")


def fig_transfer_function_depth(out: Path) -> None:
    """The transfer function of the Atlantic model for several ocean depths (C2022 Figure 7)."""
    f = np.logspace(-6, -2, 200)
    fig, (a, p) = plt.subplots(2, 1, figsize=(6.4, 6.0), sharex=True)
    styles = [("-", BLUE), ("--", RED), (":", GREY), ("-.", "#B38F00")]
    for (ls, col), d in zip(styles, [0.01, 0.1, 1.0, 4.0], strict=True):
        tf = tf_of(atlantic_site(d), f)
        a.loglog(f, abs(tf), ls=ls, color=col, label=f"d = {d * 1000:g} m")
        p.semilogx(f, np.degrees(np.angle(tf)), ls=ls, color=col)
    a.set_ylabel("|T_x|  (mV/km/nT)")
    a.legend(loc="upper left")
    p.set(xlabel="Frequency (Hz)", ylabel="phase of T_x (°)")
    a.set_title("Transfer function of the Atlantic model for four ocean depths")
    fig.tight_layout()
    save(fig, out, "transfer_function_depth.png")


def fig_synthetic_test(out: Path) -> None:
    """The six-sine test of C2022 Section 2.4.2: the record, the analytic field, and the field SCUBAS returns."""
    site = atlantic_site(4.0)
    n = 72 * 60
    t = np.arange(n) * 60.0
    b = synthetic_field(t, SYNTH_A, SYNTH_T, SYNTH_PHI)
    freqs = [1.0 / (p * 60.0) for p in SYNTH_T]
    tf = [tf_of(site, np.array([x]))[0] for x in freqs]
    e_an = analytic_efield(t, SYNTH_A, SYNTH_T, SYNTH_PHI, tf)
    frame = pd.DataFrame(
        {"X": 0.0 * b, "Y": b},
        index=pd.date_range("2020-01-01", periods=n, freq="60s", name="Time"),
    )
    model = OceanModel(site, flim=[1e-6, 1.0])
    model.to_Efields(Bfield=frame, components=["X", "Y"], p=0.1)
    e_num = np.asarray(model.Efield["X"])
    keep = slice(n // 10, -n // 10)  # the tapered ends are left out, as in C2022
    slope = np.polyfit(e_an[keep], e_num[keep], 1)[0]
    corr = np.corrcoef(e_an[keep], e_num[keep])[0, 1]
    fig = plt.figure(figsize=(9.6, 6.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.2, 1])
    ax_b = fig.add_subplot(gs[0, 0])
    ax_e = fig.add_subplot(gs[1, 0], sharex=ax_b)
    ax_s = fig.add_subplot(gs[:, 1])
    hrs = t / 3600
    ax_b.plot(hrs, b, color=BLUE)
    ax_b.set_ylabel("B (nT)")
    ax_b.set_title("Synthetic magnetic record, C2022 Eq 24")
    ax_e.plot(hrs, e_an, color=RED, label="analytic, C2022 Eq 25")
    ax_e.plot(hrs, e_num, color=BLUE, lw=1.2, label="OceanModel.to_Efields")
    ax_e.set(xlabel="Time (h)", ylabel="E (mV/km)")
    ax_e.legend(loc="lower center", bbox_to_anchor=(0.5, -0.38), ncol=2)
    ax_s.plot(e_an[keep], e_num[keep], ".", color=BLUE, ms=2)
    lim = np.array([e_an[keep].min(), e_an[keep].max()])
    ax_s.plot(lim, lim, color=RED, label="slope 1")
    ax_s.plot(lim, slope * lim, color=GREY, ls="--", label=f"fit: slope {slope:.2f}")
    ax_s.set(
        xlabel="analytic E (mV/km)",
        ylabel="SCUBAS E (mV/km)",
        title=f"correlation {corr:.5f}",
    )
    ax_s.legend(loc="upper left")
    fig.tight_layout()
    save(fig, out, "synthetic_test.png")


# --------------------------------------------------------------------------------------------------------------------------------------
# DSTL page
# --------------------------------------------------------------------------------------------------------------------------------------


def fig_line_parameters_depth(out: Path) -> None:
    """Z0 and the adjustment distance against ocean depth for the Atlantic model (C2022 Figure 8)."""
    depths = np.linspace(0.1, 5.0, 50)
    z0 = np.empty_like(depths)
    adj = np.empty_like(depths)
    for i, d in enumerate(depths):
        line = single_section(atlantic_site(d), 100.0, 1.0)
        z0[i], adj[i] = float(line.Z0), 1.0 / float(line.gma) / 1e3
    fig, (a, b) = plt.subplots(2, 1, figsize=(5.8, 5.6), sharex=True)
    a.plot(depths, z0, color=RED)
    a.set_ylabel("Z₀ (Ω)")
    b.plot(depths, adj, color=BLUE)
    b.set(xlabel="Ocean depth (km)", ylabel="adjustment distance 1/γ (km)")
    a.set_title("Line constants of the Atlantic model against ocean depth")
    fig.tight_layout()
    save(fig, out, "line_parameters_depth.png")


def fig_long_short_limits(out: Path) -> None:
    """The end potential in units of E/γ against γL: one curve that is γL/2 for short and 1 for long sections."""
    x = np.logspace(-2, 1.2, 300)
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.loglog(x, np.tanh(x / 2), color=BLUE, label="|U_A| γ/E = tanh(γL/2)")
    ax.loglog(x, x / 2, color=RED, ls="--", label="short section: γL/2  (|U_A| = LE/2)")
    ax.axhline(1.0, color=GREY, ls="--", label="long section: 1  (|U_A| = E/γ)")
    ax.axvline(1.0, color=GREY, ls=":", lw=1)
    ax.axvline(4.0, color=GREY, ls=":", lw=1)
    ax.text(1.06, 1.15, "L = 1/γ", color=GREY, fontsize=8)
    ax.text(4.3, 1.15, "L = 4/γ", color=GREY, fontsize=8)
    ax.set(
        xlabel="γ L  (section length / adjustment distance)",
        ylabel="|U_A| / (E/γ)",
        ylim=(0.005, 2.0),
    )
    ax.legend(loc="lower right")
    save(fig, out, "long_short_limits.png")


def fig_section_potentials(out: Path) -> None:
    """The potential along an electrically short and an electrically long section, E = 1 V/km (C2022 Figure 6)."""
    site = atlantic_site(4.0)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8))
    for ax, length, title in [
        (axes[0], 100.0, "electrically short, L = 100 km"),
        (axes[1], 10000.0, "electrically long, L = 10 000 km"),
    ]:
        line = single_section(site, length, 1000.0)
        cable = network_case([(site, length, 1000.0)], land_e_mv_per_km=None)
        v, d = potentials(cable)
        gamma = float(line.gma) * 1e3
        u_a = end_potential(1.0, length, gamma)
        ax.plot(d, v, color=BLUE, label="SCUBAS (Cable)")
        theory = section_potential(
            d, length, gamma, u_a, -u_a
        )  # C2022 Eq 10 with U_B = -U_A
        ax.plot(d[::6], theory[::6], "+", color=RED, label="C2022 Eqs 10 to 12")
        ax.set(
            xlabel="Distance along the section (km)",
            ylabel="Potential (V)",
            title=title,
        )
        ax.legend(loc="upper left")
    fig.tight_layout()
    save(fig, out, "section_potentials.png")


def three_models():
    """The three example models of B2023 Figure 9: (name, sections) with fields of 0.1 V/km in the deep ocean and 0.3 V/km elsewhere."""
    cs, so, do = PROFILES.CS, PROFILES.SO, PROFILES.DO
    return [
        ("Shallow sea, 600 km, 100 m deep", [(cs, 600.0, 300.0)]),
        (
            "Shallow ocean, 2000 km, 1 km deep",
            [(cs, 100.0, 300.0), (so, 1800.0, 100.0), (cs, 100.0, 300.0)],
        ),
        (
            "Deep ocean, 8000 km, 4 km deep",
            [(cs, 100.0, 300.0), (do, 7800.0, 100.0), (cs, 100.0, 300.0)],
        ),
    ]


def fig_network_models(out: Path) -> None:
    """The earth potential along the cable for the three example models of B2023 (its Figure 10)."""
    fig, axes = plt.subplots(3, 1, figsize=(6.4, 8.2))
    for ax, (title, sections) in zip(axes, three_models(), strict=True):
        v, d = potentials(network_case(sections))
        ax.plot(d, v, color=BLUE)
        ax.set(ylabel="Earth potential (V)", title=title)
        ax.annotate(
            f"{v[0]:.1f} V",
            (d[0], v[0]),
            xytext=(8, 4),
            textcoords="offset points",
            color=RED,
        )
        ax.annotate(
            f"{v[-1]:.1f} V",
            (d[-1], v[-1]),
            xytext=(-48, -14),
            textcoords="offset points",
            color=RED,
        )
    axes[-1].set_xlabel("Distance along the cable (km)")
    fig.tight_layout()
    save(fig, out, "network_models.png")


def fig_coast_potentials(out: Path) -> None:
    """Coast potential against the width of a shallow sea and of the continental shelf (B2023 Figures 11 and 12)."""
    cs, do = PROFILES.CS, PROFILES.DO
    widths_sea = np.arange(10.0, 2001.0, 70.0)
    uc_sea = [abs(potentials(network_case([(cs, w, 300.0)]))[0][0]) for w in widths_sea]
    widths_shelf = np.arange(5.0, 2001.0, 70.0)
    uc, ue = [], []
    for w in widths_shelf:
        v, d = potentials(
            network_case([(cs, w, 300.0), (do, 10000.0, 100.0), (cs, w, 300.0)])
        )
        uc.append(abs(v[0]))
        ue.append(abs(v[np.argmin(abs(d - w))]))
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))
    axes[0].plot(widths_sea, uc_sea, color=BLUE)
    axes[0].axhline(77.46, color=RED, ls="--", label="limit 77.46 V (B2023 Eq 57)")
    axes[0].set(
        xlabel="Width of the shallow sea (km)",
        ylabel="Coast potential (V)",
        title="Shallow sea, 0.3 V/km everywhere",
    )
    axes[0].legend(loc="lower right")
    axes[1].plot(widths_shelf, uc, color=BLUE, label="at the coast, U_C")
    axes[1].plot(widths_shelf, ue, color=GREY, label="at the shelf edge, U_E")
    axes[1].axhline(77.46, color=RED, ls="--", lw=1.2)
    axes[1].axhline(69.38, color=RED, ls=":", lw=1.2)
    axes[1].text(40, 77.9, "77.46 V (Eq 57)", color=RED, fontsize=8)
    axes[1].text(40, 69.9, "69.38 V (Eq 59)", color=RED, fontsize=8)
    axes[1].set(
        xlabel="Width of the continental shelf (km)",
        ylabel="Earth potential (V)",
        title="Deep ocean 0.1 V/km, shelf and land 0.3 V/km",
    )
    axes[1].legend(loc="upper right")
    fig.tight_layout()
    save(fig, out, "coast_potentials.png")


def fig_line_current(out: Path) -> None:
    """The magnetic field at the Earth's surface of a line current at 110 km height, and its induced part (the source model of ``Current``)."""
    x_km = np.linspace(-1000.0, 1000.0, 401)
    current = Current(
        h=110.0,
        a=0.0,
        I=1.0,
        omega=np.array([2 * np.pi / 300.0]),
        site=PROFILES.Quebec,
        current_unit="MA",
    )
    current.compute_B_field(
        x=x_km.copy()
    )  # the method rescales the array it is given, a known defect: pass a copy
    bx, bz = current.B_fields
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8), sharex=True)
    mu0_i_over_2pi = 4e-7 * np.pi * 1e6 / (2 * np.pi)
    h_m, x_m = 110e3, x_km * 1e3
    for ax, comp, ampere, title in [
        (axes[0], bx, mu0_i_over_2pi * h_m / (h_m**2 + x_m**2), "north component B_x"),
        (
            axes[1],
            bz,
            mu0_i_over_2pi * x_m / (h_m**2 + x_m**2),
            "vertical component B_z",
        ),
    ]:
        ax.plot(
            x_km, abs(comp["external"]) * 1e9, color=BLUE, label="primary (the current)"
        )
        ax.plot(
            x_km,
            abs(comp["internal"]) * 1e9,
            color=GREY,
            ls="--",
            label="induced (the image)",
        )
        ax.plot(x_km, abs(comp["total"]) * 1e9, color=RED, label="total")
        ax.plot(
            x_km[::25],
            abs(ampere[::25]) * 1e9,
            "o",
            color=YELLOW,
            mec=GREY,
            ms=6,
            label="Ampère's law, primary" if comp is bx else None,
        )
        ax.set(
            xlabel="Distance x from the current (km)",
            ylabel="Amplitude (nT)",
            title=title,
        )
    axes[0].legend(loc="center right")
    fig.suptitle(
        "A line current of 1 MA at 110 km height over the Quebec model, period 300 s (SCUBAS `Current`)",
        y=1.0,
        fontsize=10,
    )
    fig.tight_layout()
    save(fig, out, "line_current.png")


FIGURES = [
    fig_chain,
    fig_voltage_assembly,
    fig_attenuation_with_depth,
    fig_seafloor_ratios,
    fig_transfer_function_depth,
    fig_synthetic_test,
    fig_line_parameters_depth,
    fig_long_short_limits,
    fig_section_potentials,
    fig_network_models,
    fig_coast_potentials,
    fig_line_current,
]


def main(argv: list[str] | None = None) -> int:
    from loguru import logger

    logger.remove()
    args = sys.argv[1:] if argv is None else argv
    out = Path(args[0]) if args else DEFAULT_OUT
    apply_style()
    for make in FIGURES:
        make(out)
        print("wrote", make.__name__)
    from model_numbers import write_numbers

    write_numbers(out / "numbers.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

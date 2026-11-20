# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Generate the SCUBAS logo (SVG) in the family style of the author's other packages (BOLTS): dark navy field, two-tone name.

Aurora (the geomagnetic disturbance) fills the sky; a submarine cable runs through the sea between two shores and is grounded to the earth at
both ends, where the induced voltage appears. Colours are the Paul Tol "vibrant" set used by the SCUBAS figure style. Text is converted to
outlines with fontTools (DejaVu Sans Bold ships with matplotlib), so the SVG needs no installed font. A placeholder logo: replace the two SVGs
in ``docs/_static`` to change it. Run: ``python scripts/make_logo.py [OUTDIR]`` (default ``docs/_static``; ``make logo``). Needs ``fonttools`` and
``matplotlib``.
"""

from __future__ import annotations

import math
import pathlib
import random
import sys

import matplotlib
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

BG, SEA_HI, SEA_LO = "#0B1B2E", "#1672A8", "#0A3558"
BLUE, CYAN, ORANGE, WHITE, MIST, YELLOW = (
    "#0077BB",
    "#33BBEE",
    "#EE7733",
    "#FFFFFF",
    "#C9D6E3",
    "#FFD23F",
)
GREEN, VIOLET, LAND = "#44BB99", "#AA3377", "#1E3B2C"
LAND_Y = 176.0  # y of the ground surface in the 256 x 256 icon


def font_path() -> pathlib.Path:
    """The DejaVu Sans Bold file that ships with matplotlib."""
    return (
        pathlib.Path(matplotlib.get_data_path())
        / "fonts"
        / "ttf"
        / "DejaVuSans-Bold.ttf"
    )


def text_path(
    txt: str, size: float, x: float, y: float, tracking: float = 0.0
) -> tuple[str, float]:
    """SVG path data of ``txt`` as outlines with its baseline at (x, y); also returns the advance width."""
    face = TTFont(font_path())
    glyphs, cmap, hmtx = face.getGlyphSet(), face.getBestCmap(), face["hmtx"]
    scale = size / face["head"].unitsPerEm
    parts, cx = [], x
    for ch in txt:
        name = cmap[ord(ch)]
        pen = SVGPathPen(glyphs)
        glyphs[name].draw(TransformPen(pen, (scale, 0, 0, -scale, cx, y)))
        parts.append(pen.getCommands())
        cx += hmtx[name][0] * scale + tracking
    return " ".join(parts), cx - x - tracking


def aurora() -> str:
    """Aurora: curtains of thin vertical rays that are brightest along a wavy lower hem and fade upward, with a soft glow and a violet fringe."""
    layers = (
        (GREEN, 118, 74, 0.0, 0.85),
        (CYAN, 96, 56, 2.1, 0.6),
        (GREEN, 138, 50, 4.0, 0.55),
    )
    defs = '<filter id="blur" x="-10%" y="-30%" width="120%" height="160%"><feGaussianBlur stdDeviation="5"/></filter>'
    body = ""
    for n, (col, base, height, ph, op) in enumerate(layers):

        def hem(x: float, base: float = base, ph: float = ph) -> float:
            return base + 9 * math.sin(x / 31 + ph) + 4 * math.sin(x / 11 + 2 * ph)

        defs += (
            f'<linearGradient id="ar{n}" gradientUnits="userSpaceOnUse" x1="0" y1="{base - height - 14}" x2="0" y2="{base + 14}">'
            f'<stop offset="0" stop-color="{col}" stop-opacity="0"/><stop offset="0.75" stop-color="{col}" stop-opacity="0.75"/>'
            f'<stop offset="1" stop-color="{WHITE}" stop-opacity="0.95"/></linearGradient>'
        )
        xs = [i * 4.0 for i in range(65)]
        hem_pts = " ".join(f"{x:.1f},{hem(x):.1f}" for x in xs)
        glow = f"{hem_pts} 256,{base - height * 0.6:.1f} 0,{base - height * 0.6:.1f}"
        body += f'<polygon points="{glow}" fill="url(#ar{n})" opacity="{op * 0.45:.2f}" filter="url(#blur)"/>'
        fringe = " ".join(f"{x:.1f},{hem(x) + 6:.1f}" for x in xs)
        body += f'<polyline points="{fringe}" fill="none" stroke="{VIOLET}" stroke-width="5" opacity="{op * 0.5:.2f}" filter="url(#blur)"/>'
        rnd = random.Random(11 + n)
        for i in range(171):
            x = i * 1.5
            wave = 0.5 + 0.5 * math.sin(x / 13 + ph * 3) * math.sin(x / 37 + ph)
            h = height * (0.25 + 0.75 * wave) * rnd.uniform(0.7, 1.0)
            y = hem(x)
            body += f'<rect x="{x:.1f}" y="{y - h:.1f}" width="1.3" height="{h + 2:.1f}" fill="url(#ar{n})" opacity="{op * rnd.uniform(0.35, 0.7):.2f}"/>'
        body += f'<polyline points="{hem_pts}" fill="none" stroke="{col}" stroke-width="1.2" stroke-opacity="0.55"/>'
    return f"<defs>{defs}</defs>{body}"


def earth(x: float, y: float) -> str:
    """The earth symbol (a stem and three shortening bars) at (x, y), the PFE's earth electrode."""
    bars = "".join(
        f'<line x1="{x - w}" y1="{y + 7 + 5 * i}" x2="{x + w}" y2="{y + 7 + 5 * i}" stroke="{WHITE}" stroke-width="2"/>'
        for i, w in enumerate((9, 6, 3))
    )
    return f'<line x1="{x}" y1="{y}" x2="{x}" y2="{y + 7}" stroke="{WHITE}" stroke-width="2"/>{bars}'


def land() -> str:
    """The ground at both ends: two shores that dip into the sea."""
    left = f"M0,{LAND_Y} L80,{LAND_Y} C90,{LAND_Y} 94,{LAND_Y + 14} 100,{LAND_Y + 34} L100,256 L0,256 Z"
    right = f"M256,{LAND_Y} L176,{LAND_Y} C166,{LAND_Y} 162,{LAND_Y + 14} 156,{LAND_Y + 34} L156,256 L256,256 Z"
    return "".join(
        f'<path d="{p}" fill="{LAND}" stroke="{CYAN}" stroke-width="1.6" stroke-opacity="0.7"/>'
        for p in (left, right)
    )


def cable() -> str:
    """The cable from the left shore, down through the sea to the right shore, grounded at both ends."""
    y0, y1 = LAND_Y - 2, LAND_Y + 52
    path = f"M30,{y0} C62,{y0} 70,{y1} 92,{y1 + 4} L164,{y1 + 4} C186,{y1} 194,{y0} 226,{y0}"
    ends = "".join(
        f'<circle cx="{x}" cy="{y0}" r="3.4" fill="{ORANGE}" stroke="{WHITE}" stroke-width="1.2"/>'
        for x in (30, 226)
    )
    return f'<path d="{path}" fill="none" stroke="{WHITE}" stroke-width="2.6" stroke-linecap="round"/>{ends}'


def icon_body(uid: str = "i") -> str:
    """The 256 x 256 icon without the outer ``svg`` element."""
    rnd = random.Random(7)
    stars = "".join(
        f'<circle cx="{rnd.uniform(8, 248):.1f}" cy="{rnd.uniform(8, 90):.1f}" r="{rnd.choice([0.7, 1.0, 1.4]):.1f}" fill="{WHITE}" opacity="{rnd.uniform(0.35, 0.9):.2f}"/>'
        for _ in range(28)
    )
    ripples = "".join(
        f'<path d="M96,{LAND_Y + 36 + dy} q16,-4 32,0 t32,0" fill="none" stroke="{CYAN}" stroke-width="1" opacity="{op}"/>'
        for dy, op in ((0, 0.45), (12, 0.28), (24, 0.16))
    )
    return (
        f'<defs><linearGradient id="{uid}s" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="{SEA_HI}"/><stop offset="1" stop-color="{SEA_LO}"/></linearGradient>'
        f'<clipPath id="{uid}c"><rect width="256" height="256" rx="30"/></clipPath></defs>'
        f'<g clip-path="url(#{uid}c)"><rect width="256" height="256" fill="{BG}"/>{stars}{aurora()}'
        f'<rect y="{LAND_Y + 34}" width="256" height="{256 - LAND_Y - 34}" fill="url(#{uid}s)"/>{ripples}{land()}'
        f"{earth(30, LAND_Y)}{earth(226, LAND_Y)}{cable()}</g>"
    )


def svg(width: float, height: float, body: str, title: str) -> str:
    """Wrap ``body`` in an ``svg`` element."""
    return f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img"><title>{title}</title>{body}</svg>\n'


def wordmark() -> str:
    """The banner: icon, the two-tone name and the subtitle on the dark field."""
    scu, w1 = text_path("SCU", 112, 296, 168)
    bas, w2 = text_path("BAS", 112, 296 + w1 + 4, 168)
    sub, ws = text_path("Submarine Cables Upset by", 24, 300, 206)
    sub2, ws2 = text_path("Auroral Streams", 24, 300, 236)
    width = max(296 + w1 + w2 + 4, 300 + ws, 300 + ws2) + 40
    body = (
        f'<rect width="{width:.0f}" height="256" rx="30" fill="{BG}"/>{icon_body("w")}'
        f'<rect x="1" y="1" width="254" height="254" rx="29" fill="none" stroke="{CYAN}" stroke-opacity="0.45" stroke-width="2"/>'
        f'<path d="{scu}" fill="{ORANGE}"/><path d="{bas}" fill="{CYAN}"/>'
        f'<path d="{sub}" fill="{MIST}"/><path d="{sub2}" fill="{MIST}"/>'
    )
    return svg(round(width), 256, body, "SCUBAS")


def main() -> None:
    """Write ``logo.svg`` (the icon) and ``logo-wordmark.svg`` (the banner)."""
    out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "docs/_static")
    out.mkdir(parents=True, exist_ok=True)
    (out / "logo.svg").write_text(svg(256, 256, icon_body("i"), "SCUBAS"))
    (out / "logo-wordmark.svg").write_text(wordmark())
    print("wrote", *sorted(p.name for p in out.glob("logo*.svg")))


if __name__ == "__main__":
    main()

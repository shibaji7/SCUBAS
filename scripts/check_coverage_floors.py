#!/usr/bin/env python
# SPDX-License-Identifier: MIT
# Copyright (c) 2023 Shibaji Chakraborty
"""Fail if any module's line coverage is below its floor (``make cov-floors``, after a coverage run that wrote ``coverage.xml``).

The overall gate (``--cov-fail-under``) lets a well-tested module hide a badly tested one; the floors in ``[tool.scubas.coverage-floors]``
of ``pyproject.toml`` do not. ``default`` applies to every module, any other key (a path such as ``scubas/models.py``) overrides it.
Floors only go up: raise them when a module's coverage improves.
"""

from __future__ import annotations

import sys
import tomllib
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def module_rates(xml_path: Path) -> dict[str, float]:
    """Line coverage per source file, as a fraction, keyed by the path relative to the repository (``scubas/cables.py``).

    coverage.py writes file names relative to the ``<source>`` directory (``cables.py``, ``sources/current_sources.py``), so the
    package name is put back in front.
    """
    rates = {}
    for cls in ET.parse(xml_path).getroot().iter("class"):
        name = cls.get("filename", "")
        if name.endswith(".py"):
            parts = Path(name).parts
            key = (
                "/".join(parts[parts.index("scubas") :])
                if "scubas" in parts
                else f"scubas/{Path(name).as_posix()}"
            )
            rates[key] = float(cls.get("line-rate", 0.0))
    return rates


def failures(rates: dict[str, float], floors: dict[str, float]) -> list[str]:
    """Modules below their floor, as readable lines."""
    default = float(floors.get("default", 85))
    return [
        f"{name}: {100 * rate:.1f}% < {floors.get(name, default):g}%"
        for name, rate in sorted(rates.items())
        if 100 * rate < float(floors.get(name, default))
    ]


def main() -> int:
    """Check ``coverage.xml`` against the floors; returns the exit status."""
    xml_path = ROOT / "coverage.xml"
    if not xml_path.exists():
        print(
            "coverage.xml not found: run `make cov-xml` (pytest --cov=scubas --cov-report=xml) first",
            file=sys.stderr,
        )
        return 2
    floors = tomllib.loads((ROOT / "pyproject.toml").read_text())["tool"]["scubas"][
        "coverage-floors"
    ]
    rates = module_rates(xml_path)
    bad = failures(rates, floors)
    for name, rate in sorted(rates.items()):
        print(
            f"  {name:36s} {100 * rate:5.1f}%  (floor {floors.get(name, floors.get('default', 85)):g}%)"
        )
    if bad:
        print("\nbelow the floor:\n  " + "\n  ".join(bad), file=sys.stderr)
        return 1
    print("every module meets its coverage floor")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

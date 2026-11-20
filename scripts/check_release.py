#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Check that a release is consistent before it is tagged (``make release-check``).

The version in the code (``scubas/_version.py`` once it exists, otherwise ``[project] version`` of ``pyproject.toml``), the
newest released section of ``CHANGELOG.md`` and ``CITATION.cff`` (when the file exists) agree, and ``CHANGELOG.md`` keeps an
``[Unreleased]`` section. With ``--tag vX.Y.Z`` the tag is compared too. This script only reads files: it never tags or publishes.
"""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEMVER = r"\d+\.\d+\.\d+"


def code_version(root: Path = ROOT) -> str:
    """The version the package is built with."""
    version_file = root / "scubas" / "_version.py"
    if version_file.exists():
        m = re.search(rf'__version__\s*=\s*"({SEMVER})"', version_file.read_text())
        return m.group(1) if m else ""
    meta = tomllib.loads((root / "pyproject.toml").read_text()).get("project", {})
    version = str(meta.get("version", ""))
    return version if re.fullmatch(SEMVER, version) else ""


def citation_version(text: str) -> str:
    """The ``version:`` of CITATION.cff."""
    m = re.search(rf"^version:\s*['\"]?({SEMVER})['\"]?\s*$", text, re.M)
    return m.group(1) if m else ""


def changelog_versions(text: str) -> list[str]:
    """The released versions in the changelog, newest first (``## [0.1.9] - date``; ``[Unreleased]`` is not one)."""
    return re.findall(rf"^## \[({SEMVER})\]", text, re.M)


def problems(
    version: str, citation: str | None, changelog: str, tag: str | None = None
) -> list[str]:
    """Everything that does not agree, as readable lines (``citation`` is None while CITATION.cff does not exist yet)."""
    out = []
    released = changelog_versions(changelog)
    if not version:
        out.append(
            "no version found (scubas/_version.py or [project] version in pyproject.toml)"
        )
    if citation is not None and citation != version:
        out.append(
            f"CITATION.cff says {citation or 'nothing'}, the code says {version}"
        )
    if "## [Unreleased]" not in changelog:
        out.append("CHANGELOG.md has no [Unreleased] section")
    if not released or released[0] != version:
        out.append(
            f"the newest released section of CHANGELOG.md is {released[0] if released else 'missing'}, the code says {version}"
        )
    if tag is not None and tag.removeprefix("v") != version:
        out.append(f"the tag {tag} is not the code's version {version}")
    return out


def main(argv: list[str]) -> int:
    """Run the checks on the repository files; returns the exit status."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", help="a tag such as v0.2.0 to compare with the version")
    a = ap.parse_args(argv)
    cff = ROOT / "CITATION.cff"
    bad = problems(
        code_version(),
        citation_version(cff.read_text()) if cff.exists() else None,
        (ROOT / "CHANGELOG.md").read_text(),
        a.tag,
    )
    for line in bad:
        print("release check:", line, file=sys.stderr)
    if not bad:
        what = (
            "version, CHANGELOG.md and CITATION.cff agree"
            if cff.exists()
            else "version and CHANGELOG.md agree (no CITATION.cff yet)"
        )
        print(f"release check: {what}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

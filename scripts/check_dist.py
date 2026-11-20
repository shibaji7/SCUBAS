#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Check what is inside the built distributions (``make check`` and the build workflow, after ``python -m build``).

The wheel is the package only: every module of ``scubas``, the license files in its metadata, and nothing from the test suite. The sdist must
be enough to read the docs, run the examples and run the tests offline (docs, examples, data, scripts, tests and the compatibility baselines)
and must not carry build output or caches.
"""

from __future__ import annotations

import pathlib
import sys
import tarfile
import zipfile

WHEEL_NEEDS = [
    "scubas/__init__.py",
    "scubas/_version.py",
    "scubas/cables.py",
    "scubas/conductivity.py",
    "scubas/datasets.py",
    "scubas/errors.py",
    "scubas/future/__init__.py",
    "scubas/models.py",
    "scubas/plotlib.py",
    "scubas/sources/current_sources.py",
    "scubas/utils.py",
    ".dist-info/licenses/LICENSE",
    ".dist-info/licenses/COPYING",
    ".dist-info/licenses/NOTICE.md",
]
SDIST_NEEDS = [
    "README.md", "LICENSE", "COPYING", "NOTICE.md", "CHANGELOG.md", "ROADMAP.md", "SECURITY.md", "CITATION.cff", "pyproject.toml", "ruff.toml",
    "docs/index.md", "docs/conf.py", "docs/requirements.txt", "docs/dev/conventions.md", "examples/network_model.py",
    "examples/datasets/1989/FRD.csv", "scripts/check_dist.py", "tests/compat/test_api_surface.py", "tests/compat/baseline/api_surface.json",
    "tests/compat/baseline/behaviour.npz", "tests/compat/data/tst20200101vmin.min", "scubas/_version.py",
]  # fmt: skip
SDIST_FORBIDDEN = ("docs/_build/", "__pycache__", ".pyc", "outputs/", ".git/")


def wheel_problems(names: list[str]) -> list[str]:
    """What is missing from, or does not belong in, the wheel's file list."""
    out = [
        f"wheel is missing {n}"
        for n in WHEEL_NEEDS
        if not any(x.endswith(n) for x in names)
    ]
    out += [
        f"wheel contains {n}"
        for n in names
        if n.startswith(("tests/", "docs/")) or "/tests/" in n
    ]
    return out


def sdist_problems(names: list[str]) -> list[str]:
    """What is missing from, or does not belong in, the sdist's file list (paths relative to its top directory)."""
    rel = [n.split("/", 1)[1] if "/" in n else n for n in names]
    out = [f"sdist is missing {n}" for n in SDIST_NEEDS if n not in rel]
    out += [
        f"sdist contains {n}" for n in rel if any(bad in n for bad in SDIST_FORBIDDEN)
    ]
    return out


def main(argv: list[str]) -> int:
    """Check the wheel and the sdist found in ``dist/`` (or the paths given); returns the exit status."""
    root = pathlib.Path(__file__).resolve().parents[1]
    paths = [pathlib.Path(a) for a in argv] or sorted((root / "dist").glob("*"))
    bad, seen = [], 0
    for p in paths:
        if p.suffix == ".whl":
            bad += wheel_problems(zipfile.ZipFile(p).namelist())
            seen += 1
        elif p.name.endswith(".tar.gz"):
            with tarfile.open(p) as tf:
                bad += sdist_problems(tf.getnames())
            seen += 1
    if not seen:
        print(
            "check_dist: no wheel or sdist found: run `make build` first",
            file=sys.stderr,
        )
        return 2
    for line in bad:
        print("check_dist:", line, file=sys.stderr)
    if not bad:
        print(f"check_dist: {seen} distribution(s) have what they should")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

#!/usr/bin/env python
# SPDX-License-Identifier: MIT
# Copyright (c) 2023 Shibaji Chakraborty
"""Every Python file under scubas/, tests/, examples/, scripts/ (and the root ``setup.py`` / ``configure.py``) starts with the SPDX header.

The header is the license line and the copyright line, in that order, after an optional shebang. SCUBAS 0.1.9 has none, so the
files that lack it today are listed in ``scripts/spdx_legacy.txt``: a file outside the list must have the header, a file in the list
that now has it is stale (``--update`` removes it). The list only shrinks; adding the headers to the legacy files is a comment-only phase.
The license is NOT changed here: ``LICENSE_ID`` is the one place that names it.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LICENSE_ID = "MIT"
SPDX = f"# SPDX-License-Identifier: {LICENSE_ID}"
COPYRIGHT = "# Copyright (c) 2023 Shibaji Chakraborty"
LEGACY = ROOT / "scripts" / "spdx_legacy.txt"
DIRS = ("scubas", "tests", "examples", "scripts")
ROOT_FILES = ("setup.py", "configure.py")


def python_files() -> list[Path]:
    """The files the rule covers."""
    files = [p for d in DIRS for p in (ROOT / d).rglob("*.py")]
    files += [ROOT / n for n in ROOT_FILES if (ROOT / n).exists()]
    return sorted(
        p for p in files if "egg-info" not in str(p) and "__pycache__" not in str(p)
    )


def has_header(path: Path) -> bool:
    """True when the first two lines after an optional shebang are the license and copyright lines."""
    lines = path.read_text().splitlines()
    head = [x for x in lines[:3] if not x.startswith("#!")][:2]
    return head == [SPDX, COPYRIGHT]


def read_legacy() -> set[str]:
    """The files allowed to lack the header for now."""
    if not LEGACY.exists():
        return set()
    return {
        line.strip()
        for line in LEGACY.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }


def main(argv: list[str]) -> int:
    """Check the headers (or rewrite the legacy list with ``--update`` / ``--init``); returns the exit status."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--update",
        action="store_true",
        help="drop files that now have the header from the legacy list",
    )
    ap.add_argument(
        "--init",
        action="store_true",
        help="write the first legacy list (only when it does not exist)",
    )
    a = ap.parse_args(argv)
    missing = {
        p.relative_to(ROOT).as_posix() for p in python_files() if not has_header(p)
    }
    listed = read_legacy()
    if a.init:
        if LEGACY.exists():
            print("check_spdx: the legacy list exists; use --update", file=sys.stderr)
            return 2
        write(missing)
        print(f"check_spdx: {len(missing)} legacy files listed")
        return 0
    new = sorted(missing - listed)
    stale = sorted(listed - missing)
    if a.update:
        if new:
            print(
                "check_spdx: refusing to list new files without the header:\n  "
                + "\n  ".join(new),
                file=sys.stderr,
            )
            return 1
        write(missing)
        print(
            f"check_spdx: {len(stale)} stale entries removed, {len(missing)} files still lack the header"
        )
        return 0
    for title, rows in (
        ("missing the SPDX/copyright header", new),
        ("listed as legacy but now have the header (run --update)", stale),
    ):
        if rows:
            print(f"check_spdx: {title}:\n  " + "\n  ".join(rows), file=sys.stderr)
    if not (new or stale):
        print(
            f"check_spdx: ok ({len(python_files()) - len(missing)} files with the header, {len(missing)} legacy files listed)"
        )
    return 1 if (new or stale) else 0


def write(files: set[str]) -> None:
    """Write the legacy list."""
    head = "# Python files that predate the SPDX header (SCUBAS 0.1.9). Only shrinks: scripts/check_spdx.py --update.\n"
    LEGACY.write_text(head + "".join(f"{f}\n" for f in sorted(files)))


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

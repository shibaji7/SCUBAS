#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Refresh the table of the legacy constants inside ``docs/dev/parameters.md`` from ``scubas/defaults/legacy_constants.toml``.

Run it after changing the registry (``tests/test_provenance.py`` fails while the page and the file disagree).
"""

from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))  # this checkout, not a copy installed in site-packages

from scubas.provenance import legacy_constants  # noqa: E402

DOC = ROOT / "docs" / "dev" / "parameters.md"
START, END = "<!-- PARAMS:START -->", "<!-- PARAMS:END -->"


def main() -> None:
    """Rewrite the block between the markers."""
    ps = legacy_constants()
    text = DOC.read_text()
    block = f"{START}\n{ps.table_md()}\n{END}"
    DOC.write_text(
        re.sub(
            re.escape(START) + r".*?" + re.escape(END),
            lambda _: block,
            text,
            flags=re.S,
        )
    )
    print(f"updated {DOC.name}: {ps.counts()}")


if __name__ == "__main__":
    main()

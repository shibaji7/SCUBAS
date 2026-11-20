#!/usr/bin/env python
# SPDX-License-Identifier: MIT
# Copyright (c) 2023 Shibaji Chakraborty
"""Pre-commit hook: a commit that changes the package (``scubas/``) must also change CHANGELOG.md, and CHANGELOG.md must keep an "Unreleased" section.

Bypass for a commit with no user-visible change: ``SKIP=changelog git commit ...``.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def staged() -> list[str]:
    """Paths staged for the next commit."""
    out = subprocess.run(
        ["git", "diff", "--cached", "--name-only"],
        capture_output=True,
        text=True,
        check=True,
        cwd=ROOT,
    ).stdout
    return [x for x in out.splitlines() if x]


def main() -> int:
    """Check the changelog against the staged files; returns the exit status."""
    if "## [Unreleased]" not in (ROOT / "CHANGELOG.md").read_text():
        print("CHANGELOG.md must keep a '## [Unreleased]' section")
        return 1
    files = staged()
    if any(f.startswith("scubas/") for f in files) and "CHANGELOG.md" not in files:
        print(
            "scubas/ changed but CHANGELOG.md did not: add a line under '## [Unreleased]' (or SKIP=changelog)"
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python
# SPDX-License-Identifier: MIT
# Copyright (c) 2023 Shibaji Chakraborty
"""Keep the legacy lint debt in ``ruff.toml`` honest: it may only shrink.

``ruff.toml`` holds the full BOLTS rule set. Violations that already existed in SCUBAS 0.1.9 are listed per file between the markers
``BEGIN LEGACY DEBT`` and ``END LEGACY DEBT``. This script recomputes the violations with that block switched off and compares:

* ``--check`` (``make lint-ratchet``): fails when a violation is not in the block (new debt) or when an entry is stale (debt that was
  paid but is still listed, so it could hide a regression).
* ``--update``: rewrites the block with exactly what is still violated; refuses to add a file or a rule that is not listed.
* ``--init``: writes the first baseline (only when the block is empty).

The ruff executable is ``python -m ruff`` unless the environment variable ``RUFF`` names one.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "ruff.toml"
BEGIN, END = "# BEGIN LEGACY DEBT", "# END LEGACY DEBT"
ENTRY = re.compile(r'^"(?P<path>[^"]+)"\s*=\s*\[(?P<codes>[^\]]*)\]\s*$')


def ruff_command() -> list[str]:
    """The ruff invocation (``RUFF`` overrides)."""
    return (
        shlex.split(os.environ["RUFF"])
        if os.environ.get("RUFF")
        else [sys.executable, "-m", "ruff"]
    )


def read_block(text: str) -> dict[str, set[str]]:
    """The debt currently listed in ``ruff.toml``."""
    inside, debt = False, {}
    for line in text.splitlines():
        if line.strip() == BEGIN:
            inside = True
        elif line.strip() == END:
            break
        elif inside and (m := ENTRY.match(line.strip())):
            debt[m["path"]] = set(re.findall(r'"([A-Z]+\d+)"', m["codes"]))
    return debt


def actual_violations() -> dict[str, set[str]]:
    """Every violation of the full rule set, with the debt block removed from a temporary copy of ``ruff.toml``."""
    text = CONFIG.read_text()
    bare = text[: text.index(BEGIN)] + text[text.index(END) + len(END) :]
    tmp = (
        ROOT / ".ruff_ratchet.toml"
    )  # next to ruff.toml, so relative paths in it mean the same thing
    tmp.write_text(bare)
    try:
        cmd = [
            *ruff_command(),
            "check",
            "--no-cache",
            "--config",
            str(tmp),
            "--output-format",
            "json",
            "--exit-zero",
            "scubas",
            "tests",
            "scripts",
            "examples",
        ]
        out = subprocess.run(
            cmd, capture_output=True, text=True, cwd=ROOT, check=True
        ).stdout
    finally:
        tmp.unlink()
    found: dict[str, set[str]] = {}
    for item in json.loads(out or "[]"):
        rel = Path(item["filename"]).resolve().relative_to(ROOT).as_posix()
        found.setdefault(rel, set()).add(item["code"] or "SYNTAX")
    return found


def render(debt: dict[str, set[str]]) -> str:
    """The block, sorted, one line per file."""
    lines = [BEGIN, "[lint.extend-per-file-ignores]"]
    lines += [
        f'"{path}" = [{", ".join(chr(34) + c + chr(34) for c in sorted(codes))}]'
        for path, codes in sorted(debt.items())
    ]
    lines.append(END)
    return "\n".join(lines)


def write_block(text: str, debt: dict[str, set[str]]) -> str:
    """``text`` with the debt block replaced."""
    start, stop = text.index(BEGIN), text.index(END) + len(END)
    return text[:start] + render(debt) + text[stop:]


def compare(
    listed: dict[str, set[str]], actual: dict[str, set[str]]
) -> tuple[list[str], list[str]]:
    """(new violations, stale entries) as readable lines."""
    new = [
        f"{p}: {sorted(c - listed.get(p, set()))}"
        for p, c in sorted(actual.items())
        if c - listed.get(p, set())
    ]
    stale = [
        f"{p}: {sorted(c - actual.get(p, set()))}"
        for p, c in sorted(listed.items())
        if c - actual.get(p, set())
    ]
    return new, stale


def main(argv: list[str]) -> int:
    """Run the requested mode; returns the exit status."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--check", action="store_true", help="fail on new violations or stale entries"
    )
    mode.add_argument(
        "--update",
        action="store_true",
        help="shrink the block to what is still violated",
    )
    mode.add_argument(
        "--init",
        action="store_true",
        help="write the first baseline (block must be empty)",
    )
    a = ap.parse_args(argv)
    text = CONFIG.read_text()
    listed, actual = read_block(text), actual_violations()
    new, stale = compare(listed, actual)
    if a.init:
        if listed:
            print("ruff_ratchet: the block is not empty; use --update", file=sys.stderr)
            return 2
        CONFIG.write_text(write_block(text, actual))
        print(f"ruff_ratchet: baseline written for {len(actual)} files")
        return 0
    if a.update:
        if new:
            print(
                "ruff_ratchet: refusing to add new debt:\n  " + "\n  ".join(new),
                file=sys.stderr,
            )
            return 1
        CONFIG.write_text(
            write_block(text, {p: c for p, c in actual.items() if p in listed})
        )
        print(
            f"ruff_ratchet: {len(stale)} stale entries removed, {len(actual)} files still carry debt"
        )
        return 0
    for title, rows in (
        ("new violations (fix them, do not list them)", new),
        ("stale debt entries (run --update)", stale),
    ):
        if rows:
            print(f"ruff_ratchet: {title}:\n  " + "\n  ".join(rows), file=sys.stderr)
    if not (new or stale):
        total = sum(len(c) for c in actual.values())
        print(
            f"ruff_ratchet: legacy debt is exactly what is listed ({len(actual)} files, {total} file-rule pairs)"
        )
    return 1 if (new or stale) else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

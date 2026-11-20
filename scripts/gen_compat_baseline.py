#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Write the compatibility baselines in ``tests/compat/baseline/`` from the code that is installed right now (``make compat-baseline``).

The baselines describe what SCUBAS 0.1.9 does. Regenerate them only on the unchanged legacy code, or when a change of a
number or of a signature has been approved on purpose; a later refactor must reproduce them, never rewrite them.
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COMPAT = ROOT / "tests" / "compat"
BASELINE = COMPAT / "baseline"
sys.path.insert(0, str(COMPAT))

from compat_cases import build_cases  # noqa: E402
from compat_surface import collect_surface  # noqa: E402


def main() -> None:
    """Write ``api_surface.json``, ``behaviour.npz`` and ``behaviour_meta.json``."""
    BASELINE.mkdir(parents=True, exist_ok=True)
    surface = collect_surface()
    (BASELINE / "api_surface.json").write_text(
        json.dumps(surface, indent=1, sort_keys=True) + "\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        arrays, meta = build_cases(Path(tmp))
    np.savez_compressed(BASELINE / "behaviour.npz", **arrays)
    (BASELINE / "behaviour_meta.json").write_text(
        json.dumps(meta, indent=1, sort_keys=True, default=str) + "\n"
    )
    print(
        f"wrote {len(surface['modules'])} module surfaces, {len(arrays)} arrays, {len(meta)} metadata entries"
    )


if __name__ == "__main__":
    main()

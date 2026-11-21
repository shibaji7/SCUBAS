# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Numerical (FDTD or FEM) modelling of the coast effect.

Planned, not built: ``run`` raises ``NotImplementedYetError``. When the capability
exists this module becomes ``scubas.coastal.numerical``. See ``ROADMAP.md``.
"""

from __future__ import annotations

from typing import Any, Never

from ...errors import NotImplementedYetError


def run(*_args: Any, **_kwargs: Any) -> Never:
    """Raise ``NotImplementedYetError``: the capability is on the roadmap only."""
    raise NotImplementedYetError(
        "FDTD/FEM modelling of coastal effects", "ROADMAP: FDTD/FEM coastal effects"
    )


__all__ = ["run"]

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Auroral electrojet source model: a driving field where no magnetometer is.

Planned, not built: ``run`` raises ``NotImplementedYetError``. When the capability
exists this module becomes ``scubas.sources.electrojet``. See ``ROADMAP.md``.
"""

from __future__ import annotations

from typing import Any, Never

from ...errors import NotImplementedYetError


def run(*_args: Any, **_kwargs: Any) -> Never:
    """Raise ``NotImplementedYetError``: the capability is on the roadmap only."""
    raise NotImplementedYetError(
        "auroral electrojet source design", "ROADMAP: auroral electrojet design"
    )


__all__ = ["run"]

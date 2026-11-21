# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Uncertainty quantification: propagate input uncertainty to the cable voltage.

Planned, not built: ``run`` raises ``NotImplementedYetError``. When the capability
exists this module becomes
``scubas.uncertainty.propagation``. See ``ROADMAP.md``.
"""

from __future__ import annotations

from typing import Any, Never

from ...errors import NotImplementedYetError


def run(*_args: Any, **_kwargs: Any) -> Never:
    """Raise ``NotImplementedYetError``: the capability is on the roadmap only."""
    raise NotImplementedYetError(
        "uncertainty quantification of cable voltages",
        "ROADMAP: uncertainty quantification",
    )


__all__ = ["run"]

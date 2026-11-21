# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Automatic segmenter: split a cable route into sections.

Planned, not built: ``run`` raises ``NotImplementedYetError``. When the capability
exists this module becomes
``scubas.segmentation.auto_segmenter``. See ``ROADMAP.md``.
"""

from __future__ import annotations

from typing import Any, Never

from ...errors import NotImplementedYetError


def run(*_args: Any, **_kwargs: Any) -> Never:
    """Raise ``NotImplementedYetError``: the capability is on the roadmap only."""
    raise NotImplementedYetError(
        "automatic segmentation of a cable route", "ROADMAP: automatic segmenter"
    )


__all__ = ["run"]

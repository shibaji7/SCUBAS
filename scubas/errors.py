# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Exceptions of the packages added after 0.1.9.

The legacy modules keep raising the built-in exceptions they always raised.
"""

from __future__ import annotations


class ScubasError(Exception):
    """Base class of the errors raised by the newer parts of SCUBAS."""


class NotImplementedYetError(ScubasError, NotImplementedError):
    """A planned capability that exists only as a documented slot.

    The slot lives in ``scubas.future`` and the message names its roadmap item.
    """

    def __init__(self, what: str, roadmap: str = "see ROADMAP.md") -> None:
        super().__init__(f"{what} is not implemented yet ({roadmap}).")


__all__ = ["NotImplementedYetError", "ScubasError"]

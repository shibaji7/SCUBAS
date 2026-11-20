#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""SCUBAS: Submarine Cables Upset by Auroral Streams.

The voltage induced on submarine cables by geomagnetic disturbances, from a
transmission-line model of layered Ocean/Earth conductivity. The package imports
nothing heavy: use the modules (``scubas.cables``, ``scubas.models``,
``scubas.datasets``, ``scubas.conductivity``, ``scubas.plotlib``, ``scubas.utils``,
``scubas.sources``) directly.
"""

from ._version import __version__

__all__ = ["__version__"]

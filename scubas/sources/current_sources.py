"""
    current_sources.py: Module is used to implement ionosphereic 
                        current sources and associated ground magnetic 
                        purturbations.
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

from enum import Enum

import numpy as np
import scipy.constants as C
from loguru import logger

from scubas.datasets import PROFILES, Site


class Direction(Enum):
    GEO_NORTH_SOUTH = 1
    GEO_EAST_WEST = 2
    GEO_VERT_DOWN = 3


class Current(object):

    def __init__(
        self,
        h: float = 110.0,
        a: float = 0.0,
        I: float = 1.0,
        p: np.array = 0.0 * np.arange(10, 200),
        omega: np.array = 2 * np.pi * np.arange(10, 200),
        site: Site = PROFILES.Quebec,
        length_unit: str = "km",
        current_unit: str = "MA",
        current_direction: Direction = Direction.GEO_EAST_WEST,
        t: int = 0,
        layer: int = 0,
    ):
        """
        Current systems for one time instance / one location
            h: Height given in unit
            a: Width of the current, Cauchy distributed
            I: Currents in Amps/M Amps
            p: List of skin depths of the wave Z/(j*omega)
            omega: List of angular frequencies in rad/s
            site: Site profile default Quebec
            length_unit: Units for h, p, and a
            current_unit: Units of I
            current_direction: Direction of the currents, in geographic coordinates
            t: Time instance of the current
            layer: Layer where calculate Z and p
        """
        self.length_unit = length_unit
        self.current_unit = current_unit
        self.length_unit_multiplier = 1e3 if length_unit == "km" else 1.0
        self.current_unit_multiplier = 1e6 if current_unit == "MA" else 1.0
        self.h = h * self.length_unit_multiplier  # in meters
        self.p = p * self.length_unit_multiplier  # in meters
        self.a = a * self.length_unit_multiplier  # in meters
        self.I = I * self.current_unit_multiplier  # in Amps
        self.omega = omega
        self.freq = omega / (2 * np.pi)
        self.site = site
        self.current_direction = current_direction
        self.t = t
        # If site is given calculate the skin depth p at top layer
        if site:
            self.p = self.site.calcP(self.freq, layer=layer)
        return

    def compute_B_field(
        self,
        x: np.array = np.zeros((1)),
        length_unit: str = "km",
    ):
        """
        For a line current along y directed, we can only expect
        magnetic field along x and z directed.
        x: List of distance along x-direction given in unit
        length_unit: Units for x
        """
        length_unit_multiplier = 1e3 if length_unit == "km" else 1.0
        x *= length_unit_multiplier  # to meters
        B_field_directions = [
            Direction.GEO_EAST_WEST,
            Direction.GEO_NORTH_SOUTH,
            Direction.GEO_VERT_DOWN,
        ]
        # Remove 'y' directions only from the magnetic fields
        B_field_directions.remove(self.current_direction)
        logger.info(f"|Current| along : {self.current_direction} is {self.I} A")
        logger.info(
            f"Direction of the produced magnetic fields: {[d for d in B_field_directions]}"
        )
        self.B_fields = [
            dict(
                direction=B_field_directions[0],
                value=(C.mu_0 * self.I / (2 * np.pi))
                * (
                    (self.h / (self.h**2 + x**2))
                    + ((self.h + 2 * self.p) / ((self.h + 2 * self.p) ** 2 + x**2))
                ),
            ),
            dict(
                direction=B_field_directions[1],
                value=-1
                * (C.mu_0 * self.I / (2 * np.pi))
                * (
                    (x / (self.h**2 + x**2)) - (x / ((self.h + 2 * self.p) ** 2 + x**2))
                ),
            ),
        ]
        for b in self.B_fields:
            logger.warning(
                f"Magnetic field B along {b['direction']} is {b['value'].min()},{b['value'].max()}"
            )
        self.__compute_E_field__(x, "m")
        self.B_field_units = "T"
        return

    def __compute_E_field__(
        self,
        x: float = 0.0,
        length_unit: str = "km",
    ):
        """
        For a line current along y directed, we can only expect
        magnetic field along y directed.
        """
        length_unit_multiplier = 1e3 if length_unit == "km" else 1.0
        x *= length_unit_multiplier  # to meters
        E_field_directions = [self.current_direction]
        logger.info(f"Direction of the E_field: {E_field_directions}")
        self.E_fields = [
            dict(
                direction=E_field_directions[0],
                value=-1j
                * (self.omega * C.mu_0 * self.I / (2 * np.pi))
                * np.log(
                    np.sqrt((self.h + 2 * self.p) ** 2 + x**2)
                    / np.sqrt(self.h**2 + x**2)
                ),
            ),
        ]
        for e in self.E_fields:
            logger.warning(
                f"Magnetic field E along {e['direction']} is {e['value'].min()},{e['value'].max()}"
            )
        self.E_field_units = "V/m"
        return

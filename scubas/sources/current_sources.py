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

from loguru import logger


class Direction(Enum):
    GEO_NORTH_SOUTH = 1
    GEO_EAST_WEST = 2
    GEO_VERT_DOWN = 3


class Current(object):

    def __init__(
        self,
        h: float = 110.0,
        p: float = 0.0,
        a: float = 0.0,
        I: float = 1.0,
        length_unit: str = "km",
        current_unit: str = "MA",
        current_direction: Direction = Direction.GEO_EAST_WEST,
    ):
        """
        Current systems for one time instance / one location
            h: Height given in unit
            p: Skin depth of the wave
            a: Width of the current, Cauchy distributed
            I: Currents in Amps/M Amps
            length_unit: Units for h, p, and a
            current_unit: Units of I
            current_direction: Direction of the currents, in geographic coordinates
        """
        self.length_unit = length_unit
        self.current_unit = current_unit
        self.length_unit_multiplier = 1e3 if length_unit == "km" else 1.0
        self.current_unit_multiplier = 1e6 if current_unit == "MA" else 1.0
        self.h = h * self.length_unit_multiplier  # in meters
        self.p = p * self.length_unit_multiplier  # in meters
        self.a = a * self.length_unit_multiplier  # in meters
        self.I = I * self.current_unit_multiplier  # in Amps
        self.current_direction = current_direction
        return

    def compute_B_field(self):
        """
        For a line current along y directed, we can only expect
        magnetic field along x and z directed.
        """
        self.mag_field_directions = [
            Direction.GEO_EAST_WEST,
            Direction.GEO_NORTH_SOUTH,
            Direction.GEO_VERT_DOWN,
        ]
        # Remove 'y' directions only from the magnetic fields
        self.mag_field_directions.pop(self.current_direction)
        logger.info(f"Direction of the current: {self.current_direction}")
        return

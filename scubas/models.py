"""
    model.py: Module is used to implement Ocean Model with Layered-Earth structures.
    Class:
    ------
    OceanModel: Dedicated for ocean - Earth B and E field transfer functions
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import numpy as np
import pandas as pd
from loguru import logger
from scipy import constants as C


class OceanModel(object):
    """
    This class is to emulate electric and magnetic fields
    in a 1D ocean. We provide the 1D Earth conductivity model and
    one layer of ocean parameters (depth and resistivity) to emulate
    1D ocean.
    """

    def __init__(
        self,
        site,
        flim=[1e-6, 1e0],
    ):
        """
        Initialize the ocean model
        """
        logger.info(f"Compile Ocean-model: {site.name} E- and B-Fields")
        self.site = site
        self.flim = flim
        self.freqs = np.linspace(flim[0], flim[1], int(flim[1] / flim[0]) + 1)
        self.functions = {
            "E2B": lambda Z, Zd, kd: Zd / (np.cosh(kd) + (Zd * np.sinh(kd) / Z)),
            "B2B": lambda Z, Zd, kd: 1.0 / (np.cosh(kd) + (Zd * np.sinh(kd) / Z)),
        }
        return

    def calcZ(self, freqs=None):
        """
        Compute the Km(f) for different layers
        """
        freqs = np.copy(self.freqs) if freqs is None else freqs
        omega = 2 * C.pi * freqs
        sigma_s = 1 / self.site.layers[0].resistivity
        k2 = 1.0j * omega * C.mu_0 * sigma_s
        k = np.sqrt(k2)
        Zo = (1.0j * omega * C.mu_0 / k) / (C.mu_0 / 1.0e-3)
        Kf = self.site.calcZ(freqs, layer=1)[1]
        self.Z = {"ocean": Zo, "floor": Kf}
        return

    def calcTF(self, kinds=["E2B", "B2B"], freqs=None):
        """
        Calculate the transfer functions.
        """
        freqs = np.copy(self.freqs) if freqs is None else freqs
        self.calcZ(freqs)
        omega = 2 * C.pi * freqs
        Zd = self.Z["floor"]
        Z = self.Z["ocean"]
        TFs = {}

        omega = 2 * C.pi * freqs
        sigma_s = 1 / self.site.layers[0].resistivity
        k2 = 1.0j * omega * C.mu_0 * sigma_s
        k = np.sqrt(k2)
        kd = k * self.site.layers[0].thickness

        for kind in kinds:
            TFs[kind] = self.functions[kind](Z, Zd, kd)
        return TFs

    def get_TFs(self, key="E2B", freqs=None):
        """
        Fetch the transfer function based on key
        and frequency.
        """
        TFs = self.calcTF(freqs=freqs)
        tf = pd.DataFrame()
        tf["freq"], tf[key] = np.copy(self.freqs) if freqs is None else freqs, TFs[key]
        return tf

#!/usr/bin/env python

"""
    datasets.py: This is a data holder module with various datatypes
    1. Layer: Contains individual layer's conductivity and thickness <1D>
    2. Site: Contains sites with individual multiple layer
    3. Holds multiple types of sites/layers
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import constants as C


class Layer(object):
    """ """

    def __init__(self, name, tickness, conductivity):
        self.name = name
        self.thickness = tickness
        self.conductivity = conductivity
        self.resistivity = 1.0 / conductivity
        return


class Site(object):
    def __init__(self, layers, desciption, name):
        self.layers = layers
        self.desciption = desciption
        self.name = name
        return

    def get_thicknesses(self, index=-1):
        th = [l.thickness for l in self.layers]
        th = th[index] if index >= 0 else th
        return th

    def get_conductivities(self, index=-1):
        cd = [l.conductivity for l in self.layers]
        cd = cd[index] if index >= 0 else cd
        return cd

    def get_resistivities(self, index=-1):
        rv = [l.resistivity for l in self.layers]
        rv = rv[index] if index >= 0 else rv
        return rv

    def get_names(self, index=-1):
        na = [l.name for l in self.layers]
        na = na[index] if index >= 0 else na
        return na

    def get(self, index=-1):
        o = pd.DataFrame()
        o["names"] = self.get_names(index)
        o["conductivities"] = self.get_conductivities(index)
        o["resistivities"] = self.get_resistivities(index)
        o["thicknesses"] = self.get_thicknesses(index)
        return o
    
    def set_thicknesses(self, t, index=0):
        self.layers[index].thickness = t
        return self

    def calcZ(self, freqs, layer=0, return_all_layers=False):
        freqs = np.asarray(freqs)
        resistivities = np.asarray(self.get_resistivities())
        thicknesses = np.asarray(self.get_thicknesses())
        n = len(resistivities)
        nfreq = len(freqs)
        omega = 2 * np.pi * freqs
        complex_factor = 1j * omega * C.mu_0

        k = np.sqrt(1j * omega[np.newaxis, :] * C.mu_0 / resistivities[:, np.newaxis])
        Z = np.zeros(shape=(n, nfreq), dtype=complex)
        with np.errstate(divide="ignore", invalid="ignore"):
            Z[-1, :] = complex_factor / k[-1, :]
            r = np.zeros(shape=(n, nfreq), dtype=complex)
            for i in range(n - 2, -1, -1):
                r[i, :] = (1 - k[i, :] * Z[i + 1, :] / complex_factor) / (
                    1 + k[i, :] * Z[i + 1, :] / complex_factor
                )
                Z[i, :] = (
                    complex_factor
                    * (1 - r[i, :] * np.exp(-2 * k[i, :] * thicknesses[i]))
                    / (k[i, :] * (1 + r[i, :] * np.exp(-2 * k[i, :] * thicknesses[i])))
                )
        if freqs[0] == 0.0:
            Z[:, 0] = 0.0
        layer = layer if layer else 0
        Z_output = np.zeros(shape=(4, nfreq), dtype=complex)
        #################################################################
        # This factor '(1.0e-3 / C.mu_0)' converts K(f) and Z(f),
        # multiplying this factor makes Z -> K, devide this factor
        # regain Z.
        #################################################################
        Z_output[1, :] = Z[layer, :] * (1.0e-3 / C.mu_0)
        Z_output[2, :] = -Z_output[1, :]
        if return_all_layers:
            return Z_output, Z
        else:
            return Z_output

    def calcP(self, freqs, layer=0, return_Z=False):
        """
        Calculate the skin depth of the model
        """
        freqs = np.asarray(freqs)
        omega = 2 * np.pi * freqs

        if return_Z:
            Z_output, Z = self.calcZ(freqs=freqs, layer=layer, return_all_layers=True)
        else:
            Z_output = self.calcZ(freqs=freqs, layer=layer, return_all_layers=False)
        #################################################################
        # This factor '(1.0e-3 / C.mu_0)' converts K(f) and Z(f),
        # multiplying this factor makes Z -> K, devide this factor
        # regain Z.
        #################################################################
        p = Z_output[1, :] / (1j * omega * C.mu_0) / (1.0e-3 / C.mu_0)
        if return_Z:
            return p, Z_output, Z
        else:
            return p

    @staticmethod
    def init(conductivities, thicknesses, names, desciption, site_name):
        layers = []
        for c, t, n in zip(conductivities, thicknesses, names):
            layers.append(Layer(n, t, c))
        return Site(layers, desciption, site_name)


PROFILES = SimpleNamespace(
    **dict(
        BM=Site.init(
            conductivities=[0.2, 0.0003333, 0.02, 0.1, 1.12201],
            thicknesses=[2000, 75000, 332000, 250000, np.inf],
            names=[
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Lower Mantle",
            ],
            desciption="This model is Ben's model",
            site_name="Ben's Model",
        ),
        OM=Site.init(
            conductivities=[3.3333333, 0.2, 0.0003333, 0.02, 0.1, 1.12201],
            thicknesses=[1000, 2000, 75000, 332000, 250000, np.inf],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Lower Mantle",
            ],
            desciption="This model is modified Ben's model. An 1 km thick ocean on top of Ben's model.",
            site_name="Ocean Model",
        ),
        DB=Site.init(
            conductivities=[0.00005, 0.005, 0.001, 0.01, 0.3333333],
            thicknesses=[15000, 10000, 125000, 200000, np.inf],
            names=[
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Lower Mantle",
            ],
            desciption="This model is suggested by David.",
            site_name="David's Model",
        ),
        Quebec=Site.init(
            conductivities=[0.00005, 0.005, 0.001, 0.01, 0.3333333],
            thicknesses=[15000, 10000, 125000, 200000, np.inf],
            names=[
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Lower Mantle",
            ],
            desciption="This model is suggested by David.",
            site_name="David's Model",
        ),
        UN=Site.init(
            conductivities=[0.2, 0.2, 0.2, 0.2, 1.12201],
            thicknesses=[2000, 75000, 332000, 250000, np.inf],
            names=[
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Lower Mantle",
            ],
            desciption="This model is suggested by David.",
            site_name="Uniform Model",
        ),
        CS=Site.init(
            conductivities=[3.3333333, 0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[100, 3000, 20000, 140000, 246900, 250000, 340000],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Used in IGS study by David",
            site_name="Continental Shelf",
        ),
        SO=Site.init(
            conductivities=[3.3333333, 0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[1000, 2000, 10000, 70000, 327000, 250000, 340000],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Used in IGS study by David",
            site_name="Shallow Ocean",
        ),
        DO=Site.init(
            conductivities=[3.3333333, 0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[4000, 2000, 10000, 70000, 327000, 250000, 340000],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Used in IGS study by David",
            site_name="Deep Ocean",
        ),
        LD=Site.init(
            conductivities=[0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[1000, 20000, 140000, 249000, 250000, 340000],
            names=[
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Used in IGS study by David",
            site_name="Land",
        ),
        CS_W=Site.init(
            conductivities=[3.3333333, 0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[100, 8000, 15000, 150000, 236900, 250000, 340000],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Space Weather Paper: Cable (TAT-8) Section 1 / Used in SW 1989 Storm study",
            site_name="Continental Shelf West",
        ),
        DO_1=Site.init(
            conductivities=[3.3333333, 0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[4000, 4000, 10000, 145000, 247000, 250000, 340000],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Space Weather Paper: Cable (TAT-8) Section 1 / Used in SW 1989 Storm study",
            site_name="Deep Ocean",
        ),
        DO_2=Site.init(
            conductivities=[3.3333333, 0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[5200, 2000, 10000, 140000, 252800, 250000, 340000],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Space Weather Paper: Cable (TAT-8) Section 1 / Used in SW 1989 Storm study",
            site_name="Deep Ocean",
        ),
        DO_3=Site.init(
            conductivities=[3.3333333, 0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[4000, 2000, 10000, 140000, 254000, 250000, 340000],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Space Weather Paper: Cable (TAT-8) Section 1 / Used in SW 1989 Storm study",
            site_name="Deep Ocean",
        ),
        DO_4=Site.init(
            conductivities=[3.3333333, 0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[4800, 1000, 10000, 70000, 324200, 250000, 340000],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Space Weather Paper: Cable (TAT-8) Section 1 / Used in SW 1989 Storm study",
            site_name="Deep Ocean",
        ),
        DO_5=Site.init(
            conductivities=[3.3333333, 0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[4000, 1000, 10000, 70000, 324200, 250000, 340000],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Space Weather Paper: Cable (TAT-8) Section 1 / Used in SW 1989 Storm study",
            site_name="Deep Ocean",
        ),
        MAR=Site.init(
            conductivities=[3.3333333, 0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[3000, 0, 10000, 25000, 372000, 250000, 340000],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Space Weather Paper: Cable (TAT-8) Section 1 / Used in SW 1989 Storm study",
            site_name="Mid-Atlantic Ridge",
        ),
        DO_6=Site.init(
            conductivities=[3.3333333, 0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[4500, 1500, 10000, 70000, 324000, 250000, 340000],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Space Weather Paper: Cable (TAT-8) Section 1 / Used in SW 1989 Storm study",
            site_name="Deep Ocean",
        ),
        CS_E=Site.init(
            conductivities=[3.3333333, 0.3333333, 0.00033333333, 0.001, 0.01, 0.1, 1],
            thicknesses=[100, 3000, 20000, 120000, 266900, 250000, 340000],
            names=[
                "Seawater",
                "Sediments",
                "Crust",
                "Lithosphere",
                "Upper Mantle",
                "Transition Zone",
                "Lower Mantle",
            ],
            desciption="Space Weather Paper: Cable (TAT-8) Section 1 / Used in SW 1989 Storm study",
            site_name="Continental Shelf East",
        ),
    )
)

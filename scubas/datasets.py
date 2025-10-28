"""
Lightweight data structures describing layered conductivity models used
throughout SCUBAS workflows.
"""

from __future__ import annotations

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

from dataclasses import dataclass
from types import SimpleNamespace
from typing import Iterable, List, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from scipy import constants as C


@dataclass
class Layer:
    """
    Single geophysical layer describing thickness and electrical properties.
    """

    name: str
    thickness: float
    conductivity: float

    @property
    def resistivity(self) -> float:
        """
        Resistivity (Ohm·m) derived as the reciprocal of conductivity.
        """
        if self.conductivity == 0.0:
            raise ValueError("Conductivity must be non-zero to compute resistivity.")
        return 1.0 / self.conductivity


@dataclass
class Site:
    """
    Collection of layers associated with a location or conductivity profile.
    """

    layers: Sequence[Layer]
    description: str
    name: str

    def __post_init__(self) -> None:
        self.layers = list(self.layers)

    def get_thicknesses(self, index: int = -1) -> Union[List[float], float]:
        """
        Thickness values (km) for all layers or the specified index.
        """
        values = [layer.thickness for layer in self.layers]
        if index >= 0:
            return values[index]
        return values

    def get_conductivities(self, index: int = -1) -> Union[List[float], float]:
        """
        Conductivity values (S/m) for all layers or the specified index.
        """
        values = [layer.conductivity for layer in self.layers]
        if index >= 0:
            return values[index]
        return values

    def get_resistivities(self, index: int = -1) -> Union[List[float], float]:
        """
        Resistivity values (Ohm·m) derived from the layer conductivities.
        """
        values = [layer.resistivity for layer in self.layers]
        if index >= 0:
            return values[index]
        return values

    def get_names(self, index: int = -1) -> Union[List[str], str]:
        """
        Layer names for all layers or the specified index.
        """
        values = [layer.name for layer in self.layers]
        if index >= 0:
            return values[index]
        return values

    def get(self, index: int = -1) -> pd.DataFrame:
        """
        Return a pandas dataframe with layer metadata.
        """
        return pd.DataFrame(
            {
                "names": np.atleast_1d(self.get_names(index)),
                "conductivities": np.atleast_1d(self.get_conductivities(index)),
                "resistivities": np.atleast_1d(self.get_resistivities(index)),
                "thicknesses": np.atleast_1d(self.get_thicknesses(index)),
            }
        )

    def set_thicknesses(self, value: float, index: int = 0) -> "Site":
        """
        Mutate the thickness of a layer at the specified index.
        """
        self.layers[index].thickness = value
        return self

    def calcZ(
        self,
        freqs: Iterable[float],
        layer: int = 0,
        return_all_layers: bool = False,
    ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """
        Compute magnetotelluric impedance for the layered model.
        """
        freqs = np.asarray(freqs, dtype=float)
        resistivities = np.asarray(self.get_resistivities(), dtype=float)
        thicknesses = np.asarray(self.get_thicknesses(), dtype=float)
        n_layers = len(resistivities)
        n_freqs = len(freqs)
        omega = 2 * np.pi * freqs
        complex_factor = 1j * omega * C.mu_0

        k = np.sqrt(1j * omega[np.newaxis, :] * C.mu_0 / resistivities[:, np.newaxis])
        Z = np.zeros(shape=(n_layers, n_freqs), dtype=complex)
        with np.errstate(divide="ignore", invalid="ignore"):
            Z[-1, :] = complex_factor / k[-1, :]
            r = np.zeros(shape=(n_layers, n_freqs), dtype=complex)
            for idx in range(n_layers - 2, -1, -1):
                ratio = k[idx, :] * Z[idx + 1, :] / complex_factor
                r[idx, :] = (1 - ratio) / (1 + ratio)
                exponent = np.exp(-2 * k[idx, :] * thicknesses[idx])
                numerator = complex_factor * (1 - r[idx, :] * exponent)
                denominator = k[idx, :] * (1 + r[idx, :] * exponent)
                Z[idx, :] = numerator / denominator
        if freqs[0] == 0.0:
            Z[:, 0] = 0.0
        Z_output = np.zeros(shape=(4, n_freqs), dtype=complex)
        Z_output[1, :] = Z[layer, :] * (1.0e-3 / C.mu_0)
        Z_output[2, :] = -Z_output[1, :]
        if return_all_layers:
            return Z_output, Z
        return Z_output

    def calcP(
        self,
        freqs: Iterable[float],
        layer: int = 0,
        return_Z: bool = False,
    ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        Calculate skin depth for the model at the specified frequencies.
        """
        freqs = np.asarray(freqs, dtype=float)
        omega = 2 * np.pi * freqs

        if return_Z:
            Z_output, Z_all = self.calcZ(freqs=freqs, layer=layer, return_all_layers=True)
        else:
            Z_output = self.calcZ(freqs=freqs, layer=layer, return_all_layers=False)

        skin_depth = Z_output[1, :] / (1j * omega * C.mu_0) / (1.0e-3 / C.mu_0)
        if return_Z:
            return skin_depth, Z_output, Z_all
        return skin_depth

    @staticmethod
    def init(
        conductivities: Sequence[float],
        thicknesses: Sequence[float],
        names: Sequence[str],
        description: str,
        site_name: str,
    ) -> "Site":
        """
        Initialise a site from parallel sequences of layer properties.
        """
        layers = [
            Layer(name, thickness, conductivity)
            for conductivity, thickness, name in zip(conductivities, thicknesses, names)
        ]
        return Site(layers=layers, description=description, name=site_name)


PROFILES = SimpleNamespace(
    **{
        key: Site.init(**definition)
        for key, definition in {
            "BM": dict(
                conductivities=[0.2, 0.0003333, 0.02, 0.1, 1.12201],
                thicknesses=[2000, 75000, 332000, 250000, np.inf],
                names=[
                    "Sediments",
                    "Crust",
                    "Lithosphere",
                    "Upper Mantle",
                    "Lower Mantle",
                ],
                description="Ben's benchmark conductivity model.",
                site_name="Ben's Model",
            ),
            "OM": dict(
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
                description="Ben's model with a 1 km oceanic layer.",
                site_name="Ocean Model",
            ),
            "DB": dict(
                conductivities=[0.00005, 0.005, 0.001, 0.01, 0.3333333],
                thicknesses=[15000, 10000, 125000, 200000, np.inf],
                names=[
                    "Sediments",
                    "Crust",
                    "Lithosphere",
                    "Upper Mantle",
                    "Lower Mantle",
                ],
                description="Conductivity profile used in David's study.",
                site_name="David's Model",
            ),
            "Quebec": dict(
                conductivities=[0.00005, 0.005, 0.001, 0.01, 0.3333333],
                thicknesses=[15000, 10000, 125000, 200000, np.inf],
                names=[
                    "Sediments",
                    "Crust",
                    "Lithosphere",
                    "Upper Mantle",
                    "Lower Mantle",
                ],
                description="Conductivity model representing Quebec.",
                site_name="David's Model",
            ),
            "UN": dict(
                conductivities=[0.2, 0.2, 0.2, 0.2, 1.12201],
                thicknesses=[2000, 75000, 332000, 250000, np.inf],
                names=[
                    "Sediments",
                    "Crust",
                    "Lithosphere",
                    "Upper Mantle",
                    "Lower Mantle",
                ],
                description="Uniform conductivity model inspired by David's work.",
                site_name="Uniform Model",
            ),
            "CS": dict(
                conductivities=[
                    3.3333333,
                    0.3333333,
                    0.00033333333,
                    0.001,
                    0.01,
                    0.1,
                    1,
                ],
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
                description="IGS continental shelf reference profile.",
                site_name="Continental Shelf",
            ),
            "SO": dict(
                conductivities=[
                    3.3333333,
                    0.3333333,
                    0.00033333333,
                    0.001,
                    0.01,
                    0.1,
                    1,
                ],
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
                description="Shallow ocean conductivity profile used in IGS analyses.",
                site_name="Shallow Ocean",
            ),
            "DO": dict(
                conductivities=[
                    3.3333333,
                    0.3333333,
                    0.00033333333,
                    0.001,
                    0.01,
                    0.1,
                    1,
                ],
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
                description="Deep ocean reference profile for IGS studies.",
                site_name="Deep Ocean",
            ),
            "LD": dict(
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
                description="Land-based conductivity profile from IGS study.",
                site_name="Land",
            ),
            "CS_W": dict(
                conductivities=[
                    3.3333333,
                    0.3333333,
                    0.00033333333,
                    0.001,
                    0.01,
                    0.1,
                    1,
                ],
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
                description="Space Weather 1989: Continental shelf (western) section.",
                site_name="Continental Shelf West",
            ),
            "DO_1": dict(
                conductivities=[
                    3.3333333,
                    0.3333333,
                    0.00033333333,
                    0.001,
                    0.01,
                    0.1,
                    1,
                ],
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
                description="Space Weather 1989 study: deep ocean variant 1.",
                site_name="Deep Ocean",
            ),
            "DO_2": dict(
                conductivities=[
                    3.3333333,
                    0.3333333,
                    0.00033333333,
                    0.001,
                    0.01,
                    0.1,
                    1,
                ],
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
                description="Space Weather 1989 study: deep ocean variant 2.",
                site_name="Deep Ocean",
            ),
            "DO_3": dict(
                conductivities=[
                    3.3333333,
                    0.3333333,
                    0.00033333333,
                    0.001,
                    0.01,
                    0.1,
                    1,
                ],
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
                description="Space Weather 1989 study: deep ocean variant 3.",
                site_name="Deep Ocean",
            ),
            "DO_4": dict(
                conductivities=[
                    3.3333333,
                    0.3333333,
                    0.00033333333,
                    0.001,
                    0.01,
                    0.1,
                    1,
                ],
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
                description="Space Weather 1989 study: deep ocean variant 4.",
                site_name="Deep Ocean",
            ),
            "DO_5": dict(
                conductivities=[
                    3.3333333,
                    0.3333333,
                    0.00033333333,
                    0.001,
                    0.01,
                    0.1,
                    1,
                ],
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
                description="Space Weather 1989 study: deep ocean variant 5.",
                site_name="Deep Ocean",
            ),
            "MAR": dict(
                conductivities=[
                    3.3333333,
                    0.3333333,
                    0.00033333333,
                    0.001,
                    0.01,
                    0.1,
                    1,
                ],
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
                description="Space Weather 1989 study: Mid-Atlantic Ridge profile.",
                site_name="Mid-Atlantic Ridge",
            ),
            "DO_6": dict(
                conductivities=[
                    3.3333333,
                    0.3333333,
                    0.00033333333,
                    0.001,
                    0.01,
                    0.1,
                    1,
                ],
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
                description="Space Weather 1989 study: deep ocean variant 6.",
                site_name="Deep Ocean",
            ),
            "CS_E": dict(
                conductivities=[
                    3.3333333,
                    0.3333333,
                    0.00033333333,
                    0.001,
                    0.01,
                    0.1,
                    1,
                ],
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
                description="Space Weather 1989: Continental shelf (eastern) section.",
                site_name="Continental Shelf East",
            ),
        }.items()
    }
)

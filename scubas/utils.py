"""
Shared utility helpers covering recursive namespaces, spectral transforms,
component mappings, and great-circle distance calculations.
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

import math
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any, Mapping, Optional, Sequence, Tuple

import numpy as np


class RecursiveNamespace(SimpleNamespace):
    """
    Nested namespace that recursively converts dictionaries and lists.
    """

    def __init__(self, **kwargs: Any) -> None:
        converted = {key: self._convert(value) for key, value in kwargs.items()}
        super().__init__(**converted)

    @classmethod
    def _convert(cls, value: Any) -> Any:
        if isinstance(value, Mapping):
            return cls(**value)
        if isinstance(value, list):
            return [cls._convert(item) for item in value]
        return value


def frexp102str(x: float, precision: int = 2) -> str:
    """
    Represent ``x`` using mantissa-exponent scientific notation.
    """
    if x == 0:
        return "0"
    sign = "-" if x < 0 else ""
    value = abs(x)
    exponent = int(math.floor(math.log10(value)))
    mantissa = value / (10**exponent)
    fmt = f"{{:{'.' + str(precision) if precision >= 0 else ''}f}}"
    return f"{sign}{fmt.format(mantissa)}$\\times 10^{{{exponent}}}$"


def fft(
    X: Sequence[float],
    dT: float,
    remove_zero_frequency: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Forward FFT with consistent scaling for real-valued time series.
    """
    data = np.asarray(X, dtype=float)
    if data.size == 0:
        raise ValueError("FFT requires at least one sample.")
    if dT <= 0:
        raise ValueError("Sampling interval 'dT' must be positive.")
    spectrum = 2.0 / data.size * np.fft.rfft(data)
    freqs = np.fft.rfftfreq(data.size, d=dT)
    if remove_zero_frequency and freqs.size > 1:
        freqs[0] = freqs[1]
    return spectrum, freqs


def ifft(Y: Sequence[complex]) -> np.ndarray:
    """
    Inverse FFT matching the scaling used in :func:`fft`.
    """
    spectrum = np.asarray(Y, dtype=complex)
    if spectrum.size == 0:
        raise ValueError("IFFT requires at least one sample.")
    return np.fft.irfft(spectrum) * 2 * spectrum.size


def component_mappings(field: str = "B2E", comp: str = "X") -> str:
    """
    Map field component names between B-field and E-field conventions.
    """
    mappings: Mapping[str, Mapping[str, str]] = {
        "B2E": {"X": "Y", "Y": "X"},
    }
    try:
        return mappings[field][comp]
    except KeyError as exc:
        raise KeyError(f"Unsupported component mapping for field '{field}' and component '{comp}'.") from exc


def component_sign_mappings(fromto: str = "BxEy") -> float:
    """
    Provide sign adjustments for specific component transformations.
    """
    mappings: Mapping[str, float] = {
        "BxEy": -1.0,
        "ByEx": 1.0,
    }
    try:
        return mappings[fromto]
    except KeyError as exc:
        raise KeyError(f"Unsupported component sign mapping '{fromto}'.") from exc


@dataclass(frozen=True)
class GreatCircle:
    """
    Compute great-circle and haversine distances (km) between locations.
    """

    initial: Any
    final: Any
    Re: float = 6371.0

    @staticmethod
    def _has_location_fields(loc: Any) -> bool:
        return hasattr(loc, "lat") and hasattr(loc, "lon")

    def _to_radians(self, loc: Any) -> Tuple[float, float]:
        if not self._has_location_fields(loc):
            raise ValueError("Location objects must expose 'lat' and 'lon' attributes.")
        return math.radians(float(loc.lat)), math.radians(float(loc.lon))

    def great_circle(self, initial: Optional[Any] = None, final: Optional[Any] = None) -> float:
        """
        Great-circle distance assuming a spherical Earth.
        """
        start = initial if initial is not None else self.initial
        end = final if final is not None else self.final
        lat1, lon1 = self._to_radians(start)
        lat2, lon2 = self._to_radians(end)
        cosine_argument = (
            math.sin(lat1) * math.sin(lat2)
            + math.cos(lat1) * math.cos(lat2) * math.cos(lon1 - lon2)
        )
        cosine_argument = max(-1.0, min(1.0, cosine_argument))
        return self.Re * math.acos(cosine_argument)

    def haversine(self, initial: Optional[Any] = None, final: Optional[Any] = None) -> float:
        """
        Haversine distance offering improved numerical stability.
        """
        start = initial if initial is not None else self.initial
        end = final if final is not None else self.final
        lat1, lon1 = self._to_radians(start)
        lat2, lon2 = self._to_radians(end)
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
        a = min(1.0, max(0.0, a))
        return 2 * self.Re * math.asin(math.sqrt(a))

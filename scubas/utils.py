"""utility.py: Module is used to implement utility methods.
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import math
from types import SimpleNamespace

import numpy as np


class RecursiveNamespace(SimpleNamespace):
    """ """

    @staticmethod
    def map_entry(entry):
        if isinstance(entry, dict):
            return RecursiveNamespace(**entry)
        return entry

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        for key, val in kwargs.items():
            if type(val) == dict:
                setattr(self, key, RecursiveNamespace(**val))
            elif type(val) == list:
                setattr(self, key, list(map(self.map_entry, val)))


def frexp102str(x: float) -> str:
    """Convert to a float to `str` in the form of mantesa-exponent.

    Arguments:
        x: `float` number to be converted.
        
    Returns:
        `str` value of the `float` data in mantesa-exponent form.
    """
    exp = int(math.log10(x))
    if exp <= 0:
        exp -= 1
    m, exp = x / 10**exp, exp
    txt = r"%.2f$\times 10^{%d}$" % (m, exp)
    return txt


def fft(X: np.array, dT: float, remove_zero_frequency:bool=True) -> tuple:
    """This method runs FFT using NumPy package of a real signal (X).
    
    Arguments:
        X: Timeseries data 
        dT: Timedelta in seconds
        remove_zero_frequency: If true replace the first 0 frequency (DC component)
        
    Returns:
        FFT values with frequency as a tuple.
    """
    n = len(X)
    Y = 2.0 / n * np.fft.rfft(X)
    f = np.fft.rfftfreq(len(X)) / dT
    if remove_zero_frequency:
        f[0] = f[1]
    return (Y, f)


def ifft(Y: np.array) -> np.array:
    """This method runs Inverse FFT using NumPy package of a complex signal (Y).
    
    Arguments:
        Y: Complex FFT dataset
        
    Returns:
        Timeseries values.
    """
    n = len(Y)
    X = np.fft.irfft(Y) * 2 * n
    return X


def component_mappings(field="B2E", comp="X"):
    """
    This method holds components mapping from
    (i) B2E
    """
    _map_ = {"B2E": {"X": "Y", "Y": "X"}}
    return _map_[field][comp]


def component_sign_mappings(fromto="BxEy"):
    """
    This method holds components mapping from
    (i) B2E sign change (Bx to Ey) and (By to Ex)
    """
    _map_ = {
        "BxEy": -1.0,
        "ByEx": 1.0,
    }
    return _map_[fromto]

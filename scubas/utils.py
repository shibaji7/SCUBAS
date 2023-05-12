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


def frexp102str(x):
    """
    Convert to mantesa exponent form in txt
    """
    exp = int(math.log10(x))
    if exp <= 0:
        exp -= 1
    m, exp = x / 10**exp, exp
    txt = r"%.2f$\times 10^{%d}$" % (m, exp)
    return txt


def fft(X, dT, remove_zero_frequency=True):
    """
    This function is responsible for FFT using
    numpy package of a real signal (X).
    """
    n = len(X)
    Y = 2.0 / n * np.fft.rfft(X)
    f = np.fft.rfftfreq(len(X)) / dT
    if remove_zero_frequency:
        f[0] = f[1]
    return (Y, f)


def ifft(Y):
    """
    This function is responsible for IFFT using
    numpy package of a complex FFT signal (Y).
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

import math
from types import SimpleNamespace

import numpy as np
import pytest

from scubas.utils import (
    GreatCircle,
    RecursiveNamespace,
    component_mappings,
    component_sign_mappings,
    fft,
    frexp102str,
    ifft,
)


def test_recursive_namespace_nested_conversion():
    data = {"a": 1, "b": {"c": 2, "d": [{"e": 3}, 4]}}
    namespace = RecursiveNamespace(**data)
    assert namespace.a == 1
    assert namespace.b.c == 2
    assert namespace.b.d[0].e == 3
    assert namespace.b.d[1] == 4


@pytest.mark.parametrize(
    "value,expected",
    [
        (0.0, "0"),
        (1.234e5, "1.23$\\times 10^{5}$"),
        (-9.876e-4, "-9.88$\\times 10^{-4}$"),
    ],
)
def test_frexp102str_handles_values(value, expected):
    assert frexp102str(value) == expected


def test_fft_ifft_roundtrip():
    t = np.linspace(0, 1, 128, endpoint=False)
    signal = np.sin(2 * math.pi * 5 * t)
    spectrum, freqs = fft(signal, dT=t[1] - t[0])
    assert freqs[0] == freqs[1]
    reconstructed = ifft(spectrum)
    np.testing.assert_allclose(reconstructed, signal, atol=1e-6)


def test_component_mappings_and_signs():
    assert component_mappings("B2E", "X") == "Y"
    assert component_mappings("B2E", "Y") == "X"
    assert component_sign_mappings("BxEy") == -1.0
    assert component_sign_mappings("ByEx") == 1.0
    with pytest.raises(KeyError):
        component_mappings("E2B", "X")
    with pytest.raises(KeyError):
        component_sign_mappings("invalid")


def test_great_circle_and_haversine():
    start = SimpleNamespace(lat=0.0, lon=0.0)
    end = SimpleNamespace(lat=0.0, lon=90.0)
    gc = GreatCircle(start, end)
    expected_distance = gc.Re * math.pi / 2
    assert math.isclose(gc.great_circle(), expected_distance, rel_tol=1e-6)
    assert math.isclose(gc.haversine(), expected_distance, rel_tol=1e-6)

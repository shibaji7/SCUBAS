import numpy as np
import pandas as pd
import pytest

from scubas.conductivity import ConductivityProfile
from scubas.datasets import Site


class DummyInterpolator:
    def __init__(self, value):
        self.value = value

    def __call__(self, pts):
        length = len(pts) if isinstance(pts, np.ndarray) else 1
        return np.full(length, self.value, dtype=float)


@pytest.fixture(autouse=True)
def patch_load_earth_model(monkeypatch):
    def _dummy_load(self):
        self.lithosphere_model = {
            "water_top_depth": DummyInterpolator(0.0),
            "water_bottom_depth": DummyInterpolator(1.0),
            "upper_sediments_top_depth": DummyInterpolator(1.0),
            "upper_sediments_bottom_depth": DummyInterpolator(2.0),
            "middle_sediments_top_depth": DummyInterpolator(2.0),
            "middle_sediments_bottom_depth": DummyInterpolator(3.0),
            "lower_sediments_top_depth": DummyInterpolator(3.0),
            "lower_sediments_bottom_depth": DummyInterpolator(4.0),
            "upper_crust_top_depth": DummyInterpolator(4.0),
            "upper_crust_bottom_depth": DummyInterpolator(5.0),
            "middle_crust_top_depth": DummyInterpolator(5.0),
            "middle_crust_bottom_depth": DummyInterpolator(6.0),
            "lower_crust_top_depth": DummyInterpolator(6.0),
            "lower_crust_bottom_depth": DummyInterpolator(8.0),
            "lithosphere_top_depth": DummyInterpolator(8.0),
            "lithosphere_bottom_depth": DummyInterpolator(12.0),
            "asthenospheric_mantle_top_depth": DummyInterpolator(12.0),
        }

    monkeypatch.setattr(
        ConductivityProfile, "load_earth_model", _dummy_load, raising=False
    )


def test_get_interpolation_points():
    profile = ConductivityProfile()
    pts = profile.get_interpolation_points([0, 0], [1, 1])
    assert pts.shape[1] == 2
    assert np.allclose(pts[0], [0, 0])
    assert np.allclose(pts[-1], [1, 1])


def test_compile_profile_to_site():
    profile = ConductivityProfile()
    site = profile.compile_profile([10.2, -70.6])
    assert isinstance(site, Site)
    assert len(site.layers) == 7
    assert pytest.approx(site.layers[0].conductivity, rel=1e-6) == pytest.approx(
        1.0 / profile.seawater_resistivity
    )


def test_compile_profile_dataframe():
    profile = ConductivityProfile()
    df = profile.compile_profile([10.2, -70.6], to_site=False)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["thickness", "resistivity", "name"]
    assert (df["thickness"] > 0).all()

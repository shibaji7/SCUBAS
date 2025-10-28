import numpy as np
import pandas as pd
import pytest

from scubas.datasets import Site
from scubas.models import OceanModel


@pytest.fixture
def simple_site():
    return Site.init(
        conductivities=[3.0, 0.5, 0.1],
        thicknesses=[1.0, 5.0, np.inf],
        names=["Seawater", "Sediment", "Mantle"],
        description="Test site",
        site_name="Test",
    )


def build_bfield_dataframe():
    times = pd.date_range("2020-01-01", periods=16, freq="S")
    data = {
        "X": np.sin(np.linspace(0, 2 * np.pi, len(times))),
        "Y": np.cos(np.linspace(0, 2 * np.pi, len(times))),
    }
    return pd.DataFrame(data, index=times)


def test_ocean_model_transfer_functions(simple_site):
    model = OceanModel(simple_site, flim=(1e-4, 1e-2))
    tf = model.get_TFs("E2B")
    assert "E2B" in tf.columns
    assert (tf["E2B"].values != 0).any()


def test_ocean_model_to_Efields(simple_site):
    model = OceanModel(simple_site, flim=(1e-4, 1e-2))
    bfield = build_bfield_dataframe()
    model.to_Efields(bfield)
    assert "X" in model.Efield.columns
    assert "Y" in model.Efield.columns
    assert model.Efield["X"].shape == bfield.index.shape
    assert model.Efield["Y"].shape == bfield.index.shape

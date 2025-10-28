import numpy as np

from scubas.datasets import Layer, Site


def build_simple_site() -> Site:
    return Site.init(
        conductivities=[1.0, 0.5],
        thicknesses=[10.0, np.inf],
        names=["Water", "Mantle"],
        description="Test profile",
        site_name="Test",
    )


def test_site_init_creates_layers():
    site = build_simple_site()
    assert len(site.layers) == 2
    assert isinstance(site.layers[0], Layer)
    assert site.get_names(0) == "Water"


def test_site_calcZ_and_calcP_shapes():
    site = build_simple_site()
    freqs = np.logspace(-4, -2, 5)
    Z = site.calcZ(freqs)
    assert Z.shape == (4, freqs.size)
    p = site.calcP(freqs)
    assert p.shape == (freqs.size,)


def test_site_calcP_with_return_Z():
    site = build_simple_site()
    freqs = np.array([1e-3, 1e-2])
    skin_depth, Z_output, all_layers = site.calcP(freqs, return_Z=True)
    assert skin_depth.shape == (2,)
    assert Z_output.shape == (4, 2)
    assert all_layers.shape[0] == len(site.layers)

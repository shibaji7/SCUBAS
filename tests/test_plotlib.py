import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt

import scubas.plotlib as plotlib

PlotArtifacts = getattr(plotlib, "PlotArtifacts", None)


class DummyTransferFunction:
    def __init__(self):
        self.freq = np.logspace(-4, -2, 5)
        self.E2B = np.exp(1j * np.linspace(0, np.pi / 2, 5))


def test_plot_transfer_function_returns_artifacts():
    artifacts = plotlib.plot_transfer_function(DummyTransferFunction())
    assert isinstance(artifacts, PlotArtifacts)
    assert isinstance(artifacts.figure, plt.Figure)
    plt.close(artifacts.figure)


def test_potential_plots_return_artifacts():
    x = np.linspace(0, 10, 20)
    V = np.sin(x)
    artifacts = plotlib.potential_along_section(
        V, x, sec=1, Vi=1.0, Vk=0.5, Z=2.0, Y=3.0, gma=0.1, Z0=5.0
    )
    assert isinstance(artifacts, PlotArtifacts)
    plt.close(artifacts.figure)

    artifacts = plotlib.cable_potential(V, x)
    assert isinstance(artifacts, PlotArtifacts)
    plt.close(artifacts.figure)

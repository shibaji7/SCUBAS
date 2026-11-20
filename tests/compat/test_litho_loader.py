# SPDX-License-Identifier: MIT
# Copyright (c) 2023 Shibaji Chakraborty
"""The LITHO1.0 loader without the network: a synthetic netCDF stands in for the real model and ``urlopen`` is mocked.

``@pytest.mark.litho`` tests need a real ``LITHO1.0.nc`` (set ``SCUBAS_LITHO1_PATH``) and are skipped otherwise.
"""

import io
import os
import urllib.error
from pathlib import Path

import numpy as np
import pytest
from compat_cases import chdir, write_synthetic_litho
from loguru import logger

from scubas.conductivity import ConductivityProfile


@pytest.fixture(autouse=True)
def _quiet():
    logger.disable("scubas")
    yield
    logger.enable("scubas")


def test_defaults_of_the_constructor_are_pinned(tmp_path):
    write_synthetic_litho(tmp_path / ".scubas_config" / "LITHO1.0.nc")
    with chdir(tmp_path):
        cp = ConductivityProfile()
    assert vars(cp.cprop) == {
        "earth_model": "LITHO1.0.nc",
        "grid_interpolation_method": "nearest",
        "seawater_resistivity": 0.3,
        "sediment_resistivity": 3.0,
        "crust_resistivity": 3000.0,
        "lithosphere_resistivity": 1000.0,
        "asthenosphere_resistivity": 100.0,
        "transition_zone_resistivity": 10.0,
        "lower_mantle_resistivity": 1.0,
        "transition_zone_top": 410.0,
        "transition_zone_bot": 660.0,
        "profile_max_depth": 1000.0,
        "uri": "http://ds.iris.edu/files/products/emc/emc-files/LITHO1.0.nc",
    }


def test_overrides_are_merged_over_the_defaults(tmp_path):
    write_synthetic_litho(tmp_path / ".scubas_config" / "LITHO1.0.nc")
    with chdir(tmp_path):
        cp = ConductivityProfile(
            {"seawater_resistivity": 0.25, "grid_interpolation_method": "linear"}
        )
    assert cp.seawater_resistivity == 0.25 and cp.grid_interpolation_method == "linear"
    assert cp.crust_resistivity == 3000.0


def test_the_model_is_downloaded_once_into_the_working_directory(tmp_path, monkeypatch):
    source = tmp_path / "remote.nc"
    write_synthetic_litho(source)
    calls = []

    def fake_urlopen(request):
        calls.append(request.full_url)
        return io.BytesIO(source.read_bytes())

    monkeypatch.setattr("scubas.conductivity.urllib.request.urlopen", fake_urlopen)
    work = tmp_path / "work"
    work.mkdir()
    with chdir(work):
        ConductivityProfile()
        assert (work / ".scubas_config" / "LITHO1.0.nc").exists()
        ConductivityProfile()  # cached: no second download
    assert calls == ["http://ds.iris.edu/files/products/emc/emc-files/LITHO1.0.nc"]


def test_a_failed_download_is_a_runtime_error(tmp_path, monkeypatch):
    def failing(request):
        raise urllib.error.URLError("offline")

    monkeypatch.setattr("scubas.conductivity.urllib.request.urlopen", failing)
    with chdir(tmp_path), pytest.raises(RuntimeError, match="download"):
        ConductivityProfile()


def test_a_profile_has_seven_layers_in_the_documented_order(tmp_path):
    write_synthetic_litho(tmp_path / ".scubas_config" / "LITHO1.0.nc")
    with chdir(tmp_path):
        site = ConductivityProfile.compile_profile([30.0, -70.0], site_name="x")
    assert site.get_names() == [
        "Seawater",
        "Sediment",
        "Crust",
        "Lithosphere",
        "Upper Mantle",
        "Transition Zone",
        "Lower Mantle",
    ]
    assert site.get_resistivities() == pytest.approx(
        [0.3, 3.0, 3000.0, 1000.0, 100.0, 10.0, 1.0]
    )
    assert site.name == "x"


def test_the_transition_zone_and_lower_mantle_thicknesses_are_fixed_by_the_defaults(
    tmp_path,
):
    write_synthetic_litho(tmp_path / ".scubas_config" / "LITHO1.0.nc")
    with chdir(tmp_path):
        cp = ConductivityProfile()
    assert cp.get_transition_zone_layer() == pytest.approx(660.0 - 410.0)
    assert cp.get_lower_mantle_layer() == pytest.approx(1000.0 - 660.0)


def test_interpolation_points_include_both_edges(tmp_path):
    write_synthetic_litho(tmp_path / ".scubas_config" / "LITHO1.0.nc")
    with chdir(tmp_path):
        pts = ConductivityProfile().get_interpolation_points(
            [30.0, -70.0], [35.0, -40.0]
        )
    np.testing.assert_allclose(pts[0], [30.0, -70.0])
    np.testing.assert_allclose(pts[-1], [35.0, -40.0])
    assert pts.shape == (round(np.hypot(5.0, 30.0)) + 2, 2)


@pytest.mark.litho
def test_the_real_model_gives_an_ocean_profile(tmp_path):
    path = os.environ.get("SCUBAS_LITHO1_PATH")
    if not path or not Path(path).exists():
        pytest.skip("set SCUBAS_LITHO1_PATH to a LITHO1.0.nc to run this test")
    cache = tmp_path / ".scubas_config"
    cache.mkdir()
    (cache / "LITHO1.0.nc").symlink_to(
        Path(path).resolve()
    )  # in place, so nothing is downloaded
    with chdir(tmp_path):
        site = ConductivityProfile.compile_profile([40.0, -60.0])
    assert site.get_thicknesses(0) > 0.0  # open ocean: seawater has a thickness

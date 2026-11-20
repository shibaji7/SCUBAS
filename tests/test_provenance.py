# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Parameters with provenance: the pattern for new code, and the read-only list of the legacy constants, which must equal what the legacy code uses."""

import dataclasses
import importlib.util
import inspect
import math
import re
from pathlib import Path

import numpy as np
import pytest

from scubas.cables import CableSection, TransmissionLine
from scubas.conductivity import ConductivityProfile
from scubas.errors import ConfigError, MissingParameterError, ScubasError
from scubas.models import OceanModel, Preprocess
from scubas.provenance import (
    TBD,
    Param,
    ParameterSet,
    Status,
    legacy_constants,
    load_parameters,
    parse_parameters,
)
from scubas.sources.current_sources import Current
from scubas.utils import GreatCircle

ROOT = Path(__file__).resolve().parents[1]
LEGACY = legacy_constants()


def _compat():
    spec = importlib.util.spec_from_file_location(
        "compat_cases", ROOT / "tests" / "compat" / "compat_cases.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---- the pattern ------------------------------------------------------------------------------------------------------------------------------
def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


DEFAULTS = """
[line]
length_km = { value = 750.8, unit = "km", status = "assumed", source = "a link of the example" }
tempco = { value = "TBD", unit = "1/K", status = "unknown", source = "needed for the temperature model; not published" }
"""


def test_a_parameter_carries_unit_status_and_source(tmp_path):
    ps = load_parameters(write(tmp_path, "d.toml", DEFAULTS))
    p = ps.param("line", "length_km")
    assert (p.key, p.value, p.unit, p.status, p.source) == (
        "line.length_km",
        750.8,
        "km",
        Status.ASSUMED,
        "a link of the example",
    )
    assert ps.get("line", "length_km") == 750.8 and ps.counts() == {
        "sourced": 0,
        "assumed": 1,
        "unknown": 1,
        "override": 0,
    }


def test_an_unknown_value_raises_instead_of_being_guessed(tmp_path):
    ps = load_parameters(write(tmp_path, "d.toml", DEFAULTS))
    assert (
        not ps.has("line", "tempco")
        and ps.has("line", "length_km")
        and not ps.has("line", "nothing")
    )
    with pytest.raises(
        MissingParameterError,
        match=r"line\.tempco is unknown \(TBD\): needed for the temperature model",
    ):
        ps.get("line", "tempco")
    assert issubclass(MissingParameterError, KeyError) and issubclass(
        MissingParameterError, ScubasError
    )


def test_shipped_parameters_must_be_complete_tables(tmp_path):
    for text, message in (
        ("[a]\nx = 1.0\n", "needs value, unit, status and source"),
        (
            '[a]\nx = { unit = "m", status = "assumed", source = "s" }\n',
            "needs a 'value'",
        ),
        ('[a]\nx = { value = 1.0, unit = "m", source = "s" }\n', "missing 'status'"),
        (
            '[a]\nx = { value = 1.0, unit = "m", status = "assumed" }\n',
            "missing 'source'",
        ),
        (
            '[a]\nx = { value = 1.0, status = "assumed", source = "s" }\n',
            "missing 'unit'",
        ),
        (
            '[a]\nx = { value = 1.0, unit = "m", status = "guess", source = "s" }\n',
            "not in",
        ),
        (
            '[a]\nx = { value = 1.0, unit = "m", status = "unknown", source = "s" }\n',
            "needs the value 'TBD'",
        ),
    ):
        with pytest.raises(ConfigError, match=message):
            load_parameters(write(tmp_path, "bad.toml", text))
    with pytest.raises(ConfigError, match="not found"):
        load_parameters(tmp_path / "missing.toml")


def test_a_user_file_may_give_bare_values_which_become_overrides(tmp_path):
    d = write(tmp_path, "d.toml", DEFAULTS)
    ps = load_parameters(
        d, write(tmp_path, "u.toml", "[line]\nlength_km = 100.0\ntempco = 0.00393\n")
    )
    p = ps.param("line", "length_km")
    assert (p.value, p.unit, p.status, p.source) == (
        100.0,
        "km",
        Status.OVERRIDE,
        "user file",
    )
    assert ps.get("line", "tempco") == 0.00393  # the unknown value is now given
    with pytest.raises(ConfigError, match=r"unknown parameter line\.typo"):
        load_parameters(d, write(tmp_path, "typo.toml", "[line]\ntypo = 1.0\n"))
    with pytest.raises(ConfigError, match=r"unknown parameter other\.x"):
        load_parameters(d, write(tmp_path, "sec.toml", "[other]\nx = 1.0\n"))


def test_overrides_return_a_copy_and_the_original_is_untouched(tmp_path):
    ps = load_parameters(write(tmp_path, "d.toml", DEFAULTS))
    new = ps.with_overrides({"line.length_km": 12.5})
    assert (
        new.get("line", "length_km") == 12.5
        and new.param("line", "length_km").status is Status.OVERRIDE
    )
    assert "override of: a link of the example" in new.param("line", "length_km").source
    assert (
        ps.get("line", "length_km") == 750.8 and ps.fingerprint() != new.fingerprint()
    )
    for bad in ("length_km", "line.nothing", "nothing.x"):
        with pytest.raises(ConfigError, match="no parameter"):
            ps.with_overrides({bad: 1.0})
    assert (
        load_parameters(
            write(tmp_path, "d2.toml", DEFAULTS), overrides={"line.length_km": 1.0}
        ).get("line", "length_km")
        == 1.0
    )


def test_a_param_is_frozen_and_the_set_reports_itself(tmp_path):
    p = Param("s", "n", 1.0, "m", Status.SOURCED, "a paper")
    with pytest.raises(dataclasses.FrozenInstanceError):
        p.value = 2.0
    assert (
        p.is_unknown is False
        and Param("s", "n", TBD, "", Status.UNKNOWN, "why").is_unknown
    )
    ps = load_parameters(write(tmp_path, "d.toml", DEFAULTS))
    assert (
        ps.sections() == ["line"]
        and len(ps.params()) == 2
        and ps.table_md().splitlines()[0].startswith("| key | value")
    )
    assert (
        ps.fingerprint() == load_parameters(tmp_path / "d.toml").fingerprint()
        and len(ps.fingerprint()) == 12
    )
    with pytest.raises(ConfigError, match="must be a table"):
        parse_parameters({"a": 1})
    assert isinstance(parse_parameters({"a": {"x": 1.0}}, from_user=True), ParameterSet)


# ---- the legacy constants equal what the code uses ----------------------------------------------------------------------------------------
def _cprop(name):
    def get():
        compat = _compat()
        import os
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            compat.write_synthetic_litho(Path(tmp) / ".scubas_config" / "LITHO1.0.nc")
            old = Path.cwd()
            os.chdir(tmp)
            try:
                return getattr(ConductivityProfile().cprop, name)
            finally:
                os.chdir(old)

    return get


def _default(obj, parameter):
    return lambda: inspect.signature(obj).parameters[parameter].default


def _source(obj, pattern, cast=float):
    def get():
        found = re.search(pattern, inspect.getsource(obj))
        assert found, pattern
        return cast(found.group(1))

    return get


def _section_geometry():
    def formula(key):
        lat0, lat1, lon0, lon1 = 10.0, 11.0, 20.0, 21.0
        lam = 0.5 * (lat0 + lat1)
        cos2 = math.cos(math.radians(2 * lam))
        if key == "north":
            return (
                LEGACY.get("geometry", "lat_km_per_degree_a")
                - LEGACY.get("geometry", "lat_km_per_degree_b") * cos2
            ) * abs(lat1 - lat0)
        east = (
            LEGACY.get("geometry", "lon_km_per_degree_a")
            - LEGACY.get("geometry", "lon_km_per_degree_b") * cos2
        )
        return east * math.cos(math.radians(lam)) * abs(lon0 - lon1)

    section = CableSection(
        "g",
        {
            "edge_locations": {
                "initial": {"lat": 10.0, "lon": 20.0},
                "final": {"lat": 11.0, "lon": 21.0},
            }
        },
    )
    return formula("north") / section.length_north, formula(
        "east"
    ) / section.length_east


CHECKS = {
    "conductivity.seawater_resistivity_ohm_m": _cprop("seawater_resistivity"),
    "conductivity.sediment_resistivity_ohm_m": _cprop("sediment_resistivity"),
    "conductivity.crust_resistivity_ohm_m": _cprop("crust_resistivity"),
    "conductivity.lithosphere_resistivity_ohm_m": _cprop("lithosphere_resistivity"),
    "conductivity.asthenosphere_resistivity_ohm_m": _cprop("asthenosphere_resistivity"),
    "conductivity.transition_zone_resistivity_ohm_m": _cprop(
        "transition_zone_resistivity"
    ),
    "conductivity.lower_mantle_resistivity_ohm_m": _cprop("lower_mantle_resistivity"),
    "conductivity.transition_zone_top_km": _cprop("transition_zone_top"),
    "conductivity.transition_zone_bottom_km": _cprop("transition_zone_bot"),
    "conductivity.profile_max_depth_km": _cprop("profile_max_depth"),
    "conductivity.grid_interpolation_method": _cprop("grid_interpolation_method"),
    "conductivity.earth_model_file": _cprop("earth_model"),
    "conductivity.earth_model_uri": _cprop("uri"),
    "preprocess.taper_fraction": _default(Preprocess, "p"),
    "preprocess.baseline_samples": _source(
        Preprocess.detrend_magnetic_field, r"min\((\d+), len", int
    ),
    "models.flim_min_Hz": lambda: _default(OceanModel.__init__, "flim")()[0],
    "models.flim_max_Hz": lambda: _default(OceanModel.__init__, "flim")()[1],
    "cables.width_m": lambda: TransmissionLine("c", {"length": 1.0}).elec_params.width,
    "cables.flim_min_Hz": lambda: TransmissionLine(
        "c", {"length": 1.0}
    ).elec_params.flim[0],
    "cables.flim_max_Hz": lambda: TransmissionLine(
        "c", {"length": 1.0}
    ).elec_params.flim[1],
    "geometry.north_east_split_divisor": lambda: (
        CableSection("d", {"length": 2.0}).length
        / CableSection("d", {"length": 2.0}).length_north
    ),
    "geometry.lat_km_per_degree_a": _source(
        CableSection.compute_lengths, r"\((111\.133) -"
    ),
    "geometry.lat_km_per_degree_b": _source(
        CableSection.compute_lengths, r"- (0\.56) \*"
    ),
    "geometry.lon_km_per_degree_a": _source(
        CableSection.compute_lengths, r"\((111\.5065) -"
    ),
    "geometry.lon_km_per_degree_b": _source(
        CableSection.compute_lengths, r"- (0\.1872) \*"
    ),
    "utils.earth_radius_km": lambda: next(
        f.default for f in dataclasses.fields(GreatCircle) if f.name == "Re"
    ),
    "sources.current_height_km": _default(Current.__init__, "h"),
    "sources.cauchy_width": _default(Current.__init__, "a"),
    "sources.current_amplitude": _default(Current.__init__, "I"),
    "readers.iaga_missing_values": lambda: [
        float(x)
        for x in re.search(
            r"na_values=\[(.*?)\]", inspect.getsource(OceanModel.read_iaga)
        )
        .group(1)
        .split(",")
    ],
    "readers.gap_spline_order": _source(
        OceanModel.read_Bfield_data, r"order=(\d+)", int
    ),
    "monte_carlo.default_draws": _default(
        ConductivityProfile.compile_mcmc_bined_profiles, "n"
    ),
    "monte_carlo.default_random_seed": _default(
        ConductivityProfile.compile_mcmc_bined_profiles, "random_seed"
    ),
    "monte_carlo.shelf_thickness_min_km": _source(
        ConductivityProfile.compile_mcmc_bined_profiles, r"\[(0\.01), 1\.0\]"
    ),
    "monte_carlo.shelf_thickness_max_km": _source(
        ConductivityProfile.compile_mcmc_bined_profiles, r"\[0\.01, (1\.0)\]"
    ),
}


def test_every_legacy_constant_has_a_check_and_every_check_has_a_constant():
    assert {p.key for p in LEGACY.params()} == set(CHECKS)


@pytest.mark.parametrize("key", sorted(CHECKS))
def test_a_legacy_constant_equals_what_the_code_uses(key):
    section, name = key.split(".")
    listed, used = LEGACY.get(section, name), CHECKS[key]()
    if isinstance(listed, list):
        assert listed == used
    elif isinstance(listed, float):
        assert used == pytest.approx(listed, rel=1e-12)
    else:
        assert listed == used


def test_the_geometry_constants_reproduce_the_edge_length_of_a_section():
    """The four constants listed for the length of a degree give exactly the north and east lengths of ``CableSection``."""
    north, east = _section_geometry()
    assert north == pytest.approx(1.0, rel=1e-12) and east == pytest.approx(
        1.0, rel=1e-12
    )
    assert LEGACY.get("geometry", "north_east_split_divisor") == pytest.approx(
        np.sqrt(2.0)
    )


def test_the_registry_is_honest_about_its_origin():
    assert LEGACY.counts()["unknown"] == 0 and LEGACY.counts()["override"] == 0
    for p in LEGACY.params():
        assert p.source.strip() and p.unit is not None, p.key
        if p.status is Status.SOURCED:
            assert re.search(r"Chakraborty|IAGA|IRIS", p.source), (
                f"{p.key} is 'sourced' but names no document"
            )
    # a number the code does not document is "assumed": it must not carry an invented citation
    assert all(
        not re.search(r"\(\d{4}\)", p.source)
        for p in LEGACY.params()
        if p.status is Status.ASSUMED
    )


def test_the_docs_page_lists_the_registry_exactly():
    page = (ROOT / "docs" / "dev" / "parameters.md").read_text()
    block = re.search(
        r"<!-- PARAMS:START -->\n(.*?)\n<!-- PARAMS:END -->", page, re.S
    ).group(1)
    assert block == LEGACY.table_md(), "run `python scripts/gen_param_table.py`"

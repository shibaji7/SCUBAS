# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Every import path that users, the tutorials and the examples rely on keeps working and keeps its side-effect profile.

If a module ever moves, its old path must stay importable (a shim that re-exports the same object and warns); the identity
checks below are written so that they pass for a shim too.
"""

import importlib
import os
import subprocess
import sys
from pathlib import Path

import pytest

import scubas

# the subprocesses must import the same scubas as this test run (not a copy installed in site-packages)
ENV = {**os.environ, "PYTHONPATH": str(Path(scubas.__file__).resolve().parents[1])}

# (module, names) as used in docs/, examples/ and tests/
IMPORTS = [
    ("scubas.cables", ["Cable", "CableSection", "TransmissionLine"]),
    ("scubas.conductivity", ["ConductivityProfile"]),
    ("scubas.datasets", ["Layer", "Site", "PROFILES"]),
    ("scubas.models", ["OceanModel", "Preprocess"]),
    (
        "scubas.plotlib",
        [
            "PlotArtifacts",
            "update_rc_params",
            "plot_transfer_function",
            "potential_along_section",
            "cable_potential",
        ],
    ),
    (
        "scubas.utils",
        [
            "RecursiveNamespace",
            "frexp102str",
            "fft",
            "ifft",
            "component_mappings",
            "component_sign_mappings",
            "GreatCircle",
        ],
    ),
    ("scubas.sources.current_sources", ["Current", "Direction"]),
]


@pytest.mark.parametrize(("module", "names"), IMPORTS, ids=[m for m, _ in IMPORTS])
def test_documented_names_import(module, names):
    mod = importlib.import_module(module)
    for name in names:
        assert hasattr(mod, name), f"{module}.{name}"


def test_the_profiles_are_one_shared_object():
    from scubas import cables, datasets
    from scubas.sources import current_sources

    assert cables.PROFILES is datasets.PROFILES is current_sources.PROFILES


def test_import_scubas_is_light_and_has_no_side_effects(tmp_path):
    """The package root imports nothing heavy and writes nothing (a version attribute may be added, but no imports)."""
    out = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys, os, scubas; print(sorted(m for m in ('matplotlib','pandas','scipy','pyproj','loguru') if m in sys.modules)); print(os.listdir('.'))",
        ],
        capture_output=True,
        text=True,
        check=True,
        cwd=tmp_path,
        env=ENV,
    ).stdout.splitlines()
    assert out == ["[]", "[]"]


def test_importing_the_conductivity_module_does_not_touch_the_disk_or_network(tmp_path):
    out = subprocess.run(
        [
            sys.executable,
            "-W",
            "ignore",
            "-c",
            "import os, scubas.conductivity; print(os.listdir('.'))",
        ],
        capture_output=True,
        text=True,
        check=True,
        cwd=tmp_path,
        env=ENV,
    ).stdout.strip()
    assert (
        out == "[]"
    )  # the cache directory appears only when a ConductivityProfile is built


def test_scubas_has_no_console_script_and_no_main_module():
    """There is no CLI today; adding one is a new feature that needs approval, not a refactor."""
    from importlib import metadata

    assert not [
        e
        for e in metadata.entry_points(group="console_scripts")
        if e.value.startswith("scubas")
    ]
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("scubas.__main__")

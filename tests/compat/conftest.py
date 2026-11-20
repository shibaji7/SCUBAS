# SPDX-License-Identifier: MIT
# Copyright (c) 2023 Shibaji Chakraborty
"""Shared setup of the compatibility suite."""

import os
from pathlib import Path

import pytest

import scubas

ROOT = Path(__file__).resolve().parents[2]


def pytest_collection_modifyitems(config, items):
    """Skip the tests that need a real LITHO1.0 model unless ``-m litho`` asks for them (they skip themselves when the file is absent)."""
    if "litho" in (config.getoption("-m") or ""):
        return
    for item in items:
        if item.get_closest_marker("litho"):
            item.add_marker(
                pytest.mark.skip(
                    reason="needs a real LITHO1.0.nc: run with -m litho and SCUBAS_LITHO1_PATH"
                )
            )


@pytest.fixture(scope="session", autouse=True)
def _tests_run_against_the_repository_code():
    """The suite must exercise this checkout, not a copy in site-packages; set SCUBAS_TEST_INSTALLED=1 to test an installed wheel on purpose."""
    if os.environ.get("SCUBAS_TEST_INSTALLED"):
        return
    assert Path(scubas.__file__).resolve().is_relative_to(ROOT), (
        f"scubas imported from {scubas.__file__}, not from {ROOT}"
    )

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""The README stays true: its quick start runs, and the facts it states agree with the repository."""

import re
from pathlib import Path

import pytest
from loguru import logger

ROOT = Path(__file__).resolve().parents[1]
README = (ROOT / "README.md").read_text()


@pytest.fixture(autouse=True)
def _quiet():
    logger.disable("scubas")
    yield
    logger.enable("scubas")


def test_the_quick_start_runs_and_gives_antisymmetric_end_potentials():
    section = README.split("## Quick start", 1)[1]
    code = re.search(r"```python\n(.*?)```", section, re.S).group(1)
    code = code.replace("print(", "result = (")
    scope: dict = {}
    exec(compile(code, "README quick start", "exec"), scope)
    row = scope["result"]
    assert row["U0"] == pytest.approx(-row["U1"], abs=1e-6) and row["U1"] > 0


def test_the_status_line_names_the_current_version():
    from scubas._version import __version__ as version

    assert f"Status: **{version} (alpha)**" in README


def test_every_file_and_folder_the_readme_points_at_exists():
    for name in re.findall(
        r"`((?:scubas|tests|docs|examples)/[A-Za-z0-9_./]*)`", README
    ):
        assert (ROOT / name.rstrip("/")).exists(), name
    for name in re.findall(r"\]\(([A-Z_]+\.md)\)", README):
        assert (ROOT / name).exists(), name


def test_the_license_text_is_the_project_license():
    assert (
        "LGPL-3.0-or-later" in README
        and (ROOT / "COPYING").exists()
        and (ROOT / "NOTICE.md").exists()
    )

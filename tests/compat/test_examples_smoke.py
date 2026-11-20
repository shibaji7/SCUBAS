# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""The example scripts that regenerate the tutorial figures still run from start to finish.

They write into ``<repo>/docs/tutorial/figures`` relative to their own location, so they are run from a temporary copy of
``examples/`` and the repository is not touched.
"""

import ast
import os
import runpy
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = sorted(
    p.name
    for p in (ROOT / "examples").glob("network_*.py")
    if ast.parse(
        p.read_text()
    ).body  # network_ls.py is an empty placeholder (a header only)
)


@pytest.fixture
def sandbox(tmp_path, monkeypatch):
    shutil.copytree(ROOT / "examples", tmp_path / "examples")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setitem(os.environ, "MPLBACKEND", "Agg")
    return tmp_path


@pytest.mark.parametrize("script", SCRIPTS)
def test_example_runs_and_writes_its_figures(script, sandbox):
    pytest.importorskip("matplotlib")
    runpy.run_path(str(sandbox / "examples" / script), run_name="__main__")
    figures = list((sandbox / "docs" / "tutorial" / "figures").glob("*.png"))
    assert figures, f"{script} wrote no figure"


def test_the_unchanged_repository_is_not_written_to(sandbox):
    before = sorted(
        p.name for p in (ROOT / "docs" / "tutorial" / "figures").glob("*.png")
    )
    runpy.run_path(str(sandbox / "examples" / "network_model.py"), run_name="__main__")
    after = sorted(
        p.name for p in (ROOT / "docs" / "tutorial" / "figures").glob("*.png")
    )
    assert before == after

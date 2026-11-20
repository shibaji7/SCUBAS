# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""The pieces of ``docs/conf.py`` that are code: the development-status banner and the colour classes of the roadmap statuses."""

import importlib.util
import re
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def conf():
    spec = importlib.util.spec_from_file_location(
        "docs_conf", ROOT / "docs" / "conf.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def app():
    seen = []
    return SimpleNamespace(env=SimpleNamespace(note_dependency=seen.append), seen=seen)


def test_the_version_comes_from_pyproject(conf):
    assert re.fullmatch(r"\d+\.\d+\.\d+", conf.release)
    assert conf.version == ".".join(conf.release.split(".")[:2])


def test_every_page_starts_with_the_banner_and_keeps_its_text(conf, app):
    source = ["# Title\n\ntext\n"]
    conf.add_banner(app, "dev/cables", source)
    assert source[0].startswith(":::{warning}\n**Development status: alpha.**")
    assert source[0].endswith("# Title\n\ntext\n")
    assert "(/user/citing.md)" in source[0]


def test_front_matter_stays_first(conf, app):
    source = ["---\norphan: true\n---\n\n# Title\n"]
    conf.add_banner(app, "user/docker", source)
    assert source[0].startswith("---\norphan: true\n---\n")
    assert ":::{warning}" in source[0] and source[0].endswith("# Title\n")


def test_the_roadmap_page_is_the_root_file_with_coloured_statuses(conf, app):
    source = ["```{include} ../ROADMAP.md\n```\n"]
    conf.add_banner(app, "roadmap", source)
    text = source[0]
    assert "# Roadmap" in text and ":::{warning}" in text
    assert (
        text.count("{.status-done}") >= 4
        and "{.status-planned}" in text
        and "{.status-proposed}" in text
    )
    assert app.seen == [
        str(ROOT / "ROADMAP.md")
    ]  # a change of the file rebuilds the page


def test_other_pages_are_not_given_status_colours(conf, app):
    source = ["The step is **done**.\n"]
    conf.add_banner(app, "glossary", source)
    assert "{.status-done}" not in source[0]


def test_the_roadmap_progress_line_matches_the_table():
    text = (ROOT / "ROADMAP.md").read_text()
    rows = [
        line
        for line in text.splitlines()
        if line.startswith("| ") and "---" not in line and "Step" not in line
    ]
    done = sum("**done**" in line for line in rows)
    planned = sum("**planned**" in line for line in rows)
    proposed = sum("**proposed**" in line for line in rows)
    assert f"**Progress: {done} of {len(rows)} steps done**" in text
    assert f"{proposed} proposed, {planned} planned" in text


def test_the_toolchain_files_agree():
    rtd = (ROOT / ".readthedocs.yaml").read_text()
    assert (
        "builder: dirhtml" in rtd
        and "configuration: docs/conf.py" in rtd
        and "fail_on_warning: true" in rtd
    )
    assert not (ROOT / "mkdocs.yml").exists()

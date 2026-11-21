# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Packaging hygiene: one source of the version, the runtime dependencies of 0.1.9 are kept, the license metadata is right, and the distribution check works."""

import importlib.util
import tomllib
from pathlib import Path

import pytest

import scubas

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = tomllib.loads((ROOT / "pyproject.toml").read_text())


@pytest.fixture(scope="module")
def dist():
    spec = importlib.util.spec_from_file_location(
        "check_dist", ROOT / "scripts" / "check_dist.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_version_has_one_source():
    from scubas._version import __version__

    assert scubas.__version__ == __version__
    assert "version" in PYPROJECT["project"]["dynamic"] and "version" not in {
        k for k in PYPROJECT["project"] if k != "dynamic"
    }
    assert PYPROJECT["tool"]["setuptools"]["dynamic"]["version"] == {
        "attr": "scubas._version.__version__"
    }


def test_the_runtime_dependencies_of_0_1_9_are_all_still_required():
    """ADR 0003: no runtime dependency is removed in 0.x (numpy was added, it was imported but not declared)."""
    deps = set(PYPROJECT["project"]["dependencies"])
    assert {
        "pandas",
        "matplotlib>=3.3.4",
        "pyproj",
        "loguru",
        "scipy",
        "SciencePlots",
    } <= deps and "numpy" in deps


def test_the_extras_exist_and_the_dev_extra_has_the_tools_the_docs_name():
    extras = PYPROJECT["project"]["optional-dependencies"]
    assert set(extras) == {"plots", "data", "dev", "docs"}
    assert {
        "pytest>=8",
        "pytest-cov>=4",
        "ruff",
        "pre-commit",
        "build",
        "twine",
    } <= set(extras["dev"])


def test_the_license_metadata_matches_the_repository():
    project = PYPROJECT["project"]
    assert project["license"] == "LGPL-3.0-or-later"
    assert set(project["license-files"]) == {"LICENSE", "COPYING", "NOTICE.md"}
    assert not any(
        c.startswith("License ::") for c in project["classifiers"]
    )  # superseded by the expression
    assert PYPROJECT["build-system"]["requires"][0] == "setuptools>=77"
    assert all((ROOT / f).exists() for f in project["license-files"])


def test_the_description_and_classifiers_are_about_scubas():
    project = PYPROJECT["project"]
    assert (
        "ionospheric radio sounding" not in project["description"]
        and "Submarine Cables" in project["description"]
    )
    versions = {
        c.rsplit(" ", 1)[1]
        for c in project["classifiers"]
        if "Programming Language :: Python :: 3." in c
    }
    assert (
        versions == {"3.11", "3.12", "3.13"} and project["requires-python"] == ">=3.11"
    )
    assert "pynasonde" not in (ROOT / "MANIFEST.in").read_text()


def test_the_manifest_asks_for_everything_the_dist_check_needs(dist):
    manifest = (ROOT / "MANIFEST.in").read_text()
    for part in (
        "recursive-include docs",
        "recursive-include examples",
        "recursive-include tests",
        "recursive-include scripts",
        "prune docs/_build",
        "NOTICE.md",
    ):
        assert part in manifest
    assert all((ROOT / n).exists() for n in dist.SDIST_NEEDS)


def test_a_good_wheel_and_sdist_pass_the_check(dist):
    wheel = [
        f"scubas/{n}"
        for n in (
            "__init__.py",
            "_version.py",
            "cables.py",
            "conductivity.py",
            "datasets.py",
            "errors.py",
            "future/__init__.py",
            "models.py",
            "provenance.py",
            "defaults/legacy_constants.toml",
            "plotlib.py",
            "sources/current_sources.py",
            "utils.py",
        )
    ]
    wheel += [
        f"scubas-0.1.9.dist-info/licenses/{n}"
        for n in ("LICENSE", "COPYING", "NOTICE.md")
    ]
    sdist = ["scubas-0.1.9/" + n for n in dist.SDIST_NEEDS]
    assert dist.wheel_problems(wheel) == [] and dist.sdist_problems(sdist) == []


def test_missing_and_forbidden_files_are_named(dist):
    assert any("errors.py" in p for p in dist.wheel_problems(["scubas/__init__.py"]))
    assert any(
        "tests/" in p
        for p in dist.wheel_problems(["tests/test_a.py", "scubas/__init__.py"])
    )
    sdist = ["scubas-0.1.9/" + n for n in dist.SDIST_NEEDS if n != "docs/index.md"] + [
        "scubas-0.1.9/docs/_build/x.html",
        "scubas-0.1.9/scubas/__pycache__/a.pyc",
    ]
    got = dist.sdist_problems(sdist)
    assert (
        "sdist is missing docs/index.md" in got
        and any("_build" in p for p in got)
        and any("__pycache__" in p for p in got)
    )

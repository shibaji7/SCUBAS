# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""The release workflow is a safety mechanism, so its safety properties are tested: nothing publishes without the switch, no token is stored, and the approval gate exists."""

import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
TEXT = (ROOT / ".github" / "workflows" / "release.yml").read_text()
WORKFLOW = yaml.safe_load(TEXT)
JOBS = WORKFLOW["jobs"]
PUBLISH_JOBS = ("testpypi", "pypi", "github-release")


def test_it_runs_only_on_version_tags_or_by_hand():
    triggers = WORKFLOW[True]  # YAML reads the key `on` as the boolean True
    assert (
        triggers["push"] == {"tags": ["v*"]}
        and "workflow_dispatch" in triggers
        and "pull_request" not in triggers
    )


def test_every_publishing_job_waits_for_the_switch_and_for_the_job_before_it():
    for name in PUBLISH_JOBS:
        assert JOBS[name]["if"] == "success() && vars.PUBLISH_ENABLED == 'true'", name
    assert (
        JOBS["testpypi"]["needs"] == "build"
        and JOBS["pypi"]["needs"] == "testpypi"
        and JOBS["github-release"]["needs"] == "pypi"
    )


def test_pypi_needs_an_environment_for_the_manual_approval_and_testpypi_comes_first():
    assert (
        JOBS["pypi"]["environment"] == "pypi"
        and JOBS["testpypi"]["environment"] == "testpypi"
    )


def test_only_the_publishing_jobs_get_more_than_read_permission():
    assert (
        WORKFLOW["permissions"] == {"contents": "read"}
        and "permissions" not in JOBS["build"]
    )
    assert JOBS["testpypi"]["permissions"] == {"id-token": "write"} and JOBS["pypi"][
        "permissions"
    ] == {"id-token": "write"}
    assert JOBS["github-release"]["permissions"] == {"contents": "write"}


def test_no_token_or_password_is_used_trusted_publishing_only():
    assert not re.search(r"(?i)password|api[_-]?token|pypi_token|secrets\.", TEXT)
    assert TEXT.count("pypa/gh-action-pypi-publish@release/v1") == 2


def test_the_build_job_refuses_a_tag_that_does_not_match_and_checks_the_distributions():
    build = "\n".join(step.get("run", "") for step in JOBS["build"]["steps"])
    assert (
        'check_release.py --tag "${GITHUB_REF_NAME}"' in build
        and "check_dist.py" in build
        and "twine check" in build
    )
    assert "scubas/_version.py" in build and "release_notes.md" in build


def test_the_makefile_has_the_switch_targets_and_the_guide_describes_the_setup():
    make = (ROOT / "Makefile").read_text()
    for target in (
        "publish-on:",
        "publish-off:",
        "publish-status:",
        "publish-tag:",
        "codecov-off:",
    ):
        assert target in make
    guide = (ROOT / "docs" / "dev" / "release.md").read_text()
    for part in (
        "PUBLISH_ENABLED",
        "trusted publisher",
        "Required reviewers",
        "git tag vX.Y.Z",
        "testpypi",
        "make publish-tag",
    ):
        assert part in guide

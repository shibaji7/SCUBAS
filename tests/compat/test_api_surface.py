# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""The public surface of SCUBAS 0.1.9 is pinned: no name disappears and nothing recorded about a name changes.

Adding a name is compatible and allowed. To change the baseline on purpose (an approved API change) regenerate it with
``make compat-baseline`` and explain the diff in ``CHANGELOG.md``.
"""

import json
from pathlib import Path

import pytest
from compat_surface import collect_surface

BASELINE = json.loads(
    (Path(__file__).parent / "baseline" / "api_surface.json").read_text()
)


@pytest.fixture(scope="module")
def current():
    return collect_surface()


def _diff(where: str, keys: tuple[str, ...], base: dict, now: dict) -> list[str]:
    """One line per recorded key whose value changed."""
    return [
        f"{where}: {key} {base.get(key)} -> {now.get(key)}"
        for key in keys
        if base.get(key) != now.get(key)
    ]


def _compare(where: str, base: dict, now: dict | None) -> list[str]:
    if now is None:
        return [f"{where}: missing"]
    if base["kind"] != now["kind"]:
        return [f"{where}: kind {base['kind']} -> {now['kind']}"]
    if base["kind"] in ("function", "assignment"):
        return _diff(where, ("signature", "type"), base, now)
    out = _diff(where, ("bases", "dataclass", "members"), base, now)
    for name, spec in base["methods"].items():
        got = now["methods"].get(name)
        if got is None:
            out.append(f"{where}.{name}: missing")
        elif got != spec:
            out.append(f"{where}.{name}: {spec} -> {got}")
    return out


@pytest.mark.parametrize("module", sorted(BASELINE["modules"]))
def test_module_surface_is_unchanged(module, current):
    problems = []
    for name, spec in BASELINE["modules"][module].items():
        problems += _compare(
            f"{module}.{name}", spec, current["modules"][module].get(name)
        )
    assert not problems, "\n".join(problems)


def test_profiles_are_unchanged(current):
    """``scubas.datasets.PROFILES``: the layered Earth models are data that published results depend on."""
    problems = []
    for key, spec in BASELINE["profiles"].items():
        if current["profiles"].get(key) != spec:
            problems.append(key)
    assert not problems, f"PROFILES entries missing or changed: {problems}"
    assert len(BASELINE["profiles"]) == 19


def test_the_known_legacy_names_are_still_there(current):
    """A few names that look like typos are public API: users call them."""
    cables = current["modules"]["scubas.cables"]
    assert "calc_trasmission_line_parameters" in cables["TransmissionLine"]["methods"]
    assert "solve_admitance_matrix" in cables["Cable"]["methods"]
    assert "_update_recursive_namespace" in cables

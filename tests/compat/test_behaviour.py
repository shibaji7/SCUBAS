# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Characterization: the numbers SCUBAS 0.1.9 produces (tutorial cables, the 1989 event, readers, sources, utilities).

The baseline is ``baseline/behaviour.npz`` written by ``scripts/gen_compat_baseline.py`` from the unchanged legacy code.
Arrays are compared with a relative tolerance of 1e-9 and an absolute tolerance of 1e-10 of the array's own scale, which is
far below any physical meaning and above BLAS and platform noise. Long arrays are stored thinned plus summary statistics.
"""

import json
from pathlib import Path

import numpy as np
import pytest
from compat_cases import build_cases

BASE = Path(__file__).parent / "baseline"
EXPECTED = np.load(BASE / "behaviour.npz")
EXPECTED_META = json.loads((BASE / "behaviour_meta.json").read_text())
RTOL = 1e-9


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    arrays, meta = build_cases(tmp_path_factory.mktemp("compat"))
    return arrays, json.loads(json.dumps(meta, default=str))


def test_the_same_arrays_are_produced(run):
    arrays, _ = run
    assert sorted(arrays) == sorted(EXPECTED.files)


def _groups():
    return sorted({name.split("|")[0] for name in EXPECTED.files})


@pytest.mark.parametrize("group", _groups())
def test_values_are_unchanged(group, run):
    arrays, _ = run
    bad = []
    for name in EXPECTED.files:
        if name.split("|")[0] != group:
            continue
        want, got = EXPECTED[name], arrays[name]
        if want.dtype.kind in "iu":
            if not np.array_equal(want, got):
                bad.append(name)
            continue
        atol = 1e-10 * float(np.nanmax(np.abs(want))) if want.size else 0.0
        if name.endswith("#stats"):
            atol = 1e-9 * abs(float(want[3]))
        try:
            np.testing.assert_allclose(got, want, rtol=RTOL, atol=atol, equal_nan=True)
        except AssertionError:
            bad.append(name)
    assert not bad, f"{len(bad)} arrays moved: {bad[:8]}"


def test_metadata_is_unchanged(run):
    _, meta = run
    assert meta == EXPECTED_META


def test_the_1989_event_keeps_its_headline_numbers(run):
    """A human-readable anchor in case the baseline file is ever regenerated: the TAT-8 potential peaks stay in the published range."""
    arrays, _ = run
    stats = arrays["event1989|tot_params|Vt(v)#stats"]
    assert stats[0] == 4320
    assert 390.0 < stats[4] < 395.0  # largest |V_TAT-8| in volts

# Conventions

The rules enforced by `make qa` and by the same checks in CI and pre-commit. SCUBAS is older code, so the existing violations are
**frozen, not fixed**; the rule is that new and touched code is clean.

## Lint, format and headers

- **ruff** (`ruff.toml`): the BOLTS rule set (pycodestyle, pyflakes, isort, bugbear, pyupgrade, simplifications, naming, NumPy-style
  docstrings, pathlib, annotations, complexity 10, unused arguments). Line width 88, the width the legacy code was formatted with.
- **Legacy debt** is a per-file list in `ruff.toml` (between `BEGIN LEGACY DEBT` and `END LEGACY DEBT`) written by
  `scripts/ruff_ratchet.py`. `make lint-ratchet` fails on a new violation and on a stale entry; fix a file's violations and run
  `python scripts/ruff_ratchet.py --update` to shrink the list. Never add an entry by hand. Two legacy files are excluded from
  `ruff format` (`scubas/cables.py`, `scubas/sources/current_sources.py`) until they are edited for another reason.
- **SPDX header** on every Python file: the license line and the copyright line (`scripts/check_spdx.py`, which holds the text).
- **Names.** The physics keeps the symbols of the transmission-line literature (`Z`, `Y`, `Z0`, `gma`, `Ye`, `Yp2`, `Ie`) and the
  names users already call, including the misspelt `calc_trasmission_line_parameters` and `solve_admitance_matrix`; a linter never renames them.
  New code puts the unit in the name (`length_km`, `e_mV_per_km`, `thickness_m`).

## Changes

- Every change under `scubas/` has a line under `## [Unreleased]` in `CHANGELOG.md` (the pre-commit hook `changelog` checks it;
  `SKIP=changelog git commit` bypasses it for a change nobody will see).
- A public name that moves or is renamed keeps its old import path through a thin shim that re-exports the same object and emits a
  `DeprecationWarning` naming the new path and the version that removes it (not before 1.0.0, and not sooner than two minor releases
  after the warning first ships). The compatibility suite has the identity checks.
- New or changed physics: a docstring with the equation and its source (paper, book and page), a test with an analytic expectation,
  a CHANGELOG line and a paragraph in these docs. The pull request template has the checklist.

## Parameters

New parameters carry a unit, a status and a source. The pattern (`Param(value, unit, status, source)` with status `sourced`,
`assumed`, `unknown` or `override`, an unknown value raising a clear error) is used for **new code only**; the constants of the
existing code (resistivities, the taper, the 120-sample baseline, the geodesic constants) get their provenance gradually, and are
documented before they are moved.

## Docstrings

NumPy style (`Parameters`, `Returns`, `Raises`), as in the existing modules.

## Plots

Figures follow the SCUBAS figure style (one palette, font and layout). `scubas.plotlib.update_rc_params` is the existing helper for matplotlib settings.

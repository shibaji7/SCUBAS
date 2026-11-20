# Changelog

All notable changes are listed here. Format: [Keep a Changelog](https://keepachangelog.com/), versioning: [SemVer](https://semver.org/).

SCUBAS is in alpha (0.x). **Compatibility rule while the engineering standards are being adopted:** every public import path, function signature, file format and numerical result of 0.1.9 is pinned by `tests/compat/`; a change to any of them needs its own approved, documented entry here.

## [Unreleased]
### Added
- A compatibility suite (`tests/compat/`) that pins the public surface and the behaviour of 0.1.9: an API-surface snapshot (every module, class, method signature, the 19 `PROFILES`), a numeric baseline (the three conductivity tutorial cables, the March 1989 TAT-8 event, transfer functions, IAGA-2002 and CSV readers, the current source, FFT and distance helpers, the LITHO1.0 profile code on a synthetic grid), file-format pins (`Cable.save` output, B-field readers and gap handling), closed-form checks (long and short cable end potentials, half-space impedance), the known defects of 0.1.9 recorded as they are, import-path checks, and a smoke test of the example scripts. Regenerate the baselines only on unchanged legacy code: `make compat-baseline`.
- Developer tooling: `Makefile` (`make help`, `make qa`), `ruff.toml` (the BOLTS rule set with the existing violations frozen as a per-file debt list that may only shrink: `scripts/ruff_ratchet.py`), a per-module coverage floor check (`scripts/check_coverage_floors.py`), an SPDX header check with a legacy list (`scripts/check_spdx.py`), a release-consistency check (`scripts/check_release.py`), a changelog guard (`scripts/check_changelog.py`) and `.pre-commit-config.yaml`.

### Changed
- `tests/test_cables.py::test_transmission_line_active_termination_populated` is fixed (test only): an active termination is a `Site`, its derived admittance and current source are stored in `TransmissionLine.term_params`, and its fields need the B-field files; the old test passed a `SimpleNamespace` and read the result from the wrong attribute, so it failed on every run. Package code is unchanged.

## [0.1.9] - 2026-03-04
Published to PyPI. No changelog was kept before this file; see the git history.

## [0.1.8] - 2025-12-01
Published to PyPI.

## [0.1.7] - 2025-12-01
Published to PyPI.

## [0.1.6] - 2025-09-22
Published to PyPI.

## [0.1.5] - 2024-01-31
Published to PyPI (GitHub tag `TLM-0.1.5`).

## [0.1.4] - 2023-05-12
Published to PyPI.

## [0.1.3] - 2023-05-12
Published to PyPI (GitHub release `TLM-0.1.3`, "Transmission Line Modeling").

## [0.1.2] - 2023-03-02
Published to PyPI.

## [0.1.1] - 2023-02-20
First release on PyPI.

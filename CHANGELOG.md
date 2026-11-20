# Changelog

All notable changes are listed here. Format: [Keep a Changelog](https://keepachangelog.com/), versioning: [SemVer](https://semver.org/).

SCUBAS is in alpha (0.x). **Compatibility rule while the engineering standards are being adopted:** every public import path, function signature, file format and numerical result of 0.1.9 is pinned by `tests/compat/`; a change to any of them needs its own approved, documented entry here.

## [Unreleased]
### Fixed
- **pandas 3 is supported.** Reading IAGA-2002 files (`OceanModel.read_iaga`, so `read_Bfield_data` for `.txt` and `.min` input and `TransmissionLine.compile_oml`) failed with pandas 3 (`read_csv() got an unexpected keyword argument 'delim_whitespace'`), and a fresh `pip install scubas` installs pandas 3. `read_iaga` now splits on whitespace with `sep=r"\s+"` and parses the date and time columns itself; results are identical to 0.1.9 on pandas 2 (the compatibility suite is unchanged and passes on pandas 2.3 and 3.1). `ConductivityProfile.compile_mcmc_bined_profiles` set the seawater thickness of the end bins through a chained assignment that pandas 3 silently ignores; it now writes into the frame, with the same random draws in the same order. The only visible difference on pandas 3 is pandas' own: datetime indexes parsed from text have microsecond instead of nanosecond resolution.

### Added
- Release control: `.github/workflows/release.yml` (a pushed tag `vX.Y.Z` is checked against `scubas/_version.py` and the changelog, built and checked, published to TestPyPI, then to PyPI after a manual approval of the `pypi` environment with PyPI trusted publishing (no token is stored), then a GitHub Release from the changelog section), switched off until the repository variable `PUBLISH_ENABLED` is `true`; `make publish-on`, `publish-off`, `publish-status`, `publish-tag` and `codecov-on`, `codecov-off`, `codecov-status`; the release guide with the owner's one-time setup (`docs/dev/release.md`); ADR 0009; `tests/test_release_workflow.py`. Nothing is published, tagged or released by this change, and no repository setting was changed.
- Packaging hygiene, additive: `scubas.__version__` (a new `scubas/_version.py`, the single source of the version, read by `pyproject.toml`, the docs and the release check); `numpy` is declared (it was imported but not listed); the extras `plots`, `data`, `dev` and `docs`, next to the unchanged runtime dependencies (pandas, matplotlib, pyproj, loguru, scipy, SciencePlots); a small `environment.yml` and `make env`; `make build` and `make check`; `scripts/check_dist.py`, which checks what is inside the wheel and the sdist and runs in the build workflow; an explicit `MANIFEST.in`.
- Governance and presentation: `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md` (with the git workflow: pull requests to `main`, squash merge, a green CI), `SECURITY.md` (including the plain-HTTP LITHO1.0 download), `CITATION.cff` (preferred citation: Front. Phys. 2022), `codecov.yml`, a README in the BOLTS layout (badges, status, quick start, layout), and a generated text logo (`scripts/make_logo.py`, `make logo`; a placeholder, replace `docs/_static/logo*.svg`) used as the README banner and the docs logo.
- `scubas.future` and `scubas.errors` (`ScubasError`, `NotImplementedYetError`): the pattern for planned capabilities, as stubs that raise `NotImplementedYetError` naming their roadmap item. The package is empty until a capability is named (ADR 0008, `docs/dev/future.md`, `tests/test_future.py`). Nothing existing changes.
- Documentation moves from mkdocs-material to **Sphinx and MyST** (the BOLTS stack: Read the Docs theme, autodoc, napoleon, viewcode, intersphinx), built with the `dirhtml` builder so every existing page keeps its address; the pages stay Markdown, the mkdocs-only parts are gone (`mkdocs.yml`, the card-grid CSS, the lightbox and MathJax shim, `mkdocstrings` blocks, which are now autodoc directives). New: a development-status banner on every page (one `source-read` handler in `docs/conf.py`), a glossary, decision records (`docs/dev/adr/`), developer pages (conventions, testing and the compatibility suite, release), `ROADMAP.md` with colour-coded statuses (done is green) and pages for the roadmap and the changelog, and a Project section on the home page. `.readthedocs.yaml` builds with `fail_on_warning: true` and has an on/off switch for Read the Docs builds (the file `.rtd-skip`, with `make rtd-on`, `make rtd-off` and `make rtd-status`; absent means builds run, as before). `docs/requirements.txt` is now the short pinned list of the Sphinx stack (it was a freeze of about eighty packages the docs did not use); the docs CI job installs the same file, and Dependabot watches it in one weekly group.
- CI gates (`.github/workflows/`): `tests.yml` runs ruff (the legacy debt in `ruff.toml` is allowed, new violations are not), the debt ratchet, `ruff format --check`, the SPDX and release-consistency checks once, and the test suite on Python 3.11, 3.12 and 3.13 with pandas 2 and pandas 3 (overall coverage gate 90%, per-module floors); `build.yml` builds the sdist and wheel, runs `twine check`, installs the wheel in a clean environment and runs the compatibility suite against it; `docs.yml` builds the documentation with `mkdocs build --strict` (tools pinned in `ci/docs-requirements.txt`); `security.yml` runs `pip-audit` on every push and every Monday. Nothing is published by these workflows. Codecov upload moves to `codecov-action@v5` from one test leg; it stays on unless the repository variable `CODECOV_UPLOAD` is set to `false`.
- `.github/dependabot.yml` (GitHub Actions weekly; no `pip` entry, because Dependabot then opened a pull request per package of the old `docs/requirements.txt` freeze that Read the Docs installs), a pull request template whose checklist includes the compatibility suite, and bug and feature issue templates.
- A compatibility suite (`tests/compat/`) that pins the public surface and the behaviour of 0.1.9: an API-surface snapshot (every module, class, method signature, the 19 `PROFILES`), a numeric baseline (the three conductivity tutorial cables, the March 1989 TAT-8 event, transfer functions, IAGA-2002 and CSV readers, the current source, FFT and distance helpers, the LITHO1.0 profile code on a synthetic grid), file-format pins (`Cable.save` output, B-field readers and gap handling), closed-form checks (long and short cable end potentials, half-space impedance), the known defects of 0.1.9 recorded as they are, import-path checks, and a smoke test of the example scripts. Regenerate the baselines only on unchanged legacy code: `make compat-baseline`.
- Developer tooling: `Makefile` (`make help`, `make qa`), `ruff.toml` (the BOLTS rule set with the existing violations frozen as a per-file debt list that may only shrink: `scripts/ruff_ratchet.py`), a per-module coverage floor check (`scripts/check_coverage_floors.py`), an SPDX header check with a legacy list (`scripts/check_spdx.py`), a release-consistency check (`scripts/check_release.py`), a changelog guard (`scripts/check_changelog.py`) and `.pre-commit-config.yaml`.

### Changed
- Package metadata: the description said "Precision ionospheric radio sounding tools" (copied from another project) and the classifiers named Python 3.6 to 3.10 although the package needs Python 3.11; both, the keywords and the project URLs are corrected, and the license is declared as the SPDX expression `LGPL-3.0-or-later` with `license-files` (build requirement `setuptools>=77`, which replaces the unused `setuptools_scm`). The wheel has the same modules as before plus `_version.py`; the only new runtime requirement is `numpy`. The sdist no longer carries the `.github` folder, `.gitignore` and `.pre-commit-config.yaml` (it used to include every tracked file) but does carry `environment.yml` and what the tests need.
- **License: from MIT to LGPL-3.0-or-later.** `LICENSE` is the LGPL text and `COPYING` the GPL text it builds on; `NOTICE.md` keeps the MIT notice for releases 0.1.1 to 0.1.9, which stay available under the MIT license. Every Python file carries an SPDX header (`# SPDX-License-Identifier: LGPL-3.0-or-later`; the text lives in `scripts/check_spdx.py` and the legacy list is gone), `pyproject.toml` declares the license and its classifier, and the README, the badges and the disclaimer comment of every documentation page say so. The only runtime-visible change is the value of the module attribute `__license__` in seven modules (from `"MIT"` to `"LGPL-3.0-or-later"`). No other code changed: the syntax tree of every Python file is identical to before apart from that string.
- The documentation no longer loads a script from `polyfill.io` (a domain that was taken over and served malicious code in 2024); MathJax 3 loads from jsDelivr.
- `tests/test_cables.py::test_transmission_line_active_termination_populated` is fixed (test only): an active termination is a `Site`, its derived admittance and current source are stored in `TransmissionLine.term_params`, and its fields need the B-field files; the old test passed a `SimpleNamespace` and read the result from the wrong attribute, so it failed on every run. This test fix changes no package code.
- Dependabot ignores minor and major updates of `sphinx` and every update of `sphinx-autodoc-typehints` (Read the Docs builds with Python 3.11; Sphinx 9.1 and later, and typehints 3.6.3 and later, need Python 3.12 or 3.13); patch updates of Sphinx 9.0 are still proposed.

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

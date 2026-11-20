# Roadmap

Status key: **done** = merged and tested, **planned** = agreed, not started, **proposed** = waiting for a decision.

**Progress: 5 of 10 steps done** (the documentation step lands with the pull request that moves the docs to Sphinx), 1 proposed, 4 planned.

SCUBAS 0.1.9 (on PyPI) is the baseline. The work below brings the project to the engineering standard of its sibling BOLTS
(tests, gates, documentation, release control) **without changing the API, the file formats or the numbers of 0.1.9**; they are
pinned by the compatibility suite (`tests/compat/`).

## Engineering

| Step | Scope | Status |
|---|---|---|
| Safety net | the compatibility suite, ruff with the legacy lint debt frozen in a list that may only shrink, SPDX, release, changelog and coverage-floor checks, `make qa`, pre-commit, `CHANGELOG.md` | **done** |
| pandas 3 | `read_iaga` and the Monte Carlo end-bin thickness work with pandas 2 and 3, results unchanged | **done** |
| CI gates | tests on Python 3.11 to 3.13 with pandas 2 and 3, build and wheel test, strict docs build, weekly `pip-audit`, Dependabot for the Actions, templates | **done** |
| Documentation | development-status banner on every page, glossary, decision records, testing, conventions and release guides, roadmap and changelog pages, strict build on CI and Read the Docs, the `.rtd-skip` switch | **done** in this change (Sphinx and MyST, like BOLTS) |
| License | relicense from MIT to LGPL-3.0-or-later; earlier releases stay MIT (`NOTICE.md`) | **proposed**: draft pull request, needs the consent of every copyright holder |
| Governance files | `CITATION.cff`, `CONTRIBUTING.md`, `SECURITY.md`, `CODE_OF_CONDUCT.md`, README corrections | **planned** (after the license decision, so the wording is right) |
| Packaging hygiene | correct the package metadata (description, classifiers, keywords), declare numpy, add the `plots`, `data`, `dev` and `docs` extras next to the unchanged runtime dependencies, fix `MANIFEST.in`, a small `environment.yml` (`env.yml` stays), `scubas.__version__` | **planned** |
| Release control | a tag-driven release workflow (TestPyPI, then PyPI after a manual approval, trusted publishing) behind an on/off repository variable; `scripts/check_dist.py` | **planned**, needs the owner to create the environments, the variable and the PyPI trusted publisher |
| Provenance for new code | `Param(value, unit, status, source)` and a clear error for unknown values, used by new code only; a read-only list of the constants of the existing code with a test that they match | **planned** |
| Per-change discipline | pull request checklist, changelog guard, lint debt and coverage floors that only ratchet in one direction | **done** (ongoing) |

Not planned for now: a `src/` layout, splitting the modules, replacing loguru, making matplotlib or pyproj optional. Each would change
how users import or install the package; the reasons are in the decision records (`docs/dev/adr/`). If something ever moves, its old
import path stays through a shim that warns, with the version that removes it (not before 1.0.0).

## Known defects of 0.1.9

Pinned as they are, not fixed (see `docs/dev/testing.md`). Each fix is a separate, approved change.

## Open decisions

The consent for the license change, the version number and date of the next release, the author list for `CITATION.cff`, and which
known defects to fix and how (a new function or option when the output changes). Changes to the model itself are decided by the
maintainers and the published studies; none is committed here.

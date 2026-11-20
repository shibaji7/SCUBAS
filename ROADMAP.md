# Roadmap

Status key: **done** = merged and tested, **planned** = agreed, not started, **proposed** = waiting for a decision.

**Progress: 11 of 11 steps done**, 0 proposed, 0 planned.

SCUBAS 0.1.9 (on PyPI) is the baseline. The work below brings the project to the engineering standard of its sibling BOLTS
(tests, gates, documentation, release control) **without changing the API, the file formats or the numbers of 0.1.9**; they are
pinned by the compatibility suite (`tests/compat/`).

## Engineering

| Step | Scope | Status |
|---|---|---|
| Safety net | the compatibility suite, ruff with the legacy lint debt frozen in a list that may only shrink, SPDX, release, changelog and coverage-floor checks, `make qa`, pre-commit, `CHANGELOG.md` | **done** |
| pandas 3 | `read_iaga` and the Monte Carlo end-bin thickness work with pandas 2 and 3, results unchanged | **done** |
| CI gates | tests on Python 3.11 to 3.13 with pandas 2 and 3, build and wheel test, strict docs build, weekly `pip-audit`, Dependabot for the Actions, templates | **done** |
| Documentation | development-status banner on every page, glossary, decision records, testing, conventions and release guides, roadmap and changelog pages, strict build on CI and Read the Docs, the `.rtd-skip` switch | **done** (Sphinx and MyST, like BOLTS) |
| License | relicensed from MIT to LGPL-3.0-or-later; earlier releases stay MIT (`NOTICE.md`) | **done** |
| Governance files | `CITATION.cff`, `CONTRIBUTING.md` (with the git workflow), `SECURITY.md`, `CODE_OF_CONDUCT.md`, `codecov.yml`, a README in the BOLTS layout, a generated text logo (a placeholder: replace `docs/_static/logo*.svg`) | **done** |
| Planned capabilities | the `scubas.future` package and `NotImplementedYetError` (stubs that say what is not built, with tests); the package is empty until a capability is named | **done** (the pattern; no capability is committed) |
| Packaging hygiene | the package metadata corrected (description, classifiers, keywords, URLs, SPDX license expression with setuptools 77), numpy declared, the `plots`, `data`, `dev` and `docs` extras next to the unchanged runtime dependencies, an explicit `MANIFEST.in`, a small `environment.yml` (`env.yml` stays), `scubas.__version__`, `scripts/check_dist.py` | **done** |
| Release control | a tag-driven release workflow (TestPyPI, then PyPI after a manual approval, trusted publishing) behind an on/off repository variable, `make publish-*` switches, and a release guide with the owner's checklist | **done** (built, and inert until the owner creates the environments, the variable and the trusted publishers) |
| Provenance for new code | `Param(value, unit, status, source)` and a clear error for unknown values, used by new code only; a read-only list of the 35 constants of the existing code with a test that they match the code | **done** |
| Per-change discipline | pull request checklist, changelog guard, lint debt and coverage floors that only ratchet in one direction | **done** (ongoing) |

Not planned for now: a `src/` layout, splitting the modules, replacing loguru, making matplotlib or pyproj optional. Each would change
how users import or install the package; the reasons are in the decision records (`docs/dev/adr/`). If something ever moves, its old
import path stays through a shim that warns, with the version that removes it (not before 1.0.0).

## Planned capabilities

None is committed. A capability the maintainers decide to build is first added here and as a stub in `scubas.future` (see `docs/dev/future.md`), in the same pull request.

## Known defects of 0.1.9

Pinned as they are, not fixed (see `docs/dev/testing.md`). Each fix is a separate, approved change.

## Open decisions

The date of the next release (0.2.0 is the agreed version) and who tags it, the author list for `CITATION.cff`, and which
known defects to fix and how (a new function or option when the output changes). Changes to the model itself are decided by the
maintainers and the published studies; none is committed here.

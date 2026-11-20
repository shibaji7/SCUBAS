# Contributing to SCUBAS

Thank you for helping. SCUBAS is a small, careful scientific package that is on PyPI and cited in published studies, so the first rule is
**compatibility**: every public import, function signature, file format and number of 0.1.9 is pinned by `tests/compat/`, and a change is accepted only if
those tests stay unchanged and green (see [docs/dev/testing.md](docs/dev/testing.md) and the decision record `docs/dev/adr/0001-legacy-numerics-are-pinned.md`).

## Set up

```bash
git clone https://github.com/shibaji7/SCUBAS && cd SCUBAS
make env && conda activate scubas   # conda environment from environment.yml (or any Python >= 3.11)
make dev                            # editable install with pytest, coverage, ruff, pre-commit, build, twine, build, twine
make pre-commit-install                                        # the cheap gates on every commit
make qa                                                        # everything that must pass before a commit
```

`make help` lists every target. `make docs` builds the documentation (Sphinx, warnings are errors) and needs the tools of `docs/requirements.txt`.

## The git workflow

- **`main` is always green and nobody pushes to it directly.** Every change is a pull request.
- One branch per change, cut from the latest `main`, named for what it does (`pandas3-compat`, `phase-4-docs`). Keep it small and focused; phases of a larger
  plan are separate pull requests merged in order, never stacked unmerged.
- Commits say **why**, not only what. A change under `scubas/` also needs a line under `## [Unreleased]` in `CHANGELOG.md` (the pre-commit hook `changelog` checks it;
  `SKIP=changelog git commit` bypasses it for a change nobody will see).
- Open the pull request against `main` and fill in the template. CI must be green (`tests`, `build`, `docs`, `security`): lint, format, SPDX headers, release consistency,
  the tests on Python 3.11 to 3.13 with pandas 2 and 3, 90% coverage and a floor per module, the wheel installed in a clean environment, a strict docs build, `pip-audit`.
- Merge with **squash**, so `main` has one commit per change, and delete the branch.
- Dependabot opens pull requests for the GitHub Actions; each is merged on its own CI result.
- Never tag, release or publish from a pull request. Releases are made by the maintainer ([docs/dev/release.md](docs/dev/release.md)).

## Principles

1. **Do not change what users get.** A moved or renamed public name keeps its old import path through a thin shim that re-exports the same object and emits a
   `DeprecationWarning` naming the new path and the version that removes it (not before 1.0.0, and not sooner than two minor releases after the warning first ships).
2. **Known defects are pinned and reported, not fixed silently** ([docs/dev/testing.md](docs/dev/testing.md#known-defects)). A fix is its own pull request; one that changes
   an output goes behind a new function or option.
3. **New physics is documented and tested analytically.** A docstring with the equation and its source (paper, book and page), a test against a closed form or a limit
   (not a snapshot), a changelog line and a paragraph in the docs.
4. **New parameters carry a unit, a status and a source** (`sourced`, `assumed`, `unknown`, `override`); an unknown value raises a clear error instead of being guessed.
   The pattern is for new code; the constants of the existing code get their provenance gradually.
5. **Say what is not implemented.** A planned capability is a stub under `scubas/future/` that raises `NotImplementedYetError` naming its roadmap item
   ([docs/dev/future.md](docs/dev/future.md)), so the package and the roadmap cannot drift apart.
6. **No network and no unseeded randomness in tests.**

## Style and checks

The rules are in [docs/dev/conventions.md](docs/dev/conventions.md). `make qa` runs `make lint` (ruff, the legacy debt in `ruff.toml` is allowed and new violations are not),
`make lint-ratchet` (that debt may only shrink), `make fmt-check` (ruff format, width 88), `make spdx`, `make release-check`, `make cov-xml` (coverage and the per-module floors)
and `make docs`. `make format` writes the formatting. Every Python file starts with the SPDX header
`# SPDX-License-Identifier: LGPL-3.0-or-later` and the copyright line (`make spdx` checks it).

## License

By contributing you agree that your work is released under the project license, LGPL-3.0-or-later (`LICENSE`, with the GPL text in `COPYING`).

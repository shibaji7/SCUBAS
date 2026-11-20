# Release

## What exists today

SCUBAS is released to PyPI by hand. There is no release workflow yet, no tag for the 0.1.x versions (the only GitHub tags are the old
`TLM-0.1.3` and `TLM-0.1.5`), and nothing publishes from CI. Treat this page as the checklist; the planned automated flow is at the end.

1. Decide the version and update `version` in `pyproject.toml`.
2. Move the `## [Unreleased]` entries of `CHANGELOG.md` into a `## [X.Y.Z] - date` section and keep an empty `## [Unreleased]`.
3. `make release-check` (add `TAG=vX.Y.Z` to compare a tag): the version, `CHANGELOG.md` and, once it exists, `CITATION.cff` agree.
4. `make clean && make qa`, then `python -m build` and `python -m twine check dist/*`.
5. Upload with `twine` (a PyPI token for the project) after checking `dist/`. This is the one step that is public and cannot be
   taken back: a version number on PyPI is never reused.

`configure.py -upip` (clean, format with black, build, upload) is the old way; it formats the whole tree before building, so prefer the steps above.

## The gates that run on every pull request

| Workflow | What it checks |
|---|---|
| `tests` | ruff (legacy debt allowed, new violations not), the debt ratchet, format, SPDX, release consistency; the tests on Python 3.11 to 3.13 with pandas 2 and 3; 90% coverage and a floor per module |
| `build` | sdist and wheel build, `twine check`, the wheel installed in a clean environment and the compatibility suite run against it |
| `docs` | the Sphinx build with warnings as errors |
| `security` | `pip-audit` on every push and every Monday |

## Switches

| Switch | Controls | Default | How to change |
|---|---|---|---|
| repository variable `CODECOV_UPLOAD` | the Codecov upload in `tests.yml` | on (anything but `false`) | `gh variable set CODECOV_UPLOAD --body false` |
| file `.rtd-skip` | Read the Docs builds (they are cancelled while it exists; the published site stays) | absent: builds run | `make rtd-off` creates it, `make rtd-on` removes it, `make rtd-status` shows which; commit and push |

Setting a repository variable, a secret or an environment is the owner's action.

## Documentation

Read the Docs builds `latest` from `main` with Sphinx and the `dirhtml` builder (`.readthedocs.yaml`, `fail_on_warning: true`), so the
pages keep the addresses they had under mkdocs. Before pushing: `make docs` builds the site with warnings as errors into
`docs/_build/dirhtml` and `make docs-serve` serves it on <http://localhost:8000>. The tools are pinned in `docs/requirements.txt`.

## Planned: a tag-driven release

Not built yet; it needs actions only the owner can take. The plan: a workflow on a pushed tag `vX.Y.Z` that refuses a tag that is not
the version in the code or has no changelog section; TestPyPI first, then PyPI after a manual approval of a `pypi` environment, using
PyPI trusted publishing (no stored token); a GitHub Release from the changelog section; and an on/off repository variable
(`PUBLISH_ENABLED`) so that nothing is published by accident. Creating the environments, the variable and the trusted publisher on PyPI
and TestPyPI is done by the owner.

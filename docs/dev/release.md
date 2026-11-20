# Release

A release is a **tag**. Pushing `vX.Y.Z` starts `.github/workflows/release.yml`, which checks the tag, builds, publishes to TestPyPI, waits for a manual
approval, publishes to PyPI with trusted publishing (no token is stored anywhere) and makes a GitHub Release from the changelog. The workflow is
**switched off until the owner turns it on**, and it needs a few one-time actions that only the owner can take. Nothing is ever published from a pull request.

## One-time setup (the owner)

None of this is done by the code; each step is a setting on pypi.org, test.pypi.org or GitHub.

1. **PyPI trusted publisher.** pypi.org, project `scubas`, Publishing, "Add a new publisher", GitHub: owner `shibaji7`, repository `SCUBAS`, workflow `release.yml`,
   environment `pypi`.
2. **TestPyPI trusted publisher.** test.pypi.org: the project does not exist there yet, so add a *pending* publisher (Publishing, "Add a new pending publisher") with the
   project name `scubas`, the same owner, repository and workflow, and environment `testpypi`.
3. **GitHub environments** (repository Settings, Environments): create `testpypi` (no reviewers) and `pypi`. On `pypi` add yourself under "Required reviewers": that approval
   is the manual gate before anything reaches PyPI. Optionally restrict both environments to the tag pattern `v*`.
4. **Turn the switch on** when you are ready to publish: `make publish-on` (it sets the repository variable `PUBLISH_ENABLED=true`; or set it in Settings, Secrets and
   variables, Variables). `make publish-off` turns it off again, `make publish-status` shows the value. While it is unset or `false`, a pushed tag is only checked and built.
5. The old `twine` token (`configure.py -upip`) is no longer needed once the trusted publishers exist; revoke it on PyPI.

## Cutting a release

1. A pull request that prepares the release: set `__version__` in `scubas/_version.py`, set `version` and `date-released` in `CITATION.cff`, and move the `## [Unreleased]` entries
   of `CHANGELOG.md` into a `## [X.Y.Z] - date` section (keep an empty `## [Unreleased]`). `make release-check TAG=vX.Y.Z` says whether they all agree.
2. `make clean && make qa && make check` (builds the sdist and wheel, runs `twine check` and `scripts/check_dist.py`). Merge the pull request after a green CI.
3. The owner tags the merge commit on `main` and pushes the tag: `git tag vX.Y.Z && git push origin vX.Y.Z`. The tag must equal the version in the code and have a changelog
   section; the workflow refuses it otherwise. (An assistant never tags, releases or approves an environment.)
4. Watch the `release` workflow: `build` always runs; with `PUBLISH_ENABLED=true` it continues to TestPyPI, then waits for your approval of the `pypi` environment, then publishes to
   PyPI and creates the GitHub Release (the changelog section as its notes, the files attached).
5. Check the PyPI page (description, license, classifiers, version) and `pip install scubas==X.Y.Z` in a clean environment.

A version number on PyPI is never reused: if something is wrong after the upload, fix it and release the next version. `skip-existing` is set for TestPyPI only, so re-running the
workflow on the same tag does not fail there.

**A tag pushed while publishing was off** is built and checked but not uploaded. After `make publish-on`, upload exactly that version with `make publish-tag TAG=vX.Y.Z` (this
starts the workflow by hand on the tag; `pypi` still waits for the approval).

**The old versions.** 0.1.6 to 0.1.9 were uploaded by hand and have no git tags (only `TLM-0.1.3` and `TLM-0.1.5` exist); they are not back-filled, because a tag would be a
guess at the commit. The tag history starts with the first release made this way.

## The gates that run on every pull request

| Workflow | What it checks |
|---|---|
| `tests` | ruff (legacy debt allowed, new violations not), the debt ratchet, format, SPDX, release consistency; the tests on Python 3.11 to 3.13 with pandas 2 and 3; 90% coverage and a floor per module |
| `build` | sdist and wheel build, `twine check`, `scripts/check_dist.py`, the wheel installed in a clean environment and the compatibility suite run against it |
| `docs` | the Sphinx build with warnings as errors |
| `security` | `pip-audit` on every push and every Monday |

`main` is protected: a pull request is required and the checks above (the ten that run per commit) must pass before it merges.

## Switches

| Switch | Controls | Default | How to change |
|---|---|---|---|
| repository variable `PUBLISH_ENABLED` | the TestPyPI, PyPI and GitHub Release jobs of `release.yml` | off (unset or anything but `true`) | `make publish-on`, `make publish-off`, `make publish-status` |
| repository variable `CODECOV_UPLOAD` | the Codecov upload in `tests.yml` | on (anything but `false`) | `make codecov-off`, `make codecov-on`, `make codecov-status` |
| file `.rtd-skip` | Read the Docs builds (they are cancelled while it exists; the published site stays) | absent: builds run | `make rtd-off` creates it, `make rtd-on` removes it, `make rtd-status` shows which; commit and push |

The `make publish-*` and `make codecov-*` targets change the repository's settings through `gh`; setting a variable, a secret or an environment is the owner's action.

## Documentation

Read the Docs builds `latest` from `main` with Sphinx and the `dirhtml` builder (`.readthedocs.yaml`, `fail_on_warning: true`), so the
pages keep the addresses they had under mkdocs. Before pushing: `make docs` builds the site with warnings as errors into
`docs/_build/dirhtml` and `make docs-serve` serves it on <http://localhost:8000>. The tools are pinned in `docs/requirements.txt`.

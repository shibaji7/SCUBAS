# ADR 0009: a release is a tag, published with trusted publishing behind a switch

**Status:** accepted (2026-11).

**Context.** SCUBAS 0.1.1 to 0.1.9 were uploaded to PyPI by hand with a stored token (`configure.py -upip`), without tags, so a version cannot be matched to a commit. The repository is
public and `main` is protected; a release is the one public, irreversible step (a version number on PyPI is never reused).

**Decision.** A release is a pushed tag `vX.Y.Z`. `.github/workflows/release.yml` refuses a tag that is not the version in `scubas/_version.py` or has no changelog section, builds and
checks the distributions, then publishes to TestPyPI, waits for a manual approval of the `pypi` environment, publishes to PyPI with PyPI trusted publishing (OIDC: no token is stored) and
makes a GitHub Release from the changelog. The publishing jobs run only while the repository variable `PUBLISH_ENABLED` is `true` (`make publish-on` and `publish-off`), so the workflow is
inert until the owner has created the environments, the variable and the trusted publishers. The old versions are not tagged after the fact. The owner pushes the tag; an assistant never tags,
releases or approves an environment.

**Consequences.** The first release needs a few settings on pypi.org, test.pypi.org and GitHub (`docs/dev/release.md`). `configure.py -upip` and its token are no longer the way to publish.

**Rejected.** Publishing on every merge to `main` (no deliberate release step); a stored API token (a secret that can leak); tagging 0.1.6 to 0.1.9 retroactively (a guess at the commit).

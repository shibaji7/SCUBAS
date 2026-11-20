# ADR 0002: the documentation stays on mkdocs-material

**Status:** superseded by [ADR 0007](0007-documentation-on-sphinx.md) (2026-11): the owner then chose to move to Sphinx like BOLTS.

**Context.** BOLTS documents itself with Sphinx and MyST. SCUBAS has eight tutorials written for mkdocs-material (admonitions,
mkdocstrings directives, HTML cards, MathJax), and its Read the Docs URLs are cited in papers.

**Decision.** Keep mkdocs-material. Everything the BOLTS documentation has is added on this toolchain: a development-status banner on
every page (one hook, `docs/hooks.py`), a glossary, these records, a release guide, roadmap and changelog pages, a strict build in CI
and on Read the Docs. Every existing URL keeps working.

**Consequences.** The API pages stay mkdocstrings pages. A move to Sphinx later is possible but would be a deliberate, separate step
with a URL redirect plan.

**Rejected.** Moving to Sphinx now (about two thousand lines of markup to rewrite and URLs that papers cite would change, for no new capability).

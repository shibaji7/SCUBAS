# ADR 0007: the documentation moves to Sphinx and MyST, built as `dirhtml`

**Status:** accepted (2026-11). Supersedes [ADR 0002](0002-stay-on-mkdocs-material.md).

**Context.** The first plan kept mkdocs-material because the pages and their addresses are cited in papers. The owner then asked to
move Read the Docs to Sphinx, as the sibling project BOLTS does, so that the two projects share one documentation stack.

**Decision.** Sphinx with MyST Markdown, the Read the Docs theme, autodoc and napoleon (the package uses NumPy-style docstrings),
viewcode, intersphinx and a copy button; the pages stay Markdown. The site is built with the **`dirhtml` builder**, which writes
`tutorial/theory/index.html` rather than `tutorial/theory.html`, so every address of the mkdocs site is unchanged and nothing needs
redirecting. The development-status banner is added to every page by one `source-read` handler in `docs/conf.py`. The build is strict
(`-W`) on CI and on Read the Docs; the tools are pinned in `docs/requirements.txt`.

**Consequences.** The mkdocs-only parts are gone (`mkdocs.yml`, the card-grid CSS, the lightbox and the MathJax shim); the API pages use
autodoc directives instead of mkdocstrings blocks; a page needs a title to be in a table of contents. The search index and the look change.

**Rejected.** Sphinx's default `html` builder (every page address would change); redirects (extra moving parts for no gain).

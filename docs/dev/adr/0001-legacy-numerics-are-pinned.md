# ADR 0001: the numbers, API and file formats of 0.1.9 are pinned, not rewritten

**Status:** accepted (2026-11).

**Context.** SCUBAS 0.1.9 is on PyPI and cited in published studies. Its science code works, but it has little structure around it (57%
coverage, no gates, a red CI). The goal is the engineering standard of the sibling project BOLTS without changing what users get: old
scripts and notebooks must keep running and old results must be reproducible.

**Decision.** Before anything else changes, `tests/compat/` records what 0.1.9 does: every public name and signature, the numbers of
the tutorial cables and the 1989 event, the file formats, and the closed forms of the theory. A change is accepted only if those tests
stay unchanged and green. New standards (lint, headers, docs, provenance) apply to new and touched code; working science code is not
rewritten to satisfy them, and existing lint violations are frozen in a list that may only shrink.

**Consequences.** Refactoring is cheap to check but slow to start. Legacy files carry a debt list in `ruff.toml` and two of them are
excluded from `ruff format`. A baseline may be regenerated only on unchanged legacy code or for an explicitly approved change.

**Rejected.** Rewriting the modules to the BOLTS layout first (no way to prove the results did not move); adding lint and type gates
that fail on legacy code (CI red for the wrong reason, so people stop reading it).

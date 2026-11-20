## What changes
<!-- One or two sentences. Say if package code (scubas/) changes. -->

## Checklist
- [ ] `make qa` passes
- [ ] `tests/compat/` is unchanged and green (the public API, file formats and numbers of 0.1.9 are pinned). If a baseline or a known-defect test had to change, the change was approved and is explained below
- [ ] CHANGELOG.md updated under `[Unreleased]` (required when `scubas/` changes)
- [ ] New or changed physics has a docstring with the equation and its source, and a test with an analytic expectation
- [ ] Docs paragraph updated when behaviour changes
- [ ] A moved or renamed public name keeps its old import path (a shim that warns, with the version that removes it)

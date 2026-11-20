# ADR 0005: known defects are pinned and reported, never fixed silently

**Status:** accepted (2026-11).

**Context.** Reading the legacy code for the compatibility suite turned up behaviour that looks wrong: an uncallable
`Direction.enumerate_key`, an input array rescaled in place, end potentials rounded to 0.01 V inside the physics path that feeds
`Vt(v)`, thicknesses documented as km but stored in metres. Published results may depend on some of them.

**Decision.** Each defect is pinned by a test in `tests/compat/test_known_defects.py` exactly as it behaves, and listed in
[testing](../testing.md#known-defects). A fix is a separate pull request, approved one by one. A fix that changes an output goes behind
a new function or option and changes the pinned test and the changelog in the same pull request; it never changes a result silently.
A defect that appears only with a newer dependency is different: the pandas 3 incompatibility of `read_iaga` was fixed because it
broke a fresh install, with identical results on pandas 2.

**Consequences.** The test suite documents wrong behaviour on purpose. Each pinned test says what is wrong.

**Rejected.** Fixing everything found (the published numbers would move without notice); leaving the defects unrecorded (they would
be rediscovered and fixed by accident).

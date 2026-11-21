# ADR 0008: a planned capability is a stub that says so

**Status:** accepted (2026-11).

**Context.** BOLTS keeps its planned capabilities as modules that raise `NotImplementedYetError` under `bolts.future` (its ADR 0005), so the tree and the roadmap stay in
step. SCUBAS has no list of planned capabilities in the repository: apart from the unused `length_method` option of `CableSection.compute_lengths`, the code names nothing as
planned.

**Decision.** `scubas.future` exists with the pattern (`NotImplementedYetError` in `scubas/errors.py`, the rules in [Planned capabilities](../future.md), tests that check the
mechanism) and is **empty** until the maintainers name a capability. No capability is invented to fill it.

**Consequences.** A new package and a new error class are added (nothing existing changes). The first stub arrives with its roadmap row.

**Rejected.** Seeding stubs from what the tutorials describe as intended (the unconfirmed list would read as a promise); leaving the pattern out (the first planned capability
would then be written as half-working code).

**Update (2026-11).** The maintainers named four capabilities (uncertainty quantification, an automatic segmenter, an auroral electrojet design, FDTD/FEM for coastal effects). Each
arrived as a stub with its roadmap row and a plan page in `docs/plans/`, as this record asks; the decision itself is unchanged. The package is no longer empty, and the tests now
tie the stubs, the roadmap table and the plan pages together.

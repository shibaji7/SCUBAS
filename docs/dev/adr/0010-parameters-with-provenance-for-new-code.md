# ADR 0010: parameters carry provenance, for new code only

**Status:** accepted (2026-11).

**Context.** BOLTS records every parameter as `value, unit, status, source` and refuses to guess an unknown one (its ADR 0002). SCUBAS 0.1.9 has its constants written into the code: layer
resistivities, a taper fraction, a 120-sample baseline, the geodesic constants of `CableSection.compute_lengths`, the Earth radius, the defaults of the line-current source. For most of them the
code names no source. Rewiring the code to read them from files would change legacy modules, which ADR 0001 forbids without an approved, output-preserving change.

**Decision.** `scubas.provenance` provides the pattern (`Param`, `Status`, `ParameterSet`, `load_parameters`, `MissingParameterError`) for **new** code. The existing constants are
**listed, not rewired**: `scubas/defaults/legacy_constants.toml` records each with its unit, a status and a source, and `tests/test_provenance.py` fails if a listed value differs from the
value the code uses. A status is honest: `sourced` only where a document can be named, otherwise `assumed` with a source that says the origin is not documented; a citation is never invented
(a test refuses a citation-like year on an `assumed` entry and a `sourced` entry with no named document).

**Consequences.** New physics gets its numbers from provenance-carrying files. The list shows at a glance which legacy numbers are unsupported (25 of 35 are `assumed`), which is the first step to
replacing them. The unit is a label for now: nothing is converted. A later, separate change can move a constant out of the code, with the compatibility suite unchanged.

**Rejected.** Reading the legacy constants from the file now (changes legacy modules); leaving them undocumented (no record of what is a placeholder); marking everything `sourced` to look complete.

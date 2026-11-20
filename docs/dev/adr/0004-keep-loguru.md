# ADR 0004: keep loguru for logging

**Status:** accepted (2026-11).

**Context.** BOLTS logs through the standard library with a `NullHandler`. SCUBAS logs through `loguru` in five modules, and users may
have configured loguru sinks around it.

**Decision.** Keep loguru in existing code and use it in new code of the same modules, so there is one logging system in the package.

**Consequences.** The package is not a "silent library" in the standard-library sense; loguru is a runtime dependency. The compat
tests disable the `scubas` logger while they run.

**Rejected.** Switching to `logging` now (silently changes what users see and break sinks they configured).

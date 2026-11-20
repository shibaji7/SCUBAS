# ADR 0003: keep the flat package layout and the runtime dependencies in 0.x

**Status:** accepted (2026-11).

**Context.** BOLTS uses a `src/` layout and a runtime dependency list of numpy and scipy. SCUBAS imports `matplotlib` in `scubas.plotlib`
and `pyproj` in `scubas.conductivity` at module level, and users call both after a plain `pip install scubas`.

**Decision.** Keep the flat `scubas/` layout and every runtime dependency in 0.x. The benefit of a `src/` layout (testing the installed
package, not the checkout) is obtained in CI: `build.yml` installs the wheel in a clean environment and runs the compatibility suite
against it. Add extras (`plots`, `data`, `dev`, `docs`) and a small `environment.yml` only additively; `env.yml` stays as the
compatibility file.

**Consequences.** `pip install scubas` stays as heavy as before. Slimming the dependencies is a decision for 1.0 with lazy imports and a
clear `ImportError`.

**Rejected.** Moving to `src/` now (changes how a cloned repository is imported); making matplotlib and pyproj optional now (breaks
`import scubas.plotlib` after an install that worked before).

# Testing and the compatibility suite

`make test` runs everything; `make qa` runs every gate before a commit (lint, format, SPDX, release check, tests with coverage,
strict docs build). The same checks run in GitHub Actions (`tests`, `build`, `docs`, `security`).

## The compatibility suite (`tests/compat/`)

SCUBAS 0.1.9 is on PyPI and cited in papers, so its public API, file formats and numbers are pinned *before* anything else changes.
A change is accepted only if these tests stay unchanged and green.

| File | What it pins |
|---|---|
| `test_api_surface.py` | every module, class, method signature and the 19 `PROFILES` (snapshot in `baseline/api_surface.json`); new names are allowed, removals and changes fail |
| `test_behaviour.py` | numbers: the three conductivity tutorial cables, the March 1989 TAT-8 event, transfer functions, IAGA and CSV readers, the current source, FFT and distance helpers, the LITHO1.0 profile code on a synthetic grid (`baseline/behaviour.npz`, relative tolerance 1e-9) |
| `test_file_formats.py` | IAGA-2002 and CSV B-field input, gap handling, the two files `Cable.save` writes |
| `test_analytic.py` | closed forms from the tutorials: long and short cable end potentials, half-space impedance, equivalent-pi branch |
| `test_known_defects.py` | the known defects, below |
| `test_import_paths.py` | the import paths of the docs and examples, and that `import scubas` stays light and writes nothing |
| `test_litho_loader.py` | the LITHO1.0 loader without the network (a synthetic netCDF, `urlopen` mocked) |
| `test_examples_smoke.py` | the example scripts that make the tutorial figures still run (from a temporary copy) |

Regenerate the baselines **only** on the unchanged legacy code or for an approved change of a number or signature:
`make compat-baseline`, then review the diff and explain it in `CHANGELOG.md`. A real `LITHO1.0.nc` is not needed; a test that uses
one is marked `litho` and runs with `SCUBAS_LITHO1_PATH=/path/to/LITHO1.0.nc pytest -m litho`.

The suite runs on Python 3.11 to 3.13 with pandas 2 and pandas 3, and `build.yml` runs it against the installed wheel
(`SCUBAS_TEST_INSTALLED=1`, so it does not test the checkout by accident).

## Known defects

Behaviour of 0.1.9 that looks wrong but that users or published results may depend on is **recorded, not fixed**. Each is pinned in
`test_known_defects.py` as it is:

- `Direction.enumerate_key` cannot be called (it compares the class with a member, so no branch matches).
- `Current.compute_B_field` rescales the caller's array in place.
- `Cable._pot_end_cable_` and the other end-potential helpers round to 0.01 V, and `Vt(v)` is computed from the rounded values.
- The `unit` argument of those helpers is `"V"` or anything else, which means millivolts.
- Layer thicknesses are in metres although the `Site` docstrings say km.
- Two guards (`add_active_termination` before the E-field exists, `_pot_alongCS_` before the nodal solution) behave differently from what their messages suggest.

A fix is a separate pull request that you approve one by one; an output-changing fix goes behind a new function or option, never
silently, and updates the one test here and the changelog.

## Tests you add

A new or changed piece of physics gets an analytic test (a closed form, a limit, a symmetry such as "splitting a uniform cable
does not change the end potentials"), not a snapshot. Tests never use the network or random numbers without a fixed seed.

## Coverage

CI fails below 90% overall and below a floor per module (`[tool.scubas.coverage-floors]` in `pyproject.toml`). A floor only goes up.

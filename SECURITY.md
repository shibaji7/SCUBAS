# Security policy

SCUBAS is a scientific library: it reads geomagnetic records (IAGA-2002 text files and CSV) and cable descriptions, computes induced voltages, and writes
result files. It runs no code from its inputs. The realistic risks are a malformed input file crashing a run, a vulnerable dependency (numpy, scipy, pandas,
matplotlib, pyproj), and one network access described below.

## What touches the network

`scubas.conductivity.ConductivityProfile` downloads the LITHO1.0 Earth model (a large netCDF file) **the first time it is used**, from the
address in its `uri` parameter (by default `http://ds.iris.edu/files/products/emc/emc-files/LITHO1.0.nc`, plain HTTP) into `.scubas_config/` in the current
working directory. Nothing checks the integrity of the file. If that matters to you, put a LITHO1.0 file you trust in `.scubas_config/LITHO1.0.nc` before the first
use (nothing is downloaded when it exists) or pass your own `uri`. This is a known limitation, recorded in `docs/dev/testing.md`; changing it is a separate change.

## Supported versions

The latest release on PyPI gets fixes. SCUBAS is in alpha (0.x): older 0.x releases are not patched.

## Reporting a vulnerability

Please do not open a public issue. Use GitHub's private vulnerability reporting for this repository (Security tab, "Report a vulnerability"). Expect an
acknowledgement within a week. Dependencies are audited with `pip-audit` on every push and every week (`.github/workflows/security.yml`), and Dependabot proposes
updates to the GitHub Actions.

## Not a security issue

A wrong physical result is a bug, not a vulnerability: open a normal issue with the cable sections and the record that reproduce it.

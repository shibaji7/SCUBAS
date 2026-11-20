<!--
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is licensed under the GNU Lesser General Public License v3.0 or later
(LGPL-3.0-or-later); see LICENSE and COPYING in the root directory. Versions up to
and including 0.1.9 were released under the MIT license (see NOTICE.md).
-->

[![License: LGPL v3+](https://img.shields.io/badge/license-LGPL--3.0--or--later-blue.svg)](https://github.com/shibaji7/SCUBAS/blob/main/LICENSE) 
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3110/) 
![GitHub Stable Release (latest by date)](https://img.shields.io/github/v/release/shibaji7/SCUBAS)
[![Documentation Status](https://img.shields.io/readthedocs/SCUBAS?logo=readthedocs&label=docs)](https://SCUBAS.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/shibaji7/SCUBAS/branch/main/graph/badge.svg)](https://codecov.io/gh/shibaji7/SCUBAS)

# SCUBAS

SCUBAS: **S**ubmarine **C**ables **U**pset **b**y **A**uroral **S**treams.

**Submarine cable geomagnetic disturbance modeling.** SCUBAS provides a physics-based workflow for estimating the voltage induced on submarine cables during geomagnetic storm conditions.

SCUBAS is an open source Python-based computational model of geomagnetic induction on submarine cables. The model is used to estimate the induced voltage in the submarine cables in response to geomagnetic disturbances. It utilizes newly acquired knowledge from magnetotelluric studies and associated investigations of geomagnetically induced currents in power systems.

## Source Code

The library source code can be found on the [SCUBAS GitHub](https://github.com/shibaji7/SCUBAS) repository.

If you have any questions or concerns please submit an **Issue** on the [SCUBAS GitHub](https://github.com/shibaji7/SCUBAS) repository.

## Documentation

| Section | What it covers |
|---|---|
| [Installation](user/install.md) | environment setup, dependencies, and package installation |
| [INTERMAGNET data access](user/intermagnet.md), [Citing](user/citing.md) | where the magnetic field records come from, and citation guidance |
| [Tutorials](tutorial/theory.md) | theory, long and short cables, active terminations, network modeling, case studies, the 1989 event, extreme value analysis, uncertainty |
| [API reference](dev/index.md) | cables, conductivity, datasets, models, plotting, current sources and utilities |
| [Project](roadmap.md) | the roadmap, the changelog and the glossary |
| [Developer guide](dev/conventions.md) | conventions, testing and the compatibility suite, planned capabilities, release, decision records |

```{toctree}
:hidden:
:maxdepth: 2
:caption: Start here

user/install
user/intermagnet
user/citing
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Tutorials

tutorial/theory
tutorial/elsc
tutorial/active
tutorial/netmodel
tutorial/conduct
tutorial/1989
tutorial/eva
tutorial/uq
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: API

dev/index
dev/cables
dev/conductivity
dev/datasets
dev/models
dev/plotlib
dev/sources
dev/utils
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Project

roadmap
changelog
glossary
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Developer guide

dev/conventions
dev/testing
dev/future
dev/release
dev/adr/index
```

<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->

[![License: MIT](https://img.shields.io/badge/License%3A-MIT-green)](https://choosealicense.com/licenses/mit/) 
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3110/) 
![GitHub Stable Release (latest by date)](https://img.shields.io/github/v/release/shibaji7/SCUBAS)
[![Documentation Status](https://img.shields.io/readthedocs/SCUBAS?logo=readthedocs&label=docs)](https://SCUBAS.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/shibaji7/SCUBAS/branch/main/graph/badge.svg)](https://codecov.io/gh/shibaji7/SCUBAS)

# SCUBAS
SUCBAS: **S**ubmarine **C**ables **U**pset **b**y **A**uroral **S**treams.

SCUBAS is an open source Python-based computational model of geomagnetic induction on submarine cables. The model is used to estimate the induced voltage in the submarine cables in response to geomagnetic disturbances. It utilizes newly acquired knowledge from magnetotelluric studies and associated investigations of geomagnetically induced currents in power systems.

## Source Code 

The library source code can be found on the [SCUBAS GitHub](https://github.com/shibaji7/SCUBAS) repository. 

If you have any questions or concerns please submit an **Issue** on the [SCUBAS GitHub](https://github.com/shibaji7/SCUBAS) repository. 

## Table of Contents 
  - [Installation](user/install.md)
  - [INTERMAGNET Data Access](user/intermagnet.md)
  - [Citing](user/citing.md)
  - Tutorials
    - [TLM Theory](tutorial/theory.md)
    - [Long and Short Cables](tutorial/elsc.md)
    - [Active Termination](tutorial/active.md)
    - [Network Modeling](tutorial/netmodel.md)
    - [Case Studies/Codes](tutorial/conduct.md)
    - [Event Study (1989 Strom)](tutorial/1989.md)
    - [Extreme Value Analysis](tutorial/eva.md)
    - [Uncertainity Quantification](tutorial/uq.md)
  - Code Documentation:
    - [Cables](dev/cables.md)
    - [Ocean/Earth Conductivity](dev/conductivity.md)
    - [Datasets](dev/datasets.md)
    - [Ocean/Earth Model](dev/models.md)
    - [Plotting Library](dev/plotlib.md)
    - [Utility Module](dev/utils.md)

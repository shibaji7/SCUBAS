<!--
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is licensed under the GNU Lesser General Public License v3.0 or later
(LGPL-3.0-or-later); see LICENSE and COPYING in the root directory. Versions up to
and including 0.1.9 were released under the MIT license (see NOTICE.md).
-->

# The models

SCUBAS estimates the voltage that a geomagnetic disturbance induces on the power feed of a submarine cable. These pages explain **how the model works**: what it assumes, which equations it solves, where each equation is in the code, and how it was checked. The [tutorials](../tutorial/theory.md) show how to use it, and the [API reference](../dev/index.md) documents the classes.

| Model | Status | What it is |
|---|---|---|
| [The SCUBAS model](scubas.md) | implemented | The whole chain: the magnetic record, the transfer function to the seafloor electric field, the emf in the cable, and the cable voltage $V_C = \mathcal{E}_C + U_A - U_B$ |
| [The DSTL model](dstl.md) | implemented | The distributed-source transmission-line model of the ocean and the Earth: line constants, equivalent-π circuits, active terminations and the nodal network that give the Earth potentials $U_A$ and $U_B$ at the ends of the cable |
| [Auroral electrojet currents](electrojet.md) | planned | A modelled ionospheric source in place of observatory records; what `scubas.sources.Current` does today and what the literature offers |
| [Numerical coast models (FDTD, FEM)](fdtd.md) | planned | A two- or three-dimensional solution of the coast effect, to test where the transmission-line model stops |

The two **planned** models are not in SCUBAS; their pages say what is known, in the order of the work plans ([plans](../plans/index.md)). The two **implemented** models are one calculation in two parts, and the [SCUBAS page](scubas.md) is the place to start.

## How the pages are kept correct

* Every equation is attributed to a paper and an equation number, and is written out again, independently of the package, in `scripts/model_reference.py`. `tests/test_docs_model.py` compares those formulas with the code, and compares the code with the numbers printed in the papers (tables, and values read off figures where there is no table).
* The numbers in the text and the tables of the pages come from `docs/model/figures/numbers.json`, which the same test keeps equal to what the code computes; the figures are regenerated from the code by `scripts/make_model_figures.py` (`make model-figures`).
* Where the code differs from a paper, the page says so ([paper and code](scubas.md#paper-and-code)); nothing is smoothed over, and nothing is changed: the behaviour of 0.1.9 is pinned ([testing](../dev/testing.md)).
* Each page ends with what was read in full, what was checked against the code, and what is only cited.

## Conventions

* Axes: $x$ north, $y$ east, $z$ down (geomagnetic convention). Time dependence $e^{i\omega t}$. Displacement currents are neglected.
* Fields: $B$ in nT and $E$ in mV/km (V/km for the potentials of the examples); the transfer function $T_x = E/B$ in mV/km/nT. The code works internally in SI and the papers quote line constants per kilometre, a factor of 1000.
* Line constants are for a strip of width $w = 1$ m, as in the papers.

| Symbol | Meaning | In the code |
|---|---|---|
| $E_N$, $E_E$ | north and east electric field on the seafloor | `Efield` columns `X`, `Y` of an `OceanModel` or `TransmissionLine` |
| $B_N$, $B_E$ | north and east magnetic field at the surface | `Bfield` columns `X`, `Y` |
| $T_x$ | transfer function, seafloor $E$ over surface $B$ | `OceanModel.calcTF`, key `E2B` |
| $d$, $\rho$ or $\sigma$ | thickness and resistivity or conductivity of a layer | `Layer.thickness`, `Layer.resistivity`, `Layer.conductivity` |
| $Z_d$ | surface impedance of the Earth below the seafloor | `Site.calcZ` |
| $L$ | length of a section | `CableSection.length` |
| $\mathcal{E}_j$, $\mathcal{E}_C$ | emf of a section, of the cable | `TransmissionLine.V`, column `V(v)` of `Cable.tot_params` |
| $Z$, $Y$ | series impedance and parallel admittance per unit length of a section | `TransmissionLine.Z`, `.Y` |
| $\gamma$, $Z_0$, $1/\gamma$ | propagation constant, characteristic impedance, adjustment distance | `.gma`, `.Z0` |
| $Y_E$, $Y'/2$, $I_E$ | series admittance, shunt admittance and current source of the equivalent-π circuit | `.Ye`, `.Yp2`, `.Ie` |
| $J_N$, $Y_N$ | Norton circuit of the land at an end | `term_params[...].Jn`, `.Yn` |
| $U_A$, $U_B$ | Earth potentials at the two ends of the cable | columns `U0`, `U1` of `Cable.tot_params` |
| $V_C$ | voltage on the power feed | column `Vt(v)` |

```{toctree}
:hidden:
:maxdepth: 1

scubas
dstl
electrojet
fdtd
```

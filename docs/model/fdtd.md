<!--
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is licensed under the GNU Lesser General Public License v3.0 or later
(LGPL-3.0-or-later); see LICENSE and COPYING in the root directory. Versions up to
and including 0.1.9 were released under the MIT license (see NOTICE.md).
-->

# Numerical models of the coast effect (FDTD and FEM)

**Status: planned.** SCUBAS solves nothing in two or three dimensions: the coast effect enters only through the end potentials of the [DSTL model](dstl.md). This page says where that model stops, what a numerical solution (finite-difference time-domain, FDTD, or finite element, FEM) would add, and which results such a solution must reproduce. The work plan and its open questions are in [the coastal-model plan](../plans/coastal.md); this page does not repeat them.

## The coast effect, and what DSTL gives

The current density in seawater is far greater than on land. At a coast the difference leaves charge on the boundary, which creates a potential and so a potential gradient on both sides; the gradient raises the electric field on the resistive side and lowers it on the conducting side, and bends the currents down into the Earth so that current is continuous across the boundary (Price, 1973, as described by Boteler et al., 2024, Section 2.3). The cable stations at the ends of a cable are on this boundary, and the Earth potentials $U_A$ and $U_B$ partly cancel the electric field in the sea near the coast.

The DSTL model computes those potentials from a one-dimensional description of each section, joined at the boundaries ([DSTL page](dstl.md)). It reproduces, for the example models of Boteler et al. (2023):

* coast potentials of about 48 V for a 600 km shallow sea and 105 V for an 8000 km deep ocean, at 0.3 V/km on the shelf and the land and 0.1 V/km in the deep ocean;
* a coast potential that tends to 77.46 V for a very wide shallow sea or shelf, and a shelf-edge potential of 69.38 V for a very wide shelf next to a deep ocean;
* a coast potential that first rises and then falls with the width of the continental shelf, with the maximum near 250 km.

## Where it stops

The assumptions that bound the DSTL model are the same ones that a numerical model would relax:

* **One dimension in depth for each section; lateral changes only between sections.** Boteler et al. (2023) state that the model is valid for the thin-sheet condition, typically met for longer periods, which penetrate the seawater and the crust; for short periods, with their smaller skin depth, other methods such as finite elements are needed (Discussion). They add that their results can serve as test cases for such software.
* **A uniform field in each section, and land as a semi-infinite line** at the two ends.
* **A segmented 1-D framework "does not explicitly resolve" transverse two-dimensional gradients or localised edge effects at abrupt bathymetric transitions**, and SCUBAS "is not intended to replace" three-dimensional models (Chakraborty et al., 2026, Sections 1 and 2.3). Boteler et al. (2023) and Chakraborty et al. (2022, Section 6) mention finite-element and 3-D models as the improvement.

## What a numerical model solves

At the periods of geomagnetic variations displacement currents are negligible, so the problem is quasi-static: the electromagnetic field diffuses into a conducting Earth (the assumption of Boteler and Pirjola, 2003, and of Dong et al., 2013). Dong et al. (2013) solve it for a line current $J_s$ along $y$, over an Earth that is uniform in $y$ (a two-dimensional problem in $x$ and $z$), with the frequency-domain equation for the magnetic vector potential $A_y$ (their Eq 1):

$$
-\frac{1}{\mu}\nabla^2 A_y + j\omega\sigma A_y = J_s ,
$$

with $\sigma$ the conductivity, zero in the air, and the Galerkin finite-element method on a rectangular mesh. In their case the electric field has only a $y$ component and follows from $A_y$ with the scalar potential zero.

Their results, with a layered Earth and with two different layered Earths side by side, are the nearest published check of what such a model does at a conductivity boundary:

* The electric and magnetic fields agree with exact analytic expressions and with the Fast Hankel Transform for a free-space line current, a uniform Earth and a layered Earth to within 1 %, and the surface impedance to within 0.5 %, when the model is large enough and the elements small enough (1 km near the surface, a model 4000 km wide and about 4000 km deep).
* FEM does not need the plane-wave assumption, and gives the fields below the surface as well as at it.
* For two layered Earths side by side (Quebec and British Columbia), the fields far from the interface tend to those of the single layered model on each side, but near the interface they are not the sum or the mean of the two: the electric field is largest on the more resistive side, and the Quebec-side fields do not recover those of the Quebec model until at least 500 km from the interface.
* Each 2-D simulation takes under 10 minutes on a laptop; a 3-D finite-element model still needs a supercomputer. The authors name the coast effect at the ocean-land interface as future work.

That paper is a 2-D finite-element study of a line current in the frequency domain, **not** an FDTD scheme and not a coast with a sea in it. No paper read for this page applies FDTD to the problem; the choice between FDTD (a time-domain diffusion scheme, which would yield time series directly) and a frequency-domain FEM is one of the open questions of the plan.

## What a numerical model would give SCUBAS

* A check, and where it fails an improvement, of the end potentials $U_A$ and $U_B$ that the [DSTL model](dstl.md) computes with active terminations, in particular for short periods and for narrow or steep shelves.
* The electric field in the sea near the coast, which the 1-D sections cannot resolve.
* A way to say how much of the error of a cable-voltage calculation is the coast, and so to judge the limits quoted by the papers.

## Benchmarks it must reproduce

The test cases are fixed by the sources already read, and several are in the repository today:

* Uniform and layered half-spaces: the surface impedance of the 1-D recursion (`Site.calcZ`), which the tests of the [SCUBAS model](scubas.md) already compare with the formulas of Boteler and Pirjola (2003).
* The two-region problem of the DSTL model: the limiting potentials 77.46 V and 69.38 V and the end potentials of the three example models of Boteler et al. (2023, Figures 10 to 12), all reproduced by `Cable` in `tests/test_docs_model.py`.
* The free-space, uniform and layered cases of Dong et al. (2013), against which that paper validated FEM.
* The vertical-boundary problem of Price (1973) (not read).

A numerical result that disagrees with the DSTL values is not wrong for that reason: the DSTL model is a thin-sheet model, and the comparison is how its range is found.

## Sources and what was checked

* **Read:** Dong et al. (2013) (pages 1689 to 1698); Boteler et al. (2023); Boteler et al. (2024); Chakraborty et al. (2022, 2026); Boteler and Pirjola (2003).
* **Checked against the code by test:** the DSTL limits and the example-model potentials quoted above.
* **Cited, not read:** Price (1973); Wang et al. (2023), which checks the DSTL theory against finite-element solutions (abstract known only from a search result); the three-dimensional and global studies named in [the plan](../plans/coastal.md).

## References

* Boteler, D. H., Chakraborty, S., Shi, X., Hartinger, M. D., and Wang, X. (2023). Transmission line modelling of geomagnetic induction in the ocean/Earth conductivity structure. *International Journal of Geosciences* 14, 767-791. doi:10.4236/ijg.2023.148041
* Boteler, D. H., Chakraborty, S., Shi, X., Hartinger, M. D., and Wang, X. (2024). An examination of geomagnetic induction in submarine cables. *Space Weather* 22, e2023SW003687. doi:10.1029/2023SW003687
* Boteler, D. H., and Pirjola, R. J. (2003). Magnetic and electric fields produced in the sea during geomagnetic disturbances. *Pure and Applied Geophysics* 160, 1695-1716. doi:10.1007/s00024-003-2372-6
* Chakraborty, S., Boteler, D. H., Shi, X., Murphy, B. S., Hartinger, M. D., Wang, X., Lucas, G., and Baker, J. B. H. (2022). Modeling geomagnetic induction in submarine cables. *Frontiers in Physics* 10, 1022475. doi:10.3389/fphy.2022.1022475
* Chakraborty, S., et al. (2026). Validating SCUBAS predictions of geomagnetically induced voltage in submarine cables using legacy superstorm observations. *Space Weather*. doi:10.1029/2025SW004843
* Dong, B., Danskin, D. W., Pirjola, R. J., Boteler, D. H., and Wang, Z. Z. (2013). Evaluating the applicability of the finite element method for modelling of geoelectric fields. *Annales Geophysicae* 31, 1689-1698. doi:10.5194/angeo-31-1689-2013
* Price, A. T. (1973). The theory of geomagnetic induction. *Physics of the Earth and Planetary Interiors* 7, 227-233 (as cited by Boteler et al., 2024; not read).

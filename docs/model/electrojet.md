<!--
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is licensed under the GNU Lesser General Public License v3.0 or later
(LGPL-3.0-or-later); see LICENSE and COPYING in the root directory. Versions up to
and including 0.1.9 were released under the MIT license (see NOTICE.md).
-->

# Auroral electrojet currents

**Status: planned.** A modelled auroral electrojet as the source of the magnetic field is not yet part of the cable model. This page says why it is needed, describes the line-current source that SCUBAS already has (`scubas.sources.Current`, which is **not** connected to `Cable`), and what the literature gives for a more realistic source. The work plan and its open questions are in [the electrojet plan](../plans/electrojet.md); this page does not repeat them.

## Why the model needs a source

Today SCUBAS is driven by magnetic records: each section of the cable takes the record of the nearest land observatory ([SCUBAS model](scubas.md)). Two things limit that.

* **There are no magnetometers on the ocean.** The authors name the missing magnetic data over the oceans and small-scale ionospheric currents as the main limits of the method, and SCUBAS does not interpolate the magnetic structure between observatories (Chakraborty et al., 2026, Sections 6 and 7; Chakraborty et al., 2022, Section 5.1).
* **A cable is long compared with a disturbance.** The TAT-8 cable spans 6300 km, which is 5 hours of local time, so the field at the two ends is often different. On the evening of 13 March 1989 Fredericksburg (west) recorded a positive $B_x$ variation from an eastward electrojet in the afternoon sector while Hartland (east) recorded a negative one from a westward electrojet in the midnight sector; the electric fields at the two ends tend to cancel and the cable voltage is small. At 01:30 UT on 14 March a similar disturbance was seen along the whole route, which indicates a westward electrojet from Europe to North America; the induced field then has one sign along the whole cable and the cable showed the largest voltage of the event (Boteler et al., 2024, Section 4.3 and Figure 10).

A model of the source would supply the field along the route where no observatory exists. Boteler et al. (2024) say that future work will examine induction by modelled ionospheric and magnetospheric current sources. At auroral latitudes the non-uniformity of the source also has to be taken into account in the relation between the surface and seafloor fields (Boteler and Pirjola, 2003, Discussion, citing Pirjola et al., 2000).

## What SCUBAS has today

`scubas.sources.current_sources.Current` models an **infinitely long line current** in the east-west direction, at height $h$ above the ground, with the amplitude $I$ and a width $a$ (a Cauchy-distributed current, which raises the effective height to $a + h$), over a layered Earth `Site`. It works at one time and in the frequency domain. The Earth enters through a single **complex skin depth** at the top layer,

$$
p = \frac{Z}{i\omega\mu_0},
$$

with $Z$ the plane-wave surface impedance of the layers (`Site.calcP`); the induced currents in the Earth then act as an **image of the current at the complex depth $h + 2p$** below the surface. This is the complex image method (Thomson and Weaver, 1975; Boteler and Pirjola, 1998, as described by Pirjola and Viljanen, 1998). The code returns, with $x$ the horizontal distance from the current and $h' = a + h + 2p$:

$$
B_x = \frac{\mu_0 I}{2\pi}\left(\frac{h}{h^2 + x^2} + \frac{h'}{h'^2 + x^2}\right), \qquad
B_z = -\frac{\mu_0 I}{2\pi}\left(\frac{x}{h^2 + x^2} - \frac{x}{h'^2 + x^2}\right),
$$

$$
E_y = -\,i\omega\,\frac{\mu_0 I}{2\pi}\,\ln\frac{\sqrt{h'^2 + x^2}}{\sqrt{h^2 + x^2}} ,
$$

where $h$ stands for $a + h$ in the first terms. The first term of each magnetic component is the primary field of the current and the second is the induced one (`external` and `internal` in the returned dictionaries); $z$ is positive downward.

```{figure} figures/line_current.png
:alt: Two panels against distance from -1000 to 1000 km. Left, the north component of the magnetic field: the primary field peaks at 1818 nT under the current, the induced field is a lower and broader curve peaking near 480 nT, and the total peaks near 2270 nT. Right, the vertical component: the primary field is zero under the current and peaks near 900 nT at 110 km on either side; the total is smaller, with a deeper dip at the centre. Yellow dots lie on the primary curve.
:width: 100%

The field of a line current of 1 MA at 110 km height over the Quebec profile for a 300 s period, from `Current`. The yellow dots are Ampère's law for the primary field, $\mu_0 I h/(2\pi(h^2+x^2))$ and $\mu_0 I x/(2\pi(h^2+x^2))$. The induced field (dashed) is that of the image at depth $h + 2p$; it adds to the horizontal field and reduces the vertical field.
```

Three limits check the formulas (`tests/test_docs_model.py`):

* the primary field is Ampère's law for a line current;
* for a perfectly conducting Earth ($p \to 0$) the image is a mirror image at depth $h$: the total horizontal field is twice the primary one and the vertical field and the electric field at the surface vanish;
* the induced field is that of an image at $h + 2p$ with $p$ from `Site.calcP`.

`Current` is not connected to `Cable`: the cable is driven only by observatory records. The module has [known defects](../dev/testing.md#known-defects) (an enumeration that cannot be called, and a method that rescales the array it is given) that are pinned and not fixed.

## What the literature gives

Pirjola and Viljanen (1998) extend the complex image method from an infinite line current to a more realistic **electrojet of finite length $L$**, a horizontal line current with vertical (field-aligned) currents at its two ends, over a layered Earth. Their results, with the equation numbers of that paper:

* The "U-shaped" current can be replaced by a purely horizontal current distribution (their Figure 2) that gives the same **total** (primary plus induced) magnetic field and the same total horizontal electric field at the surface. The result for the electric field is new there, and it makes the earlier formulation of Thomson and Weaver (1975) directly applicable.
* The complex skin depth is $p = Z/(i\omega\mu_0)$ (their Eq 1) and the image of the horizontal current lies at $z = h + 2p$. A purely image-current treatment of the vertical currents does not work (it requires image charges as well), which is why the equivalence is used.
* In a numerical test (a 1000 km electrojet at 110 km height with 1 MA at a period of 300 s, over layers of 100, 2000 and 50 Ω m with thicknesses 20 and 30 km), the method agrees very closely with an exact calculation (the integral formulas of Häkkinen and Pirjola, 1986) for the three magnetic and the two horizontal electric components, at distances of 100 and 500 km from the electrojet. The exact calculation takes 1 to 2 hours on the grid used; the complex image method takes a few seconds.
* It works well for typical geomagnetic-induction values: periods of 10 to 1000 s, resistivities from 1 to $10^5$ Ω m, source height 100 km and horizontal distances of 0 to 1000 km. The method is valid when the modulus of the complex depth is small compared with the horizontal scale of the field, or very much larger.
* Near the end of the electrojet $B_y$ and $E_x$ become important; an infinitely long line current, which neglects them, may lead to wrong conclusions, for example about geomagnetically induced currents.

Pirjola et al. (2000), cited by Boteler and Pirjola (2003) for the seafloor electric field produced by a two-dimensional ionospheric current, was **not read**; the equations for the seafloor field of a non-uniform source are therefore not given here.

## What the model will be

The plan fixes the questions to settle first (geometry, parameters and their provenance, the response of the Earth, validation and the relation to `Current`). This page records only what the sources above already show about the shape of the model:

* a source described by a few parameters (position, height, width or length, strength, time dependence), each a [parameter with a status](../dev/parameters.md);
* its field at the Earth's surface, by the complex image method, for each section of the cable;
* the seafloor field from the surface field by the transfer function of the [SCUBAS model](scubas.md), which assumes a plane wave, or by a relation that accounts for the non-uniform source (the question that Pirjola et al. (2000) address);
* the existing `Cable` and the [DSTL model](dstl.md) unchanged downstream: the fields become the per-section electric fields.

The first validation case is the 1989 event with the electrojet positions of Boteler et al. (2024, Figure 10).

## Sources and what was checked

* **Read:** Pirjola and Viljanen (1998) (the text and figures of pages 1434 to 1439 and the reference list; the proof of the equivalence in the Appendix was not checked); Boteler et al. (2024); Chakraborty et al. (2026); Chakraborty et al. (2022); Boteler and Pirjola (2003).
* **Checked against the code by test:** the formulas of `Current` against Ampère's law, the perfect-conductor limit and the image depth $h + 2p$.
* **Cited, not read:** Thomson and Weaver (1975); Boteler and Pirjola (1998); Häkkinen and Pirjola (1986); Pirjola et al. (2000); Lanzerotti et al. (2001).

## References

* Boteler, D. H., and Pirjola, R. J. (1998). The complex-image method for calculating the magnetic and electric fields produced at the surface of the Earth by the auroral electrojet. *Geophysical Journal International* 132, 31-40 (as listed by Pirjola and Viljanen, 1998; not read).
* Boteler, D. H., Chakraborty, S., Shi, X., Hartinger, M. D., and Wang, X. (2024). An examination of geomagnetic induction in submarine cables. *Space Weather* 22, e2023SW003687. doi:10.1029/2023SW003687
* Boteler, D. H., and Pirjola, R. J. (2003). Magnetic and electric fields produced in the sea during geomagnetic disturbances. *Pure and Applied Geophysics* 160, 1695-1716. doi:10.1007/s00024-003-2372-6
* Chakraborty, S., et al. (2026). Validating SCUBAS predictions of geomagnetically induced voltage in submarine cables using legacy superstorm observations. *Space Weather*. doi:10.1029/2025SW004843
* Pirjola, R., and Viljanen, A. (1998). Complex image method for calculating electric and magnetic fields produced by an auroral electrojet of finite length. *Annales Geophysicae* 16, 1434-1444. doi:10.1007/s00585-998-1434-6
* Pirjola, R., Viljanen, A., and Boteler, D. H. (2000). Electric field at the seafloor due to a two-dimensional ionospheric current. *Geophysical Journal International* 140, 286-294.
* Thomson, D. J., and Weaver, J. T. (1975). The complex image approximation for induction in a multilayered Earth. *Journal of Geophysical Research* 80, 123-129 (as listed by Pirjola and Viljanen, 1998; not read).
* Häkkinen, L., and Pirjola, R. (1986). Calculation of electric and magnetic fields due to an electrojet current system above a layered earth. *Geophysica* 22, 31-44 (as listed by Pirjola and Viljanen, 1998; not read).

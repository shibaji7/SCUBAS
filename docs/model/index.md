# The SCUBAS model

SCUBAS (Submarine Cables Upset by Auroral Streams) estimates the voltage that a geomagnetic disturbance induces on the power-feed conductor of a submarine cable. This page says what
the model is, which equations it solves, what it assumes, where it stops being reliable, and how it relates to other ways of attacking the same problem. The pages for the
individual [components](components/source.md) are drafts for the authors to complete.

:::{note}
Every statement here is tied to a paper that was read, with the section or equation; what was found by search but not read is marked. See [Sources and what was checked](#sources-and-what-was-checked).
Items marked **to confirm** are the writer's reading of the papers and need the authors' check.
:::

## What the model is for

A fibre-optic cable carries light, but a copper conductor along it feeds power to the repeaters, and that conductor is exposed to geomagnetic induction. A storm can change the
voltage the power-feed equipment must supply: the TAT-8 cable showed excursions of up to 700 V in March 1989 (Boteler et al., 2024; Chakraborty et al., 2026). SCUBAS computes that
voltage from a magnetic record, the route of the cable and a model of the Earth beneath it. It is a practical, scalable estimate; its authors say it is "not intended to replace"
three-dimensional electromagnetic models (Chakraborty et al., 2026, Section 1).

## The chain, from field to voltage

```
magnetic field B at an observatory
  -> transfer function T_X of the section (seawater depth d, Earth below)       [Transfer function]
  -> electric field E at the seafloor of the section                             [Transfer function]
  -> emf of the section = E dot (north and east lengths)                         [Cable line and network]
  -> transmission-line sections joined in a nodal network -> Earth potentials    [Cable line and network]
  -> cable voltage = sum of emfs + potential at one end - potential at the other
```

**1. The field at the seafloor.** The Earth is taken as horizontal layers (1-D) under each section. The relation between the surface electric and magnetic field is the surface
impedance of the layers (recursive formulas of Weaver, 1994, and Boteler and Pirjola, 2019). The magnetic field is known at the sea surface, not on the seafloor, so the attenuation
in the seawater of depth $d$ is included (Boteler and Pirjola, 2003), which gives the transfer function $T_X(f)$ of each section. In the frequency domain, with the electric field
north $E_N$ and east $E_E$ and the magnetic field east $B_E$ and north $B_N$ (Boteler et al., 2024, Equation 16):

$$
E_N(f) = T_X(f)\,B_E(f), \qquad -E_E(f) = T_X(f)\,B_N(f).
$$

The spectra of the record are multiplied by $T_X$ and transformed back to give $E(t)$ (Equations 17 and 18 there). Deeper water means a smaller transfer function, so the deep
ocean sees weaker fields than the continental shelf (Boteler et al., 2024, Figure 6; Chakraborty et al., 2026, Figures 3 and 4).

**2. The field along the cable is the field on the seafloor.** The loop formed by the cable and the seafloor beside it has effectively no area, so Faraday's law gives no flux
change through it, and the electric field along the cable equals the seafloor field (Boteler et al., 2024, Equations 8 to 10). The emf of section $i$ is
$\varepsilon_i = E^N_i L^N_i + E^E_i L^E_i$ (Equation 19, lengths from the WGS84 formulas of Horton et al., 2012) and the total is $\varepsilon_C = \sum_i \varepsilon_i$ (Equation 20).

**3. The Earth potentials at the ends.** The currents the field drives in sea and land charge the boundary and give potentials at the cable ends (the geoelectric coast effect).
They come from a transmission-line model of the seawater and the layers below, in the distributed-source transmission-line (DSTL) form, which is equivalent to generalised thin-sheet
modelling (Ranganayaki and Madden, 1980; Wang et al., 2023). Per unit length the series impedance $Z$ comes from the conductance of the seawater and sediment, the shunt admittance $Y$
from the resistive crust and mantle lithosphere; then

$$
\gamma = \sqrt{ZY}, \qquad Z_0 = \sqrt{Z/Y}, \qquad \text{adjustment distance } 1/\gamma .
$$

Each section becomes an equivalent-$\pi$ circuit, the land at each end an active termination, and the circuits are solved as one admittance network (Boteler et al., 2023;
Boteler et al., 2024, Section 3.2). Worked equations are in the [theory tutorial](../tutorial/theory.md).

**4. The voltage.** With $U_A$ and $U_B$ the Earth potentials at the two ends (Boteler et al., 2024, Equation 12):

$$
V_C = \varepsilon_C + U_A - U_B .
$$

In the code this is `Vt(v) = U0 - U1 + V(v)` of `Cable.tot_params`, with the section sum $\varepsilon_C$ held in `V(v)` (see the [glossary](../glossary.md)). A cable whose sections are
all electrically long ($L > 4/\gamma$) has end potentials near $-E/\gamma$; a short one ($L < 1/\gamma$) near $-LE/2$ ([tutorial](../tutorial/elsc.md)).

Why the end potentials matter: when a westward electrojet covers the whole TAT-8 route the induced field has one sign along it and the voltage is large, but when the two ends see
opposite disturbances the fields cancel (Boteler et al., 2024, Section 5). The end potentials reduce the effect of the emf further, because the potential difference between the ends opposes the induced field.

## Components

| Component | In the code | Page |
|---|---|---|
| Source (the auroral electrojet and the magnetic record) | `scubas.models.OceanModel.read_Bfield_data`, `scubas.sources` | [Source](components/source.md) |
| Earth conductivity | `scubas.datasets`, `scubas.conductivity` | [Earth conductivity](components/conductivity.md) |
| Transfer function | `scubas.models.OceanModel`, `Preprocess` | [Transfer function](components/transfer-function.md) |
| Cable line and network | `scubas.cables` | [Cable line and network](components/cable-network.md) |

```{toctree}
:hidden:
:maxdepth: 1

components/source
components/conductivity
components/transfer-function
components/cable-network
```

## Assumptions and where the model stops

| Assumption | Consequence | Where stated |
|---|---|---|
| Each section has a 1-D layered Earth; lateral change only between sections (the piecewise approach of Marti et al., 2014) | Transverse 2-D gradients and edge effects at abrupt bathymetry are not resolved | Boteler et al., 2024, Section 3.1; Chakraborty et al., 2026, Section 2.3 |
| The seafloor of a section is flat | Detail of the bathymetry inside a section is lost; sections come from depth and length by hand | Boteler et al., 2024, Section 4; Chakraborty et al., 2026, Section 2.1 |
| The magnetic field of a section is that of the nearest observatory; no spatial interpolation | Works best when the disturbance is coherent along the route; small-scale currents and the lack of ocean magnetometers limit accuracy | Boteler et al., 2024, Section 5; Chakraborty et al., 2026, Sections 6 and 7 |
| DSTL is a thin-sheet, thin-wire theory | Section length should not fall far below the adjustment distance scale; the papers warn against sections shorter than a few hundred kilometres | Chakraborty et al., 2026, Sections 2.3 and 7 |
| Uniform seawater conductivity | See the note on [Earth conductivity](components/conductivity.md) | Boteler et al., 2024, Section 3.1 |
| The electric and magnetic field are related by the surface impedance of the layers | Implicitly a source field that is uniform over the scale of a section (**to confirm**) | Boteler et al., 2024, Section 3.1 |
| The result is the induced voltage | It is not the response of the power-feed equipment, which holds the current constant | Boteler et al., 2024, Section 4.4 |

## How well it works

- **TAT-8, March 1989.** Calculated voltage variations 67, 249, 406 and 763 V against about 75, 300, 450 and 700 V reported by Medford et al. (1989) for four events (Boteler et al., 2024,
  Table 4). The authors caution that this may exaggerate the general skill, because the largest events are those with uniform disturbances, which land observatories represent best.
- **TAT-1 (February 1958) and TAT-8 (March 1989), with digitised data.** A median percentage error below 30 % and a skill score of 0.85 (Chakraborty et al., 2026, abstract and Section 6);
  the 1958 case has one magnetometer for the whole route and a larger error. The model underperforms for fast wave-like variations and when the structure of the current systems is not
  resolved. The data are digitised from printed figures, which adds error (Sections 3 and 7).
- Validation data are scarce because operators hold voltage records privately (Chakraborty et al., 2026, Section 7).

## Similar models and how SCUBAS differs

| Approach | What it does | Relation to SCUBAS | Source (status) |
|---|---|---|---|
| Induction in a loop (Faraday) | Takes the voltage as the rate of change of flux through a loop of cable and sea or Earth | The return path is ambiguous and neither choice satisfied early workers; SCUBAS uses the field along the cable and the end potentials instead | Winckler et al., 1959, and others as discussed in Boteler et al., 2024, Section 2 (cited, not read) |
| Seafloor field from the surface impedance | Computes $E$ at the seafloor from $B$ at the sea surface for a layered Earth | The first step of SCUBAS, section by section | Cox et al., 1970; Boteler and Pirjola, 2003; Goto, 2015 (cited in Boteler et al., 2024; not read) |
| Analytic cable-voltage theory | Voltage proportional to the magnitude of the field change (not its time derivative), from the attenuation in seawater and reflection at the bottom; design criteria for new systems | A closed-form treatment, where SCUBAS is a numerical section-by-section model with an Earth profile and active terminations | Mecozzi, 2022, arXiv:2204.00560 (abstract and introduction read) |
| Power-system GIC networks | Induced emfs in lines plus the nodal network of the grid | SCUBAS carries the same structure to cables: transmission-line sections, equivalent-$\pi$ circuits, admittance matrix | Boteler and Pirjola, 1998 and 2017 (cited in Boteler et al., 2024; not read) |
| Observed voltages on cables | Records of the power-feed voltage during storms | The data SCUBAS is compared with | Medford et al., 1981, 1989; Lanzerotti et al., 2001; Axe, 1968 (cited; not read) |
| 3-D electromagnetic models (FEM, FDTD, integral equation) | Solve Maxwell's equations on a 3-D conductivity model, including the coast effect | More complete and far heavier; SCUBAS does not replace them. A numerical coast model is on the [roadmap](../plans/coastal.md) | Liu et al., 2018 and other studies (found by search; not read) |

## Sources and what was checked

Read in full or in the sections used: Chakraborty et al. (2022), Boteler et al. (2024), Chakraborty et al. (2026). Read in part: Mecozzi (2022). Cited by those papers and not read here:
Boteler et al. (2023), Wang et al. (2023), Ranganayaki and Madden (1980), Marti et al. (2014), Medford et al. (1981, 1989), Lanzerotti et al. (2001), Price (1973), Cox et al. (1970),
Pasyanos et al. (2014), Weaver (1994), Winckler et al. (1959), Goto (2015). Found by web search and not read: Pirjola and Viljanen (1998), Liu et al. (2018), Liu et al. (2020), the USGS
hazard studies. Open the original before using any of these in a paper.

- Chakraborty, S., Boteler, D. H., Shi, X., Murphy, B. S., Hartinger, M. D., Wang, X., Lucas, G., & Baker, J. B. H. (2022). Modeling geomagnetic induction in submarine cables. *Frontiers in Physics*, 10, 1022475. doi:10.3389/fphy.2022.1022475
- Boteler, D. H., Chakraborty, S., Shi, X., Hartinger, M. D., & Wang, X. (2023). Transmission line modelling of geomagnetic induction in the Ocean/Earth conductivity structure. *International Journal of Geosciences*, 14(8), 767-791. doi:10.4236/ijg.2023.148041
- Boteler, D. H., Chakraborty, S., Shi, X., Hartinger, M. D., & Wang, X. (2024). An examination of geomagnetic induction in submarine cables. *Space Weather*, 22(2), e2023SW003687. doi:10.1029/2023SW003687
- Chakraborty, S., Hartinger, M. D., Shi, X., Boteler, D. H., Lawrence, E., Boteler, C., Baker, J. B. H., & Macalester, M. (2026). Validating SCUBAS predictions of geomagnetically induced voltage in submarine cables using legacy superstorm observations. *Space Weather*, 24, e2025SW004843. doi:10.1029/2025SW004843
- Wang, X., Boteler, D. H., & Pirjola, R. (2023). Distributed-source transmission line theory for modeling the coast effect on geoelectric fields. *IEEE Transactions on Power Delivery*, 38(5), 3541-3550. doi:10.1109/tpwrd.2023.3279462
- Mecozzi, A. (2022). Effects of geomagnetic field perturbations on the power supply of transoceanic fiber optic cables. arXiv:2204.00560.
- Medford, L. V., Lanzerotti, L. J., Kraus, J. S., & Maclennan, C. G. (1989). Transatlantic earth potential variations during the March 1989 magnetic storms. *Geophysical Research Letters*, 16(10), 1145-1148. doi:10.1029/gl016i010p01145 (as listed in Boteler et al., 2024; not read)
- Lanzerotti, L. J., Medford, L. V., Maclennan, C. G., Kraus, J. S., Kappenman, J., & Radasky, W. (2001). Trans-Atlantic geopotentials during the July 2000 solar event and geomagnetic storm. *Solar Physics*, 204, 351-359. doi:10.1023/a:1014289410205
- Pirjola, R., & Viljanen, A. (1998). Complex image method for calculating electric and magnetic fields produced by an auroral electrojet of finite length. *Annales Geophysicae*, 16, 1434. doi:10.1007/s00585-998-1434-6 (found by search; not read)

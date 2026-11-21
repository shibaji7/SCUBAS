<!--
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is licensed under the GNU Lesser General Public License v3.0 or later
(LGPL-3.0-or-later); see LICENSE and COPYING in the root directory. Versions up to
and including 0.1.9 were released under the MIT license (see NOTICE.md).
-->

# The DSTL model

**Status: implemented.** The distributed-source transmission-line (DSTL) model gives the **Earth potentials $U_A$ and $U_B$ at the two ends of a cable** and the potential along each section. It is the second half of the [SCUBAS model](scubas.md): the first half gives the electric field on the seafloor and the emf in the cable, and this page gives the term $U_A - U_B$ of $V_C = \mathcal{E}_C + U_A - U_B$.

:::{note}
Sources: B2023 is Boteler et al. (2023), the paper this page follows; C2022 is Chakraborty et al. (2022); B2024 is Boteler et al. (2024). Equation numbers are those of the papers. Every formula is written again, independently of the package, in `scripts/model_reference.py` and compared with the code by `tests/test_docs_model.py`; the numbers in the text and the tables come from `docs/model/figures/numbers.json`, which that test keeps equal to what the code computes. The derivation of the equivalence with the generalised thin sheet is in Wang et al. (2023), which was **not read** for this page (see [sources](#sources-and-what-was-checked)).
:::

## Why the coast matters

Electric fields drive currents in the sea and the land. The current density in seawater is far larger than on land, so charge accumulates at the coast and the electric field is modified on both sides of it: the geoelectric coast effect (B2023, Introduction; B2024, Section 2.3, after Price, 1973). A cable whose power feed is grounded at the two coasts feels the potential of the Earth there. For a cable across an ocean these coast potentials act against the emf in the cable and reduce its effect on the power feed (B2023, Discussion). To compute them one cannot treat each section alone: the voltage in a section depends on the fields in its neighbours. The DSTL model treats the whole ocean and the land at both ends as one circuit.

## Why a transmission line

The conducting surface layers (seawater and sediments) lie over resistive layers (crust and mantle lithosphere), which lie over the conducting mantle. Ranganayaki and Madden (1980) showed that for frequencies that penetrate both the conducting and the resistive layers this "thin double layer" can be treated as a *generalised thin sheet*: the current flows horizontally in the conducting layers and leaks down through the resistive ones. Wang et al. (2023) showed that the generalised thin sheet is equivalent to a **transmission line** with a series impedance $Z$ set by the conducting layers and a parallel admittance $Y$ set by the resistive layers (B2023, Section 3), so that standard circuit theory solves it.

![A layered ocean-Earth model and its transmission-line equivalent](../figures/Thinsheet-Model.png "(a) the four layers of the double layer, with their thicknesses d and conductivities σ or resistivities ρ; (b) a short element of the equivalent line, with the series impedance Z dx, the source E dx and the parallel admittance Y dx. After Boteler et al. (2023, Figure 2), CC BY 4.0, modified from Wang et al.")

## The line constants

The horizontal currents in the seawater (layer 1) and the sediments (layer 2) flow in parallel, so their conductances add. For a strip of width $w$ (B2023 Eqs 1 and 2):

$$
C = \sigma_1 d_1 w + \sigma_2 d_2 w , \qquad Z = \frac{1}{C} .
$$

The vertical currents through the crust (layer 3) and the mantle lithosphere (layer 4) flow in series, so their resistances add (B2023 Eqs 3 and 4):

$$
R = \frac{\rho_3 d_3}{w} + \frac{\rho_4 d_4}{w} , \qquad Y = \frac{1}{R} .
$$

$Z$ and $Y$ are per unit length: $Z$ in Ω/m and $Y$ in S/m (the papers quote them per kilometre). The propagation constant, the characteristic impedance and the **adjustment distance** are (B2023 Eq 5)

$$
\gamma = \sqrt{ZY}, \qquad Z_0 = \sqrt{\frac{Z}{Y}}, \qquad \text{adjustment distance } \frac{1}{\gamma}.
$$

In the code these are `TransmissionLine.calc_trasmission_line_parameters`, with the width $w$ in `elec_params.width` (default and in the tables 1 m). The seawater and the sediments (layers 0 and 1 of the `Site`) give $C$, and the crust and the mantle lithosphere (layers 2 and 3) give $R$. The width scales $Z \propto 1/w$, $Y \propto w$ and $Z_0 \propto 1/w$ and leaves $\gamma$ unchanged; the papers use $w = 1$ m, and the potentials do not depend on $w$ when every section and every termination uses the same one. For **land** (a `Site` named `Land`, whose first layer is the sediments and which has no water layer) the code uses the sediments for $C$ and the crust and mantle lithosphere for $R$.

The line constants of the example models of B2023 (its Table 2) and of the Atlantic model of C2022 (its Table 3), from the code:

| Section | Depth (m) | $Z$ (Ω/km) | $Y$ (S/km) | $\gamma$ (km⁻¹) | $Z_0$ (Ω) | $1/\gamma$ (km) |
|---|---|---|---|---|---|---|
| Shelf (`CS`) | 100 | {{ dstl_cs_Z }} | {{ dstl_cs_Y }} | {{ dstl_cs_gamma }} | {{ dstl_cs_Z0 }} | {{ dstl_cs_adj }} |
| Shallow ocean (`SO`) | 1000 | {{ dstl_so_Z }} | {{ dstl_so_Y }} | {{ dstl_so_gamma }} | {{ dstl_so_Z0 }} | {{ dstl_so_adj }} |
| Deep ocean (`DO`) | 4000 | {{ dstl_do_Z }} | {{ dstl_do_Y }} | {{ dstl_do_gamma }} | {{ dstl_do_Z0 }} | {{ dstl_do_adj }} |
| Land (`LD`) | 0 | {{ dstl_ld_Z }} | {{ dstl_ld_Y }} | {{ dstl_ld_gamma }} | {{ dstl_ld_Z0 }} | {{ dstl_ld_adj }} |
| Atlantic (C2022) | 4000 | {{ dstl_atl_Z }} | {{ dstl_atl_Y }} | {{ dstl_atl_gamma }} | {{ dstl_atl_Z0 }} | {{ dstl_atl_adj }} |

These agree with the printed tables (B2023 Table 2: 516.4, 632.5, 1183.2 and 258.2 km; C2022 Table 3: $Z = 7.14\times10^{-2}$ Ω/km, $Y = 5.88\times10^{-6}$ S/km, $Z_0 = 110$ Ω, $\gamma = 6.48\times10^{-4}$ km⁻¹).

A deeper ocean carries more current, which lowers $Z$, so $Z_0$ falls and $\gamma$ falls (C2022, Section 4.4): the adjustment distance grows with depth.

```{figure} figures/line_parameters_depth.png
:alt: Two panels against ocean depth from 0.1 to 5 km. The characteristic impedance falls from about 410 ohm to 100 ohm. The adjustment distance rises from about 410 km to 1700 km.
:width: 60%

The characteristic impedance $Z_0$ and the adjustment distance $1/\gamma$ of the Atlantic model against the depth of the ocean (C2022 Figure 8), from `TransmissionLine`.
```

## The line equations

An element of the line (the lower part of the figure above) has a series impedance $Z\,dx$, a source $E\,dx$ for the electric field induced along it, and a parallel admittance $Y\,dx$ to ground. With $V$ and $I$ the voltage and current (B2023 Eqs 6 and 7):

$$
\frac{dV}{dx} + Z I = E , \qquad \frac{dI}{dx} + Y V = 0 .
$$

Differentiating and substituting (B2023 Eqs 8 to 10), for a field $E$ uniform in the section ($dE/dx = 0$):

$$
\frac{d^2V}{dx^2} - \gamma^2 V = 0 .
$$

The solution has the form (B2023 Eqs 11 and 12)

$$
V = A e^{\gamma x} + B e^{-\gamma x}, \qquad
I = \frac{E}{\gamma Z_0} - \frac{A}{Z_0}e^{\gamma x} + \frac{B}{Z_0}e^{-\gamma x},
$$

and fixing $A$ and $B$ by the voltages $V_i$ at $x = 0$ and $V_k$ at $x = L$ gives the voltage along a section of length $L$ (B2023 Eq 19; C2022 Eq 10):

$$
V(x) = \frac{V_k e^{\gamma L} - V_i}{e^{\gamma L} - e^{-\gamma L}}\,e^{-\gamma(L-x)}
     + \frac{V_i e^{\gamma L} - V_k}{e^{\gamma L} - e^{-\gamma L}}\,e^{-\gamma x} .
$$

The voltage is the sum of two terms that fall off exponentially from the two ends with the characteristic length $1/\gamma$, the adjustment distance (B2023, after Eq 20). In the code this is `Cable._pot_along_cable_` and `TransmissionLine._pot_alongCS_`; $V_i$ and $V_k$ are the end potentials from the network below.

## Electrically long and short sections

Whether the two exponential terms overlap decides how a section behaves.

* **Electrically long**, $L > 4/\gamma$: $e^{\gamma L} \gg 1 \gg e^{-\gamma L}$, the effects of the two ends do not overlap, and the section looks from either end like an infinite line.
* **Electrically short**, $L \ll 1/\gamma$: $e^{\pm\gamma L} \approx 1 \pm \gamma L$, the line is a simple resistor with the source $EL$, and nothing leaks to the ground.

For an isolated section with open ends (zero admittance to ground at each end) the end potentials are (C2022 Eqs 11 and 12)

$$
U_A = -\frac{\sinh\gamma L}{\cosh\gamma L + 1}\,\frac{E}{\gamma} , \qquad U_B = -U_A .
$$

Because $\sinh x/(\cosh x + 1) = \tanh(x/2)$, this is one formula,

$$
U_A = -\frac{E}{\gamma}\tanh\frac{\gamma L}{2} ,
$$

which contains both limits of C2022 (Eqs 17 and 22): for a long section $\tanh \to 1$ and $U_A = -E/\gamma$, the field times the adjustment distance; for a short section $\tanh(x/2) \to x/2$ and $U_A = -LE/2$, half the emf of the section. The check of the identity and of both limits is in the tests.

```{figure} figures/long_short_limits.png
:alt: Log-log plot of the end potential in units of E over gamma against gamma L from 0.01 to 15. The curve follows the straight line gamma L over 2 at small values, then bends and approaches 1 above gamma L of about 4. Vertical dotted lines mark gamma L of 1 and 4.
:width: 60%

The end potential of an isolated section against its length in adjustment distances. For $L < 1/\gamma$ it is $LE/2$; for $L > 4/\gamma$ it is within 4 % of $E/\gamma$ ($\tanh 2 = 0.964$; at $L = 10/\gamma$ the difference is 0.01 %). The "electrically long" criterion is therefore approximate: it does not mean that the potential has reached $E/\gamma$ to better than a few percent.
```

```{figure} figures/section_potentials.png
:alt: Two panels. Left, the potential along an electrically short section of 100 km rises linearly from about -50 V to +50 V. Right, along a long section of 10000 km it follows an S shape from about -1540 V to +1540 V. Red crosses lie on the blue curves.
:width: 100%

The potential along an electrically short (100 km) and an electrically long (10 000 km) section of the Atlantic model for $E = 1$ V/km (C2022 Figure 6): the potential of `Cable` (blue) and the closed form C2022 Eqs 10 to 12 (red crosses). The paper's panel title reads "$E_f = 1$ mV/km", but its values ($\pm50$ V for 100 km, $\pm1543$ V for the long section) are those of 1 V/km, as here.
```

## The equivalent-π circuit

To join sections, each is replaced by an **equivalent-π circuit** (Boteler, 2013, as cited by B2023 Section 3.2): a series element and an admittance to ground at each end, with the same terminal behaviour as the line. With a voltage source (B2023 Eqs 29 to 31):

$$
Z' = Z_0\sinh\gamma L , \qquad
\frac{Y'}{2} = \left(\cosh\gamma L - 1\right)\frac{1}{Z_0\sinh\gamma L} , \qquad
E' = \frac{E}{\gamma}\sinh\gamma L ,
$$

or, equivalently, with a **current source** $I_E$ in parallel with a series admittance $Y_E$ (the Norton form; B2023 Eqs 32 to 34), which is the form the nodal analysis needs:

$$
Y_E = \frac{1}{Z_0\sinh\gamma L} , \qquad
\frac{Y'}{2} = \left(\cosh\gamma L - 1\right) Y_E , \qquad
I_E = \frac{E}{Z} .
$$

The two forms are the same circuit: $Y_E = 1/Z'$ and $I_E Z' = E'$. In the code these are `Ye`, `Yp2` and `Ie` of `TransmissionLine` (in S and A, with $E$ in V/m), computed by `compute_eqv_pi_circuit`. For an electrically short section $Z' \to ZL$, $Y'/2 \to 0$ and $E' \to EL$ (B2023 Eqs 35 to 37), the simple resistor with a source.

![The two forms of the equivalent-pi circuit](../figures/Circuit-Equivalent.png "Equivalent-π circuit of a section: (a) with a voltage source E′ and series impedance Z′; (b) with a current source I_E and series admittance Y_E. After Boteler et al. (2023, Figure 4), CC BY 4.0")

## The land at the ends

The land at each end of a cable is a region with its own $Z$, $Y$ and field, so long that it is an electrically long line seen from its end. Looking into the end of such a line, the open-circuit voltage is $-E/\gamma$ and the short-circuit current is $E/(\gamma Z_0)$, so the line is equivalent to a **Thevenin** circuit with (B2023 Eqs 23 to 26)

$$
V_{Th} = -\frac{E}{\gamma} , \qquad Z_{Th} = Z_0 ,
$$

or to a **Norton** circuit, a current source in parallel with an admittance (B2023 Eqs 27 and 28):

$$
J_N = \frac{V_{Th}}{Z_{Th}} = \frac{E}{Z}, \qquad Y_N = \frac{1}{Z_0} .
$$

(B2023 write the Norton current as $E/Z$ and note that the Thevenin voltage at the other end of a line is $+E/\gamma$: the sign with which the source enters a node depends on whether the field points towards or away from that end, which the nodal equations below take into account, $J_1 = j_W - j_{12}$ and $J_4 = j_{34} - j_E$.) This is the **active termination**: the land is not a passive ground but a source with the field of the land. `TransmissionLine.add_active_termination` computes $J_N$ and $Y_N$ from the land profile and from the land magnetic field, and stores them in `term_params` (`Jn`, `Yn`).

Two limiting cases are given in B2023 (Section 3.3). A section short enough to be electrically short, between two active terminations $(V_1, Z_1)$ and $(V_2, Z_2)$, carries the current and has the end potentials

$$
I = \frac{V_1 + EL + V_2}{Z_1 + ZL + Z_2}, \qquad U_1 = V_1 - I Z_1 , \qquad U_2 = I Z_2 - V_2
$$

(Eqs 38 to 40; the [active-termination tutorial](../tutorial/active.md) uses these). The boundary between two regions that are both electrically long has (Eqs 41 and 42)

$$
I = \frac{V_{Th1} + V_{Th2}}{Z_{Th1} + Z_{Th2}} , \qquad U_{boundary} = V_{Th1} - I Z_{Th1} = I Z_{Th2} - V_{Th2} .
$$

For a very wide continental shelf of 0.3 V/km against land of 0.3 V/km, the Thevenin circuits are $V_{Th}$ = {{ dstl_vth_cs }} V, $Z_{Th}$ = {{ dstl_zth_cs }} Ω (shelf) and {{ dstl_vth_ld }} V, {{ dstl_zth_ld }} Ω (land), so $I$ = {{ dstl_i_cs_ld }} A and the coast potential is the limit **{{ dstl_uc_limit }} V** (B2023 Eqs 56 and 57). For a very wide shelf against a deep ocean of 0.1 V/km ($V_{Th}$ = {{ dstl_vth_do }} V, $Z_{Th}$ = {{ dstl_zth_do }} Ω), $I$ = {{ dstl_i_do_cs }} A and the potential at the shelf edge is the limit **{{ dstl_ue_limit }} V** (Eqs 58 and 59).

## The network

All the sections and the land are combined in one circuit, because the potential of one section depends on its neighbours. Each section is its equivalent-π circuit and each land end its Norton circuit; adding the admittances to ground that meet at each node gives a **nodal admittance network**.

![The network of a deep ocean with a shallow sea on each side and land at the ends](../figures/Network-Model.png "(a) Equivalent-π circuits of three sections (shallow sea A, deep ocean B, shallow sea C) and the Norton circuits of the land at the two ends; (b) the nodal admittance network. After Boteler et al. (2023, Figure 8), CC BY 4.0")

For the three sections of the figure (B2023 Eqs 43 to 46), with $y'_A/2$ the shunt admittance $Y'/2$ of section A and so on, and $y_W$, $y_E$ the Norton admittances of the land at the west and east ends:

$$
y_1 = \frac{y'_A}{2} + y_W, \quad y_2 = \frac{y'_A}{2} + \frac{y'_B}{2}, \quad y_3 = \frac{y'_B}{2} + \frac{y'_C}{2}, \quad y_4 = \frac{y'_C}{2} + y_E ,
$$

and the series admittance $Y_E$ and the current source $I_E$ of each section (Eqs 32 and 34, labelled $Y_A, Y_B, Y_C$ and $J_A, J_B, J_C$ in the figure) become those of the network, $y_{12} = Y_A$, $y_{23} = Y_B$, $y_{34} = Y_C$ and $j_{12} = J_A$, $j_{23} = J_B$, $j_{34} = J_C$. Kirchhoff's current law at each node, with the current in each line $i_{12} = (v_1 - v_2)y_{12} + j_{12}$ and so on, gives (B2023 Eqs 47 to 53; C2022 Eqs 6 to 9)

$$
\begin{bmatrix} J_1\\ J_2\\ J_3\\ J_4 \end{bmatrix} =
\begin{bmatrix}
y_1 + y_{12} & -y_{12} & 0 & 0\\
-y_{12} & y_{12} + y_{23} + y_2 & -y_{23} & 0\\
0 & -y_{23} & y_{23} + y_{34} + y_3 & -y_{34}\\
0 & 0 & -y_{34} & y_{34} + y_4
\end{bmatrix}
\begin{bmatrix} v_1\\ v_2\\ v_3\\ v_4 \end{bmatrix},
$$

with the currents entering the nodes

$$
J_1 = j_W - j_{12}, \quad J_2 = j_{12} - j_{23}, \quad J_3 = j_{23} - j_{34}, \quad J_4 = j_{34} - j_E ,
$$

that is $[J] = [Y][V]$, and the potentials at the ends of the sections are (B2023 Eq 55)

$$
[V] = [Y]^{-1}[J] .
$$

The potentials $V_1$ and $V_4$ at the two outer nodes are the Earth potentials $U_A$ and $U_B$ at the two coasts that enter $V_C$; substituting the node potentials in the section formula above gives the potential along every section. The code is `Cable.run_nodal_analysis`, which builds the matrix for each field component, and `Cable.solve_admitance_matrix`, which inverts it; the potentials of the cable ends are then `U0` and `U1` of `Cable.tot_params`.

### A worked example

The deep ocean of B2023 Figure 9c: shallow sea A of 100 km (0.3 V/km), a deep ocean B of 7800 km (0.1 V/km), shallow sea C of 100 km (0.3 V/km), and land at both ends (0.3 V/km), so the total width from coast to coast is 8000 km.

| | Section A (shelf) | Section B (deep ocean) | Section C (shelf) |
|---|---|---|---|
| Length (km) | {{ net_a_L }} | {{ net_b_L }} | {{ net_c_L }} |
| $\gamma L$ | {{ net_a_gL }} | {{ net_b_gL }} | {{ net_c_gL }} |
| $Y_E$ (S) | {{ net_a_Ye }} | {{ net_b_Ye }} | {{ net_c_Ye }} |
| $Y'/2$ (S) | {{ net_a_Yp2 }} | {{ net_b_Yp2 }} | {{ net_c_Yp2 }} |
| $I_E$ (A) | {{ net_a_Ie }} | {{ net_b_Ie }} | {{ net_c_Ie }} |

The land at each end has $Y_N$ = {{ net_land_Yn }} S and $J_N$ = {{ net_land_Jn }} A. The nodal equations are then (in S and A):

| | $v_1$ | $v_2$ | $v_3$ | $v_4$ | | $J$ (A) |
|---|---|---|---|---|---|---|
| node 1 | {{ net_Y11 }} | {{ net_Y12 }} | 0 | 0 | | {{ net_J1 }} |
| node 2 | {{ net_Y21 }} | {{ net_Y22 }} | {{ net_Y23 }} | 0 | | {{ net_J2 }} |
| node 3 | 0 | {{ net_Y32 }} | {{ net_Y33 }} | {{ net_Y34 }} | | {{ net_J3 }} |
| node 4 | 0 | 0 | {{ net_Y43 }} | {{ net_Y44 }} | | {{ net_J4 }} |

and the solution is $v_1, \dots, v_4$ = {{ net_U1 }}, {{ net_U2 }}, {{ net_U3 }}, {{ net_U4 }} V. The Earth potentials at the two coasts are therefore $\mp${{ net_deep_end }} V. (The potential varies much more at the ends than in the middle of the ocean: the middle section is electrically long, $\gamma L$ = {{ net_b_gL }}, and the 100 km shelves are electrically short.)

```{figure} figures/network_models.png
:alt: Three panels of Earth potential in volts against distance along the cable. Shallow sea of 600 km: a straight line from about -48 V to +48 V. Shallow ocean of 2000 km: an S shape from -48 V to +48 V with kinks where the shelves join. Deep ocean of 8000 km: a steeper S shape from -104.5 V to +104.5 V.
:width: 60%

The Earth potential along the cable for the three example models of B2023 Figure 9 (its Figure 10): a shallow sea of 600 km and 100 m depth, a shallow ocean of 2000 km and 1 km depth, and a deep ocean of 8000 km and 4 km depth, with 0.1 V/km in the deep ocean and 0.3 V/km elsewhere, from `Cable`. The potential is linear for the electrically short sea and S-shaped for the electrically long ocean; it is largest at the coasts.
```

## Validation

The model is validated against closed forms (C2022 Section 4.2; B2023 Sections 5 and 6), in the tests of `tests/test_docs_model.py`:

* **Line constants:** the three rows of B2023 Table 2, the land row and C2022 Table 3 (the table above).
* **One section:** the end potentials against $-(E/\gamma)\tanh(\gamma L/2)$ for lengths from 10 to 10 000 km, and the potential along the section against B2023 Eq 19. C2022 reports a correlation of $1 - 10^{-9}$ for these.
* **The network:** the potentials of the code against a separate solution of B2023 Eqs 43 to 55 for the three example models, and the end potentials against the figures of B2023: about ±48 V for the shallow sea and the shallow ocean ({{ net_shallow_end }} and {{ net_so_end }} V) and about ±105 V for the deep ocean ({{ net_deep_end }} V), and 99 V for a deep ocean with no shelf (B2023, Discussion).
* **The limits:** the potential at the coast of a very wide shelf tends to the {{ dstl_uc_limit }} V of Eq 57.

```{figure} figures/coast_potentials.png
:alt: Two panels. Left, the coast potential of a shallow sea rises from 0 at zero width to about 75 V at 2000 km, approaching a dashed limit of 77.46 V. Right, for a deep ocean, the coast potential rises from 99 V to a maximum of about 107 V at 220 km of shelf and falls to about 79 V at 2000 km; the shelf-edge potential falls from 99 V to about 70 V. Dashed lines mark the limits 77.46 V and 69.38 V.
:width: 100%

The coast potential against the width of a shallow sea (left, B2023 Figure 11) and of the continental shelf of a deep ocean (right, B2023 Figure 12), from `Cable`; the dashed lines are the limits of the preceding section. The coast potential of the deep ocean first rises and then falls with the width of the shelf: for a narrow shelf the shelf adds its emf to the current across the boundary, and for a shelf wider than about 250 km (the adjustment distance of the shelf is {{ dstl_cs_adj }} km) current leaks to the ground through the resistive layers (B2023, Section 6).
```

## Assumptions and limits

* **Thin-sheet and thin-wire.** The line is the model of a conducting double layer over a resistive one, valid for periods long enough to penetrate the water and the crust; at short periods the skin depth is smaller and other methods, such as finite elements, are needed. B2023 present its results as test cases for such methods (Discussion). Chakraborty et al. (2026, Sections 2.3 and 7) caution against sections much shorter than the adjustment distance or than a few hundred kilometres, and against the large lateral changes of a coast inside one section.
* **Uniform field in a section** ($dE/dx = 0$, B2023 Eq 10): a section is chosen to follow the large-scale changes of depth, and the field is that of the nearest observatory.
* **Land as semi-infinite** at both ends: the land beyond the shore is an electrically long line (B2023, Section 5 uses the land parameters at both ends of every model).
* **The "electrically long" criterion** $L > 4/\gamma$ is approximate (the figure above): the potential is then 96 % of $E/\gamma$.
* **The result is a voltage,** not the response of the power-feed equipment, which holds the current constant (B2024, Section 4.4).

A numerical model of the coast effect, to test where these assumptions fail, is planned ([numerical coast model](fdtd.md)).

## In the code

```python
import pandas as pd

from scubas.cables import Cable, TransmissionLine
from scubas.datasets import PROFILES

# One section of continental shelf, 600 km, a uniform north field of 300 mV/km (0.3 V/km), no land at the ends
line = TransmissionLine("shelf", directed_length={"length_north": 600.0}, elec_params={"site": PROFILES.CS})
line.compute_eqv_pi_circuit(Efield=pd.DataFrame({"X": [300.0]}, index=pd.RangeIndex(1, name="Time")), components=["X"])

print(line.Z, line.Y, line.gma, line.Z0)  # per metre: the constants of the table above are 1000 times larger
print(line.Ye["X"], line.Yp2["X"], line.Ie["X"])  # Y_E, Y'/2 and I_E of the pi circuit

cable = Cable([line], components=["X"])  # the network of one section
print(cable.tot_params[["U0", "U1", "V(v)", "Vt(v)"]])  # the end potentials, the emf and V_C
```

For the normal use with the land at the ends, give `active_termination={"left": PROFILES.LD, "right": PROFILES.LD}` and the magnetic record of the land (`bfield_data_files`); the [tutorials](../tutorial/active.md) show it. `scripts/model_cases.py` attaches the same Norton circuit from a given land field instead, which is how the papers specify the examples above.

## Sources and what was checked

* **Read in full for this page:** Boteler et al. (2023) (the equations above are those of the PDF, not of its web page); Chakraborty et al. (2022); Boteler and Pirjola (2003); Boteler et al. (2024); Chakraborty et al. (2026).
* **Checked against the code by test:** B2023 Eqs 1 to 5, 19, 23 to 42, 43 to 55, 56 to 60, Tables 2 and 3, Figures 10 to 12; C2022 Eqs 10 to 12, 14 to 23 (limits), Table 3, Figure 6.
* **Cited, not read:** Wang et al. (2023), the equivalence of the generalised thin sheet and the transmission line, with a finite-element check (IEEE Transactions on Power Delivery, doi 10.1109/TPWRD.2023.3279462); Ranganayaki and Madden (1980); Price (1973); Boteler (1997) and Boteler (2013), which B2023 cites for the distributed-source line and the equivalent-π circuit (not read; B2023 Sections 3 and 4 give the equations used here); the Supplementary Document of C2022 (Equations S1 to S13).

## References

* Boteler, D. H. (1997). Distributed-source transmission line theory for electromagnetic induction studies. *Proceedings of the 1997 Zurich EMC Symposium*, 401-408.
* Boteler, D. H. (2013). A new versatile model of geomagnetic induction of telluric currents in pipelines. *Geophysical Journal International* 193, 98-109. doi:10.1093/gji/ggs113
* Boteler, D. H., Chakraborty, S., Shi, X., Hartinger, M. D., and Wang, X. (2023). Transmission line modelling of geomagnetic induction in the ocean/Earth conductivity structure. *International Journal of Geosciences* 14, 767-791. doi:10.4236/ijg.2023.148041
* Boteler, D. H., Chakraborty, S., Shi, X., Hartinger, M. D., and Wang, X. (2024). An examination of geomagnetic induction in submarine cables. *Space Weather* 22, e2023SW003687. doi:10.1029/2023SW003687
* Chakraborty, S., Boteler, D. H., Shi, X., Murphy, B. S., Hartinger, M. D., Wang, X., Lucas, G., and Baker, J. B. H. (2022). Modeling geomagnetic induction in submarine cables. *Frontiers in Physics* 10, 1022475. doi:10.3389/fphy.2022.1022475
* Chakraborty, S., et al. (2026). Validating SCUBAS predictions of geomagnetically induced voltage in submarine cables using legacy superstorm observations. *Space Weather*. doi:10.1029/2025SW004843
* Price, A. T. (1973). The theory of geomagnetic induction. *Physics of the Earth and Planetary Interiors* 7, 227-233.
* Ranganayaki, R. P., and Madden, T. R. (1980). Generalized thin sheet analysis in magnetotellurics: an extension of Price's analysis. *Geophysical Journal International* 60, 445-457. doi:10.1111/j.1365-246X.1980.tb04820.x
* Wang, X., Boteler, D. H., and Pirjola, R. (2023). Distributed-source transmission line theory for modeling the coast effect on geoelectric fields. *IEEE Transactions on Power Delivery*. doi:10.1109/TPWRD.2023.3279462

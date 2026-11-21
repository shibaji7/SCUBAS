<!--
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is licensed under the GNU Lesser General Public License v3.0 or later
(LGPL-3.0-or-later); see LICENSE and COPYING in the root directory. Versions up to
and including 0.1.9 were released under the MIT license (see NOTICE.md).
-->

# The SCUBAS model

**Status: implemented.** This page follows the model in the order the code runs it. The transmission-line part (the Earth potentials at the ends of the cable) has its own page, [the DSTL model](dstl.md).

A fibre-optic submarine cable carries light, but a copper conductor along it feeds power to the repeaters, and that conductor is exposed to geomagnetic induction. SCUBAS (Submarine Cables Upset by Auroral Streams) computes the **voltage that a geomagnetic disturbance adds to or subtracts from the power feed**, from three inputs: a record of the magnetic field, the route of the cable, and a model of the Earth beneath it (Chakraborty et al., 2022; Boteler et al., 2024).

:::{note}
Sources are cited as BP2003 (Boteler and Pirjola, 2003), C2022 (Chakraborty et al., 2022), B2023 (Boteler et al., 2023) and B2024 (Boteler et al., 2024), with the equation number. Every formula on this page is written out again, independently of the package, in `scripts/model_reference.py` and compared with the code by `tests/test_docs_model.py`; every number that is quoted in the text or in a table comes from `docs/model/figures/numbers.json`, which that test keeps equal to what the code computes. The [last section](#sources-and-what-was-checked) says what was read and what was not.
:::

## 1. The chain

```{figure} figures/chain.png
:alt: Flow diagram. Top row, blue: magnetic record, detrend taper and FFT, multiply by the transfer function, inverse FFT, emf of the section. Bottom row, yellow: layered Earth, line constants, equivalent-pi circuit, nodal network, end potentials. Both rows lead to the cable voltage.
:width: 100%

From the magnetic record to the cable voltage, with the class or function that does each step. The blue row is the field along the cable; the yellow row is the Earth under it (the [DSTL model](dstl.md)). The electric field also drives the network, as the current sources $I_E = E/Z$.
```

The cable is cut into **sections**. In each section the Earth is a stack of horizontal layers (a 1-D model), and the magnetic field is that of the nearest observatory. The model does two things for every section:

1. It converts the magnetic record into the electric field on the seafloor, through a transfer function that depends on the depth of the water and on the Earth below it (sections 3 and 4). The electric field integrated along the section is its **emf** (section 5).
2. It models the conducting sea and Earth as a transmission line, and combines the sections and the land at the two ends in one network, which gives the **Earth potentials** $U_A$ and $U_B$ at the ends of the cable ([DSTL page](dstl.md)).

The voltage on the power feed is then (B2023 Eq 60; C2022 Eq 27; B2024 Eq 12)

$$
V_C = \mathcal{E}_C + U_A - U_B , \qquad \mathcal{E}_C = \sum_j \mathcal{E}_j .
$$

## 2. The inputs

**The magnetic record.** `OceanModel.read_Bfield_data` reads IAGA-2002 files (`.txt`, `.min`) or CSV files with a `Date` column and the components `X` (north), `Y` (east), `Z` and `F`; files are concatenated in time order and gaps are filled by a second-order spline. The records come from INTERMAGNET observatories ([data access](../user/intermagnet.md)) or are synthetic (section 7). The model takes the field at the sea surface to be that measured at the land observatory (BP2003, Seafloor Power Spectra). That is exact for a plane-wave source over a layered Earth, whose total horizontal field at the surface is twice the incident field whatever the layering (BP2003, Appendix), but it is an approximation near a coast and for a source of small scale ([limits](#assumptions-and-limits)).

**The route.** A section is given by its north and east lengths, by its total length, or by its end points. For end points the code uses the WGS84 series of C2022 Eqs 4 and 5 (the Appendix of Horton et al., 2012), with $\bar\lambda$ the mean latitude of the two ends:

$$
L^N = (111.133 - 0.56\cos 2\bar\lambda)\,\delta\lambda , \qquad
L^E = (111.5065 - 0.1872\cos 2\bar\lambda)\cos\bar\lambda\;\delta\phi ,
$$

in kilometres, with $\delta\lambda$ and $\delta\phi$ in degrees (`CableSection.compute_lengths`). The code takes the absolute value of the differences, so the lengths carry no sign (see [paper and code](#paper-and-code)).

**The Earth.** Each section has a `Site`, a list of `Layer` objects with a thickness (the package stores metres) and a conductivity. The built-in `PROFILES` hold published models; `ConductivityProfile` builds a profile from LITHO1.0 at the position of a section (Pasyanos et al., 2014; C2022 Section 2.1.1) and assigns a resistivity to each layer. The Atlantic model of C2022 Table 2, used for the examples on this page:

| Layer | Thickness (km) | Resistivity (Ω m) |
|---|---|---|
| Seawater | 4 | 0.3 |
| Sediments | 2 | 3 |
| Crust | 10 | 3000 |
| Mantle lithosphere | 140 | 1000 |
| Upper mantle | 254 | 100 |
| Transition zone | 250 | 10 |
| Lower mantle | 340 | 1 |

![The layered Ocean-Earth model, with shallow seas and land on either side](../figures/Ocean-Earth-Layer.png "Ocean/Earth conductivity structure. After Boteler et al. (2023, Figure 1), CC BY 4.0")

## 3. The seafloor field

The magnetic field is measured on land or at the sea surface, but a cable lies on the seafloor, below water that attenuates the fields. BP2003 derive the relation between the two.

**A plane wave in the sea.** Take a plane wave at vertical incidence on an ocean of depth $d$ over a horizontally layered Earth, with $x$ north, $y$ east, $z$ down and the time dependence $e^{i\omega t}$. Displacement currents are neglected, which is justified at the low frequencies of geomagnetic variations. The field in the ocean is a downgoing and a reflected wave (BP2003 Eqs 1 to 4):

$$
E_x = S e^{-kz} + R e^{kz}, \qquad
H_y = \frac{1}{Z}\left(S e^{-kz} - R e^{kz}\right), \qquad
k = (i\omega\mu_0\sigma)^{1/2}, \quad
Z = \left(\frac{i\omega\mu_0}{\sigma}\right)^{1/2},
$$

with $\sigma$ the conductivity of the seawater and $\mu_0 = 4\pi\times10^{-7}$ H/m. The ratio of the reflected to the downgoing wave is fixed by the impedance $Z_d = E_d/H_d$ of the Earth below the seafloor (BP2003 Eq 8):

$$
\alpha = \frac{R}{S} = \frac{Z_d - Z}{Z_d + Z}\,e^{-2kd}.
$$

$Z_d$ is the surface impedance of the layers under the seawater. For a layered Earth it follows from a recursion from the half-space at the bottom up to the seafloor, which is the standard technique of magnetotellurics (BP2003, Theory; `Site.calcZ`).

**The seafloor over the sea surface.** Setting $z = d$ and using $\alpha$ gives the electric and magnetic fields at the seafloor in terms of the magnetic field at the sea surface (BP2003 Eqs 10 and 13; the same relations as C2022 Eqs 1 and 2 with $B = \mu_0 H$):

$$
\frac{E_d}{H_0} = \frac{2Z}{\left(1 + Z/Z_d\right)e^{kd} - \left(1 - Z/Z_d\right)e^{-kd}}
= \frac{Z_d}{\cosh kd + (Z_d/Z)\sinh kd},
$$

$$
\frac{B_d}{B_0} = \frac{1}{\cosh kd + (Z_d/Z)\sinh kd}.
$$

The **transfer function** of SCUBAS is the first ratio divided by $\mu_0$,

$$
T_x(f) = \frac{E_d}{B_0} = \frac{1}{\mu_0}\,\frac{E_d}{H_0},
$$

computed by `OceanModel.calcTF` with the key `E2B` (the key `B2B` is $B_d/B_0$). It is in mV/km per nT, and complex: its phase is the phase delay of the electric field behind the magnetic field. Its default frequency grid is $10^{-6}$ to $1$ Hz.

**What the water does.** The propagation constant is $k = (1+i)/\delta$, with the skin depth $\delta = \sqrt{2\rho_s/(\omega\mu_0)}$ (BP2003 Eq 16), so $kd = (1+i)\,d/\delta$ and everything depends on the depth of the ocean in skin depths. For seawater of 0.25 Ω m and a 5 min period, $\delta = 4.36$ km. Two limits are useful (BP2003 Eqs 18 to 29):

* **Shallow water, $|kd| \ll 1$** (the thin sheet). The sea changes the magnetic field a little, $B_d/B_0 = 1/(1 + (Z_d/Z)kd)$, and the electric field hardly at all, $E_d/E_0 \approx 1$: the seafloor field is the field the Earth below would have without the sea.
* **Deep water, depth more than about two skin depths.** Both fields decay as $e^{-kd}$, and $E_d/H_0 \approx 2Z e^{-kd}$ when also $Z \ll Z_d$ (the factor 2 is the reflection at the seafloor).

```{figure} figures/attenuation_with_depth.png
:alt: Three panels of depth against amplitude ratio for oceans of 100 m, 1 km and 5 km. The magnetic field falls roughly linearly to 0.87, 0.35 and 0.085 at the seafloor; the electric field falls to 1.0, 0.96 and 0.57.
:width: 70%

The electric field (blue) and the magnetic field (red, dashed) against depth for a 5 min period (seawater 0.25 Ω m over a half-space of 10 Ω m): the conditions of BP2003 Figure 1, from the formulas of that paper. The yellow dots are the seafloor values of the magnetic field from `OceanModel` (the key `B2B`). In the 100 m and 1 km oceans the electric field is nearly constant with depth while the magnetic field is not; in the 5 km ocean, which is thicker than the skin depth of 4.36 km, both decay.
```

```{figure} figures/seafloor_ratios.png
:alt: Two log-log panels against frequency from 0.0001 to 0.01 Hz. Left, the magnetic ratio for 100 m stays near 1 and for 5 km falls from 0.4 to 0.04. Right, the electric ratio rises from 0.07 to 0.55 for 100 m and stays near 0.03 for 5 km; a dashed line is the response with no sea.
:width: 100%

The seafloor field over the surface magnetic field for oceans of 100 m and 5 km (seawater 0.25 Ω m over 10 Ω m), as in BP2003 Figures 4 and 5. Lines: `OceanModel`; markers: BP2003 Eq 13. The dashed line is the electric response of the Earth with no sea (depth zero). At 100 m the seafloor electric field is close to it; at 5 km the electric field is far smaller and nearly flat in frequency.
```

**The transfer function of a layered ocean-Earth.** For the Atlantic model the transfer function changes by an order of magnitude with the depth of the water (C2022 Figure 7). Deeper water means a smaller electric field at the seafloor at every frequency above about $10^{-5}$ Hz, and a phase that falls with frequency; the Earth below and the sea together act as a band-pass filter on the disturbance (C2022 Section 5.1).

```{figure} figures/transfer_function_depth.png
:alt: Amplitude on a log scale and phase in degrees of the transfer function against frequency from 1e-6 to 1e-2 Hz for ocean depths of 10 m, 100 m, 1 km and 4 km. The curves coincide at low frequency and separate above 1e-5 Hz, the 4 km curve flattening near 0.05.
:width: 70%

The transfer function $T_x$ of the Atlantic model of C2022 Table 2 for four ocean depths (C2022 Figure 7), from `OceanModel.calcTF`. All curves merge at low frequency, where the skin depth in the sea is much larger than the depth; the deep ocean shields the seafloor at higher frequency.
```

## 4. From the record to the field in time

The transfer function is applied in the frequency domain, one component at a time. With $E_N$, $E_E$ the north and east electric field and $B_N$, $B_E$ the north and east magnetic field (B2024 Eq 16; C2022 Section 2.2):

$$
E_N(f) = T_x(f)\,B_E(f), \qquad -E_E(f) = T_x(f)\,B_N(f).
$$

So the north electric field comes from the east magnetic field, and the east electric field comes from the north magnetic field with the opposite sign. `OceanModel.to_Efields` does, for each component:

1. Subtracts a baseline from the record and tapers its two ends (`Preprocess.detrend_magnetic_field`); the taper reduces the spurious frequencies that a record with abrupt ends would produce (C2022 Section 2.2).
2. Takes the real Fourier transform (`scubas.utils.fft`).
3. Multiplies by $T_x(f)$ at each frequency of the transform, with the sign and component swap of the equations above (`component_mappings`, `component_sign_mappings`).
4. Transforms back (`scubas.utils.ifft`).

The six-sine test of C2022 Section 2.4.2 checks this step against a closed form. The magnetic record is a sum of six sines with the amplitudes, periods and phases of C2022 Table 4,

$$
B(t) = \sum_{m=1}^{6} A_m \sin\left(2\pi f_m t + \Phi_m\right),
$$

and since $T_x$ is known at each of the six frequencies the electric field is (C2022 Eq 25)

$$
E_a(t) = \sum_{m=1}^{6} |T_x(f_m)|\,A_m \sin\left(2\pi f_m t + \Phi_m + \theta_m\right), \qquad \theta_m = \arg T_x(f_m).
$$

| $m$ | $A_m$ (nT) | period (min) | $\Phi_m$ (°) | amplitude of $E_m$ (mV/km) | $\Phi_m + \theta_m$ (°) |
|---|---|---|---|---|---|
| 1 | 200 | 180 | 10 | {{ synth_amp_1 }} | {{ synth_phase_1 }} |
| 2 | 90 | 80 | 20 | {{ synth_amp_2 }} | {{ synth_phase_2 }} |
| 3 | 30 | 36 | 30 | {{ synth_amp_3 }} | {{ synth_phase_3 }} |
| 4 | 17 | 15 | 40 | {{ synth_amp_4 }} | {{ synth_phase_4 }} |
| 5 | 8 | 8 | 50 | {{ synth_amp_5 }} | {{ synth_phase_5 }} |
| 6 | 3.5 | 3 | 60 | {{ synth_amp_6 }} | {{ synth_phase_6 }} |

The last two columns are computed with `OceanModel` for the Atlantic model; they reproduce the values printed in C2022 Table 4 (10.35, 4.86, 1.65, 0.95, 0.45, 0.19 mV/km and 20.66°, 25.0°, 31.22°, 36.88°, 41.93°, 36.37°).

```{figure} figures/synthetic_test.png
:alt: Top left, the synthetic magnetic record over 72 hours between -330 and +310 nT. Bottom left, the analytic electric field in red between about -17 and +17 mV/km and the SCUBAS field in blue about twice as large, with ringing at the ends. Right, a scatter of the two with a fitted slope of 2.00 and a correlation of 0.99992.
:width: 100%

The six-sine test of C2022 Section 2.4.2 on 72 hours of one-minute data (the tapered 10 % at each end is left out of the comparison, as in C2022). The correlation is {{ synth_corr }}, as in the paper, but **the field that `to_Efields` returns is {{ synth_slope }} times the analytic field of C2022 Eq 25** (see [paper and code](#paper-and-code)).
```

## 5. The emf of a section and of the cable

The loop formed by the cable and the seafloor next to it encloses no area, so Faraday's law gives no flux change through it, and the electric field along the cable is the electric field on the seafloor (B2024 Eqs 8 to 10; the approach of Cox et al., 1970). The emf of section $j$, with the field taken as uniform within the section, is (C2022 Eq 3; B2024 Eq 19)

$$
\mathcal{E}_j = E_j^N L_j^N + E_j^E L_j^E ,
$$

and for a cable of $S$ sections (C2022 Eq 26; B2024 Eq 20)

$$
\mathcal{E}_C = \sum_{j=1}^{S}\mathcal{E}_j .
$$

In the code this sum is the column `V(v)` of `Cable.tot_params`, in volts (the field is in mV/km and the length in km). Each section has the field of the observatory chosen for it, so a cable that crosses several time zones is driven by different records along its length.

## 6. The voltage on the power feed

```{figure} figures/voltage_assembly.png
:alt: Sketch of a cable between coast A and coast B under the sea. A red arrow along the cable shows the induced field and emf; dotted lines mark the Earth potentials U_A and U_B at the two coasts. The title is V_C equals emf plus U_A minus U_B.
:width: 80%

The two parts of the cable voltage: the emf induced in the cable by the electric field along it, and the difference of the Earth potentials at the two ends, where the power-feed equipment is grounded.
```

The power-feed equipment at each end is connected to the local Earth, which has the potential $U_A$ or $U_B$ produced by the same disturbance (the geoelectric coast effect; [the DSTL page](dstl.md)). The cable therefore sees (B2023 Eq 60)

$$
V_C = \mathcal{E}_C + U_A - U_B ,
$$

which is `Vt(v) = U0 - U1 + V(v)` in `Cable.tot_params`, with `U0` and `U1` the end potentials $U_A$ and $U_B$. The two terms oppose each other: for a cable across an ocean the potential difference between the coasts acts against the emf (B2023, Discussion), which is why the end potentials reduce the effect of the disturbance on the power feed. A worked example is in the [DSTL page](dstl.md#the-network): for the 100 m deep, 600 km wide shallow sea with 0.3 V/km everywhere, the emf is 180 V and the two Earth potentials are {{ net_shallow_end }} V in magnitude and of opposite sign, so $V_C = 180 - 2 \times$ {{ net_shallow_end }} $= 83.6$ V.

## 7. Verification

The model is checked against closed forms, at three levels.

1. **The transfer function** against BP2003 (Eqs 10, 11 and 13 and their figures) and against C2022 Table 4 (above).
2. **The potentials** of a single section against the long and short limits, and a three-section network against the examples of B2023 (its Figures 10 to 12, Table 2 and Eqs 56 to 59; [DSTL page](dstl.md#validation)).
3. **The chain**: the six-sine test above, and the 1989 event ([tutorial](../tutorial/1989.md)).

`tests/test_docs_model.py` runs these comparisons: formulas, the tables of the papers, and spot values read off the figures of the papers, with the tolerance of the figure's axes. C2022 reports a correlation of $1 - 10^{-6}$ for the field and $1 - 10^{-9}$ for the potentials; the code reproduces the shape of the field to {{ synth_corr }} and the potentials to the precision of the printed tables.

### How well it works

* **TAT-8, March 1989.** Calculated voltage variations of 67, 249, 406 and 763 V against about 75, 300, 450 and 700 V reported by Medford et al. (1989) for four events (B2024, Table 4); the authors caution that the largest events are those with uniform disturbances, which land observatories represent best.
* **TAT-1 (1958) and TAT-8 (1989), digitised data.** A median percentage error below 30 % and a skill score above 0.85 (Chakraborty et al., 2026, abstract). The model does less well for fast, wave-like variations and when the structure of the current systems is not resolved; the data are digitised from printed figures, which adds error, and records of cable voltages are scarce because operators hold them privately.

## Paper and code

Where the code does something the papers do not say, or differs from them, the difference is listed here and not smoothed over. None has been changed: the 0.1.9 behaviour is pinned by the [compatibility suite](../dev/testing.md) and a change is a separate, approved step ([known defects](../dev/testing.md#known-defects)).

| What | The paper | The code (0.1.9, unchanged) |
|---|---|---|
| Scale of the field | C2022 Eq 25 and Figure 5: the numerical field equals the analytic field (correlation $1 - 10^{-6}$, slope 1 in the plot of the two) | `OceanModel.to_Efields` returns **twice** the analytic field (slope {{ synth_slope }} in the test above). `scubas.utils.ifft` multiplies by `2 * n`; the factor was `n` until the commit of 2023-02-20, so `ifft(fft(x))` is `2 x`. Pinned in `tests/test_docs_model.py`. Whether the factor is intended is for the authors to decide |
| Taper | "Remove mean and trend and taper the ends" (C2022 Section 2.2) | The baseline is the median of the first 120 samples, not a mean or a trend. The cosine ramp is computed from the elapsed time in seconds divided by the number of samples, so it is a smooth ramp only for one-second data; for one-minute data it oscillates between 0 and 1 over the first and last tenth of the record. The default taper is `p = 0.1` |
| Seawater conductivity | 0.25 Ω m (BP2003); 0.3 Ω m (C2022 Table 2, B2023 Table 1); 4 S/m, that is 0.25 Ω m (B2024, Section 4.1); 3 S/m (Chakraborty et al., 2026, Section 2.3) | The built-in profiles use 0.3 Ω m. The papers' values differ by a factor of up to 1.3 and the transfer function depends on them, so a result should state which was used |
| Section lengths | $\delta\lambda$ and $\delta\phi$ are differences (C2022 Eqs 4 and 5) | Absolute differences: a section laid north-east and the same section laid south-west have the same lengths, and the sign of the emf is set only by the field |
| Giving only a length | not described | `length` is split equally between north and east (45°) |
| End potentials | exact | Rounded to 0.01 V, and `Vt(v)` is computed from the rounded values (known defect) |
| Layer thickness | km in the tables | stored in metres in `Layer` although the docstrings say km (known defect); `scripts/model_cases.py` converts |

## Assumptions and limits

| Assumption | Consequence | Where stated |
|---|---|---|
| Each section has a layered Earth, 1-D in depth; lateral change is only between sections (the piecewise approach of Marti et al., 2014) | Transverse gradients and edge effects at abrupt bathymetry are not resolved | C2022 Section 5.1; B2024 Section 3.1 |
| The seafloor of a section is flat | Detail of the bathymetry inside a section is lost; sections are chosen by depth and length | C2022 Section 5.1 |
| The field of a section is that of the nearest observatory, with no interpolation | Best when the disturbance is coherent along the route; small-scale currents and the lack of magnetometers over the ocean limit accuracy | C2022 Section 5.1; B2024 Section 5 |
| A vertical plane wave (a source uniform over the section) | The method is for variations whose source is large compared with the water depth and the skin depth | BP2003, Theory |
| The thin-sheet condition: periods long enough to penetrate the water and the crust | Short periods have a smaller skin depth, and other methods (finite elements) are needed | B2023, Discussion |
| The result is the induced voltage | It is not the response of the power-feed equipment, which holds the current constant | B2024, Section 4.4 |

## Related models

| Approach | What it does | Relation to SCUBAS |
|---|---|---|
| Seafloor field from the surface impedance | $E$ at the seafloor from $B$ at the sea surface for a layered Earth | The first step of SCUBAS, section by section (BP2003) |
| Analytic cable-voltage theory (Mecozzi, 2022, arXiv:2204.00560) | Voltage proportional to the magnitude of the field change, from the attenuation in the sea and the reflection at the bottom; design criteria for new systems | A closed-form treatment; SCUBAS is numerical, section by section, with an Earth profile and the Earth potentials at the ends |
| Power-system GIC networks | Induced emfs in lines and the nodal network of the grid | SCUBAS carries the same structure to cables: transmission-line sections, equivalent-$\pi$ circuits, an admittance matrix |
| Three-dimensional models (finite element, finite difference) | Solve Maxwell's equations on a 3-D conductivity model, including the coast effect | More complete and far heavier; SCUBAS is not meant to replace them ([the numerical coast model](fdtd.md)) |

## Where each step is in the code

| Step | Code |
|---|---|
| Magnetic record | `OceanModel.read_Bfield_data`, `OceanModel.read_iaga` |
| Earth model | `scubas.datasets` (`Layer`, `Site`, `PROFILES`), `scubas.conductivity.ConductivityProfile` |
| Surface impedance of the layers | `Site.calcZ` |
| Transfer function | `OceanModel.calcZ`, `OceanModel.calcTF`, `OceanModel.get_TFs` |
| Detrend and taper | `Preprocess` |
| Spectra | `scubas.utils.fft`, `scubas.utils.ifft` |
| Electric field | `OceanModel.to_Efields` |
| Section geometry and emf | `CableSection`, `TransmissionLine.compute_Vj` |
| Line constants, π circuit, terminations, network | `TransmissionLine`, `Cable` ([DSTL page](dstl.md)) |
| Cable voltage | `Cable.tot_params` (`V(v)`, `U0`, `U1`, `Vt(v)`) |

## Sources and what was checked

* **Read in full for this page:** Boteler and Pirjola (2003); Chakraborty et al. (2022); Boteler et al. (2023); Boteler et al. (2024); Chakraborty et al. (2026); Mecozzi (2022, abstract and introduction).
* **Checked against the code by test:** BP2003 Eqs 3 to 5, 8, 10, 11, 13, 15, 16, and the values read off its Figures 1, 4 and 5; C2022 Eqs 3 to 5, 10 to 12, 24, 25, 27 and Tables 3 and 4; B2023 Eqs 1 to 5, 19, 27 to 34, 38 to 42, 43 to 55, 56 to 60, Tables 2 and 3, and Figures 10 to 12.
* **Cited by the papers, not read:** Horton et al. (2012), Weaver (1994), Marti et al. (2014), Cox et al. (1970), Medford et al. (1989), Pasyanos et al. (2014).

## References

* Boteler, D. H., and Pirjola, R. J. (2003). Magnetic and electric fields produced in the sea during geomagnetic disturbances. *Pure and Applied Geophysics* 160, 1695-1716. doi:10.1007/s00024-003-2372-6
* Boteler, D. H., Chakraborty, S., Shi, X., Hartinger, M. D., and Wang, X. (2023). Transmission line modelling of geomagnetic induction in the ocean/Earth conductivity structure. *International Journal of Geosciences* 14, 767-791. doi:10.4236/ijg.2023.148041
* Boteler, D. H., Chakraborty, S., Shi, X., Hartinger, M. D., and Wang, X. (2024). An examination of geomagnetic induction in submarine cables. *Space Weather* 22, e2023SW003687. doi:10.1029/2023SW003687
* Chakraborty, S., Boteler, D. H., Shi, X., Murphy, B. S., Hartinger, M. D., Wang, X., Lucas, G., and Baker, J. B. H. (2022). Modeling geomagnetic induction in submarine cables. *Frontiers in Physics* 10, 1022475. doi:10.3389/fphy.2022.1022475
* Chakraborty, S., et al. (2026). Validating SCUBAS predictions of geomagnetically induced voltage in submarine cables using legacy superstorm observations. *Space Weather*. doi:10.1029/2025SW004843
* Horton, R., Boteler, D., Overbye, T. J., Pirjola, R., and Dugan, R. C. (2012). A test case for the calculation of geomagnetically induced currents. *IEEE Transactions on Power Delivery* 27, 2368-2373.
* Mecozzi, A. (2022). arXiv:2204.00560.
* Pasyanos, M. E., Masters, T. G., Laske, G., and Ma, Z. (2014). LITHO1.0: An updated crust and lithosphere model of the Earth. *Journal of Geophysical Research: Solid Earth* 119, 2153-2173.
